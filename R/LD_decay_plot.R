#' Plot LD vs distance
#'
#' Plot LD vs distance using an asymptotic regression model.
#'
#' An asymptotic regression model is fit using \link[stats]{SSasymp}.
#' The distance where r² reaches a specified threshold (e.g., 0.1 or 0.2) is calculated and returned.
#' This function calculates r² by chromosome and  follows the general logic of [MapRtools::LD.plot](https://github.com/jendelman/MapRtools/blob/master/R/plot_LD.R) but uses the function
#' \link[stats]{SSasymp} to model the non-linear LD decay. This function will give you a
#' quick and simple estimate of LD. For more sophisticated LD functions see [David Gerard](https://github.com/dcgerard) [ldfast()](https://cran.r-project.org/web/packages/ldsep/vignettes/fast.html)
#'
#' @param data Genotype matrix (markers as columns, individuals as rows). Rownames have the marker names that matcho those of the marker column in the `map`.
#' @param map Data frame with columns 'chrom', 'marker', and 'position'
#' @param max.pair Maximum number of r² pairs for the model (default: 1e4)
#' @param max.loci Maximum number of markers to use per chromosome (default: NULL)
#' @param position "bp" or "Mb" (default: "bp")
#' @param xlim_range Zoom in range on the x axis, default is c(0, 20). Set to NULL to see entire plot.
#' @param r2_threshold r² value for decay distance calculation (default: 0.2),
#' @param show_vline show a vertical line intersecting the x-axis at the the half-decay distance (default, TRUE)
#' @param show_hline show a horizontal line intersecting the y-axis at the selected r² threshold (default, TRUE)
#' @param seed Optional integer seed to make the random sampling reproducible. Default 123.
#'
#' @return A list containing:
#' \item{plot}{ggplot2 object with points and fitted curve (or NULL if fit fails)}
#' \item{half_decay_dist}{Distance where r² reaches the specified threshold (or NA if not estimable)}
#' \item{model}{Fitted nls object (or NULL if fitting failed)}
#'
#' @details
#' It processes each chromosome separately to compute all pairwise r²
#'  values within chromosomes, then it combines the results across all
#'   chromosomes into a single dataset, and fits one global LD decay model
#'   to the pooled data. No cross-chromosome marker pairs are considered.
#'
#' If the user is interested in LD per chromosome, provide a map and data file with only
#' data pertaining to the chromosome of interest.
#'
#' @note
#'  If the number of rows (number of marker pairs) is larger than max.pair (default is 1e4), we proceed to random sampling.
#'
#' @importFrom stats cor dist nls nls.control predict coef
#'
#' @export
LD_decay_plot <- function(data,
                          map,
                          max.pair = 1e4,
                          max.loci = NULL,
                          position = "bp",
                          r2_threshold = 0.2,
                          xlim_range = c(0, 20),
                          show_vline = TRUE,
                          show_hline = TRUE,
                          seed = 123)
{


  # Input validation
  if (!is.matrix(data)) stop("data must be a matrix")
  if (!all(c("chrom", "marker", "position") %in% colnames(map))) {
    stop("map must contain 'chrom', 'marker', and 'position' columns")
  }
  if (!is.numeric(r2_threshold) || r2_threshold <= 0 || r2_threshold >= 1) {
    stop("r2_threshold must be a number between 0 and 1")
  }

  chroms <- unique(map$chrom)
  result <- NULL

  for (chr in chroms) {
    chr_map <- map[map$chrom == chr, ]
    ix <- match(chr_map$marker, colnames(data))
    ix <- ix[!is.na(ix)]

    m <- length(ix)
    if (!is.null(max.loci) && m > max.loci) {
      ix <- sample(ix, max.loci)
      m <- max.loci
    }
    if (m < 2) next

    tmp <- expand.grid(col = 1:m, row = 1:m)
    tmp <- tmp[tmp$row >= tmp$col, ]  # lower triangle

    r2 <- cor(data[, ix, drop = FALSE], use = "pairwise.complete.obs")^2
    r2.vec <- r2[cbind(tmp$row, tmp$col)]

    pos <- chr_map$position[match(colnames(data)[ix], chr_map$marker)]
    d <- as.matrix(dist(matrix(pos, ncol = 1)))
    d.vec <- d[cbind(tmp$row, tmp$col)]

    if (position == "bp") {
      d.vec <- d.vec / 1e6
    }

    result <- rbind(result, data.frame(d = d.vec, r2 = r2.vec))
  }


  # Clean and  sample
  result <- result[!is.na(result$d) & !is.na(result$r2) & result$d > 0, ]
  if (nrow(result) == 0) {
    warning("No valid marker pairs with distance > 0") # Mesagge
    return(list(plot = NULL, half_decay_dist = NA, model = NULL))
  }

  if (nrow(result) > max.pair) {
    if (!is.null(seed)) set.seed(seed)
    sample_data <- result[sample(nrow(result), max.pair), ]
  } else {
    sample_data <- result
  }

  # Try fitting SSasymp model
  nls_fit <- tryCatch(
    stats::nls(r2 ~ SSasymp(d, Asym, R0, lrc),
               data = sample_data,
               control = stats::nls.control(maxiter = 200, warnOnly = TRUE)),
    error = function(e) {
      warning("nls fit failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(nls_fit)) {
    return(list(plot = NULL, half_decay_dist = NA, model = NULL))
  }

  # Predict and calculate half-decay
  dmax <- max(result$d)
  pred_df <- data.frame(d = seq(0, dmax, length.out = 500))
  pred_df$r2 <- stats::predict(nls_fit, newdata = pred_df)

  coefs <- coef(nls_fit)
  half_decay_dist <- tryCatch({
    -log((r2_threshold - coefs["Asym"]) / (coefs["R0"] - coefs["Asym"])) / exp(coefs["lrc"])
  }, error = function(e) NA_real_)

  if (is.nan(half_decay_dist) || half_decay_dist < 0 || half_decay_dist > dmax) {
    half_decay_dist <- NA
    warning("Half-decay distance could not be reliably estimated")
  }

  # Create plot
  unit_label <- switch(position,
                       bp = "Distance (Mb)",
                       Mb = "Distance (Mb)", "Distance")



  p <- ggplot(data = sample_data, aes(x = d, y = r2)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_line(data = pred_df, aes(x = d, y = r2), color = "red", linewidth = 1)

  if (!is.na(half_decay_dist) && show_vline) {
    p <- p + geom_vline(xintercept = half_decay_dist, linetype = 2, color = "red", linewidth = 0.5)
  }
  if (show_hline) {
    p <- p + geom_hline(yintercept = r2_threshold, linetype = 2, color = "red", linewidth = 0.5)
  }

  p <- p +
    labs(x = unit_label, y = expression(r^2)) +
    theme_bw()

  if (!is.null(xlim_range)) {
    p <- p + coord_cartesian(xlim = xlim_range)
  }


  return(list(plot = p,
              half_decay_dist = half_decay_dist,
              model = nls_fit))
}

