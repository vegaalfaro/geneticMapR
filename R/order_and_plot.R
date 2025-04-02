#' Optimize Marker Order
#'
#'@description
#' `order_and_plot` optimizes the marker order in a (trimmed) genotype matrix using
#' recombination frequency (RF) matrices. It performs multiple iterations of
#' marker ordering, selects the best order based on the smallest Sum of Adjusted
#' Recombination Fractions (SARF), and generates visualizations of marker order. The function uses
#' `MLEL()` and `order_markers`() from Professor Jeffrey Endelman's [MapRtools](https://github.com/jendelman/MapRtools) and streamlines the process.
#' This workflow should be easier on the user when dealing with multiple chromosomes.
#' [MapRtools](https://github.com/jendelman/MapRtools).
#'
#' @param trimmed_geno A genotype matrix where:
#'   - Rows represent markers.
#'   - Columns represent individuals.
#'   It could be the output of `trim_LG` or any genotype matrix.
#' @param pop.type Character. The population type used for ordering markers.
#'   Default is `"F2"`. `order_and_plot` would work with any of the following "DH","BC","F2","S1","RIL.self","RIL.sib". Based on `MapRtools::MLEL()`
#' @param CHR Character (optional). The chromosome identifier for labeling plots.
#' @param n.iter Integer. The number of iterations for marker ordering. Default is `6`.
#' @param prop Numeric. The proportion of individuals to include in the genotype plot.
#'   Must be between `0` and `1`. Default is `0.20` (20% of individuals). If all individuals are included the visualization quality drops.
#'
#' @return A list containing:
#'   - `"original_geno"`: The original genotype matrix.
#'   - `"ordered_geno"`: The optimized genotype matrix with markers reordered.
#'   - `"order_plot"`: A `ggplot2` object displaying the original vs. optimized marker order.
#'   - `"haplotype_plot"`: A `ggplot2` object showing the genotype haplotype plot.
#'
#' @details
#' - Computes an RF matrix using `MLEL()`, then orders markers iteratively.
#' - Runs `order_markers()` `n.iter` times and selects the order with the smallest SARF.
#' - `order_markers()` Order markers by solving the traveling salesperson problem.
#' - Plots the original vs. optimized marker order.
#' - Generates a genotype haplotype plot for a subset of individuals (`prop`).
#' - If `CHR` is provided, adds a chromosome label to the plots.
#'
#' @note This function was refined with assistance from ChatGPT to improve clarity,
#' efficiency, and visualization formatting.
#'
#' @importFrom ggplot2 ggplot geom_point facet_wrap theme_bw xlab ylab theme element_text element_blank ggtitle
#' @importFrom dplyr mutate arrange select
#' @export
order_and_plot <- function(trimmed_geno,
                           pop.type = "F2",
                           CHR = NULL,
                           n.iter = 6,
                           prop = NULL) {

  # Set default value for prop if not specified
  if (is.null(prop)) {
    prop <- 0.20
  }

  # Remind user that the functions rely on MapRtools
  if (!requireNamespace("MapRtools", quietly = TRUE)) {
    stop("The 'MapRtools' package is required but not installed.\nInstall it with:\n  remotes::install_github('jendelman/MapRtools')", call. = FALSE)
  }

  # Ensure prop is numeric and within valid range
  if (!is.numeric(prop) || prop <= 0 || prop > 1) {
    stop("Error: prop must be a numeric value between 0 and 1.")
  }

  # Store original order
  original_geno <- as.matrix(trimmed_geno)

  # Estimate RF matrix
  RFmat <- MapRtools::MLEL(trimmed_geno, pop.type = pop.type, LOD = FALSE)

  # Number of markers and individuals
  m <- nrow(trimmed_geno)
  n <- ncol(trimmed_geno)

  # Vector to store results
  vec <- vector("list", n.iter)

  # Loop that runs the ordering function n.iter times
  for (i in 1:n.iter) {
    vec[[i]] <- MapRtools::order_markers(RFmat)
  }

  # Extract the SARF values
  SARF <- sapply(vec, "[[", "distance")
  y <- sapply(lapply(vec, "[[", "order"), function(x) { match(1:m, x) })

  # Convert SARF values to a unique identifier
  SARF_labels <- paste0("SARF=", round(SARF, 3), "_iter", 1:n.iter)

  # Create dataframe for plotting
  plot.data <- data.frame(
    x = rep(1:m, n.iter),
    y = as.vector(y),
    SARF = factor(rep(SARF_labels, each = m))
  )

  # Generate the order plot
  order_plot <- ggplot(data = plot.data, aes(x = x, y = y)) +
    geom_point(size = 1.2, color = "black", alpha = 0.6) +  # Subtle points
    facet_wrap(~SARF, ncol = 2) +  # Keep facets organized
    theme_bw(base_size = 14) +  # Classic black n' white theme
    xlab("Original Marker Order") +
    ylab("Optimized Marker Order") +
    theme(
      strip.text = element_text(size = 11),  # Facet titles stand out
      panel.grid.major = element_blank(),  # Remove grid for a cleaner look
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12)
    )

  # Add chromosome title if provided
  if (!is.null(CHR)) {
    order_plot <- order_plot + ggtitle(paste("Chromosome:", CHR))
  }

  # Select the best order based on minimum SARF
  k <- which.min(SARF)
  ordered_geno <- as.matrix(trimmed_geno[vec[[k]]$order, ])

  # Extract map and plot haplotypes of ordered genotype matrix
  map <- extract_map(ordered_geno, markers = TRUE)

  # Determine number of individuals to plot
  num_individuals <- max(1, round(n * prop))  # Ensures at least 1 individual

  # Print the order plot
  print(order_plot)

  # Generate the genotype plot
  haplotype_plot <- MapRtools::plot_geno(ordered_geno[, 1:num_individuals], map = map)

  # Add chromosome title
  if (!is.null(CHR)) {
    haplotype_plot <- haplotype_plot + ggtitle(paste("Chromosome:", CHR))
  }

  # Return results, including the original genotype matrix
  return(list(original_geno = as.data.frame(original_geno),
              ordered_geno = as.data.frame(ordered_geno),
              order_plot = order_plot,
              haplotype_plot = haplotype_plot))
}



