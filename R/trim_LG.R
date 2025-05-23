#' Trim Linkage Groups Based on LOD Thresholds
#'
#' `trim_LG` is an interactive or scripted function that helps the user filter markers within a specified chromosome
#' based on linkage group (LG) assignment using LOD score thresholds. The user can interactively choose LOD thresholds
#' or provide them via arguments, visualize haplotype frequency before and after filtering, and remove outliers.
#'
#' Inspired by functions by Professor Jeffrey B. Endelman’s [MapRtools](https://github.com/jendelman/MapRtools).
#'
#' @param chromosome Character. The chromosome ID to be processed.
#' @param map A data frame with at least columns `"marker"`, `"chrom"`, and `"position"`.
#' @param geno A genotype matrix (rows = markers, columns = individuals). Preferably binned.
#' @param pop_type Character. Population type for LOD estimation. Default is `"F2"`.
#' @param drop_outliers Logical. Whether to remove outliers based on haplotype frequency. Default is `TRUE`.
#' @param n_cores Integer. Number of cores to use. If `NULL`, uses maximum available minus one.
#' @param interactive Logical. If `TRUE`, asks user for LOD thresholds. If `FALSE`, requires threshold arguments.
#' @param initial_LOD Numeric. Initial LOD value (required if `interactive = FALSE`).
#' @param end_LOD Numeric. Final LOD value (required if `interactive = FALSE`).
#' @param step_LOD Numeric. Step size for LOD sequence (required if `interactive = FALSE`).
#' @param selected_thresh Numeric. Final LOD threshold to use (required if `interactive = FALSE`).
#'
#' @return A list with filtered genotypes, final map, plots, and LOD settings.
#' @export
#'
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggtitle
#' @importFrom parallel detectCores
#'
trim_LG <- function(chromosome,
                    map,
                    geno,
                    pop_type = "F2",
                    drop_outliers = TRUE,
                    n_cores = NULL,
                    interactive = TRUE,
                    initial_LOD = NULL,
                    end_LOD = NULL,
                    step_LOD = NULL,
                    selected_thresh = NULL) {

  if (!requireNamespace("MapRtools", quietly = TRUE)) {
    stop("The 'MapRtools' package is required but not installed.\nInstall it with:\n  remotes::install_github('jendelman/MapRtools')", call. = FALSE)
  }

  # Detect CPU cores
  max_cores <- parallel::detectCores(logical = FALSE)
  usable_cores <- max_cores - 1
  if (is.null(n_cores) || n_cores > max_cores) {
    n_cores <- usable_cores
  }

  map <- dplyr::arrange(map, chrom, position)
  CHR <- as.character(chromosome)
  message("Processing chromosome: ", CHR)

  chr_markers <- map$marker[map$chrom == CHR]
  if (length(chr_markers) == 0) stop("No markers found for the specified chromosome.")

  chr_geno <- geno[chr_markers, ]

  # Initial haplotype frequency plot
  message("Plotting initial haplotype frequency...")
  initial_haplo <- MapRtools::plot_haplofreq(chr_geno)
  print(initial_haplo$plot + ggtitle(paste0("Chromosome: ", CHR)))

  # Estimate LOD matrix
  message("Estimating LOD matrix using ", n_cores, " cores...")
  LODmat_chr <- MapRtools::MLEL(geno = chr_geno, pop.type = pop_type, LOD = TRUE, n.core = n_cores)

  # Interactive or scripted input
  if (interactive) {
    initial_LOD <- as.numeric(readline("Enter the initial LOD threshold (e.g., 15): "))
    end_LOD <- as.numeric(readline("Enter the end LOD threshold (e.g., 45): "))
    step_LOD <- as.numeric(readline("Enter the step size for LOD sequence (e.g., 1): "))
    selected_thresh <- as.numeric(readline("Enter the LOD threshold to use for filtering markers: "))
  }

  # Validate inputs
  if (any(is.null(c(initial_LOD, end_LOD, step_LOD, selected_thresh))) ||
      initial_LOD >= end_LOD || step_LOD <= 0 ||
      selected_thresh < initial_LOD || selected_thresh > end_LOD) {
    stop("Invalid or missing LOD thresholds. Provide valid values or use interactive = TRUE.")
  }

  # Generate candy stripe plot
  message("Generating linkage group plot...")
  candy_plot <- MapRtools::LG(LODmat_chr, thresh = seq(initial_LOD, end_LOD, by = step_LOD))
  print(candy_plot + ggtitle(paste0("Chromosome: ", CHR)))

  # Filter based on selected threshold
  message("Filtering markers at LOD threshold = ", selected_thresh)
  LG_result <- MapRtools::LG(LODmat = LODmat_chr, thresh = selected_thresh)
  filtered_LG <- LG_result[LG_result$LG == 1, ]

  map_chr_filtered <- merge(map, filtered_LG, by = "marker")
  map_chr_filtered <- map_chr_filtered[order(map_chr_filtered$position), ]

  # Remove outliers if required
  outliers <- character(0)
  if (drop_outliers) {
    message("Removing outliers...")
    freq_plot <- MapRtools::plot_haplofreq(chr_geno)
    outliers <- freq_plot$outliers
    chr_geno <- chr_geno[!rownames(chr_geno) %in% outliers, ]
    map_chr_filtered <- map_chr_filtered[!map_chr_filtered$marker %in% outliers, ]
  }

  return(list(
    trimmed_genotype = chr_geno,
    final_map = map_chr_filtered,
    initial_haplo_plot = initial_haplo$plot,
    filtered_freq_plot = MapRtools::plot_haplofreq(chr_geno)$plot,
    selected_LOD_thresh = selected_thresh,
    initial_LOD = initial_LOD,
    ending_LOD = end_LOD,
    step = step_LOD,
    outliers = outliers
  ))
}
