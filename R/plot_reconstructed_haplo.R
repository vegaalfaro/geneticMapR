#' Plot Haplotype Reconstruction Results Before and After HMM
#'
#' `plot_reconstructed_haplo` visualizes the genotype data before and after HMM-based haplotype reconstruction
#' for a specified set of individuals. It uses `MapRtools::plot_geno()` for visualization
#' and appends informative titles based on the chromosome name.
#'
#' @param original A genotype matrix before haplotype reconstruction. Rows represent markers,
#'        and columns represent individuals.
#' @param processed A genotype matrix after HMM-based haplotype reconstruction. It must have
#'        the same dimensions as `original`.
#' @param individuals A numeric vector specifying the indices of individuals (columns) to plot.
#'        Default is `1:50`. All values must be within the column range of the input matrices.
#'
#' @return A list with two ggplot objects:
#' \describe{
#'   \item{p1}{Genotype plot before HMM, titled with the chromosome name and "Before HMM".}
#'   \item{p2}{Genotype plot after HMM, titled with the chromosome name and "After HMM".}
#' }
#'
#' @details The function uses `extract_map()` from `geneticMapR` to infer the chromosome name from the original
#'         genotype matrix, which is then used to annotate the plots. Both matrices are expected
#'         to be compatible with `MapRtools::plot_geno()`.
#'
#' @examples
#' \dontrun{
#' original_geno <- example_original
#' processed_geno <- example_processed
#' plots <- plot_reconstructed_haplo(original = original_geno,
#'                                   processed = processed_geno,
#'                                   individuals = 1:30)
#' plots$p1  # View plot before HMM
#' plots$p2  # View plot after HMM
#' }
#'
#' @importFrom ggplot2 ggtitle
#' @export
# Function to plot haplotype reconstruction results
plot_reconstructed_haplo <- function(original, processed, individuals = 1:50) {

  # Remind users that the functions rely on MapRtools
  if (!requireNamespace("MapRtools", quietly = TRUE)) {
    stop("The 'MapRtools' package is required but not installed.\nInstall it with:\n  remotes::install_github('jendelman/MapRtools')", call. = FALSE)
  }

  # Extract chromosome name
  title_prefix <- extract_map(original)
  title_prefix <- unique(title_prefix$chrom)

  # Ensure individuals provided are within the valid range
  n_individuals <- ncol(original)
  if (any(individuals > n_individuals) || any(individuals < 1)) {
    stop("Invalid individual range. Ensure selected individuals are within the available columns.")
  }

  # Generate plots
  p1 <- MapRtools::plot_geno(geno = as.matrix(original[, individuals])) +
    ggtitle(paste(title_prefix, "Before HMM"))

  p2 <- MapRtools::plot_geno(geno = as.matrix(processed[, individuals])) +
    ggtitle(paste(title_prefix, "After HMM"))

  return(list(p1 = p1, p2 = p2))
}
