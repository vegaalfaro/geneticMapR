#' Estimate Genetic Map Positions
#'
#' Calculates genetic map positions for a single chromosome using recombination frequency (RF)
#' and LOD scores from a genotype matrix. The output is standardized for compatibility with Marey plots.
#'
#' @param geno_matrix A genotype matrix (markers as rows, individuals as columns) for one chromosome or linkage group.
#' @param model Character. Mapping function to use: `"Kosambi"` (default) or `"Haldane"`.
#' @param n_point Integer. The `n.point` value for interpolation when estimating genetic positions. Default is `20`.
#' @param pop.type Character. Population type: `"DH"`, `"BC"`, `"F2"` (default), `"S1"`, `"RIL.self"`, `"RIL.sib"`.
#'
#' @return A data frame with columns:
#'   - `chrom`: Chromosome identifier
#'   - `position_Mb`: Physical position in megabases
#'   - `position_cM`: Genetic distance in centimorgans (estimated)
#'
#' @details
#' This function uses `MapRtools::MLEL()` and `MapRtools::genetic_map()` for estimating genetic positions.
#' The output is ready to be used with `plot_marey()`.
#'
#' @examples
#' \dontrun{
#' map_chr1 <- estimate_map(geno_matrix = chr1_geno_matrix, n_point = 20)
#' plot_marey(map_chr1, chrom = "CHR1")
#' }
#' @export
estimate_map <- function(geno_matrix,
                         model = "Kosambi",
                         n_point = 20,
                         pop.type = "F2") {

  if (!requireNamespace("MapRtools", quietly = TRUE)) {
    stop("The 'MapRtools' package is required. Install it from: https://github.com/jendelman/MapRtools")
  }

  # Validate inputs
  if (!(model %in% c("Kosambi", "Haldane"))) {
    stop("Invalid model. Choose either 'Kosambi' or 'Haldane'.")
  }

  if (!(pop.type %in% c("DH", "BC", "F2", "S1", "RIL.self", "RIL.sib"))) {
    stop("Invalid pop.type. Choose from: 'DH', 'BC', 'F2', 'S1', 'RIL.self', 'RIL.sib'.")
  }

  # Extract map (physical positions)
  map <- extract_map(geno_matrix)

  # Compute recombination frequency (RF) and LOD matrices
  RFmat <- MapRtools::MLEL(geno = geno_matrix, pop.type = pop.type, LOD = FALSE)
  LODmat <- MapRtools::MLEL(geno = geno_matrix, pop.type = pop.type, LOD = TRUE)

  # Convert RF to genetic distances
  x <- apply(RFmat, c(1, 2), FUN = MapRtools::map_fn, model = model)

  # Estimate genetic positions
  gen_map <- MapRtools::genetic_map(x = x, LOD = LODmat, n.point = n_point)

  # Final output
  out <- map
  out$position_Mb <- out$position / 1e6
  out$position_cM <- gen_map$position

  return(out[, c("chrom", "position_Mb", "position_cM")])
}
