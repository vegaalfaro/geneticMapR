#' Estimate Genetic Map Positions
#'
#' Designed to work with one chromosome at a time. `estimate_map` calculates genetic map positions using recombination frequency (RF)
#' and LOD scores calculated from a genotype matrix. It applies a selected mapping function
#' (Kosambi or Haldane) and outputs map positions (in centi Morgans, cM) at multiple resolutions defined by the user.
#'
#' @param geno_matrix A genotype matrix where rows represent markers and columns represent individuals.
#'        This matrix must be in a format compatible with `MapRtools::MLEL` and with markers for only one chromosome or Linkage Group.
#' @param model Character string specifying the mapping function to use. Options are `"Kosambi"` (default)
#'        or `"Haldane"`.
#' @param n_points A numeric vector indicating the different values of `n.point` to be used
#'        for estimating genetic positions. Default is `seq(10, 20, by = 5)`.
#' @param pop.type Character string specifying the population type. Supported options are
#'        `"DH"`, `"BC"`, `"F2"` (default), `"S1"`, `"RIL.self"`, and `"RIL.sib"`.
#'
#' @return A data frame containing the extracted genetic map with additional columns (`pN`)
#'         corresponding to the estimated positions for each specified `n.point` value.
#'
#' @details The function uses `MapRtools::MLEL()` to compute RF and LOD matrices and then
#'         applies the specified mapping function via `MapRtools::map_fn()`. For each `n.point`
#'         value, the function appends a new column (`pN`) to the map with the estimated genetic positions.
#'
#' @examples
#' \dontrun{
#' # Example genotype matrix
#'
#' geno_matrix <- chr1_geno_matrix
#'
#' map <- estimate_map(chr1_geno_matrix, model = "Kosambi", n_points = c(10, 15, 20), pop.type = "F2")
#'
#' # Simply
#' genetic_map_chr1 <- estimate_map(chr1_geno_matrix)
#'}
#' @export

estimate_map <- function(geno_matrix,
                         model = "Kosambi",
                         n_points = seq(10, 20, by = 5),
                         pop.type = "F2") {

# Validate model selection
if (!(model %in% c("Kosambi", "Haldane"))) {
  stop("Invalid model. Choose either 'Kosambi' or 'Haldane'.")
}

# Validate population selection
if (!(pop.type %in% c("DH", "BC", "F2", "S1", "RIL.self", "RIL.sib"))) {
  stop("Invalid population type. Choose one of the following: 'DH', 'BC', 'F2', 'S1', 'RIL.self', 'RIL.sib'.")
}

# Extract genetic map from genotype matrix
map <- extract_map(geno_matrix)

# Compute recombination frequency (RF) and LOD scores
RFmat <- MapRtools::MLEL(geno = geno_matrix, pop.type = pop.type, LOD = FALSE)
LODmat <- MapRtools::MLEL(geno = geno_matrix, pop.type = pop.type, LOD = TRUE)

# Estimate genetic positions using the chosen mapping function
x <- apply(RFmat, c(1, 2), FUN = MapRtools::map_fn, model = model)

# Generate genetic positions for different n.point values
for (n in n_points) {
  column_name <- paste0("p", n)
  map[[column_name]] <- MapRtools::genetic_map(x = x, LOD = LODmat, n.point = n)$position
}

return(map)
}

