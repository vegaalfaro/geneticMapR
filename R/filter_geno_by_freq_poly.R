
#' Filter Genotype Matrix by Dosage Frequency (Polyploid-Compatible)
#'
#' @description
#' Filters a polyploid genotype matrix based on dosage frequencies across markers.
#' Allows users to define what dosages count as heterozygous and to apply constraints
#' on maximum/ minimum genotype frequencies and heterozygous frequency ranges.
#'
#' @param geno_matrix A numeric genotype matrix or data frame where:
#'   - Rows are markers
#'   - Columns are individuals
#'   - Values are dosage values (e.g., 0 to 4 for tetraploids)
#' @param max_geno_freq Numeric. If provided, removes markers where the most frequent dosage exceeds this value.
#' @param het_freq_range Numeric vector of length 2. Keeps markers where heterozygosity frequency is within range.
#' @param min_geno_freq Numeric. If provided, removes markers where the least frequent dosage falls below this value.
#' @param het_dosages Integer vector of dosage values to treat as "heterozygous" (e.g., c(1, 2, 3)).
#'
#' @return A filtered genotype matrix containing only markers that meet all criteria.
#'
#' @details
#' - Designed for polyploid dosage matrices (not raw genotype strings).
#' - Computes dosage frequencies per marker (row-wise).
#' - `het_dosages` is used to define heterozygosity in a generalizable way.
#'
#' In diploids, it's simple:
#'
#' 0 = homozygous reference (e.g., "AA")
#'
#' 1 = heterozygous (e.g., "AB")
#'
#' 2 = homozygous alternate (e.g., "BB")
#'
#' So "1" is always the heterozygous state.
#'
#' But in polyploids (e.g., tetraploids, hexaploids), we can have more intermediate states. For example, in a tetraploid:
#'
#'  0 = "AAAA" → **homozygous ref**
#'
#' 1 = "AAAB"
#'
#' 2 = "AABB"
#'
#' 3 = "ABBB"
#'
#' 4 = "BBBB" → **homozygous alt**
#'
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' test_geno <- matrix(sample(0:4, 1000, replace = TRUE), nrow = 100, ncol = 10)
#' rownames(test_geno) <- paste0("Marker", 1:100)
#' colnames(test_geno) <- paste0("Ind", 1:10)
#'
#' # Filter markers with max dosage freq < 0.7 and heterozygote freq between 0.2 and 0.8
#' filtered <- filter_geno_by_freq_poly(
#'   test_geno,
#'   max_geno_freq = 0.7,
#'   het_freq_range = c(0.2, 0.8),
#'   min_geno_freq = 0.05,
#'   het_dosages = c(1, 2, 3)
#' )
#'
#' dim(filtered)  # Number of markers retained
#'}
#' @importFrom stats na.omit
#' @export
filter_geno_by_freq_poly <- function(geno_matrix,
                                     max_geno_freq = NULL,
                                     het_freq_range = NULL,
                                     min_geno_freq = NULL,
                                     het_dosages = c(1, 2, 3)) {

  # If no filtering is requested, return original
  if (is.null(max_geno_freq) & is.null(het_freq_range) & is.null(min_geno_freq)) {
    return(geno_matrix)
  }

  # Compute dosage frequencies per marker
  dosage_levels <- sort(unique(as.vector(na.omit(geno_matrix))))
  freq_table <- t(apply(geno_matrix, 1, function(row) {
    prop.table(table(factor(row, levels = dosage_levels)))
  }))
  colnames(freq_table) <- paste0("dosage_", dosage_levels)

  # Compute max and min genotype frequency
  max_freq <- apply(freq_table, 1, max)
  min_freq <- apply(freq_table, 1, min)

  # Compute heterozygote frequency
  het_cols <- paste0("dosage_", het_dosages)
  het_freq <- rowSums(freq_table[, colnames(freq_table) %in% het_cols, drop = FALSE])

  # Initialize filter mask
  keep <- rep(TRUE, nrow(geno_matrix))

  if (!is.null(max_geno_freq)) {
    keep <- keep & (max_freq < max_geno_freq)
  }

  if (!is.null(min_geno_freq)) {
    keep <- keep & (min_freq >= min_geno_freq)
  }

  if (!is.null(het_freq_range)) {
    keep <- keep & (het_freq >= het_freq_range[1]) & (het_freq <= het_freq_range[2])
  }

  # Filter matrix
  filtered_matrix <- geno_matrix[keep, , drop = FALSE]
  return(filtered_matrix)
}
