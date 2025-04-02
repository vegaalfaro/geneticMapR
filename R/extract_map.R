#' Create a Physical Map from Chromosome and Position  from Row Names
#'
#' @description
#'
#' This function extracts chromosome and position information from the row names
#' of a genotype matrix. The row names are expected to follow a format
#' where components (for example chromosome, position) are separated by an underscore (`"_"`).
#' The user can specify which part of the row name corresponds to the chromosome and position
#' using chrom_index and pos_index.
#'
#' @param genotype_matrix A matrix where:
#'   - Rows correspond to genetic markers.
#'   - Row names contain marker identifiers with multiple components separated by underscores
#'     (e.g., `"A_1_200"`, `"B_2_300"`, `"C_3_400"`).
#'   - Columns correspond to individuals (not used in the function but part of the input structure).
#' @param chrom_index Integer. The index of the chromosome identifier in the split row name.
#'   Default is `1` (first element).
#' @param pos_index Integer. The index of the position identifier in the split row name.
#'   Default is `2` (second element).
#' @param markers Logical. If `TRUE`, includes marker names in the output data frame.
#'   Default is `FALSE`.
#'
#' @return A data frame with columns:
#'   - `"chrom"`: Chromosome identifier extracted from row names.
#'   - `"position"`: Numeric genomic position extracted from row names.
#'   - `"marker"` (optional): Marker names (if `markers = TRUE`).
#'
#' @details
#' - Splits row names at the underscore (`"_"`) to extract chromosome and position.
#' - Uses indexes specified by user to determine which part of the row name corresponds to the chromosome and position.
#' - Ensures the position column is numeric.
#' - If `markers = TRUE`, the output includes an additional `"marker"` column, constructed as the combination of Chrom and Position.
#'
#' @examples
#' \dontrun{
#' # Example genotype matrix with row names in "Component1_Component2_Component3" format
#' geno_matrix <- matrix(nrow = 3, ncol = 2)
#' rownames(geno_matrix) <- c("A_1_200", "B_2_300", "C_3_400")
#'
#' # Default behavior (chromosome = first, position = second)
#' map1 < -extract_map(geno_matrix)
#' print(map1)
#'
#' # Custom indexes (chromosome = second, position = third)
#' map2 <- extract_map(geno_matrix, chrom_index = 2, pos_index = 3)
#' print(map2)
#'}
#' @importFrom stats setNames
#' @export
extract_map <- function(genotype_matrix, chrom_index = 1, pos_index = 2, markers = FALSE) {
  # Extract marker names from row names
  marker.name <- rownames(genotype_matrix)

  # Split marker names into chromosome and position
  chr_pos <- strsplit(marker.name, split = "_", fixed = TRUE)

  # Extract chromosome using user-specified index
  chrom <- sapply(chr_pos, function(x) x[chrom_index])

  # Extract position using user-specified index and convert to numeric
  pos <- sapply(chr_pos, function(x) as.numeric(x[pos_index]))

  # Create dataframe with chrom and position
  if (markers) {
    map <- data.frame(marker = marker.name, chrom = chrom, position = pos, stringsAsFactors = FALSE)
  } else {
    map <- data.frame(chrom = chrom, position = pos, stringsAsFactors = FALSE)
  }

  return(map)
}
