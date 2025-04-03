#' Create a Physical Map from Chromosome and Position from Row Names
#'
#' @description
#' This function extracts chromosome and position information from the row names
#' of a genotype matrix. The row names are expected to follow a format
#' where components (for example chromosome, position) are separated by a symbol (e.g., `"_"`, `"-"`, `"."`).
#' The user can specify which part of the row name corresponds to the chromosome and position
#' using `chrom_index` and `pos_index`, and can customize the delimiter using `split_symbol`.
#'
#' @param genotype_matrix A matrix where:
#'   - Rows correspond to genetic markers.
#'   - Row names contain marker identifiers with multiple components separated by a delimiter
#'     (e.g., `"A_1_200"`, `"B-2-300"`, `"C.3.400"`).
#'   - Columns correspond to individuals (not used in the function but part of the input structure).
#' @param chrom_index Integer. The index of the chromosome identifier in the split row name.
#'   Default is `1` (first element).
#' @param pos_index Integer. The index of the position identifier in the split row name.
#'   Default is `2` (second element).
#' @param markers Logical. If `TRUE`,  includes original marker names from input in the output data frame.
#'   Default is `FALSE`.
#' @param split_symbol Character. The delimiter used to split the row names.
#'   Default is `"_"`. Can be changed to `"-"`, `"."`, or other valid delimiters.
#'
#' @return A data frame with columns:
#'   - `"chrom"`: Chromosome identifier extracted from row names.
#'   - `"position"`: Numeric genomic position extracted from row names.
#'   - `"marker"` (optional): Marker names (if `markers = TRUE`).
#'
#' @details
#' - Row names are split using the specified `split_symbol` to extract chromosome and position.
#' - The user defines which components to extract using `chrom_index` and `pos_index`.
#' - Ensures the position column is numeric.
#' - If `markers = TRUE`, includes marker names in the output.
#'
#' @examples
#' \dontrun{
#' # Example genotype matrix with different delimiters
#' geno_matrix <- matrix(nrow = 3, ncol = 2)
#' rownames(geno_matrix) <- c("A_1_200", "B_2_300", "C_3_400")
#'
#' # Default delimiter "_"
#' extract_map(geno_matrix)
#'
#' # Using "-" as delimiter
#' rownames(geno_matrix) <- c("A-1-200", "B-2-300", "C-3-400")
#' extract_map(geno_matrix, chrom_index = 1, pos_index = 3, split_symbol = "-")
#' }
#'
#' @importFrom stats setNames
#' @export
extract_map <- function(genotype_matrix,
                        chrom_index = 1,
                        pos_index = 2,
                        markers = FALSE,
                        split_symbol = "_") {
  # Extract marker names from row names
  marker.name <- rownames(genotype_matrix)

  # Split marker names using user-defined symbol
  chr_pos <- strsplit(marker.name, split = split_symbol, fixed = TRUE)

  # Extract chromosome using user-specified index
  chrom <- sapply(chr_pos, function(x) x[chrom_index])

  # Extract position using user-specified index and convert to numeric
  pos <- sapply(chr_pos, function(x) as.numeric(x[pos_index]))

  # Create dataframe
  if (markers) {
    map <- data.frame(marker = marker.name, chrom = chrom, position = pos, stringsAsFactors = FALSE)
  } else {
    map <- data.frame(chrom = chrom, position = pos, stringsAsFactors = FALSE)
  }

  return(map)
}





