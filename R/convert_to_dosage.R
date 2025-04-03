#' Convert VCF Genotype Format to Dosage for Diploids
#'
#'@description
#' The function converts a VCF
#' genotype matrix (with format `"0/0"`, `"0/1"`, `"1/1"`)
#' to allele dosage values (`0`, `1`, `2`) representing the count of the alternative allele.
#' Mostly used within the package. Check out the more flexible `convert_to_dosage_flex` which works for polyploids.
#'
#' @param GT A character matrix where each entry represents a genotype in VCF format
#'   (`"0/0"`, `"0/1"`, `"1/1"`).
#'
#' @return A numeric matrix of the same dimensions as `GT`, where each value represents
#'   the dosage of the alternate allele (`0`, `1`, or `2`).
#'
#' @details
#' - Converts the alleles to numeric values and sums them to compute the alternate allele dosage.
#' - Retains the original row and column names from the input matrix.
#'
#' @examples
#' # Example genotype matrix in VCF format
#' vcf_matrix <- matrix(c("0/0", "0/1", "1/1",
#'                        "0/1", "0/0", "1/1"),
#'                      nrow = 2, ncol = 3,
#'                      dimnames = list(c("Marker1", "Marker2"),
#'                                      c("Ind1", "Ind2", "Ind3")))
#'
#' # Convert to dosage
#' dosage_matrix <- convert_to_dosage(vcf_matrix)
#' print(dosage_matrix)
#'
#' @export
convert_to_dosage <- function(GT) {
  dims <- dim(GT)  # Store original matrix dimensions

  first_allele <- as.numeric(substr(GT, 1, 1))
  second_allele <- as.numeric(substr(GT, 3, 3))

  # Convert back to matrix
  geno_matrix <- matrix(first_allele + second_allele, nrow = dims[1], ncol = dims[2])

  # Preserve row and column names
  rownames(geno_matrix) <- rownames(GT)
  colnames(geno_matrix) <- colnames(GT)

  return(geno_matrix)
}




