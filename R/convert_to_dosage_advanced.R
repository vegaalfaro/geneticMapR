#' Convert VCF Genotype Format to Dosage (Advanced, Polyploid-Compatible)
#'
#' @description
#
#' Converts a matrix of genotype calls from VCF-style format (e.g., "0/1", "1/1", "0/0/1/1")
#' into numeric dosage values. Supports variable ploidy, multi-allelic variants (e.g., "2", "3"),
#' and includes optional outputs for ploidy level, usable allele counts, and normalized dosage.
#'
#' @param GT A character matrix with genotypes in VCF format (e.g., "0/1/1/1", "1/1/1/1").
#' @param alt_alleles A character vector indicating which allele values should be counted as alternate (e.g., c("1") or c("1", "2")).
#' @param strict_missing Logical. If TRUE (default), any missing allele (i.e., ".") causes the entire dosage value to be set as NA.
#' @param normalize Logical. If TRUE, dosage values are normalized by ploidy level (i.e., scaled to 0-1). Default is FALSE.
#'
#' @return A list with the following elements:
#'   - dosage: Matrix of alternate allele dosage values.
#'   - ploidy: Matrix of the number of alleles per genotype (excluding NAs).
#'   - usable_alleles: Matrix of the number of non-missing alleles used in dosage calculation.
#'
#' @note
#' I found ways of generalizing the function to work with more complex situations from previous versions of the function.
#' This function has been tested but not rigorously. Contact author with any issues.
#'
#'
#' @examples
#' \dontrun{
#' vcf <- matrix(c("0/0/0/0", "1/1/1/1", "0/1/1/1",
#'                 "1/1/1/1", NA, "0/1/./1"),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("Marker1", "Marker2"),
#'                               c("Ind1", "Ind2", "Ind3")))
#'
#' convert_to_dosage_advanced(vcf, alt_alleles = c("1"))
#' convert_to_dosage_advanced(vcf, alt_alleles = c("1"), strict_missing = FALSE, normalize = TRUE)
#'
#'
#'
#'set.seed(123)  # for reproducibility
#'
#' # Function to generate one genotype (VCF-style)
#' generate_genotype <- function(ploidy = 4, missing_rate = 0.05) {
#'   alleles <- sample(c(0, 1, 2, "."), size = ploidy,
#'   replace = TRUE,
#'   prob = c(0.40, 0.40, 0.10, missing_rate))
#' paste(alleles, collapse = "/")
#' }
#'
#' # Parameters
#' n_markers <- 100
#' n_individuals <- 10
#' ploidy <- 4  # tetraploid
#'
# Simulate matrix
#' genotype_matrix <- matrix(
#'   data = replicate(n_markers * n_individuals, generate_genotype(ploidy = ploidy)),
#'   nrow = n_markers,
#'   ncol = n_individuals,
#'   dimnames = list(
#'     paste0("Marker", seq_len(n_markers)),
#'     paste0("Ind", seq_len(n_individuals))
#'   )
#' )
#'
#' # Preview
#' head(genotype_matrix)
#'
#'
#' result <- convert_to_dosage_advanced(genotype_matrix,
#'                                      alt_alleles = c("1", "2"),
#'                                      strict_missing = TRUE,
#'                                      normalize = FALSE)
#'
#' # View dosage matrix
#' head(result$dosage)
#'
#' # View ploidy matrix
#' head(result$ploidy)
#'
#' # View number of usable alleles
#' head(result$usable_alleles)
#'
#'}
#' @export
convert_to_dosage_advanced <- function(GT, alt_alleles = c("1"), strict_missing = TRUE, normalize = FALSE) {
  dims <- dim(GT)
  dosage <- matrix(NA_real_, nrow = dims[1], ncol = dims[2])
  ploidy <- matrix(NA_integer_, nrow = dims[1], ncol = dims[2])
  usable_alleles <- matrix(NA_integer_, nrow = dims[1], ncol = dims[2])

  for (i in seq_len(dims[1])) {
    for (j in seq_len(dims[2])) {
      gt <- GT[i, j]
      if (is.na(gt) || gt %in% c(".", "./.")) {
        dosage[i, j] <- NA_real_
        ploidy[i, j] <- NA_integer_
        usable_alleles[i, j] <- NA_integer_
      } else {
        alleles <- unlist(strsplit(gt, split = "[/|]"))
        alleles_num <- suppressWarnings(as.numeric(alleles))

        if (strict_missing && any(is.na(alleles_num))) {
          dosage[i, j] <- NA_real_
          ploidy[i, j] <- length(alleles)
          usable_alleles[i, j] <- sum(!is.na(alleles_num))
        } else {
          alt_count <- sum(alleles_num %in% as.numeric(alt_alleles), na.rm = TRUE)
          non_missing_count <- sum(!is.na(alleles_num))
          dosage[i, j] <- alt_count
          ploidy[i, j] <- length(alleles)
          usable_alleles[i, j] <- non_missing_count
        }
      }
    }
  }

  if (normalize) {
    dosage <- dosage / ploidy
  }

  rownames(dosage) <- rownames(GT)
  colnames(dosage) <- colnames(GT)
  rownames(ploidy) <- rownames(GT)
  colnames(ploidy) <- colnames(GT)
  rownames(usable_alleles) <- rownames(GT)
  colnames(usable_alleles) <- colnames(GT)

  return(list(
    dosage = dosage,
    ploidy = ploidy,
    usable_alleles = usable_alleles
  ))
}
