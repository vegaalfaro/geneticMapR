#' Convert VCF Genotype Format to Dosage (Supports Polyploids and Missing Alleles)
#'
#' @description
#' Converts a VCF genotype matrix (in the GT format e.g., "0/0", or for poplyplids "0/0/1/1") into numeric dosage values representing
#' the count of alternate alleles. Works with any ploidy (*but tested only from haploid to hexaploid*) and allows control over how missing
#' alleles (".") are handled.
#'
#' @param GT A character matrix with genotypes in VCF format (e.g., "0/1/1/1", "1/1/1/1").
#' @param strict_missing Logical. If TRUE (default), any missing allele (e.g., ".") causes the entire genotype to be set as NA.
#'   If FALSE, missing alleles are ignored and dosage is calculated from the available alleles.
#'
#' @details
#' - This function is flexible for any ploidy level — it simply counts how many `"1"` alleles exist in each genotype.
#' - Alleles are split using either `/` or `|`, so phased or unphased VCF data are supported.
#' - Missing alleles (`"."`) are handled based on the `strict_missing` argument:
#'   - `TRUE`: If *any* allele is missing in a genotype, the entire dosage is returned as `NA`.
#'   - `FALSE`: Missing alleles are ignored and dosage is calculated from known alleles.
#'
#'
#' @return A numeric matrix with the same dimensions and names as `GT`, where each value is the
#'   dosage of the alternate allele (assumed to be "1").
#'
#' @examples
#' poly_vcf <- matrix(c("0/0/0/0", "1/1/1/1", "0/1/1/1",
#'                      "1/1/1/1", NA, "0/1/./1"),
#'                    nrow = 2, byrow = TRUE,
#'                    dimnames = list(c("Marker1", "Marker2"),
#'                                    c("Ind1", "Ind2", "Ind3")))
#'
#' # Strict handling: missing allele causes full NA
#' convert_to_dosage_flex(poly_vcf, strict_missing = TRUE)   # Set whole cell to NA if "." present
#'
#' # Permissive handling: ignore missing and sum known alleles
#' convert_to_dosage_flex(poly_vcf, strict_missing = FALSE)  # Ignore "." and sum what’s there
#'
#' @export
convert_to_dosage_flex <- function(GT, strict_missing = TRUE) {
  dims <- dim(GT)
  dosage <- matrix(NA_real_, nrow = dims[1], ncol = dims[2])

  for (i in seq_len(dims[1])) {
    for (j in seq_len(dims[2])) {
      gt <- GT[i, j]
      if (is.na(gt) || gt %in% c(".", "./.")) {
        dosage[i, j] <- NA_real_
      } else {
        alleles <- unlist(strsplit(gt, split = "[/|]"))
        alleles_num <- suppressWarnings(as.numeric(alleles))

        if (strict_missing && any(is.na(alleles_num))) {
          dosage[i, j] <- NA_real_
        } else {
          dosage[i, j] <- sum(alleles_num[!is.na(alleles_num)])
        }
      }
    }
  }

  rownames(dosage) <- rownames(GT)
  colnames(dosage) <- colnames(GT)

  return(dosage)
}
