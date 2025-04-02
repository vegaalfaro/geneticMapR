#' Flexible Genotype Format Converter (Supports Polyploid Dosage)
#'
#' @description
#' Converts genotype data between character string formats (e.g., `"AABB"`, `"ABBB"`) and numeric dosage values
#' (e.g., 0, 1, 2, 3) based on user-defined reference and alternative alleles. This function is designed to handle
#' diploid and polyploid genotypes by interpreting the number of occurrences of the alternative allele.
#'
#' This is especially useful for genotype encodings that represent dosage via string repetition (e.g., `"AAAB"` → 1 alt allele),
#' or to reconstruct genotype strings from numeric dosage values (e.g., 0 → `"AAAA"`).
#'
#' @param geno A data frame containing genotype values. Values must be either character strings representing genotypes
#'   (e.g., `"AAAA"`, `"AAAB"`, `"BBBB"`), or numeric dosage values (e.g., 0 to `ploidy`).
#' @param to_numeric Logical. If `TRUE`, converts character genotypes to numeric dosage by counting the number of
#'   `alt_allele` occurrences. If `FALSE`, reconstructs genotype strings from numeric dosage values using the specified
#'   `ref_allele` and `alt_allele`.
#' @param ref_allele Character. The reference allele. This is used to define the "zero dosage" baseline and for
#'   reconstructing character genotypes from numeric dosage (e.g., "A" in "AAAB").
#' @param alt_allele Character. The alternative allele. Dosage is computed as the count of this allele in the genotype
#'   string (e.g., "B" in "AAAB" → dosage = 1).
#' @param ploidy Integer. Ploidy level of the organism. Only used when `to_numeric = FALSE` to determine how many
#'   reference and alternative alleles to paste together when reconstructing genotype strings.
#'
#' @details
#' - When `to_numeric = TRUE`, all values in `geno` must be character strings, and the function will count the
#'   number of occurrences of the `alt_allele` using `stringr::str_count()`. Missing values (`NA`) are preserved.
#' - When `to_numeric = FALSE`, all values in `geno` must be numeric dosages between 0 and `ploidy`, inclusive.
#'   The function will reconstruct character genotype strings by repeating the `ref_allele` and `alt_allele`
#'   the appropriate number of times. For example, with `ploidy = 4`, a dosage of 2 becomes `"AABB"`.
#' - Values outside the expected range (e.g., dosage > ploidy or unknown characters) are converted to `NA`.
#'
#' @return A data frame of the same structure as `geno`, with genotypes converted according to the `to_numeric` setting.
#'
#' @examples
#' # Example 1: Character to numeric (tetraploid)
#' \dontrun{
#' geno <- data.frame(
#'   Marker1 = c("AAAA", "AAAB", "AABB", "ABBB", "BBBB"),
#'   Marker2 = c("AAAB", "AABB", "ABBB", "BBBB", NA)
#' )
#' formater_flex(geno, to_numeric = TRUE, ref_allele = "A", alt_allele = "B")
#' }
#' # Example 2: Numeric to character (tetraploid)
#'  \dontrun{
#' dosage <- data.frame(
#'   Marker1 = c(0, 1, 2, 3, 4),
#'   Marker2 = c(1, 2, 3, 4, NA)
#' )
#' formater_flex(dosage, to_numeric = FALSE, ref_allele = "A", alt_allele = "B", ploidy = 4)
#' }
#' @importFrom dplyr mutate across case_when
#' @importFrom stringr str_count
#' @export
formater_flex <- function(geno, to_numeric = TRUE, ref_allele = "A", alt_allele = "B", ploidy = 2) {
  geno <- as.data.frame(geno)

  if (to_numeric) {
    geno <- geno %>%
      dplyr::mutate(dplyr::across(everything(), ~ dplyr::case_when(
        is.na(.) ~ NA_real_,
        TRUE ~ stringr::str_count(., fixed(alt_allele))  # Count alt allele in string
      )))
  } else {
    geno <- geno %>%
      dplyr::mutate(dplyr::across(everything(), ~ dplyr::case_when(
        is.na(.) ~ NA_character_,
        . >= 0 & . <= ploidy ~ paste0(strrep(ref_allele, ploidy - .), strrep(alt_allele, .)),
        TRUE ~ NA_character_  # Handle invalid dosages
      )))
  }

  return(geno)
}
