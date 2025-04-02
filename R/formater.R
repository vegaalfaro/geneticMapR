#' Convert Between Character (A, H, B) and Numeric dosage (0, 1, 2) Formats
#'
#' @description
#'
#' This function converts genotype data between character (`"A"`, `"H"`, `"B"`)
#' and numeric (`0`, `1`, `2`) formats. The conversion direction is controlled by the
#' `numeric_output` argument.
#'
#' @param geno A data frame containing genotype data, where values are either
#'   `"A"`, `"H"`, `"B"` or `0`, `1`, `2`. The function assumes the dataset consists
#'   entirely of genotypic or numeric values.
#' @param numeric_output Logical. If `TRUE`, converts `"A"` → `0`, `"H"` → `1`, `"B"` → `2`.
#'   If `FALSE`, converts `0` → `"A"`, `1` → `"H"`, `2` → `"B"`. Default is `TRUE`.
#'
#' @return A data frame with the same structure as `geno`, but with values converted
#'   based on the `numeric_output` parameter.
#'
#' @details
#' - Preserves `NA` values during conversion.
#' - Unexpected values (anything other than `"A"`, `"H"`, `"B"`, `0`, `1`, `2`)
#'   are converted to `NA`.
#'
#' @examples
#' # Example genotype data in character format
#' geno_char <- data.frame(
#'   Marker1 = c("A", "H", "B"),
#'   Marker2 = c("B", "A", "H")
#' )
#'
#' # Convert character to numeric format
#' geno_numeric <- formater(geno_char, numeric_output = TRUE)
#' print(geno_numeric)
#'
#' # Convert numeric format back to character format
#' geno_char_reversed <- formater(geno_numeric, numeric_output = FALSE)
#' print(geno_char_reversed)
#'
#' @importFrom dplyr mutate across case_when everything
#' @export
formater <- function(geno, numeric_output = TRUE) {
  # Ensure geno is a dataframe
  geno <- as.data.frame(geno)

  if (numeric_output) {
    # Convert A, H, B to 0, 1, 2
    geno <- geno %>%
      mutate(across(everything(), ~ case_when(
        . == "A" ~ 0,
        . == "H" ~ 1,
        . == "B" ~ 2,
        is.na(.) ~ NA_real_,  # Preserve NAs
        TRUE ~ as.numeric(NA)  # Convert unexpected values to NA
      )))
  } else {
    # Convert 0, 1, 2 to A, H, B
    geno <- geno %>%
      mutate(across(everything(), ~ case_when(
        . == 0 ~ "A",
        . == 1 ~ "H",
        . == 2 ~ "B",
        is.na(.) ~ NA_character_,  # Preserve NAs
        TRUE ~ as.character(NA)  # Convert unexpected values to NA
      )))
  }

  return(geno)
}
