#' Phase Genotype Marker Data Based on Parental References
#'
#' @description
#' `recode` is the heart of the package. This powerful function phases genotype marker
#' data based on two parental references (`parent1` and `parent2`).
#'
#' It phases markers according to parental allele inheritance.
#' Optionally, it can code phased markers into numeric (`0`, `1`, `2`) or character (`"A"`, `"B"`, `"H"`) formats.
#' Numeric coding is recommended for downstream analyses.
#'
#' @param geno A genotype matrix or data frame where markers are rows and individuals are columns.
#' @param parent1 Character. The name of the column representing the first parent.
#' @param parent2 Character. The name of the column representing the second parent.
#' @param numeric_output Logical. If `TRUE`, converts phased markers to numeric dosage values
#'   (`A = 0`, `H = 1`, `B = 2`). Default is `TRUE`.
#' @param handle_het_markers Logical. If `TRUE`, allows heterozygous parent markers to be included.
#'   Default is `FALSE`.
#' @param het_marker_types Character vector. Specifies which heterozygous markers to keep
#'   when `handle_het_markers = TRUE`. Options include `"AxH"`, `"HxB"`, `"HxA"`, `"BxH"`.
#'   Default is `NULL`, meaning all homozygous markers are kept. See details.
#'
#' @return A data frame containing phased genotype markers where:
#'   - `"0"` represents alleles inherited from `parent1`.
#'   - `"2"` represents alleles inherited from `parent2`.
#'   - `"1"` represents heterozygous alleles.
#'   - If `numeric_output = FALSE`,   `0`, `2`, `1`, are replaced by `"A"`, `"B"`, `"H"` respectively.
#'   - `numeric_output = TRUE` is recommended.
#'
#' @details
#' - Drops markers where either parent has `NA`.
#' - Removes non-polymorphic markers (markers where both parents have the same genotype).
#' - If `handle_het_markers = FALSE`, retains only homozygous marker where dosages are as follows: (`P1 = 0` & `P2 = 2` or `P1 = 2` & `P2 = 0`)
#' - If `handle_het_markers = TRUE`, allows heterozygous markers to be kept.
#' - Ensures that `parent1` is always `0` and `parent2` is always `2` for standardization.
#' - Returns a numeric matrix if `numeric_output = TRUE`, otherwise returns phased `"A"`, `"B"`, `"H"` values.
#' - More details on the heterozygous F2 marker types ("AxH", "HxB", "HxA", "BxH") are in [Braun et al. (2017).](https://acsess.onlinelibrary.wiley.com/doi/10.3835/plantgenome2016.10.0110)
#'
#' @examples
#' # Example genotype data
#' geno_data <- data.frame(
#'   Marker1 = c(0, 1, 2, 0, 2),
#'   Marker2 = c(2, 0, 2, 1, 0),
#'   Parent1 = c(0, 2, 2, 0, 2),
#'   Parent2 = c(2, 0, 0, 2, 0)
#' )
#'
#' # Recode genotype markers (default: numeric output)
#' phased_geno <- recode(geno_data, "Parent1", "Parent2")
#' print(phased_geno)
#'
#' # Recode genotype markers with heterozygous marker handling
#' phased_geno_het <- recode(geno_data, "Parent1", "Parent2",
#'                           handle_het_markers = TRUE,
#'                           het_marker_types = c("AxH", "HxB"))
#' print(phased_geno_het)
#'
#'
#'
#'
#' @importFrom dplyr mutate across all_of case_when
#' @export
recode <- function(geno, parent1,
                   parent2,
                   numeric_output = TRUE,
                   handle_het_markers = FALSE,
                   het_marker_types = NULL) {

  # Ensure geno is a dataframe
  geno <- as.data.frame(geno)

  # Identify and drop markers where either parent1 or parent2 has NA
  na_markers <- which(is.na(geno[[parent1]]) | is.na(geno[[parent2]]))
  if (length(na_markers) > 0) {
    geno <- geno[-na_markers, ]
    message("Dropping markers with NA in at least one of the parents. Your genotype matrix contains missing parental genotypes.")
  }

  # Identify and drop non-polymorphic markers (both parents == 0, 1, or 2)
  non_poly_markers <- which(  geno[[parent1]] == 0 & geno[[parent2]] == 0 |
                                geno[[parent1]] == 1 & geno[[parent2]] == 1 |
                                geno[[parent1]] == 2 & geno[[parent2]] == 2)

  if (length(non_poly_markers) > 0) {
    geno <- geno[-non_poly_markers, ]
    message("Dropping non-polymorphic markers. Your genotype matrix contains markers that are identical in the two parents.")
  }

  # Define the heterozygous marker types
  het_types <- list(
    "AxH" = which(geno[[parent1]] == 0 & geno[[parent2]] == 1),
    "HxB" = which(geno[[parent1]] == 1 & geno[[parent2]] == 2),
    "HxA" = which(geno[[parent1]] == 1 & geno[[parent2]] == 0),
    "BxH" = which(geno[[parent1]] == 2 & geno[[parent2]] == 1)
  )

  # Keep only homozygous markers if handle_het_markers is FALSE
  if (!handle_het_markers) {
    homozygous_markers <- which(geno[[parent1]] %in% c(0,2) & geno[[parent2]] %in% c(0,2))
    geno <- geno[homozygous_markers, ]
  } else if (!is.null(het_marker_types)) {
    keep_markers <- unlist(het_types[het_marker_types], use.names = FALSE)
    homozygous_markers <- which(geno[[parent1]] %in% c(0,2) & geno[[parent2]] %in% c(0,2))
    keep_markers <- unique(c(keep_markers, homozygous_markers))
    geno <- geno[keep_markers, ]
  }

  # Apply phasing logic across all marker columns (excluding parent columns)
  phased_geno <- geno %>%
    dplyr::mutate(dplyr::across(-dplyr::all_of(c(parent1, parent2)), ~ dplyr::case_when(
      (.data[[parent1]] == 0 & .data[[parent2]] == 2 & . == 0) ~ "A",
      (.data[[parent1]] == 0 & .data[[parent2]] == 2 & . == 1) ~ "H",
      (.data[[parent1]] == 0 & .data[[parent2]] == 2 & . == 2) ~ "B",

      (.data[[parent1]] == 2 & .data[[parent2]] == 0 & . == 0) ~ "B",
      (.data[[parent1]] == 2 & .data[[parent2]] == 0 & . == 1) ~ "H",
      (.data[[parent1]] == 2 & .data[[parent2]] == 0 & . == 2) ~ "A",

      TRUE ~ as.character(.)
    )))

  # Set parent1 = 0 and parent2 = 2
  phased_geno <- phased_geno %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(parent1), ~ dplyr::case_when(. == 1 ~ 1, TRUE ~ 0)),
      dplyr::across(dplyr::all_of(parent2), ~ dplyr::case_when(. == 1 ~ 1, TRUE ~ 2))
    )

  # Convert phased markers to numeric output
  if (numeric_output) {
    phased_geno <- phased_geno %>%
      dplyr::mutate(dplyr::across(-dplyr::all_of(c(parent1, parent2)), ~ dplyr::case_when(
        . == "A" ~ 0,
        . == "B" ~ 2,
        . == "H" ~ 1,
        TRUE ~ suppressWarnings(as.numeric(.))
      )))
  }

  return(phased_geno)
}





