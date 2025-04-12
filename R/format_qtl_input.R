#' Format Input Data for R/qtl package
#'
#' This function is part of the heart of the package. `format_qtl_input` formats genotype, map, and phenotype data for QTL analysis in [R/qtl](https://cran.r-project.org/web/packages/qtl/index.html).
#' It allows the user to specify whether the genotype data should be converted from dosage to ABH format or used as is.
#'
#' @param geno A dataframe containing genotype data with markers as rows and individuals as columns.
#' @param map A dataframe containing the genetic map with columns: marker, chrom, and position.
#' @param pheno A dataframe containing phenotype data with an ID column corresponding to individual/sample names.
#' @param numeric A logical value (TRUE or FALSE) indicating whether to convert dosages to R/qtl ABH format (default = TRUE).
#'
#' @return A dataframe combining genotype, map, and phenotype data in [R/qtl](https://cran.r-project.org/web/packages/qtl/index.html) format  and ready for \link[qtl]{read.cross} after saving it as csv.
#' @seealso [write.csv()], [read.cross()]
#' @note
#' Formats input data for the read.cross function of [R/qtl](https://rqtl.org/).
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # result <- format_qtl_input(geno2, map2, pheno2, numeric = TRUE)
#'}
#'
#' @export
format_qtl_input <- function(geno, map, pheno, numeric = TRUE) {

  # Check if map has the required columns
  required_map_cols <- c("marker", "chrom", "position")
  if (!all(required_map_cols %in% colnames(map))) {
    stop("Error: The map must contain the columns: marker, chrom, position.")
  }

  # Check if pheno has the ID column
  if (!"ID" %in% colnames(pheno)) {
    stop("Error: The phenotype dataset must include an 'ID' column.")
  }

  # Check if ID column in pheno matches colnames of geno
  if (!all(pheno$ID[-c(1,2)] %in% colnames(geno))) {
    warning("Warning: Some IDs in pheno do not match column names in geno.")
  }

  # Convert genotype to ABH format if numeric = TRUE
  if (numeric) {
    geno <- formater(geno = geno, numeric_output = FALSE)
  }

  # Convert all phenotype columns to character
  pheno <- pheno %>% dplyr::mutate(across(everything(), as.character))

  # Create empty rows for phenotypic data
  empty_rows <- as.data.frame(matrix("", nrow = 2,
                                     ncol = ncol(pheno)))

  colnames(empty_rows) <- colnames(pheno)

  # Insert empty rows at the beginning of phenotype data
  pheno <- dplyr::bind_rows(empty_rows, pheno)

  # Add placeholders for ID column
  pheno <- pheno %>% dplyr::mutate(ID = replace(ID, c(1, 2), c("chrom", "position")))

  # Transpose genotype data
  geno_t <- t(geno) %>% as.data.frame()

  # Transpose map data
  map_t <- t(map) %>% as.data.frame()
  colnames(map_t) <- colnames(geno_t)

  # Insert map as first two rows in genotype matrix
  geno_t <- dplyr::bind_rows(map_t, geno_t)

  # Remove marker row
  geno_t <- geno_t[-1, ]

  # Convert row names to column 'ID'
  geno_t <- tibble::rownames_to_column(geno_t, var = "ID")

  # Merge phenotype with genotype
  cross_df <- dplyr::full_join(pheno, geno_t, by = "ID")

  # Remove placeholder names in the first two rows
  cross_df$ID[1:2] <- ""

  # Ensure first two rows are empty
  cross_df[1:2, ] <- dplyr::mutate_all(cross_df[1:2, ], ~ifelse(. == "NA", "", .))

  return(as.data.frame(cross_df))
}

