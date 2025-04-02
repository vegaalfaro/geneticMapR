#' Drop Parents and F1 Individuals from a Data Frame
#'
#' @description
#'
#'
#' This helper function removes specified parental and F1 columns from a genotype
#' dataset.
#'
#' @param y A data frame where individuals (including parents and F1s) are columns.
#' @param parent1 Character. The name of the first parent column to be removed. Default is `"P1"`.
#' @param parent2 Character. The name of the second parent column to be removed. Default is `"P2"`.
#' @param F1 (Optional) Character. The name of the F1 hybrid column to be removed. Default is `NULL`.
#'
#' @return A data frame with the specified parent and F1 columns removed.
#'
#' @details
#' -  `y` must be coercible to a data frame.
#' - Checks if the specified columns exist before attempting to drop them.
#' - Only removes columns that are present in the data frame.
#' - Preserves all other columns and their structure.
#' - If multiple F1s are present they can be declared using the `c()` function
#'
#' @examples
#' # Example dataset
#' geno_data <- data.frame(
#'   Marker1 = c(0, 1, 2),
#'   P1 = c(0, 0, 0),
#'   P2 = c(2, 2, 2),
#'   F1 = c(1, 1, 1),
#'   Ind1 = c(0, 1, 2),
#'   Ind2 = c(2, 0, 1)
#' )
#'
#' # Drop parents P1 and P2
#' filtered_data <- drop_parents(geno_data)
#' print(filtered_data)
#'
#' # Drop parents and F1
#' filtered_data_f1 <- drop_parents(geno_data, F1 = "F1")
#' print(filtered_data_f1)
#'
#' # Example dataset2
#' geno_data2 <- data.frame(
#'   Marker1 = c(0, 1, 2),
#'   P1 = c(0, 0, 0),
#'   P2 = c(2, 2, 2),
#'   F1a = c(1, 1, 1),
#'   F1b = c(1, 1, 1),
#'   Ind1 = c(0, 1, 2),
#'   Ind2 = c(2, 0, 1)
#' )
#'
#'
#' # Drop parents and multiple F1s
#' filtered_data_f1 <- drop_parents(geno_data, F1 = c("F1a", "F1b"))
#' print(filtered_data_f1)
#'
#' @export
drop_parents <- function(y, parent1 = "P1", parent2 = "P2", F1 = NULL) {
  y <- as.data.frame(y)  # Ensure y is a data frame
  cols_to_drop <- c(parent1, parent2, F1)
  cols_to_drop <- cols_to_drop[cols_to_drop %in% colnames(y)]  # Only drop existing columns
  y <- y[, !colnames(y) %in% cols_to_drop, drop = FALSE]  # Drop selected columns
  return(y)
}
