#' Filter Genotype Matrix Based on Parental Genotype
#'
#'
#'@description
#'
#' `filter_geno_by_parents` filters a genotype matrix to retain only homozygous polymorphic markers
#' between two specified parents. It removes markers that do not meet the homozygous
#' polymorphism criteria (`P1 = 0` & `P2 = 2` or `P1 = 2` & `P2 = 0`).
#'
#' @param geno A genotype matrix or data frame where markers are rows and individuals are columns.
#'   Must be coercible to a data frame.
#' @param parent1 Character. The name of the column representing the first parent.
#' @param parent2 Character. The name of the column representing the second parent.
#'
#' @return A filtered data frame containing only markers that are homozygous polymorphic
#'   between the two parents.
#'
#' @details
#' - Retains markers where `parent1` is `0` and `parent2` is `2` (`A × B`)
#'   or `parent1` is `2` and `parent2` is `0` (`B × A`). See the methods section
#'   of Braun et al. 2017 for  more info on marker types. [The Plant Genome Vol. 10 No. 3.](https://acsess.onlinelibrary.wiley.com/doi/pdf/10.3835/plantgenome2016.10.0110)
#'
#' @examples
#' # Example genotype matrix
#' geno_data <- data.frame(
#'   Marker1 = c(0, 1, 2, 0, 2),
#'   Marker2 = c(2, 0, 2, 1, 0),
#'   Parent1 = c(0, 2, 2, 0, 2),
#'   Parent2 = c(2, 0, 0, 2, 0)
#' )
#'
#' # Filter markers based on parents
#' filtered_geno <- filter_geno_by_parents(geno_data, "Parent1", "Parent2")
#' print(filtered_geno)
#'
#' @importFrom dplyr filter mutate across where
#' @export
filter_geno_by_parents <- function(geno, parent1, parent2) {
  # Ensure geno is a dataframe
  geno <- as.data.frame(geno)

  # Convert character columns to numeric
  geno <- geno %>% mutate(across(where(is.character), as.numeric))

  # Filter based on specified parent columns
  ans <- geno %>%
    filter((.data[[parent1]] == 0 & .data[[parent2]] == 2) |  # A × B
             (.data[[parent1]] == 2 & .data[[parent2]] == 0))   # B × A (y = 2 - x)

  return(ans)
}



