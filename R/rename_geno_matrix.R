#' Rename Parents, F1, and F2 Individuals in a Genotype Matrix
#'
#' This function standardizes the column names of a genotype matrix by renaming
#' parental genotypes (`P1`, `P2`), F1 individuals (`F1.1`, `F1.2`, ...), and
#' F2 individuals (`F2.1`, `F2.2`, ...). This renaming helps reduce the name size
#' and exchange it with a consistent labeling system for dendrogram visualization
#' constrains
#'
#' @param geno A genotype matrix or data frame where:
#'   - Rows represent genetic markers.
#'   - Columns represent individuals (parents, F1, and F2).
#' @param parent1 Character. The column name corresponding to the first parent.
#' @param parent2 Character. The column name corresponding to the second parent.
#' @param f1 Character vector. The column names corresponding to F1 individuals.
#'
#' @return A genotype matrix with updated column names:
#'   - `"P1"` for `parent1`.
#'   - `"P2"` for `parent2`.
#'   - `"F1.1"`, `"F1.2"`, ... for F1 individuals.
#'   - `"F2.1"`, `"F2.2"`, ... for all other individuals.
#'
#' @examples
#' # Example genotype matrix
#' geno_matrix <- matrix(sample(0:2, 30, replace = TRUE),
#'                       nrow = 5, ncol = 6,
#'                       dimnames = list(
#'                         paste0("Marker", 1:5),
#'                         c("ParentA-Plate1-WellAH", "ParentB-Plate2-WellAJ",
#'                         "F1abc", "F1bcd", "Ind1", "Ind2")
#'                       ))
#'
#' # Rename genotype matrix
#' renamed_geno <- rename_geno_matrix(geno_matrix,
#'                                     parent1 = "ParentA",
#'                                     parent2 = "ParentB",
#'                                     f1 = c("F1abc", "F1bcd"))
#'
#' # Print renamed genotype matrix
#' print(colnames(renamed_geno))
#'
#' @export
rename_geno_matrix <- function(geno, parent1, parent2, f1) {
  # Get column names
  col_names <- colnames(geno)

  # Create a vector for new names
  new_names <- col_names

  # Rename Parents
  new_names[new_names == parent1] <- "P1"
  new_names[new_names == parent2] <- "P2"

  # Rename F1 individuals
  f1_indices <- col_names %in% f1
  new_names[f1_indices] <- paste0("F1.", seq_len(sum(f1_indices)))

  # Identify F2 individuals (those that are neither P1, P2, nor F1)
  f2_indices <- !(col_names %in% c(parent1, parent2, f1))

  # Rename F2 individuals sequentially
  new_names[f2_indices] <- paste0("F2.", seq_len(sum(f2_indices)))

  # Assign the new column names to the geno matrix
  colnames(geno) <- new_names

  return(geno)
}

