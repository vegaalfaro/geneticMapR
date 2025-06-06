#'
#' F2 diploid population Haplotype Reconstruction
#'
#' This function performs haplotype reconstruction using a Hidden Markov Model (HMM).
#' It applies the Viterbi algorithm to infer the most likely sequence of true genotypic
#' states, accounting for genotyping errors and missing data.
#'
#' @param geno_matrix A numeric genotype matrix where:
#'   - Rows represent genetic markers.
#'   - Columns represent progeny (individuals).
#'   - Values are expected to be `0`, `1`, `2`, or `NA` for missing data.
#' @param error_rate Numeric. The assumed genotyping error rate. Default is `0.05` (5% error rate).
#' @param r Numeric. The recombination rate used in the transition probability matrix (`T.mat`).
#'   Default is `0.01`. The user can use the `MLEL` function from [MapRtools](https://github.com/jendelman/MapRtools)
#'   to get an average estimate of recombination frequency with the adjacent argument set to `TRUE`
#'   `result <- MLEL(geno = geno_matrix, pop.type = "f2", LOD = FALSE, adjacent = TRUE)`
#'   In this case the geno_matrix should be of only one chromosome. `r` could be estimated as `mean(result$value, na.rm = TRUE)`
#'
#' @return A genotype matrix with inferred haplotypes based on HMM correction.
#'
#' @details
#' - Estimates the missing data rate in the input genotype matrix.
#' - Initializes an HMM using:
#'   - Three hidden states (`"0"`, `"1"`, `"2"`).
#'   - Four observable symbols (`"0"`, `"1"`, `"2"`, `"NA"` for missing data).
#'   - Transition probabilities generated via `T.mat(r = r)`, using the user-defined recombination rate.
#'   - Emission probabilities computed using `E.mat(error = error_rate, missing = missing)`.
#'
#' @note
#' Works only for F2 populations on experimental crosses. Based on the 615 Genetic Mapping Class notes by
#' Prof. Jeffrey Endelman. Spring 2021.
#'
#' @examples
#' \dontrun{
#' # Example genotype matrix with missing values
#' geno_data <- matrix(c(0, 1, NA, 2, 0, 1, NA, 2, NA, NA, 0, 1),
#'                     nrow = 4, ncol = 3,
#'                     dimnames = list(c("Marker1", "Marker2", "Marker3", "Marker4"),
#'                                     c("Ind1", "Ind2", "Ind3")))
#'
#' # Perform haplotype reconstruction with default parameters
#' reconstructed_geno <- haplotype_reconstruction(geno_data)
#' print(reconstructed_geno)
#'
#' # Perform haplotype reconstruction with a modified recombination rate
#' reconstructed_geno_custom_r <- haplotype_reconstruction(geno_data, error_rate = 0.05, r = 0.02)
#' print(reconstructed_geno_custom_r)
#'}
#' @importFrom HMM initHMM viterbi
#' @export
haplotype_reconstruction <- function(geno_matrix, error_rate = 0.05, r = 0.01) {
  n <- ncol(geno_matrix)  # Number of progeny
  m <- nrow(geno_matrix)  # Number of markers

  # Transition matrix for an F2
  # r is recombination frequency
  # s = 1-r
  T.mat <- function(r){
    s <- 1 - r
    mat <- rbind(
      c(s^2, 2*r* s, r^2),
      c(r*s, s^2 + r^2, r*s),
      c(r^2, 2*r*s, s^2)
    )
    rownames(mat) <- colnames(mat) <- c("0", "1", "2")
    return(mat)
  }

  # Emission probability matrix: depends on error rate and missing data
  E.mat <- function(error, missing = 0) {
    mat <- rbind(c(1-error, error/2, error/2),
                 c(error/2, 1-error, error/2),
                 c(error/2, error/2, 1-error)
    )
    mat <- cbind(mat - missing/3, rep(missing, 3))
    rownames(mat) <- c("0", "1", "2")
    colnames(mat) <- c("0", "1", "2", "NA")
    return(mat)
  }

  # Estimate missing rate
  missing <- sum(is.na(geno_matrix)) / (m * n)

  # Define HMM
  hmm <- initHMM(States = c("0", "1", "2"),
                 Symbols = c("0", "1", "2", "NA"),
                 startProbs = c(1, 2, 1) / 4,
                 transProbs = T.mat(r = r),  # Use the user-specified or default r
                 emissionProbs = E.mat(error = error_rate, missing = missing))

  # Apply Viterbi algorithm
  geno_processed <- geno_matrix
  for (i in 1:n) {
    input <- as.character(geno_matrix[, i])
    input[is.na(input)] <- "NA"
    out <- viterbi(hmm, observation = input)
    geno_processed[, i] <- as.integer(out)
  }

  return(geno_processed)
}
