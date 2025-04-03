#' Convert (Melt) Genotype Frequency Data to Long Format
#'
#' @description
#'
#' Transforms a genotype frequency data frame from wide format
#' (where genotype categories are columns) to long format, making it easier
#' to use with `ggplot2` and other tidyverse functions. Works well with `freq`
#' for data visualization.
#'
#' @param freq_df Usually the output from `freq`. A data frame where rows represent
#' markers or individuals,  and columns represent genotype categories with
#'  their frequencies.
#'
#' @return A data frame in long format with three columns:
#'   - `"Marker"`: The marker or individual ID.
#'   - `"Dosage"`: The genotype category.
#'   - `"Frequency"`: The relative frequency of the genotype category.
#'
#' @details
#' - Uses `tibble::rownames_to_column()` to preserve marker or individual names.
#' - Reshapes the data using `pivot_longer()` from `tidyverse`.
#' - Ideal for visualization with `ggplot2` or further data analysis on genotype
#' frequencies by marker or by individual.
#'
#' @examples
#' # Example frequency data frame
#' freq_data <- data.frame(
#'   `0` = c(0.2, 0.3, 0.4),
#'   `1` = c(0.5, 0.4, 0.4),
#'   `2` = c(0.3, 0.3, 0.2),
#'   row.names = c("Marker1", "Marker2", "Marker3")
#' )
#'
#' # Convert to long format
#' melt_freq(freq_data)
#'
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#'
#' @export
melt_freq <- function(freq_df) {
  Marker <- NULL
  # Convert to long format
  melted_df <- freq_df %>%
    tibble::rownames_to_column("Marker") %>%
    pivot_longer(cols = -Marker, names_to = "Dosage", values_to = "Frequency")

  return(melted_df)
}


