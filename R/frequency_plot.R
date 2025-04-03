#' Generate a Histogram of Genotype Frequencies
#'
#'
#'
#' @description
#' The function creates a histogram of genotype frequencies for each genotype category
#' (e.g., `"0", "1", "2"`) based on a frequency data frame. The function also
#' returns the processed long-format data as an attribute. It will not work if your data
#' is in A, H, B format. Use `formater` to change to dosage before using `frequency_plot.`
#'
#' @param freq_df Usually the output of `freq`. A data frame where rows represent markers or individuals,
#'   and columns represent genotype categories with their frequencies.
#'
#' @return A `ggplot2` histogram visualizing the distribution of genotype frequencies.
#'   The processed long-format data is attached as an attribute (`attr(output, "data")`).
#'
#' @details
#' - Converts the input data frame to long format using `pivot_longer()`.
#' - Ensures correct ordering of genotype categories.
#' - Generates a faceted histogram where each panel represents a genotype category.
#' - Stores both the generated plot and the processed data but returns only the plot by default.
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
#' # Generate the frequency histogram
#' p <- frequency_plot(freq_data)
#' print(p)  # Display the plot
#'
#' # Access the processed long-format data
#' attr(p, "data")
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap labs theme
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @export
frequency_plot <- function(freq_df) {

  Marker <- NULL
  Dosage <- NULL
  Frequency <- NULL

  # Convert to long format
  melted_df <- freq_df %>%
    tibble::rownames_to_column("Marker") %>%
    pivot_longer(cols = -Marker, names_to = "Dosage", values_to = "Frequency") %>%
    mutate(Dosage = factor(Dosage, levels = c("0", "1", "2")))  # Ensure correct ordering

  # Create the histogram plot
  p <- ggplot(melted_df, aes(x=Frequency)) +
    geom_histogram(fill="black") +
    facet_wrap(~Dosage) +
    labs(x = "Frequency", y = "Count") +
    theme(
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.border = element_blank()       # Remove panel border
    )

  # Store both plot and data, but return only the plot by default
  output <- structure(p, data = melted_df)

  return(output)
}


