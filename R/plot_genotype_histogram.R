#' Plot a Histogram of Genotypic Values
#'
#' Generates a histogram of genotypic values (`0`, `1`, `2`)
#' from a genotype matrix and shows the distribution of genotypic classes.
#'
#' @param ans A genotype matrix or data frame where:
#'   - Rows correspond to markers.
#'   - Columns correspond to individuals.
#'   - Values are numeric genotypes (`0`, `1`, `2`).
#'
#' @return A `ggplot2` histogram object.
#'
#' @examples
#' # Example genotype matrix
#' geno_matrix <- matrix(
#'   sample(0:2, 30, replace = TRUE),
#'   nrow = 10, ncol = 3,
#'   dimnames = list(paste0("Marker", 1:10), paste0("Ind", 1:3))
#' )
#'
#' # Generate and display the histogram
#' hist_plot <- plot_genotype_histogram(geno_matrix)
#' print(hist_plot)
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom scales comma
#' @export
plot_genotype_histogram <- function(ans) {
  Marker <- NULL
  Value <- NULL
  # Convert genotype matrix to long format
  genotype_data_long2 <- as.data.frame(ans) %>%
    tibble::rownames_to_column(var = "Marker") %>%
    tidyr::pivot_longer(-Marker, names_to = "Genotype", values_to = "Value")

  # Generate histogram
  geno_hist_plot <- ggplot(genotype_data_long2, aes(x = Value)) +
    geom_histogram(binwidth = 0.2,   # Adjust bin width
                   fill = "#A680B8", # Use purple color for bars
                   color = "#2C3E50", # Dark navy border for contrast
                   alpha = 0.8) +    # Slight transparency for visual appeal
    labs(
      x = "Genotype",
      y = "Count") +
    theme_minimal(base_size = 14) +  # Change theme
    scale_x_continuous(breaks = c(0, 1, 2)) +
    scale_y_continuous(labels = scales::comma) +  # Use comma format for Y-axis
    theme(
      plot.title = element_text(size = 13),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 12.5),
      axis.title.x = element_text(size = 12.5),
      panel.background = element_rect(fill = "#FAF3E0", color = NA), # Light beige background
      panel.grid.major = element_line(color = "#D5DBDB", linewidth = 0.5),
      panel.grid.minor = element_blank(), # Remove minor grid lines
      legend.position = "none" # No legend needed for a histogram
    )

  return(geno_hist_plot)
}
