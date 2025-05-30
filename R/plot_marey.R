#' Generate a Marey Plot for one Chromosome
#'
#' This function creates a **Marey plot**, which visualizes the relationship between
#' physical position (in megabases) and genetic distance (in centimorgans) for a
#' single chromosome.
#'
#' @param map A data frame containing marker mapping information with at least the following columns:
#'   - `"position_Mb"`: Physical position of markers (best in Mb).
#'   - `"position_cM"`: Genetic distance in centimorgans (cM).
#' @param chrom Character or numeric. The chromosome identifier to be displayed in the plot title.
#'
#' @return A `ggplot2` scatter plot showing the Marey plot.
#'
#' @examples
#' \dontrun{
#' # Example map data
#' example_map <- data.frame(
#'   position = c(1000000, 5000000, 10000000, 20000000, 30000000),
#'   p20 = c(0, 5, 10, 20, 30)
#' )
#'
#' # Generate Marey plot for Chromosome 1
#' plot_marey(example_map, chrom = 1)
#'}
#' @importFrom ggplot2 ggplot aes geom_point labs
#' @export
plot_marey <- function(map, chrom) {
  p <- ggplot(map, aes(position_Mb, position_cM)) +
    geom_point(shape = 21, fill = "black", color = "black", size = 1.8) +
    labs(
      x = "Physical position (Mb)",
      y = "Genetic distance (cM)",
      title = paste("Chromosome", chrom)
    )
  print(p)
}

