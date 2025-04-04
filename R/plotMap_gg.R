#' Plot Genetic Map Using ggplot2 in R/qtl data environment
#'
#' This function generates a **genetic linkage map** visualization using [ggplot2](https://ggplot2.tidyverse.org/).
#' It plots chromosomes as vertical or horizontal lines and displays markers as colored
#' segments along each chromosome using the data structure built in [R/qtl](https://rqtl.org/). Replaces \link[qtl]{plotMap} which uses base R.
#'
#' @param map  `cross` object representing a genetic map. The output of R/qtl's `read.cross()`, `jittermap()` or `calc.genoprob()`
#'   - If a `cross` object, the function extracts the map using `pull.map()`.
#' @param horizontal Logical. If `TRUE`, the map is plotted with chromosomes arranged horizontally.
#'   Default is `FALSE` (vertical layout).
#' @param color Character. The color for the marker segments. Default is `"red"` if `NULL`.
#'
#' @return A `ggplot2` object displaying the genetic map.
#'

#' @examples
#'
#' # Default plot (vertical)
#' # plotMap_gg(example_map)
#'
#' # Horizontal plot with custom color
#' # plotMap_gg(example_map, horizontal = TRUE, color = "blue")
#'
#' @import ggplot2
#' @import dplyr
#' @export
plotMap_gg <- function(map, horizontal = FALSE, color = NULL) {

  # Convert map into a tidy format if necessary
  if (inherits(map, "cross")) {
    map <- pull.map(map)
  }

  # Convert list to dataframe
  df <- bind_rows(
    lapply(names(map), function(chr) {
      data.frame(
        Chromosome = chr,
        Position = map[[chr]],
        Marker = names(map[[chr]])
      )
    })
  )

  # Ensure Chromosome is a factor for correct ordering
  df$Chromosome <- factor(df$Chromosome, levels = unique(df$Chromosome))

  # Get chromosome segment positions (min/max per chromosome)
  chrom_data <- df %>%
    group_by(Chromosome) %>%
    summarise(y_start = min(Position), y_end = max(Position), .groups = "drop")

  # Define plot using geom_segment
  p <- ggplot() +
    # Chromosomes as vertical lines spanning their actual marker range
    geom_segment(data = chrom_data,
                 aes(x = Chromosome, xend = Chromosome,
                     y = y_start, yend = y_end),
                 linewidth = 1, color = "black")

  # Add markers with color handling
  if (!is.null(color)) {
    p <- p + geom_segment(data = df,
                          aes(x = Chromosome, xend = Chromosome,
                              y = Position, yend = Position + 0.20),
                          linewidth = 5, color = color)
  } else {
    p <- p + geom_segment(data = df,
                          aes(x = Chromosome, xend = Chromosome,
                              y = Position, yend = Position + 0.20),
                          linewidth = 5, color = "red")  # Default color
  }

  # Apply common plot settings
  p <- p + theme_minimal() +
    labs(title = "Genetic Map",
         x = "Chromosome",
         y = "Position (cM)") +
    theme(legend.position = "none")

  # Adjust for horizontal layout
  if (horizontal) {
    p <- p + coord_flip() +
      labs(y = "Position (cM)", x = "Chromosome")
  } else {
    p <- p + scale_y_reverse() # Reverse the y-axis when not horizontal
  }

  print(p)
}
