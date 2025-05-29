#' Plot Genome-wide Coverage
#'
#' @description
#' Generates a genome-wide coverage plot, displaying the positions
#' of markers across chromosomes. It is a customized version of `plot_cover`
#' from **MapRtools** and allows additional aesthetic modifications.
#'
#' @param map A data frame with two columns:
#'   - `"chrom"`: Chromosome identifier (e.g., `"CHR1", "CHR2", ..."`).
#'   - `"position"`: The genomic position of markers (in base pairs).
#' @param limits (Optional) A data frame specifying the maximum position for each chromosome.
#'   If `NULL` (default), the function computes chromosome limits from `map`.
#' @param customize Logical. If `TRUE`, applies additional visual customizations.
#'   Default is `TRUE`.
#'
#' @return A `ggplot2` object representing genome-wide coverage.
#'
#' @details
#' - Converts genomic positions from base pairs to megabases (Mb).
#' - If `limits` are not provided, the function calculates the maximum position
#'   for each chromosome from `map`.
#' - Orders chromosomes and positions correctly for visualization.
#' - Uses `geom_segment()` to generate a chromosome-wide coverage plot.
#' - When `customize = TRUE`, applies a minimalistic theme for enhanced visualization.
#'
#' @note
#' Inspired by functions of Professor Jeffrey Endelman's MapRtools
#'
#' @examples
#' # Example dataset
#' map_data <- data.frame(
#'   chrom = c("CHR1", "CHR1", "CHR2", "CHR2", "CHR3"),
#'   position = c(500000, 1200000, 800000, 1600000, 2000000)
#' )
#'
#' # Basic coverage plot
#' plot_cover(map_data)
#'
#' # Coverage plot with custom aesthetics
#' plot_cover(map_data, customize = TRUE)
#'
#' @import ggplot2
#' @export
plot_cover <- function(map, limits = NULL, customize = TRUE) {
  colnames(map) <- c("chrom", "position")
  map$position <- map$position / 1e6

  # If limits are not provided, calculate the chromosome limits from the map
  if (is.null(limits)) {
    chrom_limits <- tapply(map$position, map$chrom, max)
    limits <- data.frame(chrom = names(chrom_limits), position = as.numeric(chrom_limits))
  } else {
    limits$position <- limits$position / 1e6
  }

  limits <- limits[order(limits$chrom, limits$position), ]
  map <- map[order(map$chrom, map$position), ]
  map$chrom <- as.character(map$chrom)

  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  k <- match(map$chrom, chroms)

  y <- c(k - 0.1, 1:n.chrom)
  yend <- c(k + 0.1, 1:n.chrom)
  x <- c(map$position, rep(0, n.chrom))
  xend <- c(map$position, limits$position)

  p <- ggplot(data = data.frame(x = x, y = y, xend = xend, yend = yend),
              aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment()

  # Apply custom aesthetics if requested
  if (customize) {
    p <- p +
      theme_minimal(base_size = 14) +
      labs(
        x = "Position (Mb)",
        y = "Chromosome"
      ) +
      scale_y_continuous(breaks = 1:n.chrom, labels = chroms, minor_breaks = NULL) +
      theme(
        plot.title = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12.5),
        axis.title.x = element_text(size = 12.5),
        panel.background = element_rect(fill = "#FAF3E0", color = NA),
        panel.grid.major = element_line(color = "#D5DBDB", size = 0.5),
        panel.grid.minor = element_blank()
      ) +
      geom_segment(color = "#2C3E50")
  } else {
    p <- p + theme_bw() + xlab("Position (Mb)") +
      scale_y_continuous(name = "Chromosome", breaks = 1:n.chrom, labels = chroms, minor_breaks = NULL)
  }

  return(p)
}





