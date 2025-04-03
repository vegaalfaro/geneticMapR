#' Plot LOD Trace
#'
#' Generates a line plot showing QTL traces (LOD scores) across chromosomes for one or more traits.
#' Supports customization of visual appearance and annotations such as significance thresholds and vertical lines.
#'
#' @param qtl_df A data frame containing QTL mapping results. Must include columns: `pos` (genetic position in cM),
#'         `lod` (LOD score), `response_var` (trait name), and `chr` (chromosome). Optional `phys.pos` (physical position in Mb).
#' @param thresholds_df A data frame with LOD score thresholds per trait. Must include columns: `chr`, `hline`, and `response_var`.
#' @param vline_df Optional. A data frame specifying positions for vertical lines. Must include `chr` and `vline` columns.
#'        Can also include a `label` column for text annotation at each vertical line.
#' @param use_physical_pos Logical. If TRUE, use physical position (`phys.pos`) for the x-axis. If FALSE (default), use genetic position (`pos`).
#' @param x_angle Integer. Angle of the x-axis text labels (default is 0).
#' @param trait_colors Optional. A named vector of colors for each trait. If NULL (default), colors are automatically assigned.
#' @param trait_labels Optional. A named vector of human-readable trait labels (can include expressions) to use in the legend.
#' @param x_label Character. Label for the x-axis (default is `"Position (cM)"`).
#' @param y_label Character. Label for the y-axis (default is `"LOD"`).
#' @param plot_title Optional. Title for the plot (default is NULL).
#' @param show_legend Logical. Whether to display the legend (defaults to TRUE).
#' @param facet_nrow Integer. Number of rows in the facet layout (defaults to 1).
#'
#' @return A `ggplot` object displaying QTL traces per trait and chromosome, with optional custom styling.
#'
#' @import ggplot2
#' @import dplyr
#' @import ggpubr
#' @importFrom scales hue_pal
#' @export
plot_qtl_trace <- function(qtl_df,
                           thresholds_df,
                           vline_df = NULL,
                           use_physical_pos = FALSE,
                           x_angle = 0,
                           trait_colors = NULL,
                           trait_labels = NULL,
                           x_label = "Position (cM)",
                           y_label = "LOD",
                           plot_title = NULL,
                           show_legend = TRUE,
                           facet_nrow = 1) {

  # Use physical or genetic position
  qtl_df <- qtl_df %>%
    mutate(chr = as.character(chr),
           plot_pos = if (use_physical_pos) phys.pos else pos)

  if (!is.null(vline_df)) {
    vline_df <- vline_df %>%
      mutate(chr = as.character(chr))
  }

  # Auto-generate trait colors if not provided
  if (is.null(trait_colors)) {
    unique_traits <- unique(qtl_df$response_var)
    trait_colors <- setNames(scales::hue_pal()(length(unique_traits)), unique_traits)
  }

  # Base plot
  p <- ggplot(qtl_df, aes(x = plot_pos, y = lod, color = response_var)) +
    geom_line(linewidth = 0.85, alpha = 0.6) +
    facet_wrap(~chr, nrow = facet_nrow) +
    labs(x = x_label, y = y_label, color = "Trait", title = plot_title) +
    theme_test() +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      axis.text.x = element_text(angle = x_angle, size = 10, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10),
      strip.text = element_text(size = 9)
    ) +
    scale_color_manual(
      values = trait_colors,
      labels = if (!is.null(trait_labels)) trait_labels[names(trait_colors)] else waiver()
    )

  # Add thresholds
  if (!is.null(thresholds_df)) {
    p <- p + geom_hline(data = thresholds_df,
                        aes(yintercept = hline, group = response_var),
                        size = 1.05, alpha = 0.75, show.legend = FALSE)
  }

  # Add vertical lines with optional labels
  if (!is.null(vline_df)) {
    p <- p +
      geom_vline(data = vline_df, aes(xintercept = vline),
                 linetype = "dotted", color = "black", alpha = 0.85, size = 0.63)

    if ("label" %in% colnames(vline_df)) {
      # Calculate max LOD per chromosome to place labels just above peaks
      lod_max_per_chr <- qtl_df %>%
        group_by(chr) %>%
        summarise(y_pos = max(lod, na.rm = TRUE), .groups = "drop")

      vline_df <- vline_df %>%
        left_join(lod_max_per_chr, by = "chr") %>%
        mutate(y_pos = y_pos * 1.15)

      p <- p +
        geom_text(data = vline_df,
                  aes(x = vline, y = y_pos, label = label),
                  inherit.aes = FALSE,
                  angle = 90, vjust = -0.5, size = 2.8)
    }
  }

  return(p)
}
