#' Powerful Effect Plots
#'
#' This function generates powerful effect plots to show the relationship between genotype classes at a a marker (QTL) on a given trait.
#' **Works on polyploids**
#' The function computes summary statistics for each genotype class and overlays them on the plot.
#' Optionally, it can flip the coordinate axes.
#' Shows the distributions and measures of spread of individuals and not simply a boxplot.
#'
#' @param effects_df A data frame containing at least the columns for marker and trait of interest.
#' Could work with the output of `format_qtl_input`with a few modifications (first two rows removed).
#' @param marker_name A character string specifying the name of the marker column in `effects_df`.
#' @param trait_name A character string specifying the name of the trait column in `effects_df`.
#' @param genotype_levels A character vector indicating the levels (genotype categories) for the marker. Default is \code{c("A", "H", "B")}.
#' @param flip Logical; if \code{TRUE}, the plot is displayed with flipped coordinates (horizontal layout). Default is \code{TRUE}.
#' @param trait_label Optional character string to use as a y-axis label. If \code{NULL}, the trait name will be used.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{[[1]]}}{A \code{ggplot2} object showing the distribution of trait values for each genotype at the given marker.}
#'   \item{\code{[[2]]}}{A data frame with summary statistics (median, max, min, count, standard deviation) by genotype.}
#'   \item{\code{[[3]]}}{A data frame with population-level summary statistics (mean, median, sd, max) for the trait.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks for presence of specified columns.
#'   \item Calculates summary statistics by genotype.
#'   \item Constructs a half-eye plot with boxplot overlays and annotations.
#' }
#' It uses `ggplot2` for plotting and `ggdist` for the half-eye visualization.
#'
#' @note
#' The half-eye visualization is inspired in work by [Cedric Scherer](https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/)
#'
#' @importFrom ggdist stat_halfeye
#' @importFrom dplyr group_by summarise
#'
#' @examples
#' \dontrun{
#' effect_plot(effects_df = my_data,
#'                 marker_name = "SNP_1",
#'                 trait_name = "PlantHeight",
#'                 genotype_levels = c("A", "H", "B"),
#'                 flip = TRUE,
#'                 trait_label = "Height (cm)")
#'
#'  #Works for polyploids if genotype levels are specified and match the `effects_df`.
#'
#' effect_plot(effects_df = my_data,
#'                 marker_name = "SNP_1",
#'                 trait_name = "PlantHeight",
#'                 genotype_levels = c("0", "1", "2", "3", "4"), # or
#'                   c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")
#'                 flip = TRUE,
#'                 trait_label = "Height (cm)")
#'
#'
#' }
#'
#' @export
effect_plot <- function(effects_df,
                            marker_name,
                            trait_name,
                            genotype_levels = c("A", "H", "B"),
                            flip = TRUE,
                            trait_label = NULL) {
  # Check marker and trait exist
  if (!marker_name %in% colnames(effects_df)) {
    stop(paste("Marker", marker_name, "not found in the dataset."))
  }
  if (!trait_name %in% colnames(effects_df)) {
    stop(paste("Trait", trait_name, "not found in the dataset."))
  }

  # Remove NAs in the marker column
  effects_df <- effects_df[!is.na(effects_df[[marker_name]]), ]

  # Convert to factor
  effects_df[[marker_name]] <- factor(effects_df[[marker_name]], levels = genotype_levels)

  # Compute medians and group stats
  medians <- effects_df %>%
    dplyr::group_by(.data[[marker_name]]) %>%
    dplyr::summarise(
      median = median(.data[[trait_name]], na.rm = TRUE),
      max = max(.data[[trait_name]], na.rm = TRUE),
      min = min(.data[[trait_name]], na.rm = TRUE),
      count = n(),
      sd = sd(.data[[trait_name]], na.rm = TRUE)
    )

  # Trait Standard Deviation (SD)
  # Compute Standard Deviation of Trait
  SD <- effects_df %>%
    dplyr::summarise(
      pop_sd = round(sd(.data[[trait_name]], na.rm = TRUE), 2),
      pop_mean = round(mean(.data[[trait_name]], na.rm = TRUE), 2),
      pop_median = round(median(.data[[trait_name]], na.rm = TRUE), 2),
      pop_max = round(max(.data[[trait_name]], na.rm = TRUE), 2)
    )

  # Label fallback
  y_label <- ifelse(is.null(trait_label), trait_name, trait_label)

  # Base plot
  p <- ggplot(effects_df, aes(x = .data[[marker_name]], y = .data[[trait_name]], fill = .data[[marker_name]])) +
    labs(x = paste("Locus", marker_name), y = y_label) +
    ggdist::stat_halfeye(
      adjust = 0.5,
      alpha = 0.85,
      width = 0.60,
      .width = 0,
      justification = -0.2,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.15,
      alpha = 0.40,
      outlier.shape = NA,
      show.legend = FALSE
    ) +
    geom_text(
      data = medians,
      aes(x = .data[[marker_name]], y = median, label = round(median, 1)),
      vjust = 0,
      hjust = 2.2,
      show.legend = FALSE,
      size = 4
    ) +
    geom_text(
      data = medians,
      aes(x = .data[[marker_name]], y = median, label = paste("n =", count)),
      vjust = 2.0,
      hjust = 1.75,
      show.legend = FALSE
    ) +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = "none", color = "none")

  # Flip coordinates optionally
  if (flip) {
    p <- ggplot(effects_df, aes(x = .data[[marker_name]], y = .data[[trait_name]], fill = .data[[marker_name]])) +
      labs(x = paste("Locus", marker_name), y = y_label) +
      ggdist::stat_halfeye(
        adjust = 0.5,
        alpha = 0.85,
        width = 0.60,
        .width = 0,
        justification = -0.2,
        point_colour = NA
      ) +
      geom_boxplot(
        width = 0.15,
        alpha = 0.40,
        outlier.shape = NA,
        show.legend = FALSE
      ) +
      geom_text(
        data = medians,
        aes(x = .data[[marker_name]], y = median, label = round(median, 1)),
        vjust = 2.2,
        show.legend = FALSE,
        size = 4
      ) +
      geom_text(
        data = medians,
        aes(x = .data[[marker_name]], y = median, label = paste("n =", count)),
        vjust = 3.5,
        show.legend = FALSE
      ) +
      theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = "none", color = "none")+
      coord_flip()
  }

  return(list(p,
              medians,
              SD))
}

