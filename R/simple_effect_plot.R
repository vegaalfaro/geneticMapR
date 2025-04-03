#' Simple Effect Plot
#'
#' Generates a clean and simple effect plot. Shows the relationship between genotype classes at a marker on a trait.
#' Optionally includes annotation of medians and sample sizes, and allows flipping of axes.
#' **Works on polyploids.**
#'
#' @param effects_df A data frame containing at least the columns for marker genotypes and trait values. Usually the output of `format_qtl_input`with a few modifications (first two rows removed).
#' @param marker_name A character string specifying the name of the marker column in `effects_df`.
#' @param trait_name A character string specifying the name of the trait column in `effects_df`.
#' @param genotype_levels A character vector indicating the expected genotype categories (factor levels). Default is \code{c("A", "H", "B")}.
#' @param annotate Logical; if \code{TRUE}, the plot includes annotations of median values and sample sizes per genotype. Default is \code{TRUE}.
#' @param flip Logical; if \code{TRUE}, the coordinates of the plot are flipped (horizontal layout). Default is \code{FALSE}.
#' @param trait_label Optional character string to use as the y-axis label. If \code{NULL}, the trait name is used.
#'
#' @return A \code{ggplot2} object displaying the distribution of trait values across genotype classes (alleles) at the specified marker.
#'
#' @details
#' The function performs the following:
#' \itemize{
#'   \item Verifies the presence of marker and trait columns in the input data frame.
#'   \item Removes rows with missing genotype data.
#'   \item Converts the marker column to a factor with specified genotype levels.
#'   \item Creates a boxplot with optional annotations of median values and counts.
#'   \item Optionally flips the plot horizontally.
#' }
#' @note
#' Inspired by  Fig. 2  [Caraza-Harter & Endelman](https://doi.org/10.1007/s00122-022-04159-z)
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise
#'
#' @examples
#' \dontrun{
#' simple_effect_plot(effects_df = my_data,
#'                         marker_name = "SNP_42",
#'                         trait_name = "RootWeight",
#'                         genotype_levels = c("A", "H", "B"),
#'                         annotate = TRUE,
#'                         flip = FALSE,
#'                         trait_label = "Weight (g)")
#'
#' # Works for polyploids if genotype levels are specified and match the `effects_df`.
#'
#' simple_effect_plot(effects_df = my_data,
#'                 marker_name = "SNP_1",
#'                 trait_name = "PlantHeight",
#'                 genotype_levels = c("0", "1", "2", "3", "4"),
#'                 # or  c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")
#'                 flip = TRUE,
#'
#' }
#'
#'
#' @export
simple_effect_plot <- function(effects_df,
                                   marker_name,
                                   trait_name,
                                   genotype_levels = c("A", "H", "B"),
                                   annotate = TRUE,
                                   flip = FALSE,
                                   trait_label = NULL) {
  # Check inputs
  if (!marker_name %in% colnames(effects_df)) {
    stop(paste("Marker", marker_name, "not found in the dataset."))
  }
  if (!trait_name %in% colnames(effects_df)) {
    stop(paste("Trait", trait_name, "not found in the dataset."))
  }

  # Clean and prep data
  effects_df <- effects_df[!is.na(effects_df[[marker_name]]), ]
  effects_df[[marker_name]] <- factor(effects_df[[marker_name]], levels = genotype_levels)
  y_label <- ifelse(is.null(trait_label), trait_name, trait_label)

  # Compute summary if needed
  medians <- effects_df %>%
    dplyr::group_by(.data[[marker_name]]) %>%
    dplyr::summarise(
      median = median(.data[[trait_name]], na.rm = TRUE),
      count = n()
    )

  # Base plot
  p <- ggplot(effects_df, aes(x = .data[[marker_name]],
                              y = .data[[trait_name]]
  )) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    labs(x = paste("Genotype at", marker_name), y = y_label) +
    theme_minimal() +
    theme(legend.position = "none")

  # Annotations
  if (annotate) {
    p <- p +
      geom_text(data = medians,
                aes(x = .data[[marker_name]], y = median, label = round(median, 1)),
                vjust = -1.2, size = 3.5) +
      geom_text(data = medians,
                aes(x = .data[[marker_name]], y = median, label = paste0("n=", count)),
                vjust = 1.8, size = 3)
  }

  # Optionally flip
  if (flip) {
    p <- p + coord_flip()
  }

  return(p)
}
