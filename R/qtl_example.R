#' QTL Dataset Example
#'
#' A list containing QTL LOD scores, thresholds, and vline (annotation) data for demonstrating the functions in `geneticMapR`. It also include a colors character vector and labels for publication ready plots.
#'
#' @format A list with the elements:
#' \describe{
#'   \item{qtl_df}{A scanone and data.frame class object containing LOD scores for each loci in the map for 7 trauts}
#'   \item{thresholds}{A scanone and dara.frame class object containing significant thresholds for each variable }
#'   \item{vline}{A data frame with information for plot annotation}
#'   \item{colors}{A vector with custom colors for plot customization}
#'   \item{labes}{A vector with publication ready labels for plot customization}
#' }
#'
#' @source Own author's QTL analysis using R/qtl on a Table beet shape mapping population (F2).
"qtl_example"

