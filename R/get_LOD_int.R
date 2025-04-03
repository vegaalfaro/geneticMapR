#' Get LOD Support Intervals for QTL Peaks (Single or Multiple Traits)
#'
#' Computes LOD support intervals for detected QTLs from \link[qtl]{scanone()}  or \link[qtl]{cim()}.
#' Works with either a single model and qtl results or a list of models and results for multiple traits.
#'
#' @param cross_obj A `cross` object from the \pkg{qtl} package. Usually the output of \link[qtl]{read.cross()} or \link[qtl]{calc.genoprob()}.
#' @param model_obj Either a scanone/CIM object (class "scanone") or a named list of such objects. Usuall the output of `cim` or `scanone` from \pkg{qtl}.
#' @param results_obj Either a data frame of QTL results (with "chr" and "pos") or a named list of such data frames.
#' @param trait Optional. Name of the trait to analyze. Required if model_obj and results_obj are lists (i.e., when working with multiple traits).
#' @param drop Numeric. LOD drop to define the support interval (default = 1.5).
#'
#' @return A named list of QTL interval summaries, flanking markers, peak marker (pseudo marker or physical marker), and usable marker (closest physical marker if peak is on a pseudo marker).
#'
#' @examples
#' \dontrun{
#' # Multiple traits:
#' # For the length-width trait
#' lw_result_pop1 <- get_LOD_int(
#'   cross_obj = M1,
#'   model_obj = cim_qtl_results1,
#'   results_obj = results_pop1,
#'   trait = "length_width_ratio",
#'   drop = 1.5
#' )
#'
#' # Single traits:
#' mod <- cim_qtl_results1[["length_width_ratio"]]
#' res <- results_pop1[["length_width_ratio"]]
#'
#' single <- get_LOD_int(
#'   cross_obj = M1,
#'   model_obj = mod,
#'   results_obj = res
#' )
#' }
#'
#' @export
get_LOD_int <- function(cross_obj,
                        model_obj,
                        results_obj,
                        trait = NULL,
                        drop = 1.5) {

  # Detect single vs. multiple trait mode
  if (inherits(model_obj, "scanone") && is.data.frame(results_obj)) {
    model <- model_obj
    results <- results_obj
  } else if (is.list(model_obj) && is.list(results_obj)) {
    if (is.null(trait)) stop("Please provide a trait name when using lists.")
    model <- model_obj[[trait]]
    results <- results_obj[[trait]]
  } else {
    stop("Invalid input: model_obj must be a 'scanone' object or a list of 'scanone' objects.")
  }

  # Check for required columns
  if (!all(c("chr", "pos") %in% names(results))) {
    stop("QTL summary must contain 'chr' and 'pos' columns.")
  }

  # Early return if no QTLs were detected
  if (nrow(results) == 0) {
    message("No QTLs found for the specified trait. Returning empty list.")
    return(list())
  }

  # Process each QTL
  qtl_results <- lapply(seq_len(nrow(results)), function(i) {
    chr <- results[i, "chr"]
    pos <- results[i, "pos"]

    flanking <- find.flanking(cross_obj, chr = chr, pos = pos)


    lod_interval <- lodint(model,
                           chr = chr,
                           expandtomarkers = TRUE,
                           drop = drop)

    lod_peak_marker <- rownames(lod_interval)[2]
    is_pseudomarker <- grepl("^c\\d+\\.loc", lod_peak_marker)
    nearest_real_marker <- if (is_pseudomarker) as.character(flanking$left) else lod_peak_marker


    list(
      chr = chr,
      pos = pos,
      flanking = flanking,
      lod_interval = lod_interval,
      lod_peak_marker = lod_peak_marker,
      usable_marker = nearest_real_marker
    )
  })

  names(qtl_results) <- paste0("QTL_", seq_len(nrow(results)))
  return(qtl_results) # This is old functionality
}


