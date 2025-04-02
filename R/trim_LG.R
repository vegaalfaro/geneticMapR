#' Trim Linkage Groups Based on LOD Thresholds
#'
#' `trim_LG`  is an interactive function that helps the user filter markers within a specified chromosome based on linkage
#' group (LG) assignment using LOD score thresholds. `trim_LG`  should make the process easy by allowing the user
#' to choose the thresholds interactively, visualize haplotype frequency before and after filtering,
#' and  remove outliers. Inspired by functions by Professor Jeffrey B. Endelman's [MapRtools](https://github.com/jendelman/MapRtools)
#' available in Github. `trim_LG` provides a simple interphase for trimming LGs easily.
#'
#' @param chromosome Character. The chromosome ID to be processed.
#' @param map A data frame containing marker map information with at least the following columns:
#'   - `"marker"`: Marker names.
#'   - `"chrom"`: Chromosome identifier.
#'   - `"position"`: Physical position of markers.
#' @param geno A genotype matrix where:
#'   - Rows represent genetic markers.
#'   - Columns represent individuals.
#'   - Values represent genotype calls.
#'   Preferably, a binned genotype matrix. Binning can be performed using `MapRtools::LDbin`.
#' @param pop_type Character. The population type used for linkage estimation. Default is `"F2"`. Would work with any of the following "DH","BC","F2","S1","RIL.self","RIL.sib". Based on `MapRtools::MLEL()`
#' @param drop_outliers Logical. If `TRUE`, removes markers identified as outliers based on haplotype frequency. Default is `TRUE`.
#' @param n_cores Integer. The number of CPU cores to use for linkage estimation.
#'   If `NULL`, the function selects the maximum available minus one.
#'
#' @return A list containing:
#'   - `"trimmed_genotype"`: The filtered genotype matrix.
#'   - `"final_map"`: The updated marker map after filtering.
#'   - `"initial_haplo_plot"`: A `ggplot2` object showing the initial haplotype plot.
#'   - `"filtered_freq_plot"`: A `ggplot2` object showing the haplotype plot after filtering.
#'   - `"final_freq_plot"`: A `ggplot2` object showing the haplotype frequency after removing outliers.
#'   - `"starting_LOD"`: A record of the user-specified initial LOD threshold.
#'   - `"ending_LOD"`: A record of the the user-specified final LOD threshold.
#'   - `"step"`: A record of the LOD sequence step size.
#'   - `"selected_LOD_thresh"`: The final LOD threshold used for filtering markers.
#'   - `"outliers"`: Markers removed as outliers saved as character vector.
#'
#' @details
#' - Computes LOD scores for marker pairs using `MapRtools::MLEL()`, parallelized with multiple cores.
#' - Displays an interactive "candy stripe" plot for visualizing linkage groups using `MapRtools::LG()`.
#' - Asks the user to define a LOD range and a final LOD threshold for marker filtering.
#' - If `drop_outliers = TRUE`, removes markers identified as outliers based on haplotype frequency.
#'
#' @note Requires a previous installation of [MapRtools](https://github.com/jendelman/MapRtools).
#' This function was refined with assistance from ChatGPT to improve readability,
#' efficiency, and visualization formatting.
#'
#' @examples
#' \dontrun{
#' # Example dataset (user should provide actual data)
#' map_data <- data.frame(
#'   marker = c("M1", "M2", "M3", "M4"),
#'   chrom = c("CHR1", "CHR1", "CHR1", "CHR1"),
#'   position = c(10, 20, 30, 40)
#' )
#' geno_matrix <- matrix(sample(0:1, 16, replace = TRUE),
#'                       nrow = 4, ncol = 4,
#'                       dimnames = list(c("M1", "M2", "M3", "M4"),
#'                                       c("Ind1", "Ind2", "Ind3", "Ind4")))
#'
#' # Run trimming function (with user input required for LOD selection)
#' result <- trim_LG(chromosome = "CHR1",
#'                   map = map_data,
#'                   geno = geno_matrix,
#'                   pop_type = "F2",
#'                   drop_outliers = TRUE,
#'                   n_cores = 2)
#'
#' # Access trimmed genotype matrix
#' print(result$trimmed_genotype)
#'}
#' @importFrom dplyr arrange select
#' @importFrom ggplot2 ggtitle
#' @importFrom parallel detectCores
#' @export
trim_LG <- function(chromosome,
                    map,
                    geno,
                    pop_type = "F2",
                    drop_outliers = TRUE,
                    n_cores = NULL) {

  if (!requireNamespace("MapRtools", quietly = TRUE)) {
    stop("The 'MapRtools' package is required but not installed.\nInstall it with:\n  remotes::install_github('jendelman/MapRtools')", call. = FALSE)
  }

  # Detect maximum available cores
  max_cores <- detectCores(logical = FALSE)  # Counts physical cores
  usable_cores <- max_cores - 1

  # Ensure n_cores does not exceed available cores
  if (!is.null(n_cores) && n_cores > max_cores) {
    message("Warning: Requested ", n_cores, " cores, but only ", max_cores,
            " are available. Using ", usable_cores, " cores instead.")
    n_cores <- usable_cores
  }

  # Order the input map
  map <- map %>%
    dplyr::arrange(chrom, position)

  # Convert chromosome to character
  CHR <- as.character(chromosome)
  message("Processing chromosome: ", CHR)

  # Extract marker names for the given chromosome
  chr_markers <- map$marker[which(map$chrom == CHR)]

  if (length(chr_markers) == 0) {
    stop("No markers found for the specified chromosome.")
  }

  # Plot initial haplotype frequency
  message("Plotting initial haplotype frequency...")
  initial_haplo <- MapRtools::plot_haplofreq(geno[chr_markers, ])
  print(initial_haplo$plot + ggtitle(paste0("Chromosome:", CHR)))

  # Subset the genotype matrix for the given chromosome
  chr_geno <- geno[chr_markers, ]

  # Estimate LOD matrix
  message("Estimating LOD matrix using ", n_cores, " cores...")
  LODmat_chr <- MapRtools::MLEL(geno = chr_geno, pop.type = pop_type, LOD = TRUE, n.core = n_cores)

  # User-defined LOD range
  initial_LOD <- as.numeric(readline(prompt = "Enter the initial LOD threshold (e.g., 15): "))
  end_LOD <- as.numeric(readline(prompt = "Enter the end LOD threshold (e.g., 45): "))
  step_LOD <- as.numeric(readline(prompt = "Enter the step size for LOD sequence (e.g., 1): "))

  # Validate LOD inputs
  if (is.na(initial_LOD) | is.na(end_LOD) | is.na(step_LOD) | initial_LOD >= end_LOD | step_LOD <= 0) {
    stop("Invalid LOD range inputs. Ensure initial < end and step > 0.")
  }

  # Generate linkage group plot
  message("Generating linkage group plot...")
  candy_plot <- MapRtools::LG(LODmat_chr, thresh = seq(initial_LOD, end_LOD, by = step_LOD))
  print(candy_plot + ggtitle(paste0("Chromosome:", CHR)))

  # User-defined LOD threshold for filtering
  selected_thresh <- as.numeric(readline(prompt = "Enter the LOD threshold to use for filtering markers: "))

  if (is.na(selected_thresh) | selected_thresh < initial_LOD | selected_thresh > end_LOD) {
    stop("Invalid LOD threshold selected. It must be within the specified range.")
  }

  # Filter markers
  message("Filtering markers based on LOD threshold ", selected_thresh, "...")
  LG_result <- LG(LODmat = LODmat_chr, thresh = selected_thresh)
  map_chr_filtered <- merge(map, LG_result[LG_result$LG == 1, ], by = "marker")
  map_chr_filtered <- map_chr_filtered[order(map_chr_filtered$position), ]

  # Remove outliers if specified
  if (drop_outliers) {
    message("Removing outliers...")
    outliers <- MapRtools::plot_haplofreq(chr_geno)$outliers
    chr_geno <- chr_geno[!rownames(chr_geno) %in% outliers, ]
    map_chr_filtered <- map_chr_filtered[!map_chr_filtered$marker %in% outliers, ]
  }

  return(list(
    trimmed_genotype = chr_geno,
    final_map = map_chr_filtered,
    initial_haplo_plot = initial_haplo$plot,
    filtered_freq_plot = MapRtools::plot_haplofreq(chr_geno)$plot,
    selected_LOD_thresh = selected_thresh,
    outliers = outliers
  ))
}
