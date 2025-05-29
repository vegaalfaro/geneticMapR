#' Plot Coverage Map with Candidate Gene Annotations
#'
#'
#' @description
#'
#' Visualizes chromosome positions, QTLs, and gene of interest annotations on a multi-chromosome physical map.
#' Useful for displaying genomic regions of interest, highlighting genetic features such as
#' QTLs and candidate genes or protein families.
#'
#' @param map A data frame with at least two columns: \code{chrom} (chromosome ID) and \code{position}
#'   (genomic coordinate in base pairs). This forms the base map.
#' @param limits Optional data frame with chromosome end positions. If \code{NULL}, limits are computed automatically.
#' @param qtls Optional data frame of QTLs. Should contain \code{chrom}, \code{position}, and optionally \code{trait}.
#' @param qtls2 Optional second QTL set (e.g., from another population). Same format as \code{qtls}.
#' @param prot1 First protein/gene annotation data frame. Should contain \code{chrom}, \code{position}, and \code{name}.
#' @param prot2 Second protein/gene annotation data frame. Same format as \code{prot1}.
#' @param prot3 Third protein/gene annotation data frame. Same format as \code{prot1}.
#' @param labels Optional vector of labels (length 1–3) for protein layers (e.g., \code{c("OFP", "IQD", "TRM")}).
#' @param protein_colors Optional vector of fill colors for proteins, matching the order in \code{labels}.
#' @param qtl_labels Character vector of length 2 defining the legend labels for \code{qtls} and \code{qtls2}.
#' @param qtl_colors Character vector of length 2 defining the fill colors used for the two QTL types.
#' @param show_labels Logical; whether to display text labels for QTL traits and protein names (default: \code{TRUE}).
#'
#' @return A \code{ggplot} object showing the annotated chromosome map.
#'
#' @examples
#'
#' # Simulated example
#' set.seed(123)
#'
#' # Create basic map for 3 chromosomes
#' example_map <- data.frame(
#'   chrom = rep(paste0("CHR", 1:3), each = 100),
#'   position = rep(seq(1e6, 100e6, length.out = 100), 3)
#' )
#'
#' # QTL sets
#' example_qtls <- data.frame(
#'   chrom = sample(paste0("CHR", 1:3), 10, replace = TRUE),
#'   position = runif(5, min = 1e6, max = 100e6),
#'   trait = paste("Trait", 1:10)
#' )
#' example_qtls2 <- data.frame(
#'   chrom = sample(paste0("CHR", 1:3), 6, replace = TRUE),
#'   position = runif(3, min = 1e6, max = 100e6),
#'   trait = paste("AltTrait", 1:6)
#' )
#'
#' # Protein annotations
#' ofp_data <- data.frame(
#'   chrom = sample(paste0("CHR", 1:3), 5, replace = TRUE),
#'   position = runif(5, 1e6, 100e6),
#'   name = paste("OFP", 1:5)
#' )
#'
#' iqd_data <- data.frame(
#'   chrom = sample(paste0("CHR", 1:3), 4, replace = TRUE),
#'   position = runif(4, 1e6, 100e6),
#'   name = paste("IQD", 1:4)
#' )
#'
#' trm_data <- data.frame(
#'   chrom = sample(paste0("CHR", 1:3), 3, replace = TRUE),
#'   position = runif(3, 1e6, 100e6),
#'   name = paste("TRM", 1:3)
#' )
#'
#' # Plot annotated coverage map
#' plot_coverage_annotate(
#'   map = example_map,
#'   qtls = example_qtls,
#'   qtls2 = example_qtls2,
#'   prot1 = ofp_data,
#'   prot2 = iqd_data,
#'   prot3 = trm_data,
#'   labels = c("OFP", "IQD", "TRM"),
#'   qtl_labels = c("QTL Set 1", "QTL Set 2"),
#'   qtl_colors = c("steelblue", "darkorange"),
#'   show_labels = TRUE
#' )
#'
#' @import ggplot2
#' @export
plot_coverage_annotate <- function(map,
                                   limits = NULL,
                                   qtls = NULL,
                                   qtls2 = NULL,
                                   prot1 = NULL,
                                   prot2 = NULL,
                                   prot3 = NULL,
                                   labels = NULL,
                                   protein_colors = NULL,
                                   qtl_labels = c("QTLs Full Panel", "QTLs Table Beet Only"),
                                   qtl_colors = c("#1E90FF", "#D95F02"),
                                   show_labels = TRUE)

{
  # Custom fallback operator. Use x if it is not NULL; otherwise use y
  `%||%` <- function(x, y) if (!is.null(x)) x else y # Returns x if not NULL, otherwise y

  # Validate and prep map
  if (!all(c("chrom", "position") %in% colnames(map))) {
    stop("Input 'map' must contain 'chrom' and 'position' columns.")
  }
  map$chrom <- as.character(map$chrom)
  map$position <- as.numeric(map$position) / 1e6
  map <- map[order(map$chrom, map$position), ]
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)

  # Process limits
  if (is.null(limits)) {
    tmp <- tapply(map$position, map$chrom, max)
    limits <- data.frame(chrom = names(tmp), position = as.numeric(tmp))
  } else {
    limits$position <- as.numeric(limits$position) / 1e6
  }
  limits$chrom <- as.character(limits$chrom)
  limits <- limits[order(limits$chrom), ]

  # Plot backbone
  map$y <- match(map$chrom, chroms)
  limits$y <- match(limits$chrom, chroms)

  p <- ggplot() +
    geom_segment(data = map,
                 aes(x = position, y = y - 0.1, xend = position, yend = y + 0.1)) +
    geom_segment(data = limits,
                 aes(x = 0, y = y, xend = position, yend = y)) +
    theme_bw() +
    xlab("Position (Mb)") +
    scale_y_continuous(name = "Chromosome", breaks = 1:n.chrom, labels = chroms, minor_breaks = NULL)

  # Helper to standardize QTL/protein datasets
  add_layer <- function(data, offset, label) {
    data$chrom <- as.character(data$chrom)
    data$position <- as.numeric(data$position) / 1e6
    data$y <- match(data$chrom, chroms) + offset
    data$type <- label
    return(data)
  }

  # QTL layers
  if (!is.null(qtls)) {
    qtls <- add_layer(qtls, 0.2, qtl_labels[1])
    p <- p + geom_point(data = qtls,
                        aes(x = position, y = y, fill = type),
                        shape = 25, size = 3,
                        color = scales::alpha("#1E90FF", 0.85),
                        inherit.aes = FALSE)

    if (show_labels && "trait" %in% colnames(qtls)) {
      p <- p + geom_text(data = qtls, aes(x = position, y = y + 0.2, label = trait),
                         size = 2.5, hjust = 0, inherit.aes = FALSE)
    }
  }

  if (!is.null(qtls2)) {
    qtls2 <- add_layer(qtls2, 0.4, qtl_labels[2])

    p <- p + geom_point(data = qtls2,
                        aes(x = position, y = y, fill = type),
                        shape = 25, size = 3,
                        color = scales::alpha("#D95F02", 0.6),
                        inherit.aes = FALSE)

    if (show_labels && "trait" %in% colnames(qtls2)) {
      p <- p + geom_text(data = qtls2, aes(x = position, y = y + 0.2, label = trait),
                         size = 2.5, hjust = 0, inherit.aes = FALSE)
    }
  }

  # Combine protein layers
  protein_list <- list(prot1, prot2, prot3)
  n_prot <- sum(!sapply(protein_list, is.null))
  if (!is.null(labels) && length(labels) != n_prot) {
    stop("Length of 'labels' must match the number of non-NULL protein datasets.")
  }

  labels <- labels %||% paste("Protein", seq_len(n_prot))
  protein_colors <- protein_colors %||% c("#CCCC00", "#E7298A", "#66A61E")[1:n_prot]
  names(protein_colors) <- labels

  for (i in seq_len(n_prot)) {
    dat <- protein_list[[i]]
    if (!is.null(dat)) {
      dat <- add_layer(dat, 0, labels[i])
      p <- p + geom_tile(data = dat,
                         aes(x = position, y = y, fill = type),
                         width = 0.35, height = 0.5, alpha = 0.9, inherit.aes = FALSE)
      if (show_labels && "name" %in% colnames(dat)) {
        p <- p + geom_text(data = dat,
                           aes(x = position, y = y + 0.4, label = name),
                           size = 2.3, hjust = 0, inherit.aes = FALSE)
      }
    }
  }

  # Build legend
  legend_colors <- c(setNames(qtl_colors, qtl_labels), protein_colors)

  p <- p + scale_fill_manual(name = "", values = legend_colors) +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12))

  return(p)
}
