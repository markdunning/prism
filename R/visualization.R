#' Plot QC Status Distribution
#'
#' Visualizes the distribution of probe QC status
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param colors Named vector of colors for each status (optional)
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal labs coord_flip
#' @examples
#' \dontrun{
#'   plot_qc_status(qc_results)
#' }
plot_qc_status <- function(qc_results, colors = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  # Count status
  status_counts <- as.data.frame(table(qc_results$qc_status))
  colnames(status_counts) <- c("Status", "Count")

  # Calculate percentages
  status_counts$Percentage <- round(100 * status_counts$Count / sum(status_counts$Count), 1)

  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "PASS" = "#2ecc71",
      "DEPRECATED" = "#e74c3c",
      "MULTI_MAPPING" = "#f39c12",
      "CROSS_HYBRIDIZATION" = "#e67e22",
      "AMBIGUOUS" = "#c0392b",
      "INTERGENIC" = "#95a5a6"
    )
  }

  # Create plot
  p <- ggplot2::ggplot(status_counts, ggplot2::aes(x = Status, y = Count, fill = Status)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Count, "\n(", Percentage, "%)")),
                       vjust = -0.5, size = 3.5) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Probe QC Status Distribution",
      subtitle = paste("Total probes:", nrow(qc_results)),
      x = "QC Status",
      y = "Number of Probes"
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  return(p)
}

#' Plot Alignment Statistics
#'
#' Visualizes probe alignment distribution
#'
#' @param alignments GRanges object from align_probes_to_genome
#' @param max_locations Maximum locations to show on x-axis (default: 10)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_alignment_stats(alignments)
#' }
plot_alignment_stats <- function(alignments, max_locations = 10) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Count locations per probe
  location_counts <- as.data.frame(table(alignments$probe_id))
  colnames(location_counts) <- c("probe_id", "n_locations")

  # Bin high values
  location_counts$n_locations_binned <- ifelse(
    location_counts$n_locations > max_locations,
    paste0(">", max_locations),
    as.character(location_counts$n_locations)
  )

  # Count probes in each bin
  binned_counts <- as.data.frame(table(location_counts$n_locations_binned))
  colnames(binned_counts) <- c("Locations", "Probes")

  # Order factor levels
  level_order <- c(as.character(1:max_locations), paste0(">", max_locations))
  binned_counts$Locations <- factor(binned_counts$Locations,
                                    levels = level_order[level_order %in% binned_counts$Locations])

  # Create plot
  p <- ggplot2::ggplot(binned_counts, ggplot2::aes(x = Locations, y = Probes)) +
    ggplot2::geom_bar(stat = "identity", fill = "#3498db") +
    ggplot2::geom_text(ggplot2::aes(label = Probes), vjust = -0.5, size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Distribution of Genomic Locations per Probe",
      subtitle = paste("Total probes:", length(unique(alignments$probe_id))),
      x = "Number of Genomic Locations",
      y = "Number of Probes"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 0)
    )

  return(p)
}

#' Plot Annotation Changes
#'
#' Visualizes how annotations changed from original to updated
#'
#' @param comparison_results Data frame from compare_annotations
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_annotation_changes(comparison)
#' }
plot_annotation_changes <- function(comparison_results) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Count each type of change
  change_counts <- data.frame(
    Category = c("Unchanged", "Changed", "Newly Annotated",
                 "Annotation Lost", "No Annotation"),
    Count = c(
      sum(comparison_results$annotation_status == "unchanged"),
      sum(comparison_results$annotation_status == "changed"),
      sum(comparison_results$annotation_status == "newly_annotated"),
      sum(comparison_results$annotation_status == "annotation_lost"),
      sum(comparison_results$annotation_status == "no_annotation")
    )
  )

  # Calculate percentages
  change_counts$Percentage <- round(100 * change_counts$Count / sum(change_counts$Count), 1)

  # Remove zero counts
  change_counts <- change_counts[change_counts$Count > 0, ]

  # Colors
  colors <- c(
    "Unchanged" = "#2ecc71",
    "Changed" = "#f39c12",
    "Newly Annotated" = "#3498db",
    "Annotation Lost" = "#e74c3c",
    "No Annotation" = "#95a5a6"
  )

  # Create plot
  p <- ggplot2::ggplot(change_counts,
                       ggplot2::aes(x = "", y = Count, fill = Category)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Percentage, "%")),
                       position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Annotation Changes: Original vs Updated",
      subtitle = paste("Total probes:", nrow(comparison_results))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )

  return(p)
}

#' Plot Gene Biotype Distribution
#'
#' Shows distribution of gene biotypes in annotations
#'
#' @param annotations Data frame from update_gene_annotations
#' @param top_n Show top N biotypes (default: 10)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_gene_biotypes(annotations)
#' }
plot_gene_biotypes <- function(annotations, top_n = 10) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (!"gene_biotype" %in% colnames(annotations)) {
    stop("annotations must contain 'gene_biotype' column")
  }

  # Count unique genes per biotype
  biotype_counts <- annotations %>%
    dplyr::group_by(gene_biotype) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(ensembl_gene_id),
      n_probes = dplyr::n_distinct(probe_id),.groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_genes))

  # Take top N
  if (nrow(biotype_counts) > top_n) {
    top_biotypes <- biotype_counts[1:top_n, ]
    other_count <- sum(biotype_counts$n_genes[(top_n + 1):nrow(biotype_counts)])
    other_probes <- sum(biotype_counts$n_probes[(top_n + 1):nrow(biotype_counts)])

    top_biotypes <- rbind(
      top_biotypes,
      data.frame(gene_biotype = "Other", n_genes = other_count, n_probes = other_probes)
    )
  } else {
    top_biotypes <- biotype_counts
  }

  # Reorder by count
  top_biotypes$gene_biotype <- factor(top_biotypes$gene_biotype,
                                      levels = rev(top_biotypes$gene_biotype))

  # Create plot
  p <- ggplot2::ggplot(top_biotypes,
                       ggplot2::aes(x = gene_biotype, y = n_genes, fill = gene_biotype)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::geom_text(ggplot2::aes(label = n_genes), hjust = -0.2, size = 3.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Gene Biotype Distribution",
      subtitle = paste("Top", top_n, "biotypes"),
      x = "Gene Biotype",
      y = "Number of Genes"
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))

  return(p)
}

#' Plot Probes per Gene Distribution
#'
#' Shows how many probes map to each gene
#'
#' @param probe_gene_map Mapping from create_probe_gene_mapping
#' @param max_probes Maximum probes to show on x-axis (default: 20)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_probes_per_gene(mapping)
#' }
plot_probes_per_gene <- function(probe_gene_map, max_probes = 20) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Count probes per gene
  probes_per_gene <- probe_gene_map %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(n_probes = dplyr::n_distinct(probe_id),.groups = "drop")

  # Bin high values
  probes_per_gene$n_probes_binned <- ifelse(
    probes_per_gene$n_probes > max_probes,
    paste0(">", max_probes),
    as.character(probes_per_gene$n_probes)
  )

  # Count genes in each bin
  binned_counts <- as.data.frame(table(probes_per_gene$n_probes_binned))
  colnames(binned_counts) <- c("Probes", "Genes")

  # Order factor levels
  level_order <- c(as.character(1:max_probes), paste0(">", max_probes))
  binned_counts$Probes <- factor(binned_counts$Probes,
                                 levels = level_order[level_order %in% binned_counts$Probes])

  # Create plot
  p <- ggplot2::ggplot(binned_counts, ggplot2::aes(x = Probes, y = Genes)) +
    ggplot2::geom_bar(stat = "identity", fill = "#9b59b6") +
    ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.5, size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Distribution of Probes per Gene",
      subtitle = paste("Total genes:", nrow(probes_per_gene)),
      x = "Number of Probes",
      y = "Number of Genes"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}

#' Plot Probe Distribution Across Chromosomes
#'
#' Shows how probes are distributed across chromosomes
#'
#' @param alignments GRanges object from align_probes_to_genome
#' @param chr_order Optional vector specifying chromosome order
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_chromosome_distribution(alignments)
#' }
plot_chromosome_distribution <- function(alignments, chr_order = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Count unique probes per chromosome
  chr_counts <- as.data.frame(table(seqnames(alignments)))
  colnames(chr_counts) <- c("Chromosome", "Alignments")

  # Order chromosomes
  if (is.null(chr_order)) {
    # Natural ordering: chr1, chr2,..., chr22, chrX, chrY, chrM
    chr_counts$Chromosome <- as.character(chr_counts$Chromosome)

    # Extract numeric part
    chr_num <- suppressWarnings(as.numeric(gsub("chr", "", chr_counts$Chromosome)))

    # Sort: numeric first, then alphabetic
    chr_counts <- chr_counts[order(chr_num, chr_counts$Chromosome, na.last = TRUE), ]

    chr_counts$Chromosome <- factor(chr_counts$Chromosome,
                                    levels = unique(chr_counts$Chromosome))
  } else {
    chr_counts$Chromosome <- factor(chr_counts$Chromosome, levels = chr_order)
  }

  # Create plot
  p <- ggplot2::ggplot(chr_counts,
                       ggplot2::aes(x = Chromosome, y = Alignments, fill = Chromosome)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Probe Alignments per Chromosome",
      subtitle = paste("Total alignments:", sum(chr_counts$Alignments)),
      x = "Chromosome",
      y = "Number of Alignments"
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}

#' Plot Probe Consistency Heatmap
#'
#' Visualizes correlation between probes for the same gene
#'
#' @param expression_matrix Probe-level expression matrix
#' @param probe_gene_map Mapping from create_probe_gene_mapping
#' @param gene Gene name to plot
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_probe_consistency(probe_expr, mapping, "TP53")
#' }
plot_probe_consistency <- function(expression_matrix, probe_gene_map, gene) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Get probes for this gene
  probes <- probe_gene_map$probe_id[probe_gene_map$gene_name == gene]
  probes <- intersect(probes, rownames(expression_matrix))

  if (length(probes) < 2) {
    stop("Gene ", gene, " has fewer than 2 probes")
  }

  # Get expression values
  probe_expr <- expression_matrix[probes, , drop = FALSE]

  # Calculate correlation matrix
  cor_mat <- cor(t(probe_expr), use = "pairwise.complete.obs")

  # Convert to long format for ggplot
  cor_long <- reshape2::melt(cor_mat)
  colnames(cor_long) <- c("Probe1", "Probe2", "Correlation")

  # Create heatmap
  p <- ggplot2::ggplot(cor_long,
                       ggplot2::aes(x = Probe1, y = Probe2, fill = Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#e74c3c", mid = "white", high = "#2ecc71",
                                  midpoint = 0, limits = c(-1, 1)) +
    ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)), size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Probe Consistency for", gene),
      subtitle = paste(length(probes), "probes"),
      x = "Probe",
      y = "Probe"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}

#' Plot Expression Level Comparison
#'
#' Compares probe-level vs gene-level expression distributions
#'
#' @param probe_expr Probe-level expression matrix
#' @param gene_expr Gene-level expression matrix
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_expression_comparison(probe_expr, gene_expr)
#' }
plot_expression_comparison <- function(probe_expr, gene_expr) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Sample data if too large
  max_values <- 10000

  probe_values <- as.vector(probe_expr)
  if (length(probe_values) > max_values) {
    probe_values <- sample(probe_values, max_values)
  }

  gene_values <- as.vector(gene_expr)
  if (length(gene_values) > max_values) {
    gene_values <- sample(gene_values, max_values)
  }

  # Create data frame
  plot_data <- data.frame(
    Expression = c(probe_values, gene_values),
    Level = rep(c("Probe-level", "Gene-level"),
                c(length(probe_values), length(gene_values)))
  )

  # Remove NAs
  plot_data <- plot_data[!is.na(plot_data$Expression), ]

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Expression, fill = Level)) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::scale_fill_manual(values = c("Probe-level" = "#3498db",
                                          "Gene-level" = "#e74c3c")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Expression Distribution: Probe-level vs Gene-level",
      x = "Expression Value",
      y = "Density",
      fill = "Level"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "top"
    )

  return(p)
}

#' Plot Annotation Coverage
#'
#' Shows what percentage of probes have annotations
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_annotation_coverage(qc_results)
#' }
plot_annotation_coverage <- function(qc_results) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Calculate coverage metrics
  coverage_data <- data.frame(
    Category = c("Has Genomic Location", "No Genomic Location",
                 "Has Gene Annotation", "No Gene Annotation",
                 "Passes QC", "Fails QC"),
    Count = c(
      sum(qc_results$n_locations > 0),
      sum(qc_results$n_locations == 0),
      sum(qc_results$n_genes > 0),
      sum(qc_results$n_genes == 0),
      sum(qc_results$qc_status == "PASS"),
      sum(qc_results$qc_status != "PASS")
    )
  )

  coverage_data$Percentage <- round(100 * coverage_data$Count / nrow(qc_results), 1)

  # Add grouping
  coverage_data$Group <- rep(c("Genomic Mapping", "Gene Annotation", "QC Status"),
                             each = 2)

  # Create plot
  p <- ggplot2::ggplot(coverage_data,
                       ggplot2::aes(x = Category, y = Percentage, fill = Category)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Percentage, "%\n(n=", Count, ")")),
                       vjust = -0.3, size = 3) +
    ggplot2::facet_wrap(~ Group, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Annotation Coverage Summary",
      subtitle = paste("Total probes:", nrow(qc_results)),
      x = "",
      y = "Percentage of Probes"
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::ylim(0, 105)

  return(p)
}

#' Create QC Dashboard
#'
#' Generates a multi-panel QC visualization
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param alignments GRanges object from align_probes_to_genome
#' @param annotations Data frame from update_gene_annotations
#' @return Combined ggplot object (requires patchwork or gridExtra)
#' @export
#' @examples
#' \dontrun{
#'   dashboard <- create_qc_dashboard(qc_results, alignments, annotations)
#'   print(dashboard)
#' }
create_qc_dashboard <- function(qc_results, alignments, annotations) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Generate individual plots
  p1 <- plot_qc_status(qc_results)
  p2 <- plot_alignment_stats(alignments)
  p3 <- plot_gene_biotypes(annotations)
  p4 <- plot_annotation_coverage(qc_results)

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p1 | p2) / (p3 | p4) +
      patchwork::plot_annotation(
        title = "PRISM Quality Control Dashboard",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
      )
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2,
                                        top = "PRISM Quality Control Dashboard")
  } else {
    warning("Install 'patchwork' or 'gridExtra' for combined plots. Returning list instead.")
    combined <- list(
      qc_status = p1,
      alignment_stats = p2,
      gene_biotypes = p3,
      coverage = p4
    )
  }

  return(combined)
}

#' Plot Gene Expression Heatmap
#'
#' Creates a heatmap of gene expression across samples
#'
#' @param gene_expr Gene expression matrix
#' @param top_genes Number of top variable genes to plot (default: 50)
#' @param scale Scale rows ("row"), columns ("column"), or none ("none")
#' @param cluster_rows Cluster rows (default: TRUE)
#' @param cluster_cols Cluster columns (default: TRUE)
#' @return Heatmap plot
#' @export
#' @examples
#' \dontrun{
#'   plot_gene_heatmap(gene_expr, top_genes = 100)
#' }
plot_gene_heatmap <- function(gene_expr,
                              top_genes = 50,
                              scale = "row",
                              cluster_rows = TRUE,
                              cluster_cols = TRUE) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with: install.packages('pheatmap')")
  }

  # Calculate variance for each gene
  gene_var <- apply(gene_expr, 1, var, na.rm = TRUE)

  # Select top variable genes
  top_var_genes <- names(sort(gene_var, decreasing = TRUE)[1:min(top_genes, length(gene_var))])

  expr_subset <- gene_expr[top_var_genes, , drop = FALSE]

  # Remove genes with NA
  expr_subset <- expr_subset[complete.cases(expr_subset), , drop = FALSE]

  # Create heatmap
  pheatmap::pheatmap(
    expr_subset,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    main = paste("Top", nrow(expr_subset), "Variable Genes"),
    fontsize_row = 8,
    fontsize_col = 8,
    color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100)
  )
}

#' Plot Probe Locations on Genome
#'
#' Visualizes probe locations along a chromosome
#'
#' @param alignments GRanges object from align_probes_to_genome
#' @param chromosome Chromosome to plot (e.g., "chr1")
#' @param start Start position (optional)
#' @param end End position (optional)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_probe_locations(alignments, "chr1", start = 1e6, end = 2e6)
#' }
plot_probe_locations <- function(alignments, chromosome, start = NULL, end = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Filter to chromosome
  chr_alignments <- alignments[seqnames(alignments) == chromosome]

  if (length(chr_alignments) == 0) {
    stop("No alignments found on chromosome ", chromosome)
  }

  # Filter to region if specified
  if (!is.null(start) && !is.null(end)) {
    chr_alignments <- chr_alignments[start(chr_alignments) >= start &
                                       end(chr_alignments) <= end]

    if (length(chr_alignments) == 0) {
      stop("No alignments found in specified region")
    }
  }

  # Convert to data frame
  plot_data <- data.frame(
    start = start(chr_alignments),
    end = end(chr_alignments),
    strand = as.character(strand(chr_alignments)),
    probe_id = chr_alignments$probe_id
  )

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = start, xend = end,
                                               y = 1, yend = 1, color = strand)) +
    ggplot2::geom_segment(size = 2, alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("+" = "#3498db", "-" = "#e74c3c")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Probe Locations on", chromosome),
      subtitle = paste(length(chr_alignments), "alignments"),
      x = "Genomic Position (bp)",
      y = "",
      color = "Strand"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma)

  return(p)
}

#' Plot Sample Correlation Heatmap
#'
#' Shows correlation between samples before and after remapping
#'
#' @param expression_matrix Expression matrix (probe or gene level)
#' @param method Correlation method: "pearson", "spearman" (default: "pearson")
#' @param title Plot title
#' @return Heatmap plot
#' @export
#' @examples
#' \dontrun{
#'   plot_sample_correlation(gene_expr, title = "Gene-level Correlations")
#' }
plot_sample_correlation <- function(expression_matrix,
                                    method = "pearson",
                                    title = "Sample Correlation") {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.")
  }

  # Calculate correlation
  cor_mat <- cor(expression_matrix, use = "pairwise.complete.obs", method = method)

  # Create heatmap
  pheatmap::pheatmap(
    cor_mat,
    main = title,
    display_numbers = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "#3498db", "#2c3e50"))(100),
    breaks = seq(min(cor_mat, na.rm = TRUE), 1, length.out = 101)
  )
}

#' Plot PCA of Expression Data
#'
#' Principal component analysis visualization
#'
#' @param expression_matrix Expression matrix (genes x samples)
#' @param sample_groups Optional vector of sample group labels
#' @param n_top Number of top variable genes to use (default: 500)
#' @param pcs Principal components to plot (default: c(1, 2))
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_pca(gene_expr, sample_groups = c(rep("Control", 5), rep("Treatment", 5)))
#' }
plot_pca <- function(expression_matrix,
                     sample_groups = NULL,
                     n_top = 500,
                     pcs = c(1, 2)) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Select top variable genes
  gene_var <- apply(expression_matrix, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_var, decreasing = TRUE)[1:min(n_top, length(gene_var))])

  expr_subset <- expression_matrix[top_genes, , drop = FALSE]
  expr_subset <- expr_subset[complete.cases(expr_subset), , drop = FALSE]

  # Transpose for PCA (samples as rows)
  expr_t <- t(expr_subset)

  # Perform PCA
  pca_result <- prcomp(expr_t, scale. = TRUE, center = TRUE)

  # Calculate variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

  # Create plot data
  plot_data <- data.frame(
    PC1 = pca_result$x[, pcs[1]],
    PC2 = pca_result$x[, pcs[2]],
    Sample = rownames(pca_result$x)
  )

  # Add groups if provided
  if (!is.null(sample_groups)) {
    if (length(sample_groups) != nrow(plot_data)) {
      warning("Length of sample_groups does not match number of samples")
    } else {
      plot_data$Group <- sample_groups
    }
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA of Gene Expression",
      subtitle = paste("Top", nrow(expr_subset), "variable genes"),
      x = paste0("PC", pcs[1], " (", var_explained[pcs[1]], "%)"),
      y = paste0("PC", pcs[2], " (", var_explained[pcs[2]], "%)")
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  # Add points with or without grouping
  if ("Group" %in% colnames(plot_data)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = Group), size = 3, alpha = 0.7) +
      ggplot2::stat_ellipse(ggplot2::aes(color = Group), type = "norm", level = 0.95)
  } else {
    p <- p + ggplot2::geom_point(size = 3, alpha = 0.7, color = "#3498db")
  }

  return(p)
}

#' Plot Variance Explained by Remapping
#'
#' Shows how much variance is retained after probe-to-gene aggregation
#'
#' @param probe_expr Probe-level expression matrix
#' @param gene_expr Gene-level expression matrix
#' @param probe_gene_map Mapping from create_probe_gene_mapping
#' @param n_genes Number of genes to sample for analysis (default: 100)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   plot_variance_explained(probe_expr, gene_expr, mapping)
#' }
plot_variance_explained <- function(probe_expr, gene_expr,
                                    probe_gene_map, n_genes = 100) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Get genes with multiple probes
  multi_probe_genes <- probe_gene_map %>%
    dplyr::filter(n_probes > 1) %>%
    dplyr::pull(gene_name) %>%
    unique()

  # Sample genes if too many
  if (length(multi_probe_genes) > n_genes) {
    multi_probe_genes <- sample(multi_probe_genes, n_genes)
  }

  message("Analyzing variance for ", length(multi_probe_genes), " genes...")

  # Calculate variance explained for each gene
  variance_data <- lapply(multi_probe_genes, function(gene) {
    # Get probes
    probes <- probe_gene_map$probe_id[probe_gene_map$gene_name == gene]
    probes <- intersect(probes, rownames(probe_expr))

    if (length(probes) < 2 || !gene %in% rownames(gene_expr)) {
      return(NULL)
    }

    # Probe-level variance (total across all probes)
    probe_values <- probe_expr[probes, , drop = FALSE]
    total_var <- sum(apply(probe_values, 1, var, na.rm = TRUE))

    # Gene-level variance
    gene_values <- gene_expr[gene, ]
    gene_var <- var(gene_values, na.rm = TRUE)

    data.frame(
      gene = gene,
      n_probes = length(probes),
      total_probe_var = total_var,
      gene_var = gene_var,
      var_ratio = gene_var / total_var
    )
  })

  variance_df <- do.call(rbind, variance_data)
  variance_df <- variance_df[!is.na(variance_df$var_ratio), ]

  # Create plot
  p <- ggplot2::ggplot(variance_df,
                       ggplot2::aes(x = n_probes, y = var_ratio)) +
    ggplot2::geom_point(alpha = 0.5, color = "#3498db") +
    ggplot2::geom_smooth(method = "loess", color = "#e74c3c") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Variance Retained After Probe Aggregation",
      subtitle = paste("Analyzed", nrow(variance_df), "genes"),
      x = "Number of Probes per Gene",
      y = "Variance Ratio (Gene / Total Probe)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    ) +
    ggplot2::ylim(0, 1)

  return(p)
}

#' Save All QC Plots to PDF
#'
#' Generates a comprehensive PDF report with all QC visualizations
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param alignments GRanges object from align_probes_to_genome
#' @param annotations Data frame from update_gene_annotations
#' @param probe_expr Optional probe expression matrix
#' @param gene_expr Optional gene expression matrix
#' @param probe_gene_map Optional probe-gene mapping
#' @param output_file Output PDF file path (default: "PRISM_QC_report.pdf")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @export
#' @examples
#' \dontrun{
#'   save_qc_plots(
#'     qc_results,
#'     alignments,
#'     annotations,
#'     output_file = "my_qc_report.pdf"
#'   )
#' }
save_qc_plots <- function(qc_results,
                          alignments,
                          annotations,
                          probe_expr = NULL,
                          gene_expr = NULL,
                          probe_gene_map = NULL,
                          output_file = "PRISM_QC_report.pdf",
                          width = 10,
                          height = 8) {

  message("Generating QC plots and saving to: ", output_file)

  pdf(output_file, width = width, height = height)

  # Page 1: QC Status
  print(plot_qc_status(qc_results))

  # Page 2: Alignment Stats
  print(plot_alignment_stats(alignments))

  # Page 3: Gene Biotypes
  if ("gene_biotype" %in% colnames(annotations)) {
    print(plot_gene_biotypes(annotations))
  }

  # Page 4: Annotation Coverage
  print(plot_annotation_coverage(qc_results))

  # Page 5: Chromosome Distribution
  print(plot_chromosome_distribution(alignments))

  # Page 6: Probes per Gene
  if (!is.null(probe_gene_map)) {
    print(plot_probes_per_gene(probe_gene_map))
  }

  # Page 7: Expression Comparison
  if (!is.null(probe_expr) && !is.null(gene_expr)) {
    print(plot_expression_comparison(probe_expr, gene_expr))
  }

  # Page 8: PCA
  if (!is.null(gene_expr)) {
    print(plot_pca(gene_expr))
  }

  dev.off()

  message("QC report saved successfully!")

  invisible(output_file)
}

#' Save All QC Plots to PDF
#'
#' Generates a comprehensive PDF report with all QC visualizations
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param alignments GRanges object from align_probes_to_genome
#' @param annotations Data frame from update_gene_annotations
#' @param probe_expr Optional probe expression matrix
#' @param gene_expr Optional gene expression matrix
#' @param probe_gene_map Optional probe-gene mapping
#' @param output_file Output PDF file path (default: "PRISM_QC_report.pdf")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @export
#' @examples
#' \dontrun{
#'   save_qc_plots(
#'     qc_results,
#'     alignments,
#'     annotations,
#'     output_file = "my_qc_report.pdf"
#'   )
#' }
save_qc_plots <- function(qc_results,
                         alignments,
                         annotations,
                         probe_expr = NULL,
                         gene_expr = NULL,
                         probe_gene_map = NULL,
                         output_file = "PRISM_QC_report.pdf",
                         width = 10,
                         height = 8) {

  message("Generating QC plots and saving to: ", output_file)

  pdf(output_file, width = width, height = height)

  # Page 1: QC Status
  print(plot_qc_status(qc_results))

  # Page 2: Alignment Stats
  print(plot_alignment_stats(alignments))

  # Page 3: Gene Biotypes
  if ("gene_biotype" %in% colnames(annotations)) {
    print(plot_gene_biotypes(annotations))
  }

  # Page 4: Annotation Coverage
  print(plot_annotation_coverage(qc_results))

  # Page 5: Chromosome Distribution
  print(plot_chromosome_distribution(alignments))

  # Page 6: Probes per Gene
  if (!is.null(probe_gene_map)) {
    print(plot_probes_per_gene(probe_gene_map))
  }

  # Page 7: Expression Comparison
  if (!is.null(probe_expr) && !is.null(gene_expr)) {
    print(plot_expression_comparison(probe_expr, gene_expr))
  }

  # Page 8: PCA
  if (!is.null(gene_expr)) {
    print(plot_pca(gene_expr))
  }

  dev.off()

  message("QC report saved successfully!")

  invisible(output_file)
}

#' Create Interactive QC Plot
#'
#' Creates an interactive plotly visualization
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @return plotly object
#' @export
#' @examples
#' \dontrun{
#'   interactive_plot <- plot_qc_interactive(qc_results)
#'   interactive_plot
#' }
plot_qc_interactive <- function(qc_results) {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  }

  # Create base ggplot
  p <- plot_qc_status(qc_results)

  # Convert to plotly
  interactive_p <- plotly::ggplotly(p, tooltip = c("x", "y", "fill"))

  return(interactive_p)
}

#' Plot Annotation Change Timeline
#'
#' Shows annotation changes over different genome builds or time points
#'
#' @param comparison_list List of comparison results from different versions
#' @param version_names Names for each version
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'   comparisons <- list(
#'     hg19_vs_hg38 = comparison1,
#'     original_vs_current = comparison2
#'   )
#'   plot_annotation_timeline(comparisons)
#' }
plot_annotation_timeline <- function(comparison_list, version_names = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (is.null(version_names)) {
    version_names <- names(comparison_list)
  }

  # Extract change statistics from each comparison
  timeline_data <- lapply(seq_along(comparison_list), function(i) {
    comp <- comparison_list[[i]]

    data.frame(
      Version = version_names[i],
      Unchanged = sum(comp$annotation_status == "unchanged"),
      Changed = sum(comp$annotation_status == "changed"),
      Gained = sum(comp$annotation_status == "newly_annotated"),
      Lost = sum(comp$annotation_status == "annotation_lost")
    )
  })

  timeline_df <- do.call(rbind, timeline_data)

  # Convert to long format
  timeline_long <- tidyr::pivot_longer(
    timeline_df,
    cols = c("Unchanged", "Changed", "Gained", "Lost"),
    names_to = "Category",
    values_to = "Count"
  )

  # Create plot
  p <- ggplot2::ggplot(timeline_long,
                       ggplot2::aes(x = Version, y = Count, fill = Category)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = c(
      "Unchanged" = "#2ecc71",
      "Changed" = "#f39c12",
      "Gained" = "#3498db",
      "Lost" = "#e74c3c"
    )) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Annotation Changes Across Versions",
      x = "Comparison",
      y = "Number of Probes",
      fill = "Status"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}

#' Plot Diagnostic for Specific Gene
#'
#' Shows detailed information about probe-gene mapping for a specific gene
#'
#' @param gene Gene name to visualize
#' @param probe_expr Probe expression matrix
#' @param gene_expr Gene expression matrix
#' @param probe_gene_map Probe-gene mapping
#' @param alignments GRanges with probe alignments
#' @return Combined plot
#' @export
#' @examples
#' \dontrun{
#'   plot_gene_diagnostic("TP53", probe_expr, gene_expr, mapping, alignments)
#' }
plot_gene_diagnostic <- function(gene, probe_expr, gene_expr,
                                 probe_gene_map, alignments) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Get probes for this gene
  probes <- probe_gene_map$probe_id[probe_gene_map$gene_name == gene]
  probes <- intersect(probes, rownames(probe_expr))

  if (length(probes) == 0) {
    stop("No probes found for gene: ", gene)
  }

  # Get probe expression
  probe_values <- probe_expr[probes, , drop = FALSE]

  # Get gene expression
  if (gene %in% rownames(gene_expr)) {
    gene_values <- gene_expr[gene, ]
  } else {
    stop("Gene not found in gene expression matrix: ", gene)
  }

  # Create long format data
  probe_long <- reshape2::melt(probe_values)
  colnames(probe_long) <- c("Probe", "Sample", "Expression")
  probe_long$Type <- "Probe"

  gene_long <- data.frame(
    Probe = gene,
    Sample = names(gene_values),
    Expression = gene_values,
    Type = "Gene"
  )

  plot_data <- rbind(probe_long, gene_long)

  # Create boxplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Probe, y = Expression, fill = Type)) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(values = c("Probe" = "#3498db", "Gene" = "#e74c3c")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Expression Profile for", gene),
      subtitle = paste(length(probes), "probes"),
      x = "",
      y = "Expression"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}

#' Get Available Visualization Functions
#'
#' Lists all available plotting functions in PRISM
#'
#' @return Character vector of function names
#' @export
list_visualization_functions <- function() {

  viz_functions <- c(
    "plot_qc_status" = "QC status distribution",
    "plot_alignment_stats" = "Alignment statistics",
    "plot_annotation_changes" = "Annotation changes pie chart",
    "plot_gene_biotypes" = "Gene biotype distribution",
    "plot_probes_per_gene" = "Probes per gene distribution",
    "plot_chromosome_distribution" = "Probe distribution across chromosomes",
    "plot_probe_consistency" = "Probe correlation heatmap for a gene",
    "plot_expression_comparison" = "Probe vs gene expression distributions",
    "plot_annotation_coverage" = "Annotation coverage summary",
    "create_qc_dashboard" = "Multi-panel QC dashboard",
    "plot_gene_heatmap" = "Gene expression heatmap",
    "plot_probe_locations" = "Probe locations on genome",
    "plot_sample_correlation" = "Sample correlation heatmap",
    "plot_pca" = "PCA plot",
    "plot_variance_explained" = "Variance retained after aggregation",
    "save_qc_plots" = "Save all plots to PDF",
    "plot_qc_interactive" = "Interactive plotly visualization",
    "plot_annotation_timeline" = "Annotation changes over versions",
    "plot_gene_diagnostic" = "Detailed gene diagnostic plot"
  )

  cat("PRISM Visualization Functions\n")
  cat("==============================\n\n")

  for (i in seq_along(viz_functions)) {
    cat(sprintf("%-30s : %s\n", names(viz_functions)[i], viz_functions[i]))
  }

  invisible(names(viz_functions))
}
