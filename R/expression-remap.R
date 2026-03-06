#' Remap Expression Matrix with Updated Annotations
#'
#' Collapses probe-level expression data to gene-level using updated annotations.
#' Handles multiple probes per gene and provides various aggregation methods.
#'
#' @param expression_matrix Matrix with probes as rows, samples as columns
#' @param updated_annotations Data frame with probe-to-gene mappings
#' @param qc_flags QC results from flag_problematic_probes (optional)
#' @param method Aggregation method: "mean", "median", "max", "iqr" (default: "median")
#' @param remove_ambiguous Remove ambiguous/multi-mapping probes (default: TRUE)
#' @param min_probes Minimum probes required per gene (default: 1)
#' @param log_transform Apply log2 transformation before aggregation (default: FALSE)
#' @return Matrix with genes as rows, samples as columns
#' @export
#' @examples
#' \dontrun{
#'   # Basic remapping with median aggregation
#'   gene_expr <- remap_expression_matrix(
#'     expression_matrix,
#'     annotations
#'   )
#'
#'   # Use mean and filter QC
#'   gene_expr <- remap_expression_matrix(
#'     expression_matrix,
#'     annotations,
#'     qc_flags = qc_results,
#'     method = "mean",
#'     remove_ambiguous = TRUE
#'   )
#' }
remap_expression_matrix <- function(expression_matrix,
                                    updated_annotations,
                                    qc_flags = NULL,
                                    method = "median",
                                    remove_ambiguous = TRUE,
                                    min_probes = 1,
                                    log_transform = FALSE) {

  method <- match.arg(method, c("mean", "median", "max", "iqr", "first"))

  message("Remapping expression matrix to gene-level...")
  message("  Input: ", nrow(expression_matrix), " probes × ",
          ncol(expression_matrix), " samples")
  message("  Method: ", method)

  # Validate inputs
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data frame")
  }

  if (nrow(updated_annotations) == 0) {
    stop("updated_annotations is empty")
  }

  # Convert to matrix if needed
  if (is.data.frame(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  # Filter by QC if provided
  if (!is.null(qc_flags)) {
    message("  Applying QC filters...")

    if (remove_ambiguous) {
      passing_probes <- qc_flags$probe_id[qc_flags$qc_status == "PASS"]
    } else {
      # Keep everything except deprecated
      passing_probes <- qc_flags$probe_id[qc_flags$qc_status != "DEPRECATED"]
    }

    # Filter annotations
    updated_annotations <- updated_annotations[
      updated_annotations$probe_id %in% passing_probes,
    ]

    message("  Retained ", nrow(updated_annotations), " probe-gene mappings after QC")
  }

  # Match probes between expression matrix and annotations
  probe_ids <- rownames(expression_matrix)

  if (is.null(probe_ids)) {
    stop("expression_matrix must have probe IDs as rownames")
  }

  # Find common probes
  common_probes <- intersect(probe_ids, updated_annotations$probe_id)

  if (length(common_probes) == 0) {
    stop("No matching probes found between expression matrix and annotations")
  }

  message("  Matched ", length(common_probes), " probes between expression and annotations")

  # Subset to common probes
  expr_subset <- expression_matrix[common_probes, , drop = FALSE]
  annot_subset <- updated_annotations[updated_annotations$probe_id %in% common_probes, ]

  # Log transform if requested
  if (log_transform) {
    message("  Applying log2 transformation...")
    # Add small constant to avoid log(0)
    expr_subset <- log2(expr_subset + 1)
  }

  # Create probe-to-gene mapping
  probe_gene_map <- annot_subset[, c("probe_id", "gene_name", "ensembl_gene_id")]

  # Remove duplicates (in case probe maps to same gene multiple times)
  probe_gene_map <- unique(probe_gene_map)

  # Count probes per gene
  probes_per_gene <- table(probe_gene_map$gene_name)

  message("  Aggregating expression data...")
  message("    Genes with 1 probe: ", sum(probes_per_gene == 1))
  message("    Genes with 2+ probes: ", sum(probes_per_gene > 1))
  message("    Max probes per gene: ", max(probes_per_gene))

  # Get unique genes
  unique_genes <- unique(probe_gene_map$gene_name)
  unique_genes <- unique_genes[!is.na(unique_genes)]

  # Initialize gene expression matrix
  gene_expr <- matrix(
    NA,
    nrow = length(unique_genes),
    ncol = ncol(expr_subset),
    dimnames = list(unique_genes, colnames(expr_subset))
  )

  # Aggregate for each gene
  for (gene in unique_genes) {
    # Get probes for this gene
    gene_probes <- probe_gene_map$probe_id[probe_gene_map$gene_name == gene]

    # Skip if too few probes
    if (length(gene_probes) < min_probes) {
      next
    }

    # Get expression values
    if (length(gene_probes) == 1) {
      gene_values <- expr_subset[gene_probes, , drop = TRUE]
    } else {
      gene_values <- expr_subset[gene_probes, , drop = FALSE]

      # Apply aggregation method
      gene_values <- switch(method,
                            "mean" = colMeans(gene_values, na.rm = TRUE),
                            "median" = apply(gene_values, 2, median, na.rm = TRUE),
                            "max" = apply(gene_values, 2, max, na.rm = TRUE),
                            "iqr" = apply(gene_values, 2, function(x) {
                              # Use probe with highest IQR (most variable)
                              iqrs <- apply(gene_values, 1, IQR, na.rm = TRUE)
                              gene_values[which.max(iqrs), ]
                            }),
                            "first" = gene_values[1, ]
      )
    }

    gene_expr[gene, ] <- gene_values
  }

  # Remove genes with all NA (didn't meet min_probes)
  genes_with_data <- rowSums(!is.na(gene_expr)) > 0
  gene_expr <- gene_expr[genes_with_data, , drop = FALSE]

  message("  Output: ", nrow(gene_expr), " genes × ",
          ncol(gene_expr), " samples")

  return(gene_expr)
}

#' Create Probe-to-Gene Mapping Table
#'
#' Generates a lookup table showing which probes map to which genes
#'
#' @param updated_annotations Data frame from update_gene_annotations
#' @param qc_flags Optional QC results to include status
#' @return Data frame with probe-gene mappings and counts
#' @export
#' @examples
#' \dontrun{
#'   mapping <- create_probe_gene_mapping(annotations, qc_results)
#' }
create_probe_gene_mapping <- function(updated_annotations, qc_flags = NULL) {

  message("Creating probe-to-gene mapping table...")

  # Basic mapping
  mapping <- updated_annotations %>%
    dplyr::select(probe_id, ensembl_gene_id, gene_name, gene_biotype) %>%
    dplyr::distinct()

  # Add QC status if available
  if (!is.null(qc_flags)) {
    qc_info <- qc_flags %>%
      dplyr::select(probe_id, qc_status, n_locations, n_genes)

    mapping <- dplyr::left_join(mapping, qc_info, by = "probe_id")
  }

  # Count genes per probe
  genes_per_probe <- mapping %>%
    dplyr::group_by(probe_id) %>%
    dplyr::summarise(n_genes_mapped = dplyr::n(),.groups = "drop")

  mapping <- dplyr::left_join(mapping, genes_per_probe, by = "probe_id")

  # Count probes per gene
  probes_per_gene <- mapping %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(n_probes = dplyr::n_distinct(probe_id),.groups = "drop")

  mapping <- dplyr::left_join(mapping, probes_per_gene, by = "gene_name")

  message("  Created mapping for ", nrow(mapping), " probe-gene pairs")
  message("  Unique probes: ", length(unique(mapping$probe_id)))
  message("  Unique genes: ", length(unique(mapping$gene_name)))

  return(as.data.frame(mapping))
}

#' Compare Probe vs Gene Expression
#'
#' Compares expression patterns before and after remapping
#'
#' @param probe_expr Matrix with probe-level expression
#' @param gene_expr Matrix with gene-level expression
#' @param probe_gene_map Mapping from create_probe_gene_mapping
#' @param sample_subset Optional vector of sample names to compare
#' @return List with comparison statistics
#' @export
#' @examples
#' \dontrun{
#'   comparison <- compare_expression_levels(
#'     probe_expression,
#'     gene_expression,
#'     mapping
#'   )
#' }
compare_expression_levels <- function(probe_expr,
                                      gene_expr,
                                      probe_gene_map,
                                      sample_subset = NULL) {

  message("Comparing probe-level vs gene-level expression...")

  # Subset samples if requested
  if (!is.null(sample_subset)) {
    probe_expr <- probe_expr[, sample_subset, drop = FALSE]
    gene_expr <- gene_expr[, sample_subset, drop = FALSE]
  }

  # Calculate summary statistics
  probe_stats <- list(
    n_features = nrow(probe_expr),
    mean_expr = mean(probe_expr, na.rm = TRUE),
    median_expr = median(probe_expr, na.rm = TRUE),
    sd_expr = sd(as.vector(probe_expr), na.rm = TRUE),
    range_expr = range(probe_expr, na.rm = TRUE)
  )

  gene_stats <- list(
    n_features = nrow(gene_expr),
    mean_expr = mean(gene_expr, na.rm = TRUE),
    median_expr = median(gene_expr, na.rm = TRUE),
    sd_expr = sd(as.vector(gene_expr), na.rm = TRUE),
    range_expr = range(gene_expr, na.rm = TRUE)
  )

  # Calculate correlation for genes with multiple probes
  genes_multi_probe <- probe_gene_map %>%
    dplyr::filter(n_probes > 1) %>%
    dplyr::pull(gene_name) %>%
    unique()

  if (length(genes_multi_probe) > 0) {
    # Sample a few genes to check correlation
    sample_genes <- sample(genes_multi_probe, min(10, length(genes_multi_probe)))

    correlations <- sapply(sample_genes, function(gene) {
      probes <- probe_gene_map$probe_id[probe_gene_map$gene_name == gene]
      probes <- intersect(probes, rownames(probe_expr))

      if (length(probes) < 2) return(NA)

      probe_vals <- probe_expr[probes, , drop = FALSE]
      gene_val <- gene_expr[gene, ]

      # Correlation between gene value and mean of probes
      cor(colMeans(probe_vals, na.rm = TRUE), gene_val, use = "complete.obs")
    })

    mean_correlation <- mean(correlations, na.rm = TRUE)
  } else {
    mean_correlation <- NA
  }

  comparison <- list(
    probe_level = probe_stats,
    gene_level = gene_stats,
    reduction_ratio = nrow(probe_expr) / nrow(gene_expr),
    mean_correlation = mean_correlation
  )

  class(comparison) <- c("prism_expr_comparison", "list")

  return(comparison)
}

#' @export
print.prism_expr_comparison <- function(x,...) {
  cat("PRISM Expression Comparison\n")
  cat("===========================\n\n")

  cat("Probe-level:\n")
  cat("  Features:", x$probe_level$n_features, "\n")
  cat("  Mean expression:", round(x$probe_level$mean_expr, 3), "\n")
  cat("  Median expression:", round(x$probe_level$median_expr, 3), "\n")
  cat("  SD:", round(x$probe_level$sd_expr, 3), "\n")
  cat("  Range:", round(x$probe_level$range_expr[1], 3), "to",
      round(x$probe_level$range_expr[2], 3), "\n\n")

  cat("Gene-level:\n")
  cat("  Features:", x$gene_level$n_features, "\n")
  cat("  Mean expression:", round(x$gene_level$mean_expr, 3), "\n")
  cat("  Median expression:", round(x$gene_level$median_expr, 3), "\n")
  cat("  SD:", round(x$gene_level$sd_expr, 3), "\n")
  cat("  Range:", round(x$gene_level$range_expr[1], 3), "to",
      round(x$gene_level$range_expr[2], 3), "\n\n")

  cat("Reduction ratio:", round(x$reduction_ratio, 2), "probes per gene\n")

  if (!is.na(x$mean_correlation)) {
    cat("Mean probe-gene correlation:", round(x$mean_correlation, 3), "\n")
  }

  invisible(x)
}

#' Batch Remap Multiple Expression Matrices
#'
#' Applies the same annotation and remapping to multiple datasets
#'
#' @param expression_list List of expression matrices
#' @param updated_annotations Data frame with probe-to-gene mappings
#' @param qc_flags QC results (optional)
#' @param method Aggregation method (default: "median")
#' @param remove_ambiguous Remove ambiguous probes (default: TRUE)
#' @param keep_common_genes Keep only genes common to all datasets (default: FALSE)
#' @return List of remapped gene expression matrices
#' @export
#' @examples
#' \dontrun{
#'   datasets <- list(
#'     study1 = expr_matrix1,
#'     study2 = expr_matrix2,
#'     study3 = expr_matrix3
#'   )
#'
#'   remapped <- batch_remap_expression(
#'     datasets,
#'     annotations,
#'     qc_flags = qc_results
#'   )
#' }
batch_remap_expression <- function(expression_list,
                                   updated_annotations,
                                   qc_flags = NULL,
                                   method = "median",
                                   remove_ambiguous = TRUE,
                                   keep_common_genes = FALSE) {

  if (!is.list(expression_list)) {
    stop("expression_list must be a list of matrices")
  }

  message("Batch remapping ", length(expression_list), " datasets...")

  # Remap each dataset
  remapped_list <- lapply(names(expression_list), function(name) {
    message("\nProcessing: ", name)

    tryCatch({
      remap_expression_matrix(
        expression_matrix = expression_list[[name]],
        updated_annotations = updated_annotations,
        qc_flags = qc_flags,
        method = method,
        remove_ambiguous = remove_ambiguous
      )
    }, error = function(e) {
      warning("Failed to remap ", name, ": ", e$message)
      return(NULL)
    })
  })

  names(remapped_list) <- names(expression_list)

  # Remove failed datasets
  remapped_list <- remapped_list[!sapply(remapped_list, is.null)]

  # Keep only common genes if requested
  if (keep_common_genes && length(remapped_list) > 1) {
    message("\nFinding common genes across datasets...")

    all_genes <- lapply(remapped_list, rownames)
    common_genes <- Reduce(intersect, all_genes)

    message("  Common genes: ", length(common_genes))

    remapped_list <- lapply(remapped_list, function(mat) {
      mat[common_genes, , drop = FALSE]
    })
  }

  message("\nBatch remapping complete!")
  message("  Successfully remapped: ", length(remapped_list), " datasets")

  return(remapped_list)
}
#' Merge Multiple Gene Expression Matrices
#'
#' Combines remapped expression matrices from different studies
#'
#' @param expression_list List of gene expression matrices
#' @param merge_method Method for handling genes: "common", "union" (default: "common")
#' @param add_batch_prefix Add dataset name prefix to sample IDs (default: TRUE)
#' @return Single merged expression matrix
#' @export
#' @examples
#' \dontrun{
#'   merged <- merge_expression_matrices(
#'     remapped_datasets,
#'     merge_method = "common"
#'   )
#' }
merge_expression_matrices <- function(expression_list,
                                      merge_method = "common",
                                      add_batch_prefix = TRUE) {

  merge_method <- match.arg(merge_method, c("common", "union"))

  if (length(expression_list) < 2) {
    stop("Need at least 2 matrices to merge")
  }

  message("Merging ", length(expression_list), " expression matrices...")
  message("  Method: ", merge_method)

  # Add batch prefixes to sample names if requested
  if (add_batch_prefix) {
    expression_list <- lapply(names(expression_list), function(name) {
      mat <- expression_list[[name]]
      colnames(mat) <- paste0(name, "_", colnames(mat))
      return(mat)
    })
    names(expression_list) <- names(expression_list)
  }

  # Get gene lists
  all_genes_list <- lapply(expression_list, rownames)

  if (merge_method == "common") {
    # Keep only common genes
    genes_to_keep <- Reduce(intersect, all_genes_list)
    message("  Common genes: ", length(genes_to_keep))

    if (length(genes_to_keep) == 0) {
      stop("No common genes found across datasets")
    }

    # Subset each matrix
    expression_list <- lapply(expression_list, function(mat) {
      mat[genes_to_keep, , drop = FALSE]
    })

  } else {
    # Union: keep all genes, fill missing with NA
    genes_to_keep <- Reduce(union, all_genes_list)
    message("  Total unique genes: ", length(genes_to_keep))

    # Expand each matrix to include all genes
    expression_list <- lapply(expression_list, function(mat) {
      missing_genes <- setdiff(genes_to_keep, rownames(mat))

      if (length(missing_genes) > 0) {
        # Create matrix of NAs for missing genes
        na_mat <- matrix(
          NA,
          nrow = length(missing_genes),
          ncol = ncol(mat),
          dimnames = list(missing_genes, colnames(mat))
        )

        # Combine
        mat <- rbind(mat, na_mat)
      }

      # Ensure same gene order
      mat[genes_to_keep, , drop = FALSE]
    })
  }

  # Merge by column binding
  merged <- do.call(cbind, expression_list)

  message("  Merged matrix: ", nrow(merged), " genes × ",
          ncol(merged), " samples")

  return(merged)
}

#' Export Remapped Expression Matrix
#'
#' Saves gene expression matrix to file
#'
#' @param gene_expr Gene expression matrix
#' @param file Output file path
#' @param format File format: "csv", "tsv", "rds" (default: "csv")
#' @param include_gene_info Include gene metadata columns (default: FALSE)
#' @param annotations Optional annotations to add metadata
#' @export
#' @examples
#' \dontrun{
#'   export_expression_matrix(gene_expr, "gene_expression.csv")
#'
#'   # With gene info
#'   export_expression_matrix(
#'     gene_expr,
#'     "gene_expression.csv",
#'     include_gene_info = TRUE,
#'     annotations = annotations
#'   )
#' }
export_expression_matrix <- function(gene_expr,
                                     file,
                                     format = "csv",
                                     include_gene_info = FALSE,
                                     annotations = NULL) {

  format <- match.arg(format, c("csv", "tsv", "rds"))

  message("Exporting expression matrix to: ", file)

  # Convert to data frame
  expr_df <- as.data.frame(gene_expr)
  expr_df$gene_name <- rownames(expr_df)

  # Add gene info if requested
  if (include_gene_info && !is.null(annotations)) {
    gene_info <- annotations %>%
      dplyr::select(gene_name, ensembl_gene_id, gene_biotype, gene_description) %>%
      dplyr::distinct()

    expr_df <- dplyr::left_join(
      expr_df,
      gene_info,
      by = "gene_name"
    )

    # Reorder columns: gene info first, then samples
    sample_cols <- setdiff(colnames(expr_df),
                           c("gene_name", "ensembl_gene_id", "gene_biotype", "gene_description"))
    expr_df <- expr_df[, c("gene_name", "ensembl_gene_id", "gene_biotype",
                           "gene_description", sample_cols)]
  }

  # Export based on format
  if (format == "csv") {
    write.csv(expr_df, file, row.names = FALSE)
  } else if (format == "tsv") {
    write.table(expr_df, file, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (format == "rds") {
    saveRDS(gene_expr, file)  # Save as matrix, not data frame
  }

  message("Expression matrix exported successfully")

  invisible(expr_df)
}

#' Calculate Probe Consistency for Multi-Probe Genes
#'
#' Evaluates how consistently probes for the same gene behave
#'
#' @param expression_matrix Probe-level expression matrix
#' @param probe_gene_map Mapping from create_probe_gene_mapping
#' @param min_probes Minimum probes per gene to evaluate (default: 2)
#' @return Data frame with consistency metrics per gene
#' @export
#' @examples
#' \dontrun{
#'   consistency <- calculate_probe_consistency(
#'     probe_expression,
#'     mapping
#'   )
#' }
calculate_probe_consistency <- function(expression_matrix,
                                        probe_gene_map,
                                        min_probes = 2) {

  message("Calculating probe consistency for multi-probe genes...")

  # Get genes with multiple probes
  multi_probe_genes <- probe_gene_map %>%
    dplyr::filter(n_probes >= min_probes) %>%
    dplyr::pull(gene_name) %>%
    unique()

  message("  Evaluating ", length(multi_probe_genes), " genes with ",
          min_probes, "+ probes")

  # Calculate consistency for each gene
  consistency_list <- lapply(multi_probe_genes, function(gene) {
    # Get probes for this gene
    probes <- probe_gene_map$probe_id[probe_gene_map$gene_name == gene]
    probes <- intersect(probes, rownames(expression_matrix))

    if (length(probes) < min_probes) return(NULL)

    # Get expression values
    probe_expr <- expression_matrix[probes, , drop = FALSE]

    # Calculate pairwise correlations
    cor_mat <- cor(t(probe_expr), use = "pairwise.complete.obs")
    cor_values <- cor_mat[upper.tri(cor_mat)]

    # Calculate coefficient of variation across probes for each sample
    cv_per_sample <- apply(probe_expr, 2, function(x) {
      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
    })

    data.frame(
      gene_name = gene,
      n_probes = length(probes),
      mean_correlation = mean(cor_values, na.rm = TRUE),
      median_correlation = median(cor_values, na.rm = TRUE),
      min_correlation = min(cor_values, na.rm = TRUE),
      mean_cv = mean(cv_per_sample, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  consistency <- do.call(rbind, consistency_list)

  # Flag low consistency genes
  consistency$low_consistency <- consistency$median_correlation < 0.5

  message("\nConsistency Summary:")
  message("  Mean correlation: ", round(mean(consistency$mean_correlation, na.rm = TRUE), 3))
  message("  Median correlation: ", round(median(consistency$median_correlation, na.rm = TRUE), 3))
  message("  Low consistency genes (<0.5): ", sum(consistency$low_consistency))

  return(consistency)
}

#' Get Remapping Summary
#'
#' Summarizes the remapping process
#'
#' @param probe_expr Original probe expression matrix
#' @param gene_expr Remapped gene expression matrix
#' @param probe_gene_map Probe-gene mapping
#' @return List with summary statistics
#' @export
get_remapping_summary <- function(probe_expr, gene_expr, probe_gene_map) {

  summary_list <- list(
    n_probes_input = nrow(probe_expr),
    n_probes_mapped = length(unique(probe_gene_map$probe_id)),
    n_genes_output = nrow(gene_expr),
    n_samples = ncol(gene_expr),
    reduction_ratio = nrow(probe_expr) / nrow(gene_expr),
    probes_per_gene = nrow(probe_expr) / nrow(gene_expr),
    mapping_rate = length(unique(probe_gene_map$probe_id)) / nrow(probe_expr)
  )

  class(summary_list) <- c("prism_remap_summary", "list")

  return(summary_list)
}

#' @export
print.prism_remap_summary <- function(x,...) {
  cat("PRISM Remapping Summary\n")
  cat("=======================\n\n")

  cat("Input:\n")
  cat("  Probes:", x$n_probes_input, "\n")
  cat("  Samples:", x$n_samples, "\n\n")

  cat("Mapping:\n")
  cat("  Probes mapped:", x$n_probes_mapped,
      sprintf("(%.1f%%)\n", 100 * x$mapping_rate))
  cat("  Genes output:", x$n_genes_output, "\n")
  cat("  Avg probes per gene:", round(x$probes_per_gene, 2), "\n")
  cat("  Reduction ratio:", round(x$reduction_ratio, 2), ":1\n")

  invisible(x)
}
