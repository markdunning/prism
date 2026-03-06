#' Flag Problematic Probes
#'
#' Identifies probes with quality issues including multiple mappings,
#' no mappings, cross-hybridization, and deprecated targets
#'
#' @param probe_alignments GRanges object from align_probes_to_genome
#' @param annotations Data frame from update_gene_annotations
#' @param original_annotations Optional: original probe annotations for comparison
#' @param max_locations Maximum allowed genomic locations (default: 1)
#' @param flag_intergenic Flag probes in intergenic regions (default: TRUE)
#' @return Data frame with probe quality flags
#' @export
#' @examples
#' \dontrun{
#'   qc_results <- flag_problematic_probes(
#'     alignments,
#'     annotations,
#'     max_locations = 1
#'   )
#' }
flag_problematic_probes <- function(probe_alignments,
                                    annotations,
                                    original_annotations = NULL,
                                    max_locations = 1,
                                    flag_intergenic = TRUE) {

  message("Running quality control on probe annotations...")

  # Get all unique probe IDs from alignments
  all_probe_ids <- unique(probe_alignments$probe_id)

  # Initialize QC data frame
  qc_report <- data.frame(
    probe_id = all_probe_ids,
    stringsAsFactors = FALSE
  )

  # Count genomic locations per probe
  message("  Checking genomic mapping locations...")
  location_counts <- as.data.frame(table(probe_alignments$probe_id))
  colnames(location_counts) <- c("probe_id", "n_locations")
  location_counts$probe_id <- as.character(location_counts$probe_id)

  qc_report <- dplyr::left_join(qc_report, location_counts, by = "probe_id")
  qc_report$n_locations[is.na(qc_report$n_locations)] <- 0
  # Flag multiple genomic locations
  qc_report$multiple_locations <- qc_report$n_locations > max_locations

  # Count gene targets per probe
  message("  Checking gene mappings...")
  if (nrow(annotations) > 0) {
    gene_counts <- annotations %>%
      dplyr::group_by(probe_id) %>%
      dplyr::summarise(
        n_genes = dplyr::n_distinct(ensembl_gene_id),
        gene_names = paste(unique(gene_name), collapse = ";"),
        gene_biotypes = paste(unique(gene_biotype), collapse = ";"),.groups = "drop"
      ) %>%
      as.data.frame()

    qc_report <- dplyr::left_join(qc_report, gene_counts, by = "probe_id")
    qc_report$n_genes[is.na(qc_report$n_genes)] <- 0
  } else {
    qc_report$n_genes <- 0
    qc_report$gene_names <- NA
    qc_report$gene_biotypes <- NA
  }

  # Flag cross-hybridization (multiple gene targets)
  qc_report$cross_hybridization <- qc_report$n_genes > 1

  # Flag no annotation
  qc_report$no_annotation <- qc_report$n_genes == 0

  # Flag intergenic probes (has genomic location but no gene)
  if (flag_intergenic) {
    qc_report$intergenic <- qc_report$n_locations > 0 & qc_report$n_genes == 0
  }

  # Check for non-protein-coding targets
  if ("gene_biotypes" %in% colnames(qc_report)) {
    qc_report$non_protein_coding <- !is.na(qc_report$gene_biotypes) &
      !grepl("protein_coding", qc_report$gene_biotypes)
  }

  # Compare with original annotations if provided
  if (!is.null(original_annotations)) {
    message("  Comparing with original annotations...")
    comparison <- compare_with_original(qc_report, original_annotations)
    qc_report <- dplyr::left_join(qc_report, comparison, by = "probe_id")
  }

  # Assign overall QC status
  message("  Assigning QC status...")
  qc_report$qc_status <- apply(qc_report, 1, function(row) {
    if (row["n_locations"] == 0) {
      return("DEPRECATED")
    } else if (row["multiple_locations"] && row["cross_hybridization"]) {
      return("AMBIGUOUS")
    } else if (row["cross_hybridization"]) {
      return("CROSS_HYBRIDIZATION")
    } else if (row["multiple_locations"]) {
      return("MULTI_MAPPING")
    } else if (row["no_annotation"]) {
      return("INTERGENIC")
    } else {
      return("PASS")
    }
  })

  # Generate summary
  message("\n=== Probe Quality Control Summary ===")
  message("Total probes evaluated: ", nrow(qc_report))
  message("PASS: ", sum(qc_report$qc_status == "PASS"),
          " (", round(100 * sum(qc_report$qc_status == "PASS") / nrow(qc_report), 1), "%)")
  message("DEPRECATED (no mapping): ", sum(qc_report$qc_status == "DEPRECATED"))
  message("MULTI_MAPPING: ", sum(qc_report$qc_status == "MULTI_MAPPING"))
  message("CROSS_HYBRIDIZATION: ", sum(qc_report$qc_status == "CROSS_HYBRIDIZATION"))
  message("AMBIGUOUS: ", sum(qc_report$qc_status == "AMBIGUOUS"))
  message("INTERGENIC: ", sum(qc_report$qc_status == "INTERGENIC"))

  class(qc_report) <- c("prism_qc", "data.frame")

  return(qc_report)
}

#' Compare with Original Annotations (Internal Helper)
#'
#' @param qc_report Current QC report
#' @param original_annotations Original probe annotations
#' @return Data frame with comparison flags
#' @keywords internal
compare_with_original <- function(qc_report, original_annotations) {

  # Extract probe IDs and gene symbols from original
  if ("ID" %in% colnames(original_annotations)) {
    orig <- data.frame(
      probe_id = original_annotations$ID,
      original_gene = original_annotations[["Gene Symbol"]] %||%
        original_annotations[["GENE_SYMBOL"]] %||%
        original_annotations[["gene_symbol"]],
      stringsAsFactors = FALSE
    )
  } else {
    warning("Cannot find probe ID column in original annotations")
    return(data.frame(probe_id = qc_report$probe_id))
  }

  # Clean up
  orig$original_gene <- trimws(as.character(orig$original_gene))
  orig$original_gene[orig$original_gene %in% c("", "---", "NA")] <- NA

  # Flag annotation changes
  comparison <- qc_report %>%
    dplyr::select(probe_id, gene_names) %>%
    dplyr::left_join(orig, by = "probe_id") %>%
    dplyr::mutate(
      annotation_changed = !is.na(original_gene) &
        !is.na(gene_names) &
        original_gene != gene_names,
      annotation_lost = !is.na(original_gene) & is.na(gene_names),
      annotation_gained = is.na(original_gene) & !is.na(gene_names)
    ) %>%
    dplyr::select(probe_id, original_gene, annotation_changed,
                  annotation_lost, annotation_gained)

  return(comparison)
}

# Helper for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Filter Probes Based on QC Results
#'
#' Removes problematic probes based on QC flags
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param remove_deprecated Remove probes with no genomic mapping (default: TRUE)
#' @param remove_multi_mapping Remove multi-mapping probes (default: TRUE)
#' @param remove_cross_hyb Remove cross-hybridizing probes (default: TRUE)
#' @param remove_intergenic Remove intergenic probes (default: FALSE)
#' @param keep_status Vector of QC status values to keep (overrides other params)
#' @return Character vector of probe IDs that pass QC
#' @export
#' @examples
#' \dontrun{
#'   # Keep only perfectly mapping probes
#'   clean_probes <- filter_probes_by_qc(
#'     qc_results,
#'     remove_deprecated = TRUE,
#'     remove_multi_mapping = TRUE,
#'     remove_cross_hyb = TRUE
#'   )
#'
#'   # Custom filtering
#'   clean_probes <- filter_probes_by_qc(
#'     qc_results,
#'     keep_status = c("PASS", "INTERGENIC")
#'   )
#' }
filter_probes_by_qc <- function(qc_results,
                                remove_deprecated = TRUE,
                                remove_multi_mapping = TRUE,
                                remove_cross_hyb = TRUE,
                                remove_intergenic = FALSE,
                                keep_status = NULL) {

  if (!is.null(keep_status)) {
    # Use explicit status list
    passing_probes <- qc_results$probe_id[qc_results$qc_status %in% keep_status]
  } else {
    # Build filter based on parameters
    keep <- rep(TRUE, nrow(qc_results))

    if (remove_deprecated) {
      keep <- keep & qc_results$qc_status != "DEPRECATED"
    }

    if (remove_multi_mapping) {
      keep <- keep & qc_results$qc_status != "MULTI_MAPPING"
    }

    if (remove_cross_hyb) {
      keep <- keep & !qc_results$qc_status %in% c("CROSS_HYBRIDIZATION", "AMBIGUOUS")
    }

    if (remove_intergenic) {
      keep <- keep & qc_results$qc_status != "INTERGENIC"
    }

    passing_probes <- qc_results$probe_id[keep]
  }

  message("QC filtering: ", nrow(qc_results), " -> ", length(passing_probes),
          " probes (", round(100 * length(passing_probes) / nrow(qc_results), 1), "% retained)")

  return(passing_probes)
}

#' Identify Specific Probe Issues
#'
#' Detailed analysis of probe problems
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param issue_type Type of issue: "deprecated", "multi_mapping", "cross_hyb", "changed"
#' @return Data frame with probes having the specified issue
#' @export
#' @examples
#' \dontrun{
#'   # Get all deprecated probes
#'   deprecated <- identify_probe_issues(qc_results, "deprecated")
#'
#'   # Get probes with changed annotations
#'   changed <- identify_probe_issues(qc_results, "changed")
#' }
identify_probe_issues <- function(qc_results, issue_type = "deprecated") {

  issue_type <- match.arg(issue_type,
                          c("deprecated", "multi_mapping", "cross_hyb",
                            "changed", "intergenic", "all"))

  if (issue_type == "deprecated") {
    issues <- qc_results[qc_results$qc_status == "DEPRECATED", ]
  } else if (issue_type == "multi_mapping") {
    issues <- qc_results[qc_results$multiple_locations == TRUE, ]
  } else if (issue_type == "cross_hyb") {
    issues <- qc_results[qc_results$cross_hybridization == TRUE, ]
  } else if (issue_type == "changed") {
    if ("annotation_changed" %in% colnames(qc_results)) {
      issues <- qc_results[qc_results$annotation_changed == TRUE, ]
    } else {
      warning("No annotation comparison available")
      return(data.frame())
    }
  } else if (issue_type == "intergenic") {
    issues <- qc_results[qc_results$qc_status == "INTERGENIC", ]
  } else if (issue_type == "all") {
    issues <- qc_results[qc_results$qc_status != "PASS", ]
  }

  message("Found ", nrow(issues), " probes with issue: ", issue_type)

  return(issues)
}

#' Generate Comprehensive QC Report
#'
#' Creates a detailed quality control report with statistics and plots
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param output_file Optional file path to save report (default: NULL)
#' @return List with QC statistics
#' @export
#' @examples
#' \dontrun{
#'   report <- generate_qc_report(qc_results)
#'
#'   # Save to file
#'   report <- generate_qc_report(qc_results, "qc_report.txt")
#' }
generate_qc_report <- function(qc_results, output_file = NULL) {

  # Calculate statistics
  total_probes <- nrow(qc_results)

  status_counts <- table(qc_results$qc_status)
  status_pct <- round(100 * status_counts / total_probes, 2)

  # Location statistics
  location_stats <- list(
    mean_locations = mean(qc_results$n_locations, na.rm = TRUE),
    median_locations = median(qc_results$n_locations, na.rm = TRUE),
    max_locations = max(qc_results$n_locations, na.rm = TRUE),
    no_location = sum(qc_results$n_locations == 0),
    unique_mapping = sum(qc_results$n_locations == 1),
    multi_mapping = sum(qc_results$n_locations > 1)
  )

  # Gene statistics
  gene_stats <- list(
    mean_genes = mean(qc_results$n_genes, na.rm = TRUE),
    median_genes = median(qc_results$n_genes, na.rm = TRUE),
    max_genes = max(qc_results$n_genes, na.rm = TRUE),
    no_genes = sum(qc_results$n_genes == 0),
    single_gene = sum(qc_results$n_genes == 1),
    multi_gene = sum(qc_results$n_genes > 1)
  )

  # Annotation change statistics (if available)
  if ("annotation_changed" %in% colnames(qc_results)) {
    change_stats <- list(
      unchanged = sum(!qc_results$annotation_changed, na.rm = TRUE),
      changed = sum(qc_results$annotation_changed, na.rm = TRUE),
      lost = sum(qc_results$annotation_lost, na.rm = TRUE),
      gained = sum(qc_results$annotation_gained, na.rm = TRUE)
    )
  } else {
    change_stats <- NULL
  }

  # Compile report
  report <- list(
    total_probes = total_probes,
    status_summary = data.frame(
      Status = names(status_counts),
      Count = as.numeric(status_counts),
      Percentage = as.numeric(status_pct)
    ),
    location_stats = location_stats,
    gene_stats = gene_stats,
    change_stats = change_stats,
    timestamp = Sys.time()
  )

  class(report) <- c("prism_qc_report", "list")

  # Save to file if requested
  if (!is.null(output_file)) {
    message("Saving QC report to: ", output_file)
    sink(output_file)
    print(report)
    sink()
  }

  return(report)
}

#' @export
print.prism_qc_report <- function(x,...) {
  cat("=======================================\n")
  cat("PRISM Quality Control Report\n")
  cat("=======================================\n")
  cat("Generated:", format(x$timestamp), "\n\n")

  cat("Total Probes:", x$total_probes, "\n\n")

  cat("QC Status Summary:\n")
  cat("------------------\n")
  print(x$status_summary, row.names = FALSE)
  cat("\n")

  cat("Genomic Location Statistics:\n")
  cat("-----------------------------\n")
  cat("  Mean locations per probe:", round(x$location_stats$mean_locations, 2), "\n")
  cat("  Median locations per probe:", x$location_stats$median_locations, "\n")
  cat("  Max locations per probe:", x$location_stats$max_locations, "\n")
  cat("  No genomic location:", x$location_stats$no_location,
      sprintf("(%.1f%%)\n", 100 * x$location_stats$no_location / x$total_probes))
  cat("  Unique mapping:", x$location_stats$unique_mapping,
      sprintf("(%.1f%%)\n", 100 * x$location_stats$unique_mapping / x$total_probes))
  cat("  Multi-mapping:", x$location_stats$multi_mapping,
      sprintf("(%.1f%%)\n", 100 * x$location_stats$multi_mapping / x$total_probes))
  cat("\n")

  cat("Gene Annotation Statistics:\n")
  cat("----------------------------\n")
  cat("  Mean genes per probe:", round(x$gene_stats$mean_genes, 2), "\n")
  cat("  Median genes per probe:", x$gene_stats$median_genes, "\n")
  cat("  Max genes per probe:", x$gene_stats$max_genes, "\n")
  cat("  No gene annotation:", x$gene_stats$no_genes,
      sprintf("(%.1f%%)\n", 100 * x$gene_stats$no_genes / x$total_probes))
  cat("  Single gene:", x$gene_stats$single_gene,
      sprintf("(%.1f%%)\n", 100 * x$gene_stats$single_gene / x$total_probes))
  cat("  Multiple genes:", x$gene_stats$multi_gene,
      sprintf("(%.1f%%)\n", 100 * x$gene_stats$multi_gene / x$total_probes))
  cat("\n")

  if (!is.null(x$change_stats)) {
    cat("Annotation Changes (vs Original):\n")
    cat("----------------------------------\n")
    cat("  Unchanged:", x$change_stats$unchanged,
        sprintf("(%.1f%%)\n", 100 * x$change_stats$unchanged / x$total_probes))
    cat("  Changed:", x$change_stats$changed,
        sprintf("(%.1f%%)\n", 100 * x$change_stats$changed / x$total_probes))
    cat("  Lost:", x$change_stats$lost,
        sprintf("(%.1f%%)\n", 100 * x$change_stats$lost / x$total_probes))
    cat("  Gained:", x$change_stats$gained,
        sprintf("(%.1f%%)\n", 100 * x$change_stats$gained / x$total_probes))
    cat("\n")
  }

  cat("=======================================\n")

  invisible(x)
}

#' Export QC Results to File
#'
#' Saves QC results in various formats
#'
#' @param qc_results Data frame from flag_problematic_probes
#' @param file Output file path
#' @param format File format: "csv", "tsv", "excel" (default: "csv")
#' @param include_pass Include probes that passed QC (default: TRUE)
#' @export
#' @examples
#' \dontrun{
#'   # Export all results to CSV
#'   export_qc_results(qc_results, "qc_results.csv")
#'
#'   # Export only problematic probes
#'   export_qc_results(qc_results, "problems.csv", include_pass = FALSE)
#' }
export_qc_results <- function(qc_results,
                              file,
                              format = "csv",
                              include_pass = TRUE) {

  format <- match.arg(format, c("csv", "tsv", "excel"))

  # Filter if needed
  if (!include_pass) {
    qc_results <- qc_results[qc_results$qc_status != "PASS", ]
    message("Exporting ", nrow(qc_results), " problematic probes")
  } else {
    message("Exporting ", nrow(qc_results), " probes")
  }

  # Export based on format
  if (format == "csv") {
    write.csv(qc_results, file, row.names = FALSE)
  } else if (format == "tsv") {
    write.table(qc_results, file, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (format == "excel") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' is required for Excel export. Install with: install.packages('writexl')")
    }
    writexl::write_xlsx(qc_results, file)
  }

  message("QC results saved to: ", file)

  invisible(qc_results)
}

#' Validate Probe Sequences
#'
#' Checks probe sequences for common issues
#'
#' @param probe_sequences DNAStringSet of probe sequences
#' @param min_length Minimum acceptable length (default: 20)
#' @param max_length Maximum acceptable length (default: 70)
#' @param check_complexity Check for low complexity sequences (default: TRUE)
#' @return Data frame with sequence validation results
#' @export
#' @importFrom Biostrings letterFrequency width
#' @examples
#' \dontrun{
#'   seq_qc <- validate_probe_sequences(probe_sequences)
#' }
validate_probe_sequences <- function(probe_sequences,
                                     min_length = 20,
                                     max_length = 70,
                                     check_complexity = TRUE) {

  message("Validating ", length(probe_sequences), " probe sequences...")

  # Get sequence lengths
  seq_lengths <- Biostrings::width(probe_sequences)

  # Check for invalid characters
  valid_chars <- Biostrings::letterFrequency(probe_sequences,
                                             letters = c("A", "C", "G", "T"))
  total_chars <- rowSums(valid_chars)
  has_invalid <- total_chars != seq_lengths

  # Calculate GC content
  gc_counts <- Biostrings::letterFrequency(probe_sequences,
                                           letters = c("G", "C"))
  gc_content <- rowSums(gc_counts) / seq_lengths

  # Check for low complexity (homopolymers)
  if (check_complexity) {
    has_homopolymer <- sapply(as.character(probe_sequences), function(seq) {
      grepl("AAAAA|TTTTT|GGGGG|CCCCC", seq)
    })
  } else {
    has_homopolymer <- rep(FALSE, length(probe_sequences))
  }

  # Compile results
  validation <- data.frame(
    probe_id = names(probe_sequences),
    length = seq_lengths,
    gc_content = round(gc_content, 3),
    too_short = seq_lengths < min_length,
    too_long = seq_lengths > max_length,
    invalid_chars = has_invalid,
    low_complexity = has_homopolymer,
    stringsAsFactors = FALSE
  )

  # Overall pass/fail
  validation$sequence_valid <- !(validation$too_short |
                                   validation$too_long |
                                   validation$invalid_chars |
                                   validation$low_complexity)

  # Summary
  message("\nSequence Validation Summary:")
  message("  Valid sequences: ", sum(validation$sequence_valid),
          " (", round(100 * sum(validation$sequence_valid) / nrow(validation), 1), "%)")
  message("  Too short: ", sum(validation$too_short))
  message("  Too long: ", sum(validation$too_long))
  message("  Invalid characters: ", sum(validation$invalid_chars))
  message("  Low complexity: ", sum(validation$low_complexity))

  return(validation)
}

#' Print QC Object
#'
#' @export
print.prism_qc <- function(x,...) {
  cat("PRISM QC Results\n")
  cat("================\n")
  cat("Total probes:", nrow(x), "\n\n")

  cat("QC Status:\n")
  print(table(x$qc_status))
  cat("\n")

  if ("annotation_changed" %in% colnames(x)) {
    cat("Annotation Changes:\n")
    cat("  Changed:", sum(x$annotation_changed, na.rm = TRUE), "\n")
    cat("  Lost:", sum(x$annotation_lost, na.rm = TRUE), "\n")
    cat("  Gained:", sum(x$annotation_gained, na.rm = TRUE), "\n")
  }

  invisible(x)
}
