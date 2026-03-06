#' Align Probe Sequences to Reference Genome
#'
#' Maps probe sequences to genome using exact and fuzzy matching. Identifies
#' probe locations on both forward and reverse strands.
#'
#' @param probe_sequences DNAStringSet of probe sequences
#' @param genome BSgenome object or DNAStringSet of chromosomes
#' @param max_mismatches Maximum number of mismatches allowed (default: 2)
#' @param min_length Minimum probe length to align (default: 20)
#' @param report_progress Show progress messages (default: TRUE)
#' @return GRanges object with probe alignments
#' @export
#' @importFrom Biostrings matchPDict PDict reverseComplement width
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @examples
#' \dontrun{
#'   library(BSgenome.Hsapiens.UCSC.hg38)
#'
#'   # Get probe sequences
#'   sequences <- get_geo_probe_sequences("GPL570")
#'
#'   # Align to genome
#'   alignments <- align_probes_to_genome(
#'     sequences,
#'     BSgenome.Hsapiens.UCSC.hg38,
#'     max_mismatches = 1
#'   )
#' }
align_probes_to_genome <- function(probe_sequences,
                                   genome,
                                   max_mismatches = 2,
                                   min_length = 20,
                                   report_progress = TRUE) {

  # Validate inputs
  if (!inherits(probe_sequences, "DNAStringSet")) {
    stop("probe_sequences must be a DNAStringSet object")
  }

  if (!inherits(genome, "BSgenome") && !inherits(genome, "DNAStringSet")) {
    stop("genome must be a BSgenome or DNAStringSet object")
  }

  # Filter short probes
  probe_lengths <- Biostrings::width(probe_sequences)
  valid_probes <- probe_lengths >= min_length

  if (sum(valid_probes) == 0) {
    stop("No probes meet minimum length requirement (", min_length, " bp)")
  }

  if (sum(!valid_probes) > 0 && report_progress) {
    message("Excluding ", sum(!valid_probes),
            " probes shorter than ", min_length, " bp")
    probe_sequences <- probe_sequences[valid_probes]
  }

  if (report_progress) {
    message("Aligning ", length(probe_sequences), " probes to genome")
    message("Parameters: max_mismatches = ", max_mismatches,
            ", min_length = ", min_length)
  }

  # Create pattern dictionary for efficient matching
  pdict <- Biostrings::PDict(probe_sequences, max.mismatch = max_mismatches)

  # Get chromosome names
  if (inherits(genome, "BSgenome")) {
    chr_names <- GenomeInfoDb::seqnames(genome)
  } else {
    chr_names <- names(genome)
  }

  # Align to each chromosome
  alignments_list <- list()

  for (i in seq_along(chr_names)) {
    chr <- chr_names[i]

    if (report_progress) {
      message("Processing ", chr, " (", i, "/", length(chr_names), ")")
    }

    # Get chromosome sequence
    chr_seq <- genome[[chr]]

    # Find matches on forward strand
    matches_fwd <- Biostrings::matchPDict(pdict, chr_seq,
                                          max.mismatch = max_mismatches)

    # Find matches on reverse strand
    chr_seq_rev <- Biostrings::reverseComplement(chr_seq)
    matches_rev <- Biostrings::matchPDict(pdict, chr_seq_rev,
                                          max.mismatch = max_mismatches)

    # Process forward strand matches
    for (j in seq_along(matches_fwd)) {
      if (length(matches_fwd[[j]]) > 0) {
        n_matches <- length(matches_fwd[[j]])
        alignments_list[[length(alignments_list) + 1]] <- data.frame(
          probe_id = rep(names(probe_sequences)[j], n_matches),
          seqnames = rep(as.character(chr), n_matches),
          start = start(matches_fwd[[j]]),
          end = end(matches_fwd[[j]]),
          strand = rep("+", n_matches),
          probe_length = rep(probe_lengths[valid_probes][j], n_matches),
          stringsAsFactors = FALSE
        )
      }
    }

    # Process reverse strand matches
    for (j in seq_along(matches_rev)) {
      if (length(matches_rev[[j]]) > 0) {
        n_matches <- length(matches_rev[[j]])
        # Convert coordinates back to forward strand
        chr_length <- length(chr_seq)
        fwd_starts <- chr_length - end(matches_rev[[j]]) + 1
        fwd_ends <- chr_length - start(matches_rev[[j]]) + 1

        alignments_list[[length(alignments_list) + 1]] <- data.frame(
          probe_id = rep(names(probe_sequences)[j], n_matches),
          seqnames = rep(as.character(chr), n_matches),
          start = fwd_starts,
          end = fwd_ends,
          strand = rep("-", n_matches),
          probe_length = rep(probe_lengths[valid_probes][j], n_matches),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Combine results
  if (length(alignments_list) == 0) {
    warning("No alignments found!")
    return(GenomicRanges::GRanges())
  }

  alignments_df <- do.call(rbind, alignments_list)

  # Calculate number of mismatches for each alignment
  alignments_df$max_mismatches <- max_mismatches

  # Convert to GRanges
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    alignments_df,
    keep.extra.columns = TRUE
  )

  if (report_progress) {
    n_probes_aligned <- length(unique(alignments_df$probe_id))
    n_total_alignments <- length(gr)
    message("\nAlignment Summary:")
    message("  Probes with alignments: ", n_probes_aligned, " / ",
            length(probe_sequences),
            " (", round(100 * n_probes_aligned / length(probe_sequences), 1), "%)")
    message("  Total alignments found: ", n_total_alignments)
    message("  Avg alignments per probe: ",
            round(n_total_alignments / n_probes_aligned, 2))
  }

  return(gr)
}

#' Count Exact Mismatches Between Probe and Genome
#'
#' Calculates the exact number of mismatches for aligned probes
#'
#' @param probe_sequences DNAStringSet of probe sequences
#' @param genome BSgenome object
#' @param alignments GRanges object from align_probes_to_genome
#' @param max_alignments Maximum alignments to check per probe (default: 100)
#' @return GRanges with added 'n_mismatches' column
#' @export
#' @importFrom Biostrings getSeq reverseComplement
count_alignment_mismatches <- function(probe_sequences,
                                       genome,
                                       alignments,
                                       max_alignments = 100) {

  if (length(alignments) == 0) {
    return(alignments)
  }

  message("Counting exact mismatches for ", length(alignments), " alignments...")

  # Limit number of alignments to check
  if (length(alignments) > max_alignments) {
    message("Sampling ", max_alignments, " alignments for mismatch counting")
    alignments <- alignments[sample(length(alignments), max_alignments)]
  }

  # Extract genomic sequences at alignment positions
  genomic_seqs <- Biostrings::getSeq(genome, alignments)

  # Get corresponding probe sequences
  probe_ids <- alignments$probe_id
  probe_seqs <- probe_sequences[probe_ids]

  # Reverse complement for minus strand
  is_minus <- as.character(strand(alignments)) == "-"
  probe_seqs[is_minus] <- Biostrings::reverseComplement(probe_seqs[is_minus])

  # Count mismatches
  n_mismatches <- sapply(seq_along(probe_seqs), function(i) {
    sum(as.character(probe_seqs[[i]]) != as.character(genomic_seqs[[i]]))
  })

  # Add to GRanges
  alignments$n_mismatches <- n_mismatches

  message("Mismatch distribution:")
  print(table(n_mismatches))

  return(alignments)
}

#' Summarize Probe Alignment Results
#'
#' Generates summary statistics for probe alignments
#'
#' @param alignments GRanges object from align_probes_to_genome
#' @param total_probes Total number of probes attempted (optional)
#' @return List with summary statistics
#' @export
summarize_alignments <- function(alignments, total_probes = NULL) {

  if (length(alignments) == 0) {
    return(list(
      total_probes = total_probes,
      aligned_probes = 0,
      total_alignments = 0,
      unmapped_probes = total_probes,
      unique_mapping = 0,
      multi_mapping = 0,
      avg_alignments_per_probe = 0
    ))
  }

  # Count alignments per probe
  alignment_counts <- table(alignments$probe_id)

  n_aligned <- length(alignment_counts)
  n_total_alignments <- length(alignments)
  n_unique <- sum(alignment_counts == 1)
  n_multi <- sum(alignment_counts > 1)

  # Calculate unmapped if total provided
  n_unmapped <- if (!is.null(total_probes)) {
    total_probes - n_aligned
  } else {
    NA
  }

  summary_list <- list(
    total_probes = total_probes,
    aligned_probes = n_aligned,
    unmapped_probes = n_unmapped,
    total_alignments = n_total_alignments,
    unique_mapping = n_unique,
    multi_mapping = n_multi,
    avg_alignments_per_probe = round(n_total_alignments / n_aligned, 2),
    max_alignments_per_probe = max(alignment_counts),
    alignment_distribution = as.data.frame(table(alignment_counts))
  )

  class(summary_list) <- c("prism_alignment_summary", "list")

  return(summary_list)
}

#' @export
print.prism_alignment_summary <- function(x,...) {
  cat("PRISM Alignment Summary\n")
  cat("=======================\n\n")

  if (!is.null(x$total_probes)) {
    cat("Total probes:", x$total_probes, "\n")
    cat("Aligned probes:", x$aligned_probes,
        sprintf("(%.1f%%)\n", 100 * x$aligned_probes / x$total_probes))
    if (!is.na(x$unmapped_probes)) {
      cat("Unmapped probes:", x$unmapped_probes,
          sprintf("(%.1f%%)\n", 100 * x$unmapped_probes / x$total_probes))
    }
  } else {
    cat("Aligned probes:", x$aligned_probes, "\n")
  }

  cat("\nAlignment Details:\n")
  cat("  Total alignments:", x$total_alignments, "\n")
  cat("  Unique mapping:", x$unique_mapping,
      sprintf("(%.1f%%)\n", 100 * x$unique_mapping / x$aligned_probes))
  cat("  Multi-mapping:", x$multi_mapping,
      sprintf("(%.1f%%)\n", 100 * x$multi_mapping / x$aligned_probes))
  cat("  Avg per probe:", x$avg_alignments_per_probe, "\n")
  cat("  Max per probe:", x$max_alignments_per_probe, "\n")

  cat("\nAlignment Distribution:\n")
  print(x$alignment_distribution, row.names = FALSE)

  invisible(x)
}

#' Filter Alignments by Quality
#'
#' Filters probe alignments based on various quality criteria
#'
#' @param alignments GRanges object from align_probes_to_genome
#' @param max_locations Maximum genomic locations per probe (default: 10)
#' @param min_length Minimum probe length (default: 20)
#' @param max_mismatches Maximum mismatches allowed (default: 2)
#' @param keep_unique_only Keep only uniquely mapping probes (default: FALSE)
#' @return Filtered GRanges object
#' @export
filter_alignments <- function(alignments,
                              max_locations = 10,
                              min_length = 20,
                              max_mismatches = 2,
                              keep_unique_only = FALSE) {

  if (length(alignments) == 0) {
    return(alignments)
  }

  n_original <- length(alignments)

  # Filter by probe length
  if ("probe_length" %in% names(mcols(alignments))) {
    alignments <- alignments[alignments$probe_length >= min_length]
  }

  # Filter by mismatches if available
  if ("n_mismatches" %in% names(mcols(alignments))) {
    alignments <- alignments[alignments$n_mismatches <= max_mismatches]
  }

  # Filter by number of locations
  if (keep_unique_only) {
    max_locations <- 1
  }

  probe_counts <- table(alignments$probe_id)
  valid_probes <- names(probe_counts)[probe_counts <= max_locations]
  alignments <- alignments[alignments$probe_id %in% valid_probes]

  n_filtered <- length(alignments)
  message("Filtered alignments: ", n_original, " -> ", n_filtered,
          " (", round(100 * n_filtered / n_original, 1), "% retained)")

  return(alignments)
}

