#' PRISM Array Data Object
#'
#' S3 class for storing microarray data with annotations
#'
#' @param expression Numeric matrix of expression values
#' @param probes Data frame with probe information
#' @param samples Data frame with sample metadata
#' @param platform Character string identifying platform
#' @param source Character string indicating data source
#' @return A prism_array object
#' @export
new_prism_array <- function(expression, probes, samples,
                            platform = NULL, source = NULL) {

  stopifnot(is.matrix(expression) || is.data.frame(expression))
  stopifnot(is.data.frame(probes))
  stopifnot(is.data.frame(samples))

  structure(
    list(
      expression = as.matrix(expression),
      probes = probes,
      samples = samples,
      platform = platform,
      source = source,
      created = Sys.time()
    ),
    class = "prism_array"
  )
}

#' @export
print.prism_array <- function(x,...) {
  cat("PRISM Array Object\n")
  cat("==================\n")
  cat("Platform:", x$platform %||% "Unknown", "\n")
  cat("Source:", x$source %||% "Unknown", "\n")
  cat("Probes:", nrow(x$expression), "\n")
  cat("Samples:", ncol(x$expression), "\n")
  cat("Created:", format(x$created), "\n")
  invisible(x)
}

#' PRISM Annotation Results
#'
#' S3 class for storing re-annotation results
#'
#' @param alignments GRanges object with probe alignments
#' @param annotations Data frame with gene annotations
#' @param qc_flags Data frame with quality control flags
#' @param summary List with summary statistics
#' @return A prism_annotation object
#' @export
new_prism_annotation <- function(alignments, annotations,
                                 qc_flags, summary) {

  structure(
    list(
      alignments = alignments,
      annotations = annotations,
      qc_flags = qc_flags,
      summary = summary,
      timestamp = Sys.time()
    ),
    class = "prism_annotation"
  )
}

#' @export
print.prism_annotation <- function(x,...) {
  cat("PRISM Annotation Results\n")
  cat("========================\n")
  cat("Total probes:", x$summary$total_probes, "\n")
  cat("Aligned probes:", x$summary$aligned_probes, "\n")
  cat("Annotated genes:", x$summary$unique_genes, "\n")
  cat("QC pass rate:",
      sprintf("%.1f%%", x$summary$qc_pass_rate * 100), "\n")
  cat("Timestamp:", format(x$timestamp), "\n")
  invisible(x)
}

# Helper function
`%||%` <- function(x, y) if (is.null(x)) y else x
