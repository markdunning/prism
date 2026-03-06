#' Utility Functions for PRISM Package
#'
#' @name utils
#' @keywords internal
NULL

#' NULL Coalescing Operator
#'
#' Returns the right-hand side if left-hand side is NULL
#'
#' @param x First value
#' @param y Second value (returned if x is NULL)
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Not In Operator
#'
#' Negation of %in%
#'
#' @param x Vector of values
#' @param y Vector to check against
#' @return Logical vector
#' @keywords internal
`%nin%` <- function(x, y) {
  !(x %in% y)
}

#' Validate Expression Matrix
#'
#' Checks if expression matrix is properly formatted
#'
#' @param expr_matrix Expression matrix to validate
#' @param require_rownames Require row names (default: TRUE)
#' @param require_colnames Require column names (default: TRUE)
#' @return TRUE if valid, stops with error otherwise
#' @export
#' @examples
#' \dontrun{
#'   validate_expression_matrix(my_matrix)
#' }
validate_expression_matrix <- function(expr_matrix,
                                       require_rownames = TRUE,
                                       require_colnames = TRUE) {

  # Check if it's a matrix or data frame
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("Expression data must be a matrix or data frame")
  }

  # Check for numeric values
  if (!is.numeric(as.matrix(expr_matrix))) {
    stop("Expression matrix must contain numeric values")
  }

  # Check dimensions
  if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    stop("Expression matrix has zero rows or columns")
  }

  # Check row names
  if (require_rownames && is.null(rownames(expr_matrix))) {
    stop("Expression matrix must have row names (probe/gene IDs)")
  }

  # Check column names
  if (require_colnames && is.null(colnames(expr_matrix))) {
    stop("Expression matrix must have column names (sample IDs)")
  }

  # Check for all NA
  if (all(is.na(expr_matrix))) {
    stop("Expression matrix contains only NA values")
  }

  # Warnings for common issues
  na_proportion <- sum(is.na(expr_matrix)) / length(expr_matrix)
  if (na_proportion > 0.5) {
    warning("Expression matrix contains >50% missing values (",
            round(na_proportion * 100, 1), "%)")
  }

  if (any(duplicated(rownames(expr_matrix)))) {
    warning("Expression matrix contains duplicated row names")
  }

  if (any(duplicated(colnames(expr_matrix)))) {
    warning("Expression matrix contains duplicated column names")
  }

  return(TRUE)
}


#' Validate Probe Sequences
#'
#' Checks if probe sequences are valid DNA sequences
#'
#' @param sequences Character vector or DNAStringSet of sequences
#' @return TRUE if valid, stops with error otherwise
#' @export
validate_probe_sequences <- function(sequences) {

  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }

  if (!is.character(sequences)) {
    stop("Sequences must be character vector or DNAStringSet")
  }

  if (length(sequences) == 0) {
    stop("No sequences provided")
  }

  # Check for valid DNA characters
  invalid_chars <- grepl("[^ACGTNacgtn]", sequences)
  if (any(invalid_chars)) {
    n_invalid <- sum(invalid_chars)
    stop("Found ", n_invalid, " sequences with invalid characters (not ACGTN)")
  }

  # Check for empty sequences
  empty_seqs <- nchar(sequences) == 0
  if (any(empty_seqs)) {
    stop("Found ", sum(empty_seqs), " empty sequences")
  }

  return(TRUE)
}


#' Validate GRanges Object
#'
#' Checks if GRanges object is properly formatted for PRISM
#'
#' @param granges GRanges object to validate
#' @param require_probe_id Require probe_id metadata column (default: TRUE)
#' @return TRUE if valid, stops with error otherwise
#' @export
validate_granges <- function(granges, require_probe_id = TRUE) {

  if (!inherits(granges, "GRanges")) {
    stop("Object must be a GRanges object")
  }

  if (length(granges) == 0) {
    stop("GRanges object is empty")
  }

  if (require_probe_id && !"probe_id" %in% names(mcols(granges))) {
    stop("GRanges object must have 'probe_id' metadata column")
  }

  return(TRUE)
}

#' Convert Expression Set to Matrix
#'
#' Extracts expression matrix from ExpressionSet or SummarizedExperiment
#'
#' @param eset ExpressionSet or SummarizedExperiment object
#' @return Expression matrix
#' @export
#' @examples
#' \dontrun{
#'   expr_matrix <- eset_to_matrix(my_eset)
#' }
eset_to_matrix <- function(eset) {

  if (inherits(eset, "ExpressionSet")) {
    return(Biobase::exprs(eset))
  } else if (inherits(eset, "SummarizedExperiment")) {
    return(SummarizedExperiment::assay(eset))
  } else {
    stop("Object must be ExpressionSet or SummarizedExperiment")
  }
}


#' Convert Data Frame to Matrix
#'
#' Safely converts data frame to numeric matrix
#'
#' @param df Data frame to convert
#' @param row_names Column to use as row names (default: 1)
#' @return Numeric matrix
#' @export
df_to_matrix <- function(df, row_names = 1) {

  if (!is.data.frame(df)) {
    stop("Input must be a data frame")
  }

  # Extract row names if specified
  if (is.numeric(row_names)) {
    rnames <- df[, row_names]
    df <- df[, -row_names, drop = FALSE]
  } else if (is.character(row_names)) {
    rnames <- df[[row_names]]
    df[[row_names]] <- NULL
  } else {
    rnames <- rownames(df)
  }

  # Convert to matrix
  mat <- as.matrix(df)

  # Check if numeric
  if (!is.numeric(mat)) {
    stop("Data frame contains non-numeric columns")
  }

  rownames(mat) <- rnames

  return(mat)
}


#' Standardize Chromosome Names
#'
#' Converts chromosome names to standard format
#'
#' @param chr_names Vector of chromosome names
#' @param style Style to use: "UCSC" (chr1), "NCBI" (1), "Ensembl" (1)
#' @return Standardized chromosome names
#' @export
#' @examples
#' standardize_chr_names(c("1", "2", "X"), style = "UCSC")
standardize_chr_names <- function(chr_names, style = "UCSC") {

  style <- match.arg(style, c("UCSC", "NCBI", "Ensembl"))

  chr_names <- as.character(chr_names)

  if (style == "UCSC") {
    # Add "chr" prefix if not present
    chr_names <- ifelse(grepl("^chr", chr_names), chr_names, paste0("chr", chr_names))
  } else {
    # Remove "chr" prefix if present
    chr_names <- gsub("^chr", "", chr_names)
  }

  return(chr_names)
}

#' Calculate Robust Mean
#'
#' Calculates mean after removing outliers
#'
#' @param x Numeric vector
#' @param trim Proportion to trim from each end (default: 0.1)
#' @param na.rm Remove NA values (default: TRUE)
#' @return Trimmed mean
#' @export
robust_mean <- function(x, trim = 0.1, na.rm = TRUE) {
  mean(x, trim = trim, na.rm = na.rm)
}


#' Calculate Coefficient of Variation
#'
#' Calculates CV (sd/mean)
#'
#' @param x Numeric vector
#' @param na.rm Remove NA values (default: TRUE)
#' @return Coefficient of variation
#' @export
cv <- function(x, na.rm = TRUE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}


#' Calculate Median Absolute Deviation
#'
#' Robust measure of variability
#'
#' @param x Numeric vector
#' @param na.rm Remove NA values (default: TRUE)
#' @return MAD value
#' @export
mad_score <- function(x, na.rm = TRUE) {
  median(abs(x - median(x, na.rm = na.rm)), na.rm = na.rm)
}


#' Normalize Expression Data
#'
#' Various normalization methods for expression data
#'
#' @param expr_matrix Expression matrix
#' @param method Normalization method: "quantile", "zscore", "log2", "scale"
#' @return Normalized expression matrix
#' @export
#' @examples
#' \dontrun{
#'   normalized <- normalize_expression(expr_matrix, method = "quantile")
#' }
normalize_expression <- function(expr_matrix, method = "quantile") {

  method <- match.arg(method, c("quantile", "zscore", "log2", "scale", "none"))

  if (method == "none") {
    return(expr_matrix)
  }

  if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package 'preprocessCore' required for quantile normalization.\n",
           "Install from Bioconductor: BiocManager::install('preprocessCore')")
    }
    return(preprocessCore::normalize.quantiles(expr_matrix))

  } else if (method == "zscore") {
    # Z-score normalization (by gene/probe)
    return(t(scale(t(expr_matrix))))

  } else if (method == "log2") {
    # Log2 transformation
    return(log2(expr_matrix + 1))

  } else if (method == "scale") {
    # Scale to 0-1 range
    return((expr_matrix - min(expr_matrix, na.rm = TRUE)) /
             (max(expr_matrix, na.rm = TRUE) - min(expr_matrix, na.rm = TRUE)))
  }
}

#' Smart Read Table
#'
#' Automatically detects file format and reads appropriately
#'
#' @param file Path to file
#' @param... Additional arguments passed to read function
#' @return Data frame
#' @export
smart_read <- function(file,...) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))

  data <- switch(ext,
                 "csv" = read.csv(file,...),
                 "tsv" = read.delim(file,...),
                 "txt" = read.delim(file,...),
                 "rds" = readRDS(file),
                 "rdata" = {
                   env <- new.env()
                   load(file, envir = env)
                   as.list(env)
                 },
                 stop("Unsupported file format: ", ext)
  )

  return(data)
}


#' Smart Write Table
#'
#' Automatically writes in appropriate format based on extension
#'
#' @param data Data to write
#' @param file Output file path
#' @param... Additional arguments
#' @export
smart_write <- function(data, file,...) {

  ext <- tolower(tools::file_ext(file))

  switch(ext,
         "csv" = write.csv(data, file, row.names = FALSE,...),
         "tsv" = write.table(data, file, sep = "\t", row.names = FALSE, quote = FALSE,...),
         "txt" = write.table(data, file, sep = "\t", row.names = FALSE, quote = FALSE,...),
         "rds" = saveRDS(data, file,...),
         stop("Unsupported file format: ", ext)
  )

  message("Data written to: ", file)
  invisible(file)
}


#' Check File Size
#'
#' Returns human-readable file size
#'
#' @param file Path to file
#' @return Character string with file size
#' @export
get_file_size <- function(file) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  size_bytes <- file.info(file)$size

  if (size_bytes < 1024) {
    return(paste(size_bytes, "bytes"))
  } else if (size_bytes < 1024^2) {
    return(paste(round(size_bytes / 1024, 2), "KB"))
  } else if (size_bytes < 1024^3) {
    return(paste(round(size_bytes / 1024^2, 2), "MB"))
  } else {
    return(paste(round(size_bytes / 1024^3, 2), "GB"))
  }
}

#' Progress Bar Wrapper
#'
#' Creates a simple progress bar
#'
#' @param n Total number of iterations
#' @param title Progress bar title
#' @return Progress bar object
#' @keywords internal
create_progress_bar <- function(n, title = "Progress") {

  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = paste0(title, " [:bar] :percent eta: :eta"),
      total = n,
      clear = FALSE,
      width = 60
    )
    return(pb)
  } else {
    # Fallback: simple counter
    list(
      tick = function() message("."),
      update = function(ratio) NULL
    )
  }
}


#' Log Message with Timestamp
#'
#' Prints message with timestamp
#'
#' @param... Message components
#' @param level Log level: "INFO", "WARNING", "ERROR"
#' @export
log_message <- function(..., level = "INFO") {

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", timestamp, "] [", level, "] ",...)

  if (level == "ERROR") {
    stop(msg, call. = FALSE)
  } else if (level == "WARNING") {
    warning(msg, call. = FALSE)
  } else {
    message(msg)
  }
}


#' Create Log File
#'
#' Initializes a log file for tracking analysis
#'
#' @param log_file Path to log file
#' @param append Append to existing file (default: FALSE)
#' @return Log file connection
#' @export
init_log_file <- function(log_file, append = FALSE) {

  if (!append && file.exists(log_file)) {
    warning("Log file already exists and will be overwritten")
  }

  # Write header
  if (!append) {
    cat("PRISM Analysis Log\n", file = log_file)
    cat("==================\n", file = log_file, append = TRUE)
    cat("Started:", format(Sys.time()), "\n\n", file = log_file, append = TRUE)
  }

  message("Log file initialized: ", log_file)
  return(log_file)
}


#' Write to Log File
#'
#' Appends message to log file
#'
#' @param log_file Path to log file
#' @param... Message components
#' @export
write_log <- function(log_file,...) {

  if (!file.exists(log_file)) {
    stop("Log file does not exist: ", log_file)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", timestamp, "] ",..., "\n")

  cat(msg, file = log_file, append = TRUE)
}
#' Check Available Memory
#'
#' Returns available system memory
#'
#' @return Character string with available memory
#' @export
check_memory <- function() {

  if (.Platform$OS.type == "unix") {
    # Linux/Mac
    mem_info <- system("free -h", intern = TRUE)
    return(mem_info[2])
  } else {
    # Windows
    mem_info <- system("wmic OS get FreePhysicalMemory", intern = TRUE)
    free_kb <- as.numeric(gsub("[^0-9]", "", mem_info[2]))
    free_gb <- round(free_kb / 1024^2, 2)
    return(paste(free_gb, "GB available"))
  }
}


#' Estimate Memory Usage
#'
#' Estimates memory required for expression matrix
#'
#' @param n_probes Number of probes/genes
#' @param n_samples Number of samples
#' @param bytes_per_value Bytes per numeric value (default: 8)
#' @return Character string with estimated memory
#' @export
#' @examples
#' estimate_memory(50000, 100)
estimate_memory <- function(n_probes, n_samples, bytes_per_value = 8) {

  total_bytes <- n_probes * n_samples * bytes_per_value

  if (total_bytes < 1024^2) {
    return(paste(round(total_bytes / 1024, 2), "KB"))
  } else if (total_bytes < 1024^3) {
    return(paste(round(total_bytes / 1024^2, 2), "MB"))
  } else {
    return(paste(round(total_bytes / 1024^3, 2), "GB"))
  }
}


#' Time Function Execution
#'
#' Measures execution time of a function
#'
#' @param expr Expression to time
#' @param label Label for the timing (optional)
#' @return Result of expression
#' @export
#' @examples
#' \dontrun{
#'   result <- time_execution(my_function(data), "My Function")
#' }
time_execution <- function(expr, label = NULL) {

  start_time <- Sys.time()
  result <- expr
  end_time <- Sys.time()

  elapsed <- difftime(end_time, start_time, units = "auto")

  if (!is.null(label)) {
    message(label, " completed in ", round(elapsed, 2), " ", attr(elapsed, "units"))
  } else {
    message("Execution time: ", round(elapsed, 2), " ", attr(elapsed, "units"))
  }

  return(result)
}


#' Batch Process with Memory Management
#'
#' Processes data in batches to manage memory
#'
#' @param data Data to process
#' @param batch_size Size of each batch
#' @param fun Function to apply to each batch
#' @param... Additional arguments to fun
#' @return Combined results
#' @export
batch_process <- function(data, batch_size, fun,...) {

  n_total <- length(data)
  n_batches <- ceiling(n_total / batch_size)

  message("Processing ", n_total, " items in ", n_batches, " batches...")

  results <- list()

  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_total)

    batch_data <- data[start_idx:end_idx]

    message("  Batch ", i, "/", n_batches, " (items ", start_idx, "-", end_idx, ")")

    results[[i]] <- fun(batch_data,...)

    # Garbage collection after each batch
    gc(verbose = FALSE)
  }

  return(results)
}

#' Clean Gene Names
#'
#' Standardizes gene names by removing special characters
#'
#' @param gene_names Vector of gene names
#' @param to_upper Convert to uppercase (default: TRUE)
#' @return Cleaned gene names
#' @export
#' @examples
#' clean_gene_names(c("tp53", "BRCA-1", "egfr "))
clean_gene_names <- function(gene_names, to_upper = TRUE) {

  # Remove leading/trailing whitespace
  gene_names <- trimws(gene_names)

  # Remove special characters (keep alphanumeric and hyphens)
  gene_names <- gsub("[^A-Za-z0-9-]", "", gene_names)

  # Convert to uppercase if requested
  if (to_upper) {
    gene_names <- toupper(gene_names)
  }

  return(gene_names)
}


#' Truncate String
#'
#' Truncates long strings with ellipsis
#'
#' @param x Character vector
#' @param max_length Maximum length (default: 50)
#' @return Truncated strings
#' @export
truncate_string <- function(x, max_length = 50) {

  ifelse(nchar(x) > max_length,
         paste0(substr(x, 1, max_length - 3), "..."),
         x)
}


#' Format Number with Commas
#'
#' Adds thousand separators to numbers
#'
#' @param x Numeric vector
#' @param digits Number of decimal places (default: 0)
#' @return Formatted character vector
#' @export
#' @examples
#' format_number(1234567.89, digits = 2)
format_number <- function(x, digits = 0) {
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}


#' Create Safe File Name
#'
#' Converts string to safe file name
#'
#' @param x Character string
#' @param replacement Replacement for invalid characters (default: "_")
#' @return Safe file name
#' @export
#' @examples
#' safe_filename("My Analysis: Results (2024)")
safe_filename <- function(x, replacement = "_") {

  # Remove or replace invalid characters
  x <- gsub("[^A-Za-z0-9._-]", replacement, x)

  # Remove multiple consecutive replacements
  x <- gsub(paste0(replacement, "+"), replacement, x)

  # Remove leading/trailing replacements
  x <- gsub(paste0("^", replacement, "|", replacement, "$"), "", x)

  return(x)
}

#' Quick Summary of Expression Matrix
#'
#' Provides quick overview of expression data
#'
#' @param expr_matrix Expression matrix
#' @return List with summary statistics
#' @export
#' @examples
#' \dontrun{
#'   summary_expr(my_matrix)
#' }
summary_expr <- function(expr_matrix) {

  summary_list <- list(
    dimensions = dim(expr_matrix),
    n_features = nrow(expr_matrix),
    n_samples = ncol(expr_matrix),
    total_values = length(expr_matrix),
    n_missing = sum(is.na(expr_matrix)),
    pct_missing = round(100 * sum(is.na(expr_matrix)) / length(expr_matrix), 2),
    range = range(expr_matrix, na.rm = TRUE),
    mean = mean(expr_matrix, na.rm = TRUE),
    median = median(expr_matrix, na.rm = TRUE),
    sd = sd(as.vector(expr_matrix), na.rm = TRUE)
  )

  class(summary_list) <- c("prism_expr_summary", "list")

  return(summary_list)
}


#' @export
print.prism_expr_summary <- function(x,...) {
  cat("Expression Matrix Summary\n")
  cat("=========================\n")
  cat("Dimensions:", x$n_features, "×", x$n_samples, "\n")
  cat("Total values:", format_number(x$total_values), "\n")
  cat("Missing values:", format_number(x$n_missing),
      "(", x$pct_missing, "%)\n")
  cat("\nValue Statistics:\n")
  cat("  Range:", round(x$range[1], 3), "to", round(x$range[2], 3), "\n")
  cat("  Mean:", round(x$mean, 3), "\n")
  cat("  Median:", round(x$median, 3), "\n")
  cat("  SD:", round(x$sd, 3), "\n")

  invisible(x)
}


#' Count Unique Values
#'
#' Counts unique values in each column
#'
#' @param data Data frame or matrix
#' @return Named vector of unique counts
#' @export
count_unique <- function(data) {

  sapply(as.data.frame(data), function(x) length(unique(x)))
}


#' Detect Outliers
#'
#' Identifies outliers using IQR method
#'
#' @param x Numeric vector
#' @param k IQR multiplier (default: 1.5)
#' @return Logical vector indicating outliers
#' @export
detect_outliers <- function(x, k = 1.5) {

  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1

  lower_bound <- q1 - k * iqr
  upper_bound <- q3 + k * iqr

  return(x < lower_bound | x > upper_bound)
}

#' Find Common Elements
#'
#' Finds elements common to all vectors in a list
#'
#' @param... Vectors or list of vectors
#' @return Vector of common elements
#' @export
#' @examples
#' find_common(c("A", "B", "C"), c("B", "C", "D"), c("C", "D", "E"))
find_common <- function(...) {

  vectors <- list(...)

  if (length(vectors) == 1 && is.list(vectors[[1]])) {
    vectors <- vectors[[1]]
  }

  Reduce(intersect, vectors)
}


#' Match with Tolerance
#'
#' Matches numeric values with tolerance
#'
#' @param x Numeric vector to match
#' @param table Numeric vector to match against
#' @param tolerance Matching tolerance (default: 0.01)
#' @return Integer vector of match positions
#' @export
match_tolerance <- function(x, table, tolerance = 0.01) {

  sapply(x, function(val) {
    diffs <- abs(table - val)
    if (min(diffs) <= tolerance) {
      which.min(diffs)
    } else {
      NA_integer_
    }
  })
}


#' Fuzzy String Match
#'
#' Finds approximate string matches
#'
#' @param x Character vector to match
#' @param table Character vector to match against
#' @param max_distance Maximum edit distance (default: 2)
#' @return Data frame with matches
#' @export
fuzzy_match <- function(x, table, max_distance = 2) {

  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("Package 'stringdist' required. Install with: install.packages('stringdist')")
  }

  matches <- lapply(x, function(query) {
    distances <- stringdist::stringdist(query, table, method = "lv")
    best_match_idx <- which.min(distances)
    best_distance <- distances[best_match_idx]

    if (best_distance <= max_distance) {
      data.frame(
        query = query,
        match = table[best_match_idx],
        distance = best_distance,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        query = query,
        match = NA_character_,
        distance = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, matches)
}

#' Get PRISM Session Info
#'
#' Returns detailed session information for reproducibility
#'
#' @return List with session details
#' @export
get_prism_session_info <- function() {

  session_info <- list(
    prism_version = packageVersion("PRISM"),
    r_version = R.version.string,
    platform = R.version$platform,
    os = Sys.info()["sysname"],
    date = Sys.Date(),
    time = Sys.time(),
    locale = Sys.getlocale(),
    packages = sessionInfo()$otherPkgs
  )

  class(session_info) <- c("prism_session_info", "list")

  return(session_info)
}


#' @export
print.prism_session_info <- function(x, ...) {
  cat("PRISM Session Information\n")
  cat("=========================\n\n")

  cat("PRISM version:", as.character(x$prism_version), "\n")
  cat("R version:", x$r_version, "\n")
  cat("Platform:", x$platform, "\n")
  cat("OS:", x$os, "\n")
  cat("Date:", format(x$date), "\n")
  cat("Time:", format(x$time), "\n\n")

  if (length(x$packages) > 0) {
    cat("Loaded packages:\n")
    pkg_info <- sapply(x$packages, function(p) {
      paste0("  ", p$Package, " (", p$Version, ")")
    })
    cat(paste(pkg_info, collapse = "\n"), "\n")
  }

  invisible(x)
}


#' Save Session Info to File
#'
#' Saves session information for reproducibility
#'
#' @param file Output file path
#' @export
save_session_info <- function(file = "PRISM_session_info.txt") {

  info <- get_prism_session_info()

  sink(file)
  print(info)
  cat("\n\nFull Session Info:\n")
  cat("==================\n")
  print(sessionInfo())
  sink()

  message("Session info saved to: ", file)
  invisible(file)
}

#' Check Package Installation
#'
#' Checks if required packages are installed
#'
#' @param packages Character vector of package names
#' @param bioc Logical indicating if packages are from Bioconductor
#' @return Logical vector indicating which packages are installed
#' @export
check_packages <- function(packages, bioc = FALSE) {

  installed <- sapply(packages, requireNamespace, quietly = TRUE)

  if (any(!installed)) {
    missing <- packages[!installed]

    if (bioc) {
      message("Missing Bioconductor packages: ", paste(missing, collapse = ", "))
      message("Install with: BiocManager::install(c('",
              paste(missing, collapse = "', '"), "'))")
    } else {
      message("Missing CRAN packages: ", paste(missing, collapse = ", "))
      message("Install with: install.packages(c('",
              paste(missing, collapse = "', '"), "'))")
    }
  }

  return(installed)
}


#' Get PRISM Example Data
#'
#' Returns path to example data files
#'
#' @param dataset Name of example dataset
#' @return Path to example data
#' @export
#' @examples
#' \dontrun{
#'   example_file <- get_example_data("sample_probes")
#' }
get_example_data <- function(dataset = NULL) {

  data_dir <- system.file("extdata", package = "PRISM")

  if (is.null(dataset)) {
    # List available datasets
    files <- list.files(data_dir)
    message("Available example datasets:")
    message(paste("  -", files, collapse = "\n"))
    return(invisible(files))
  }

  file_path <- file.path(data_dir, dataset)

  if (!file.exists(file_path)) {
    stop("Example dataset not found: ", dataset)
  }

  return(file_path)
}


#' Create PRISM Directory Structure
#'
#' Sets up organized directory structure for analysis
#'
#' @param base_dir Base directory path (default: current directory)
#' @return List of created directories
#' @export
create_prism_dirs <- function(base_dir = ".") {

  dirs <- c(
    "data/raw",
    "data/processed",
    "results/alignments",
    "results/annotations",
    "results/qc",
    "results/expression",
    "plots",
    "logs"
  )

  full_paths <- file.path(base_dir, dirs)

  for (dir_path in full_paths) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Created: ", dir_path)
    }
  }

  message("\nPRISM directory structure created in: ", base_dir)

  return(invisible(full_paths))
}

