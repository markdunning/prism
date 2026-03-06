#' Import Expression Matrix from File
#'
#' Reads expression data from various file formats
#'
#' @param file Path to expression file (CSV, TSV, or RDS)
#' @param probe_col Column name or index for probe IDs
#' @param sep Field separator (auto-detected if NULL)
#' @return Matrix with probes as rows, samples as columns
#' @export
#' @examples
#' \dontrun{
#'   expr <- import_expression_file("expression.csv", probe_col = 1)
#' }
import_expression_file <- function(file, probe_col = 1, sep = NULL) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # Detect file type
  ext <- tools::file_ext(file)

  data <- switch(ext,
                 "rds" = readRDS(file),
                 "RDS" = readRDS(file),
                 "csv" = read.csv(file, row.names = probe_col, check.names = FALSE),
                 "tsv" = read.delim(file, row.names = probe_col, check.names = FALSE),
                 "txt" = {
                   if (is.null(sep)) sep <- "\t"
                   read.delim(file, sep = sep, row.names = probe_col, check.names = FALSE)
                 },
                 stop("Unsupported file format: ", ext)
  )

  # Convert to matrix if needed
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  message("Imported expression data: ", nrow(data), " probes × ",
          ncol(data), " samples")

  return(data)
}

#' Import Probe Sequences from FASTA
#'
#' Reads probe sequences from FASTA file
#'
#' @param file Path to FASTA file
#' @return DNAStringSet with probe sequences
#' @export
#' @importFrom Biostrings readDNAStringSet
#' @examples
#' \dontrun{
#'   sequences <- import_probe_fasta("probes.fasta")
#' }
import_probe_fasta <- function(file) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  sequences <- Biostrings::readDNAStringSet(file)

  message("Imported ", length(sequences), " probe sequences")

  return(sequences)
}
