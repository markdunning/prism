#' Import Microarray Data from GEO
#'
#' Downloads and processes microarray data from NCBI GEO database
#'
#' @param geo_id GEO accession number (e.g., "GSE12345")
#' @param platform Platform GPL ID (optional, auto-detected if NULL)
#' @param destdir Directory to store downloaded files
#' @return A prism_array object
#' @export
#' @importFrom GEOquery getGEO
#' @examples
#' \dontrun{
#'   geo_data <- import_geo_data("GSE12345")
#' }
import_geo_data <- function(geo_id, platform = NULL,
                            destdir = tempdir()) {

  # Validate GEO ID format
  if (!grepl("^GSE\\d+$", geo_id)) {
    stop("Invalid GEO ID format. Expected format: GSE12345")
  }

  message("Downloading GEO dataset: ", geo_id)

  # Download data
  gse <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE, destdir = destdir)

  # Handle multiple platforms
  if (length(gse) > 1 && is.null(platform)) {
    platforms <- sapply(gse, function(x) annotation(x))
    stop("Multiple platforms detected: ", paste(platforms, collapse = ", "),
         "\nPlease specify platform parameter.")
  }

  # Extract expression set
  eset <- if (length(gse) == 1) gse[[1]] else {
    idx <- which(sapply(gse, function(x) annotation(x) == platform))
    if (length(idx) == 0) {
      stop("Platform ", platform, " not found in dataset")
    }
    gse[[idx]]
  }

  # Extract components
  expr_matrix <- Biobase::exprs(eset)
  probe_info <- Biobase::fData(eset)
  sample_info <- Biobase::pData(eset)
  platform_id <- Biobase::annotation(eset)

  # Create PRISM object
  prism_obj <- new_prism_array(
    expression = expr_matrix,
    probes = probe_info,
    samples = sample_info,
    platform = platform_id,
    source = paste0("GEO:", geo_id)
  )

  message("Successfully imported ", nrow(expr_matrix), " probes and ",
          ncol(expr_matrix), " samples")

  return(prism_obj)
}

#' Extract Probe Sequences from GEO Platform
#'
#' Retrieves probe sequences from GEO platform annotation
#'
#' @param platform GPL platform ID (e.g., "GPL570")
#' @param destdir Directory to store downloaded files
#' @return DNAStringSet with probe sequences
#' @export
#' @importFrom Biostrings DNAStringSet
#' @examples
#' \dontrun{
#'   sequences <- get_geo_probe_sequences("GPL570")
#' }
get_geo_probe_sequences <- function(platform, destdir = tempdir()) {

  message("Downloading platform annotation: ", platform)

  gpl <- GEOquery::getGEO(platform, destdir = destdir)
  platform_table <- GEOquery::Table(gpl)

  # Look for sequence column (different names on different platforms)
  seq_cols <- c("SEQUENCE", "Probe_Sequence", "probe_sequence",
                "Target_Sequence", "Sequence")
  seq_col <- seq_cols[seq_cols %in% colnames(platform_table)]

  if (length(seq_col) == 0) {
    stop("No probe sequence column found in platform annotation.\n",
         "Available columns: ", paste(colnames(platform_table), collapse = ", "))
  }

  seq_col <- seq_col[1]  # Use first match

  # Extract sequences
  probe_ids <- platform_table$ID
  sequences <- platform_table[[seq_col]]

  # Remove empty sequences
  valid <- !is.na(sequences) & sequences != "" & sequences != "---"

  if (sum(valid) == 0) {
    stop("No valid probe sequences found")
  }

  message("Found ", sum(valid), " valid probe sequences")

  # Create DNAStringSet
  dna_seqs <- Biostrings::DNAStringSet(sequences[valid])
  names(dna_seqs) <- probe_ids[valid]

  return(dna_seqs)
}
