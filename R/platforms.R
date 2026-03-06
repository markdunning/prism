#' Platform-Specific Functions for PRISM
#'
#' Handles different microarray platforms (Affymetrix, Illumina, Agilent, etc.)
#'
#' @name platforms
#' @keywords internal
NULL


#' Detect Microarray Platform
#'
#' Automatically detects the microarray platform from probe IDs or platform annotation
#'
#' @param probe_ids Vector of probe IDs
#' @param platform_id Optional platform ID (e.g., "GPL570")
#' @return Character string with detected platform
#' @export
#' @examples
#' detect_platform(c("1007_s_at", "1053_at", "117_at"))
#' detect_platform(platform_id = "GPL570")
detect_platform <- function(probe_ids = NULL, platform_id = NULL) {

  # Check platform ID first
  if (!is.null(platform_id)) {
    platform <- detect_from_gpl(platform_id)
    if (!is.null(platform)) {
      message("Detected platform from GPL ID: ", platform)
      return(platform)
    }
  }

  # Check probe ID patterns
  if (!is.null(probe_ids)) {
    platform <- detect_from_probe_ids(probe_ids)
    if (!is.null(platform)) {
      message("Detected platform from probe IDs: ", platform)
      return(platform)
    }
  }

  warning("Could not detect platform automatically")
  return("unknown")
}


#' Detect Platform from GPL ID
#'
#' @param gpl_id GEO platform ID
#' @return Platform name or NULL
#' @keywords internal
detect_from_gpl <- function(gpl_id) {

  # Common Affymetrix platforms
  affymetrix_platforms <- c(
    "GPL96" = "Affymetrix HG-U133A",
    "GPL97" = "Affymetrix HG-U133B",
    "GPL570" = "Affymetrix HG-U133 Plus 2.0",
    "GPL571" = "Affymetrix HG-U133A 2.0",
    "GPL1261" = "Affymetrix Mouse Genome 430 2.0",
    "GPL6244" = "Affymetrix Human Gene 1.0 ST",
    "GPL6246" = "Affymetrix Mouse Gene 1.0 ST",
    "GPL16686" = "Affymetrix Human Gene 2.0 ST",
    "GPL17586" = "Affymetrix Human Transcriptome Array 2.0"
  )

  # Common Illumina platforms
  illumina_platforms <- c(
    "GPL6883" = "Illumina HumanRef-8 v3",
    "GPL6884" = "Illumina HumanWG-6 v3",
    "GPL6947" = "Illumina HumanHT-12 V3",
    "GPL10558" = "Illumina HumanHT-12 V4",
    "GPL6885" = "Illumina MouseRef-8 v2",
    "GPL6887" = "Illumina MouseWG-6 v2"
  )

  # Common Agilent platforms
  agilent_platforms <- c(
    "GPL4133" = "Agilent Human 1A",
    "GPL6480" = "Agilent Human 4x44K v2",
    "GPL13497" = "Agilent SurePrint G3 Human GE 8x60K",
    "GPL14550" = "Agilent Human 8x60K v2"
  )

  if (gpl_id %in% names(affymetrix_platforms)) {
    return(affymetrix_platforms[gpl_id])
  } else if (gpl_id %in% names(illumina_platforms)) {
    return(illumina_platforms[gpl_id])
  } else if (gpl_id %in% names(agilent_platforms)) {
    return(agilent_platforms[gpl_id])
  }

  return(NULL)
}


#' Detect Platform from Probe IDs
#'
#' @param probe_ids Vector of probe IDs
#' @return Platform name or NULL
#' @keywords internal
detect_from_probe_ids <- function(probe_ids) {

  # Sample probe IDs if too many
  if (length(probe_ids) > 100) {
    probe_ids <- sample(probe_ids, 100)
  }

  # Affymetrix patterns
  if (any(grepl("_at$|_s_at$|_x_at$", probe_ids))) {
    return("Affymetrix")
  }

  # Illumina patterns (typically numeric with optional prefix)
  if (any(grepl("^ILMN_\\d+$", probe_ids))) {
    return("Illumina")
  }

  # Agilent patterns (typically alphanumeric)
  if (any(grepl("^A_\\d+_P\\d+$", probe_ids))) {
    return("Agilent")
  }

  # CodeLink patterns
  if (any(grepl("^GE\\d+$", probe_ids))) {
    return("CodeLink")
  }

  return(NULL)
}

#' Parse Affymetrix Probe IDs
#'
#' Extracts information from Affymetrix probe IDs
#'
#' @param probe_ids Vector of Affymetrix probe IDs
#' @return Data frame with parsed information
#' @export
#' @examples
#' parse_affymetrix_ids(c("1007_s_at", "1053_at", "117_at"))
parse_affymetrix_ids <- function(probe_ids) {

  # Extract probe set ID and suffix
  parsed <- data.frame(
    probe_id = probe_ids,
    probeset = gsub("_.*$", "", probe_ids),
    suffix = gsub("^.*_", "", probe_ids),
    stringsAsFactors = FALSE
  )

  # Classify probe type based on suffix
  parsed$probe_type <- dplyr::case_when(
    parsed$suffix == "at" ~ "antisense_transcript",
    parsed$suffix == "s_at" ~ "sense_transcript",
    parsed$suffix == "x_at" ~ "cross_hybridizing",
    parsed$suffix == "a_at" ~ "alternative_transcript",
    parsed$suffix == "st" ~ "sense_target",
    TRUE ~ "unknown"
  )

  return(parsed)
}


#' Get Affymetrix Annotation Package
#'
#' Returns the appropriate annotation package for Affymetrix platform
#'
#' @param platform_id GPL platform ID
#' @return Annotation package name
#' @export
get_affymetrix_annotation <- function(platform_id) {

  annotation_packages <- c(
    "GPL96" = "hgu133a.db",
    "GPL97" = "hgu133b.db",
    "GPL570" = "hgu133plus2.db",
    "GPL571" = "hgu133a2.db",
    "GPL1261" = "mouse4302.db",
    "GPL6244" = "hugene10sttranscriptcluster.db",
    "GPL6246" = "mogene10sttranscriptcluster.db"
  )

  if (platform_id %in% names(annotation_packages)) {
    return(annotation_packages[platform_id])
  } else {
    warning("No annotation package found for platform: ", platform_id)
    return(NULL)
  }
}


#' Load Affymetrix Probe Sequences
#'
#' Loads probe sequences for Affymetrix platforms
#'
#' @param platform_id GPL platform ID
#' @return DNAStringSet with probe sequences
#' @export
#' @examples
#' \dontrun{
#'   sequences <- load_affymetrix_sequences("GPL570")
#' }
load_affymetrix_sequences <- function(platform_id) {

  message("Loading Affymetrix probe sequences for ", platform_id)

  # Try to get from GEO
  sequences <- tryCatch({
    get_geo_probe_sequences(platform_id)
  }, error = function(e) {
    warning("Could not load sequences from GEO: ", e$message)
    return(NULL)
  })

  if (is.null(sequences)) {
    stop("Failed to load probe sequences for platform: ", platform_id)
  }

  return(sequences)
}

#' Parse Illumina Probe IDs
#'
#' Extracts information from Illumina probe IDs
#'
#' @param probe_ids Vector of Illumina probe IDs
#' @return Data frame with parsed information
#' @export
#' @examples
#' parse_illumina_ids(c("ILMN_1234567", "ILMN_2345678"))
parse_illumina_ids <- function(probe_ids) {

  parsed <- data.frame(
    probe_id = probe_ids,
    numeric_id = gsub("^ILMN_", "", probe_ids),
    platform = "Illumina",
    stringsAsFactors = FALSE
  )

  return(parsed)
}


#' Get Illumina Annotation Package
#'
#' Returns the appropriate annotation package for Illumina platform
#'
#' @param platform_id GPL platform ID
#' @return Annotation package name
#' @export
get_illumina_annotation <- function(platform_id) {

  annotation_packages <- c(
    "GPL6883" = "illuminaHumanv3.db",
    "GPL6884" = "illuminaHumanv3.db",
    "GPL6947" = "illuminaHumanv3.db",
    "GPL10558" = "illuminaHumanv4.db",
    "GPL6885" = "illuminaMousev2.db",
    "GPL6887" = "illuminaMousev2.db"
  )

  if (platform_id %in% names(annotation_packages)) {
    return(annotation_packages[platform_id])
  } else {
    warning("No annotation package found for platform: ", platform_id)
    return(NULL)
  }
}


#' Load Illumina Probe Sequences
#'
#' Loads probe sequences for Illumina platforms
#'
#' @param platform_id GPL platform ID
#' @return DNAStringSet with probe sequences
#' @export
load_illumina_sequences <- function(platform_id) {

  message("Loading Illumina probe sequences for ", platform_id)

  # Try to get from GEO
  sequences <- tryCatch({
    get_geo_probe_sequences(platform_id)
  }, error = function(e) {
    warning("Could not load sequences from GEO: ", e$message)
    return(NULL)
  })

  if (is.null(sequences)) {
    stop("Failed to load probe sequences for platform: ", platform_id)
  }

  return(sequences)
}


#' Handle Illumina Control Probes
#'
#' Identifies and filters Illumina control probes
#'
#' @param probe_ids Vector of probe IDs
#' @return Logical vector indicating control probes
#' @export
is_illumina_control <- function(probe_ids) {

  # Illumina control probes typically start with specific prefixes
  control_patterns <- c(
    "^NEGATIVE",
    "^ERCC",
    "^Housekeeping",
    "^Biotin",
    "^Labeling",
    "^Hybridization",
    "^Low_Stringency_Hyb"
  )

  is_control <- rep(FALSE, length(probe_ids))

  for (pattern in control_patterns) {
    is_control <- is_control | grepl(pattern, probe_ids, ignore.case = TRUE)
  }

  return(is_control)
}

#' Parse Agilent Probe IDs
#'
#' Extracts information from Agilent probe IDs
#'
#' @param probe_ids Vector of Agilent probe IDs
#' @return Data frame with parsed information
#' @export
#' @examples
#' parse_agilent_ids(c("A_23_P123456", "A_24_P234567"))
parse_agilent_ids <- function(probe_ids) {

  parsed <- data.frame(
    probe_id = probe_ids,
    platform = "Agilent",
    stringsAsFactors = FALSE
  )

  # Extract array design and probe number if possible
  if (any(grepl("^A_\\d+_P\\d+$", probe_ids))) {
    parsed$array_design <- gsub("^A_(\\d+)_P\\d+$", "\\1", probe_ids)
    parsed$probe_number <- gsub("^A_\\d+_P(\\d+)$", "\\1", probe_ids)
  }

  return(parsed)
}


#' Load Agilent Probe Sequences
#'
#' Loads probe sequences for Agilent platforms
#'
#' @param platform_id GPL platform ID
#' @return DNAStringSet with probe sequences
#' @export
load_agilent_sequences <- function(platform_id) {

  message("Loading Agilent probe sequences for ", platform_id)

  # Try to get from GEO
  sequences <- tryCatch({
    get_geo_probe_sequences(platform_id)
  }, error = function(e) {
    warning("Could not load sequences from GEO: ", e$message)
    return(NULL)
  })

  if (is.null(sequences)) {
    stop("Failed to load probe sequences for platform: ", platform_id)
  }

  return(sequences)
}


#' Handle Agilent Control Probes
#'
#' Identifies and filters Agilent control probes
#'
#' @param probe_ids Vector of probe IDs
#' @param annotations Optional probe annotation data frame
#' @return Logical vector indicating control probes
#' @export
is_agilent_control <- function(probe_ids, annotations = NULL) {

  # Agilent control probes
  control_patterns <- c(
    "^DarkCorner",
    "^BrightCorner",
    "^GE_BrightCorner",
    "^NegativeControl",
    "^Positive",
    "^ERCC"
  )

  is_control <- rep(FALSE, length(probe_ids))

  for (pattern in control_patterns) {
    is_control <- is_control | grepl(pattern, probe_ids, ignore.case = TRUE)
  }

  # Check annotations if provided
  if (!is.null(annotations) && "ControlType" %in% colnames(annotations)) {
    is_control <- is_control | (!is.na(annotations$ControlType) &
                                  annotations$ControlType != "0")
  }

  return(is_control)
}

#' Load Platform Probe Sequences
#'
#' Generic function to load probe sequences for any platform
#'
#' @param platform_id GPL platform ID or platform name
#' @param platform_type Optional: "Affymetrix", "Illumina", "Agilent"
#' @return DNAStringSet with probe sequences
#' @export
#' @examples
#' \dontrun{
#'   sequences <- load_platform_sequences("GPL570")
#' }
load_platform_sequences <- function(platform_id, platform_type = NULL) {

  # Detect platform type if not provided
  if (is.null(platform_type)) {
    platform_type <- detect_platform(platform_id = platform_id)
  }

  # Load sequences based on platform
  if (grepl("Affymetrix", platform_type, ignore.case = TRUE)) {
    return(load_affymetrix_sequences(platform_id))
  } else if (grepl("Illumina", platform_type, ignore.case = TRUE)) {
    return(load_illumina_sequences(platform_id))
  } else if (grepl("Agilent", platform_type, ignore.case = TRUE)) {
    return(load_agilent_sequences(platform_id))
  } else {
    # Try generic GEO approach
    message("Using generic sequence loading for platform: ", platform_type)
    return(get_geo_probe_sequences(platform_id))
  }
}


#' Get Platform Annotation Package
#'
#' Returns appropriate annotation package for any platform
#'
#' @param platform_id GPL platform ID
#' @return Annotation package name
#' @export
get_platform_annotation <- function(platform_id) {

  platform_type <- detect_platform(platform_id = platform_id)

  if (grepl("Affymetrix", platform_type, ignore.case = TRUE)) {
    return(get_affymetrix_annotation(platform_id))
  } else if (grepl("Illumina", platform_type, ignore.case = TRUE)) {
    return(get_illumina_annotation(platform_id))
  } else {
    warning("No annotation package available for platform: ", platform_id)
    return(NULL)
  }
}


#' Filter Control Probes
#'
#' Removes control probes from any platform
#'
#' @param probe_ids Vector of probe IDs
#' @param platform_type Platform type
#' @param annotations Optional probe annotations
#' @return Logical vector indicating non-control probes
#' @export
filter_control_probes <- function(probe_ids, platform_type = NULL, annotations = NULL) {

  if (is.null(platform_type)) {
    platform_type <- detect_platform(probe_ids = probe_ids)
  }

  if (grepl("Illumina", platform_type, ignore.case = TRUE)) {
    is_control <- is_illumina_control(probe_ids)
  } else if (grepl("Agilent", platform_type, ignore.case = TRUE)) {
    is_control <- is_agilent_control(probe_ids, annotations)
  } else {
    # Generic control detection
    is_control <- grepl("control|negative|positive|ercc|spike",
                        probe_ids, ignore.case = TRUE)
  }

  n_control <- sum(is_control)
  if (n_control > 0) {
    message("Identified ", n_control, " control probes")
  }

  return(!is_control)
}

#' Compare Platform Coverage
#'
#' Compares gene coverage across different platforms
#'
#' @param platform_ids Vector of GPL platform IDs
#' @param species Species name (default: "hsapiens")
#' @return Data frame with coverage comparison
#' @export
#' @examples
#' \dontrun{
#'   comparison <- compare_platform_coverage(c("GPL570", "GPL6244"))
#' }
compare_platform_coverage <- function(platform_ids, species = "hsapiens") {

  message("Comparing coverage for ", length(platform_ids), " platforms...")

  coverage_list <- lapply(platform_ids, function(platform_id) {

    message("  Processing ", platform_id)

    # Get platform info
    platform_name <- detect_platform(platform_id = platform_id)

    # Try to get probe count
    tryCatch({
      gpl <- GEOquery::getGEO(platform_id)
      probe_table <- GEOquery::Table(gpl)

      data.frame(
        platform_id = platform_id,
        platform_name = platform_name,
        n_probes = nrow(probe_table),
        n_genes = length(unique(probe_table$GENE_SYMBOL)),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        platform_id = platform_id,
        platform_name = platform_name,
        n_probes = NA,
        n_genes = NA,
        stringsAsFactors = FALSE
      )
    })
  })

  coverage_df <- do.call(rbind, coverage_list)

  return(coverage_df)
}


#' Get Platform Specifications
#'
#' Returns detailed specifications for a platform
#'
#' @param platform_id GPL platform ID
#' @return List with platform specifications
#' @export
get_platform_specs <- function(platform_id) {

  message("Retrieving specifications for ", platform_id)

  # Get platform data from GEO
  gpl <- tryCatch({
    GEOquery::getGEO(platform_id)
  }, error = function(e) {
    stop("Could not retrieve platform information: ", e$message)
  })

  # Extract metadata
  specs <- list(
    platform_id = platform_id,
    title = GEOquery::Meta(gpl)$title,
    organism = GEOquery::Meta(gpl)$organism,
    manufacturer = GEOquery::Meta(gpl)$manufacturer,
    technology = GEOquery::Meta(gpl)$technology,
    distribution = GEOquery::Meta(gpl)$distribution,
    n_probes = nrow(GEOquery::Table(gpl)),
    submission_date = GEOquery::Meta(gpl)$submission_date,
    last_update = GEOquery::Meta(gpl)$last_update_date
  )

  class(specs) <- c("prism_platform_specs", "list")

  return(specs)
}


#' @export
print.prism_platform_specs <- function(x,...) {
  cat("Platform Specifications\n")
  cat("=======================\n\n")
  cat("Platform ID:", x$platform_id, "\n")
  cat("Title:", x$title, "\n")
  cat("Organism:", x$organism, "\n")
  cat("Manufacturer:", x$manufacturer, "\n")
  cat("Technology:", x$technology, "\n")
  cat("Number of probes:", format_number(x$n_probes), "\n")
  cat("Submission date:", x$submission_date, "\n")
  cat("Last update:", x$last_update, "\n")

  invisible(x)
}

#' Preprocess Platform Data
#'
#' Applies platform-specific preprocessing steps
#'
#' @param expression_matrix Expression matrix
#' @param platform_type Platform type
#' @param normalize Apply normalization (default: TRUE)
#' @param log_transform Apply log transformation (default: TRUE)
#' @param remove_controls Remove control probes (default: TRUE)
#' @return Preprocessed expression matrix
#' @export
preprocess_platform_data <- function(expression_matrix,
                                     platform_type,
                                     normalize = TRUE,
                                     log_transform = TRUE,
                                     remove_controls = TRUE) {

  message("Preprocessing ", platform_type, " data...")

  # Remove control probes
  if (remove_controls) {
    probe_ids <- rownames(expression_matrix)
    keep_probes <- filter_control_probes(probe_ids, platform_type)
    expression_matrix <- expression_matrix[keep_probes, , drop = FALSE]
    message("  Retained ", nrow(expression_matrix), " non-control probes")
  }

  # Log transformation
  if (log_transform) {
    # Check if already log-transformed
    max_val <- max(expression_matrix, na.rm = TRUE)

    if (max_val > 100) {
      message("  Applying log2 transformation...")
      expression_matrix <- log2(expression_matrix + 1)
    } else {
      message("  Data appears to be already log-transformed")
    }
  }

  # Normalization
  if (normalize) {
    message("  Applying quantile normalization...")

    if (grepl("Affymetrix", platform_type, ignore.case = TRUE)) {
      # RMA-style normalization for Affymetrix
      expression_matrix <- normalize_expression(expression_matrix, method = "quantile")

    } else if (grepl("Illumina", platform_type, ignore.case = TRUE)) {
      # Quantile normalization for Illumina
      expression_matrix <- normalize_expression(expression_matrix, method = "quantile")

    } else if (grepl("Agilent", platform_type, ignore.case = TRUE)) {
      # Quantile normalization for Agilent
      expression_matrix <- normalize_expression(expression_matrix, method = "quantile")
    }
  }

  message("Preprocessing complete!")

  return(expression_matrix)
}

#' List Supported Platforms
#'
#' Returns list of platforms with built-in support
#'
#' @return Data frame with supported platforms
#' @export
list_supported_platforms <- function() {

  platforms <- data.frame(
    GPL_ID = c(
      "GPL96", "GPL97", "GPL570", "GPL571", "GPL1261",
      "GPL6244", "GPL6246", "GPL16686", "GPL17586",
      "GPL6883", "GPL6884", "GPL6947", "GPL10558",
      "GPL6885", "GPL6887",
      "GPL4133", "GPL6480", "GPL13497", "GPL14550"
    ),
    Platform_Name = c(
      "Affymetrix HG-U133A", "Affymetrix HG-U133B",
      "Affymetrix HG-U133 Plus 2.0", "Affymetrix HG-U133A 2.0",
      "Affymetrix Mouse Genome 430 2.0",
      "Affymetrix Human Gene 1.0 ST", "Affymetrix Mouse Gene 1.0 ST",
      "Affymetrix Human Gene 2.0 ST", "Affymetrix Human Transcriptome Array 2.0",
      "Illumina HumanRef-8 v3", "Illumina HumanWG-6 v3",
      "Illumina HumanHT-12 V3", "Illumina HumanHT-12 V4",
      "Illumina MouseRef-8 v2", "Illumina MouseWG-6 v2",
      "Agilent Human 1A", "Agilent Human 4x44K v2",
      "Agilent SurePrint G3 Human GE 8x60K", "Agilent Human 8x60K v2"
    ),
    Manufacturer = c(
      rep("Affymetrix", 9),
      rep("Illumina", 6),
      rep("Agilent", 4)
    ),
    Organism = c(
      "Human", "Human", "Human", "Human", "Mouse",
      "Human", "Mouse", "Human", "Human",
      "Human", "Human", "Human", "Human",
      "Mouse", "Mouse",
      "Human", "Human", "Human", "Human"
    ),
    stringsAsFactors = FALSE
  )

  return(platforms)
}

#' List Supported Platforms
#'
#' Returns list of platforms with built-in support
#'
#' @return Data frame with supported platforms
#' @export
list_supported_platforms <- function() {

  platforms <- data.frame(
    GPL_ID = c(
      "GPL96", "GPL97", "GPL570", "GPL571", "GPL1261",
      "GPL6244", "GPL6246", "GPL16686", "GPL17586",
      "GPL6883", "GPL6884", "GPL6947", "GPL10558",
      "GPL6885", "GPL6887",
      "GPL4133", "GPL6480", "GPL13497", "GPL14550"
    ),
    Platform_Name = c(
      "Affymetrix HG-U133A", "Affymetrix HG-U133B",
      "Affymetrix HG-U133 Plus 2.0", "Affymetrix HG-U133A 2.0",
      "Affymetrix Mouse Genome 430 2.0",
      "Affymetrix Human Gene 1.0 ST", "Affymetrix Mouse Gene 1.0 ST",
      "Affymetrix Human Gene 2.0 ST", "Affymetrix Human Transcriptome Array 2.0",
      "Illumina HumanRef-8 v3", "Illumina HumanWG-6 v3",
      "Illumina HumanHT-12 V3", "Illumina HumanHT-12 V4",
      "Illumina MouseRef-8 v2", "Illumina MouseWG-6 v2",
      "Agilent Human 1A", "Agilent Human 4x44K v2",
      "Agilent SurePrint G3 Human GE 8x60K", "Agilent Human 8x60K v2"
    ),
    Manufacturer = c(
      rep("Affymetrix", 9),
      rep("Illumina", 6),
      rep("Agilent", 4)
    ),
    Organism = c(
      "Human", "Human", "Human", "Human", "Mouse",
      "Human", "Mouse", "Human", "Human",
      "Human", "Human", "Human", "Human",
      "Mouse", "Mouse",
      "Human", "Human", "Human", "Human"
    ),
    stringsAsFactors = FALSE
  )

  return(platforms)
}


#' Check Platform Compatibility
#'
#' Checks if platform is compatible with PRISM workflow
#'
#' @param platform_id GPL platform ID
#' @return Logical indicating compatibility
#' @export
check_platform_compatibility <- function(platform_id) {

  supported <- list_supported_platforms()

  if (platform_id %in% supported$GPL_ID) {
    message("Platform ", platform_id, " is fully supported")
    return(TRUE)
  } else {
    message("Platform ", platform_id, " may work but is not officially supported")
    message("PRISM will attempt to process using generic methods")
    return(FALSE)
  }
}


#' Get Platform Genome Build
#'
#' Returns the recommended genome build for a platform
#'
#' @param platform_id GPL platform ID
#' @return Character string with genome build
#' @export
get_platform_genome <- function(platform_id) {

  # Get platform specs
  specs <- tryCatch({
    get_platform_specs(platform_id)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(specs)) {
    warning("Could not determine genome build")
    return(NULL)
  }

  organism <- specs$organism

  # Map organism to genome build
  genome_map <- list(
    "Homo sapiens" = "hg38",
    "Mus musculus" = "mm10",
    "Rattus norvegicus" = "rn6",
    "Danio rerio" = "danRer11",
    "Drosophila melanogaster" = "dm6",
    "Caenorhabditis elegans" = "ce11"
  )

  if (organism %in% names(genome_map)) {
    genome <- genome_map[[organism]]
    message("Recommended genome build for ", organism, ": ", genome)
    return(genome)
  } else {
    warning("Unknown organism: ", organism)
    return(NULL)
  }
}

#' Register Custom Platform
#'
#' Adds support for a custom microarray platform
#'
#' @param platform_id Platform identifier
#' @param platform_name Platform name
#' @param manufacturer Manufacturer name
#' @param organism Organism
#' @param probe_id_pattern Regex pattern for probe IDs
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#'   register_custom_platform(
#'     platform_id = "GPL_CUSTOM",
#'     platform_name = "Custom Array v1",
#'     manufacturer = "CustomCo",
#'     organism = "Homo sapiens",
#'     probe_id_pattern = "^CUST_\\d+$"
#'   )
#' }
register_custom_platform <- function(platform_id,
                                     platform_name,
                                     manufacturer,
                                     organism,
                                     probe_id_pattern = NULL) {

  # Store in package environment
  if (!exists("custom_platforms", envir =.GlobalEnv)) {
    assign("custom_platforms", list(), envir =.GlobalEnv)
  }

  custom_platforms <- get("custom_platforms", envir =.GlobalEnv)

  custom_platforms[[platform_id]] <- list(
    platform_id = platform_id,
    platform_name = platform_name,
    manufacturer = manufacturer,
    organism = organism,
    probe_id_pattern = probe_id_pattern
  )

  assign("custom_platforms", custom_platforms, envir =.GlobalEnv)

  message("Custom platform registered: ", platform_id)

  invisible(NULL)
}


#' Load Custom Platform Sequences
#'
#' Loads probe sequences from a custom file
#'
#' @param sequence_file Path to FASTA file with probe sequences
#' @param platform_id Platform identifier
#' @return DNAStringSet with probe sequences
#' @export
load_custom_sequences <- function(sequence_file, platform_id = NULL) {

  if (!file.exists(sequence_file)) {
    stop("Sequence file not found: ", sequence_file)
  }

  message("Loading custom probe sequences from: ", sequence_file)

  sequences <- Biostrings::readDNAStringSet(sequence_file)

  message("Loaded ", length(sequences), " probe sequences")

  if (!is.null(platform_id)) {
    attr(sequences, "platform_id") <- platform_id
  }

  return(sequences)
}


#' Create Platform Annotation from File
#'
#' Creates platform annotation from custom annotation file
#'
#' @param annotation_file Path to annotation file (CSV/TSV)
#' @param probe_col Column name for probe IDs
#' @param gene_col Column name for gene symbols
#' @param sequence_col Optional column name for probe sequences
#' @return Data frame with platform annotations
#' @export
create_custom_annotation <- function(annotation_file,
                                     probe_col = "Probe_ID",
                                     gene_col = "Gene_Symbol",
                                     sequence_col = NULL) {

  if (!file.exists(annotation_file)) {
    stop("Annotation file not found: ", annotation_file)
  }

  message("Loading custom platform annotation from: ", annotation_file)

  # Read file
  annotations <- smart_read(annotation_file)

  # Validate required columns
  if (!probe_col %in% colnames(annotations)) {
    stop("Probe ID column not found: ", probe_col)
  }

  if (!gene_col %in% colnames(annotations)) {
    stop("Gene symbol column not found: ", gene_col)
  }

  # Standardize column names
  annotations$probe_id <- annotations[[probe_col]]
  annotations$gene_symbol <- annotations[[gene_col]]

  if (!is.null(sequence_col) && sequence_col %in% colnames(annotations)) {
    annotations$sequence <- annotations[[sequence_col]]
  }

  message("Loaded annotations for ", nrow(annotations), " probes")

  return(annotations)
}

#' Convert Between Platform Formats
#'
#' Converts probe IDs between different platform formats
#'
#' @param probe_ids Vector of probe IDs
#' @param from_platform Source platform
#' @param to_platform Target platform
#' @return Vector of converted probe IDs
#' @export
convert_probe_ids <- function(probe_ids, from_platform, to_platform) {

  message("Converting probe IDs from ", from_platform, " to ", to_platform)

  # This is a placeholder - actual implementation would require
  # platform-specific mapping tables
  warning("Probe ID conversion not yet implemented")

  return(probe_ids)
}


#' Map Probes Across Platforms
#'
#' Finds equivalent probes across different platforms
#'
#' @param probe_ids Vector of probe IDs from source platform
#' @param source_platform Source platform ID
#' @param target_platform Target platform ID
#' @param method Mapping method: "sequence", "gene", "both"
#' @return Data frame with probe mappings
#' @export
map_probes_across_platforms <- function(probe_ids,
                                        source_platform,
                                        target_platform,
                                        method = "both") {

  method <- match.arg(method, c("sequence", "gene", "both"))

  message("Mapping probes from ", source_platform, " to ", target_platform)
  message("Method: ", method)

  # Load sequences for both platforms
  source_seqs <- load_platform_sequences(source_platform)
  target_seqs <- load_platform_sequences(target_platform)

  # Filter to requested probes
  source_seqs <- source_seqs[names(source_seqs) %in% probe_ids]

  if (method == "sequence" || method == "both") {
    # Map by sequence similarity
    message("Mapping by sequence similarity...")

    # This would require sequence alignment - placeholder
    warning("Sequence-based mapping not yet fully implemented")
  }

  if (method == "gene" || method == "both") {
    # Map by gene annotation
    message("Mapping by gene annotation...")

    # This would require gene mappings - placeholder
    warning("Gene-based mapping not yet fully implemented")
  }

  # Return empty data frame for now
  return(data.frame(
    source_probe = character(),
    target_probe = character(),
    mapping_method = character(),
    confidence = numeric()
  ))
}

#' Summarize Platform Information
#'
#' Creates comprehensive summary of platform characteristics
#'
#' @param platform_id GPL platform ID
#' @param include_sequences Include sequence statistics (default: FALSE)
#' @return List with platform summary
#' @export
summarize_platform <- function(platform_id, include_sequences = FALSE) {

  message("Summarizing platform: ", platform_id)

  # Get basic specs
  specs <- get_platform_specs(platform_id)

  # Get platform type
  platform_type <- detect_platform(platform_id = platform_id)

  # Get annotation package
  annot_pkg <- get_platform_annotation(platform_id)

  # Get genome build
  genome <- get_platform_genome(platform_id)

  summary_list <- list(
    platform_id = platform_id,
    platform_name = specs$title,
    platform_type = platform_type,
    organism = specs$organism,
    manufacturer = specs$manufacturer,
    n_probes = specs$n_probes,
    annotation_package = annot_pkg,
    genome_build = genome,
    submission_date = specs$submission_date,
    last_update = specs$last_update
  )

  # Add sequence statistics if requested
  if (include_sequences) {
    message("  Loading and analyzing sequences...")

    sequences <- tryCatch({
      load_platform_sequences(platform_id)
    }, error = function(e) {
      NULL
    })

    if (!is.null(sequences)) {
      summary_list$sequence_stats <- list(
        n_sequences = length(sequences),
        mean_length = mean(Biostrings::width(sequences)),
        median_length = median(Biostrings::width(sequences)),
        min_length = min(Biostrings::width(sequences)),
        max_length = max(Biostrings::width(sequences))
      )
    }
  }

  class(summary_list) <- c("prism_platform_summary", "list")

  return(summary_list)
}


#' @export
print.prism_platform_summary <- function(x,...) {
  cat("PRISM Platform Summary\n")
  cat("======================\n\n")

  cat("Platform ID:", x$platform_id, "\n")
  cat("Name:", x$platform_name, "\n")
  cat("Type:", x$platform_type, "\n")
  cat("Organism:", x$organism, "\n")
  cat("Manufacturer:", x$manufacturer, "\n")
  cat("Number of probes:", format_number(x$n_probes), "\n")

  if (!is.null(x$annotation_package)) {
    cat("Annotation package:", x$annotation_package, "\n")
  }

  if (!is.null(x$genome_build)) {
    cat("Recommended genome:", x$genome_build, "\n")
  }

  cat("\nDates:\n")
  cat("  Submitted:", x$submission_date, "\n")
  cat("  Last updated:", x$last_update, "\n")

  if (!is.null(x$sequence_stats)) {
    cat("\nSequence Statistics:\n")
    cat("  Available sequences:", format_number(x$sequence_stats$n_sequences), "\n")
    cat("  Mean length:", round(x$sequence_stats$mean_length, 1), "bp\n")
    cat("  Length range:", x$sequence_stats$min_length, "-",
        x$sequence_stats$max_length, "bp\n")
  }

  invisible(x)
}


#' Export Platform Information
#'
#' Exports platform information to file
#'
#' @param platform_id GPL platform ID
#' @param output_file Output file path
#' @param format Output format: "txt", "json", "yaml"
#' @export
export_platform_info <- function(platform_id,
                                 output_file,
                                 format = "txt") {

  format <- match.arg(format, c("txt", "json", "yaml"))

  summary <- summarize_platform(platform_id, include_sequences = TRUE)

  if (format == "txt") {
    sink(output_file)
    print(summary)
    sink()
  } else if (format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' required for JSON export")
    }
    jsonlite::write_json(summary, output_file, pretty = TRUE, auto_unbox = TRUE)
  } else if (format == "yaml") {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' required for YAML export")
    }
    yaml::write_yaml(summary, output_file)
  }

  message("Platform information exported to: ", output_file)

  invisible(output_file)
}
