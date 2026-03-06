#' Package Startup and Attachment Functions
#'
#' Functions that run when PRISM is loaded or attached
#'
#' @name zzz
#' @keywords internal
NULL


#' Package Load Hook
#'
#' Runs when the package namespace is loaded
#'
#' @param libname Library name
#' @param pkgname Package name
#' @keywords internal
#'

.onLoad <- function(libname, pkgname) {

# Set package options with defaults
op <- options()
op.prism <- list(
  prism.verbose = TRUE,
  prism.max_mismatches = 2,
  prism.min_probe_length = 20,
  prism.default_species = "hsapiens",
  prism.default_genome = "hg38",
  prism.batch_size = 500,
  prism.parallel_cores = 1,
  prism.cache_dir = tempdir()
)

# Only set options that are not already set
toset <- !(names(op.prism) %in% names(op))
if (any(toset)) options(op.prism[toset])

invisible()
}


#' Package Attach Hook
#'
#' Runs when the package is attached via library() or require()
#'
#' @param libname Library name
#' @param pkgname Package name
#' @keywords internal
#'
.onAttach <- function(libname, pkgname) {

# Get package version
version <- package_version("prism")

# Create startup message
msg <- paste0(
  "\n",
  "╔═══════════════════════════════════════════════════════╗\n",
  "║                     PRISM v", version, "                      ║\n",
  "║   Probe Re-Identification and Sequence Mapping       ║\n",
  "╚═══════════════════════════════════════════════════════╝\n",
  "\n",
  "Re-annotate legacy microarray datasets with current\n",
  "genome builds and gene annotations.\n",
  "\n",
  "Getting started:\n",
  "  • Tutorial: vignette('introduction', package = 'PRISM')\n",
  "  • Help: ?PRISM\n",
  "  • Report issues: https://github.com/yourusername/PRISM/issues\n",
  "\n"
)

packageStartupMessage(msg)

# Check for important dependencies
check_critical_dependencies()

invisible()
}


#' Check Critical Dependencies
#'
#' Verifies that critical Bioconductor packages are available
#'
#' @keywords internal
check_critical_dependencies <- function() {

  critical_packages <- c(
    "Biostrings",
    "GenomicRanges",
    "IRanges"
  )

  missing <- character()

  for (pkg in critical_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing) > 0) {
    packageStartupMessage(
      "\n⚠ Warning: Missing critical Bioconductor packages:\n",
      "  ", paste(missing, collapse = ", "), "\n",
      "\nInstall with:\n",
      "  BiocManager::install(c('", paste(missing, collapse = "', '"), "'))\n"
    )
  }

  invisible()
}


#' Package Unload Hook
#'
#' Runs when the package is unloaded
#'
#' @param libpath Library path
#' @keywords internal
#'

.onUnload <- function(libpath) {

# Clean up any temporary files in cache
cache_dir <- getOption("prism.cache_dir")

if (!is.null(cache_dir) && dir.exists(cache_dir)) {
  # Only clean PRISM-specific temp files
  prism_files <- list.files(cache_dir, pattern = "^prism_", full.names = TRUE)

  if (length(prism_files) > 0) {
    unlink(prism_files, recursive = TRUE)
  }
}

invisible()
}


#' Package Detach Hook
#'
#' Runs when the package is detached
#'
#' @param libpath Library path
#' @keywords internal
#'
.onDetach <- function(libpath) {

packageStartupMessage("PRISM detached. Thank you for using PRISM!")

invisible()
}


#' Get PRISM Options
#'
#' Retrieves PRISM package options
#'
#' @param option Option name (if NULL, returns all PRISM options)
#' @return Option value or list of all options
#' @export
#' @examples
#' get_prism_option("verbose")
#' get_prism_option()  # Get all options
get_prism_option <- function(option = NULL) {

  prism_options <- c(
    "prism.verbose",
    "prism.max_mismatches",
    "prism.min_probe_length",
    "prism.default_species",
    "prism.default_genome",
    "prism.batch_size",
    "prism.parallel_cores",
    "prism.cache_dir"
  )

  if (is.null(option)) {
    # Return all PRISM options
    opts <- options()[prism_options]
    names(opts) <- gsub("^prism\\.", "", names(opts))
    return(opts)
  } else {
    # Return specific option
    opt_name <- if (grepl("^prism\\.", option)) option else paste0("prism.", option)
    return(getOption(opt_name))
  }
}


#' Set PRISM Options
#'
#' Sets PRISM package options
#'
#' @param... Named arguments for options to set
#' @return Invisible list of previous option values
#' @export
#' @examples
#' set_prism_option(verbose = FALSE, max_mismatches = 3)
set_prism_option <- function(...) {

  opts <- list(...)

  if (length(opts) == 0) {
    stop("No options specified")
  }

  # Add prism. prefix if not present
  names(opts) <- ifelse(
    grepl("^prism\\.", names(opts)),
    names(opts),
    paste0("prism.", names(opts))
  )

  # Store old values
  old_opts <- options()[names(opts)]

  # Set new values
  do.call(options, opts)

  message("PRISM options updated:")
  for (name in names(opts)) {
    message("  ", gsub("^prism\\.", "", name), " = ", opts[[name]])
  }

  invisible(old_opts)
}


#' Reset PRISM Options to Defaults
#'
#' Resets all PRISM options to their default values
#'
#' @export
reset_prism_options <- function() {

  defaults <- list(
    prism.verbose = TRUE,
    prism.max_mismatches = 2,
    prism.min_probe_length = 20,
    prism.default_species = "hsapiens",
    prism.default_genome = "hg38",
    prism.batch_size = 500,
    prism.parallel_cores = 1,
    prism.cache_dir = tempdir()
  )

  options(defaults)

  message("PRISM options reset to defaults")

  invisible(defaults)
}


#' Show PRISM Configuration
#'
#' Displays current PRISM configuration
#'
#' @export
show_prism_config <- function() {

  cat("PRISM Configuration\n")
  cat("===================\n\n")

  cat("Package Version:", as.character(packageVersion("PRISM")), "\n\n")

  cat("Options:\n")
  opts <- get_prism_option()
  for (name in names(opts)) {
    cat("  ", name, ": ", opts[[name]], "\n", sep = "")
  }

  cat("\nSession Info:\n")
  cat("  R version:", R.version.string, "\n")
  cat("  Platform:", R.version$platform, "\n")

  cat("\nCritical Dependencies:\n")
  deps <- c("Biostrings", "GenomicRanges", "IRanges", "GEOquery", "biomaRt")
  for (pkg in deps) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      ver <- packageVersion(pkg)
      cat("  ✓", pkg, "(", as.character(ver), ")\n")
    } else {
      cat("  ✗", pkg, "(not installed)\n")
    }
  }

  invisible()
}


#' Check PRISM Installation
#'
#' Verifies that PRISM is properly installed and configured
#'
#' @return Logical indicating if installation is complete
#' @export
check_prism_installation <- function() {

  cat("Checking PRISM Installation\n")
  cat("============================\n\n")

  all_ok <- TRUE

  # Check R version
  cat("1. Checking R version... ")
  r_version <- as.numeric(paste0(R.version$major, ".", R.version$minor))
  if (r_version >= 4.0) {
    cat("✓ (", R.version.string, ")\n", sep = "")
  } else {
    cat("✗ R >= 4.0.0 required\n")
    all_ok <- FALSE
  }

  # Check critical packages
  cat("\n2. Checking critical packages:\n")
  critical <- c("Biostrings", "GenomicRanges", "IRanges")
  for (pkg in critical) {
    cat("   ", pkg, "... ", sep = "")
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("✓\n")
    } else {
      cat("✗ (not installed)\n")
      all_ok <- FALSE
    }
  }

  # Check recommended packages
  cat("\n3. Checking recommended packages:\n")
  recommended <- c("GEOquery", "biomaRt", "ggplot2", "dplyr")
  for (pkg in recommended) {
    cat("   ", pkg, "... ", sep = "")
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("✓\n")
    } else {
      cat("⚠ (recommended but not required)\n")
    }
  }

  # Check write permissions
  cat("\n4. Checking cache directory... ")
  cache_dir <- getOption("prism.cache_dir")
  if (dir.exists(cache_dir) && file.access(cache_dir, 2) == 0) {
    cat("✓\n")
  } else {
    cat("⚠ (may have permission issues)\n")
  }

  # Summary
  cat("\n")
  if (all_ok) {
    cat("✓ PRISM installation is complete and ready to use!\n")
  } else {
    cat("✗ PRISM installation is incomplete. Please install missing packages.\n")
    cat("\nInstall missing Bioconductor packages with:\n")
    cat("  BiocManager::install(c('Biostrings', 'GenomicRanges', 'IRanges'))\n")
  }

  invisible(all_ok)
}
