#' Update Gene Annotations from Ensembl
#'
#' Retrieves current gene annotations for probe genomic locations using biomaRt.
#' Maps probes to genes based on genomic overlap.
#'
#' @param probe_ranges GRanges object with probe locations (from align_probes_to_genome)
#' @param species Species name: "hsapiens", "mmusculus", "rnorvegicus", etc. (default: "hsapiens")
#' @param ensembl_version Ensembl version to use (default: NULL for current)
#' @param host Ensembl host URL (default: NULL for main site)
#' @param gene_biotypes Gene biotypes to include (default: NULL for all)
#' @param expand_range Expand probe range by N bp for nearby genes (default: 0)
#' @return Data frame with updated gene annotations
#' @export
#' @importFrom biomaRt useMart getBM listAttributes
#' @importFrom GenomicRanges findOverlaps
#' @examples
#' \dontrun{
#'   # Update annotations for human probes
#'   annotations <- update_gene_annotations(
#'     probe_ranges,
#'     species = "hsapiens"
#'   )
#'
#'   # Include nearby genes within 5kb
#'   annotations <- update_gene_annotations(
#'     probe_ranges,
#'     species = "hsapiens",
#'     expand_range = 5000
#'   )
#' }
update_gene_annotations <- function(probe_ranges,
                                    species = "hsapiens",
                                    ensembl_version = NULL,
                                    host = NULL,
                                    gene_biotypes = NULL,
                                    expand_range = 0) {

  # Validate inputs
  if (!inherits(probe_ranges, "GRanges")) {
    stop("probe_ranges must be a GRanges object")
  }

  if (length(probe_ranges) == 0) {
    stop("probe_ranges is empty")
  }

  message("Updating gene annotations for ", length(probe_ranges), " probe locations")
  message("Species: ", species)

  # Expand ranges if requested
  if (expand_range > 0) {
    message("Expanding probe ranges by +/- ", expand_range, " bp")
    probe_ranges <- GenomicRanges::resize(probe_ranges,
                                          width = width(probe_ranges) + 2 * expand_range,
                                          fix = "center")
  }

  # Connect to Ensembl
  message("Connecting to Ensembl BioMart...")

  dataset_name <- paste0(species, "_gene_ensembl")

  mart <- tryCatch({
    if (!is.null(host)) {
      biomaRt::useMart("ensembl", dataset = dataset_name, host = host)
    } else if (!is.null(ensembl_version)) {
      biomaRt::useMart("ensembl", dataset = dataset_name, version = ensembl_version)
    } else {
      biomaRt::useMart("ensembl", dataset = dataset_name)
    }
  }, error = function(e) {
    stop("Failed to connect to Ensembl: ", e$message,
         "\nCheck species name and internet connection.")
  })

  message("Connected to: ", mart@dataset)

  # Get unique chromosomes and ranges
  chromosomes <- unique(as.character(seqnames(probe_ranges)))

  # Remove non-standard chromosomes if needed
  chromosomes <- chromosomes[!grepl("_|\\.", chromosomes)]

  message("Querying ", length(chromosomes), " chromosomes...")

  # Define attributes to retrieve
  attributes <- c(
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "description",
    "entrezgene_id"
  )

  # Query BioMart in batches by chromosome
  all_genes <- list()

  for (chr in chromosomes) {
    message("  Querying chromosome ", chr, "...")

    # Get probes on this chromosome
    chr_probes <- probe_ranges[seqnames(probe_ranges) == chr]

    if (length(chr_probes) == 0) next

    # Get range for this chromosome
    chr_start <- min(start(chr_probes))
    chr_end <- max(end(chr_probes))

    # Query biomaRt
    filters <- c("chromosome_name", "start", "end")
    values <- list(
      chromosome_name = chr,
      start = chr_start,
      end = chr_end
    )

    genes <- tryCatch({
      biomaRt::getBM(
        attributes = attributes,
        filters = filters,
        values = values,
        mart = mart
      )
    }, error = function(e) {
      warning("Failed to query chromosome ", chr, ": ", e$message)
      return(NULL)
    })

    if (!is.null(genes) && nrow(genes) > 0) {
      all_genes[[chr]] <- genes
    }
  }

  # Combine all chromosomes
  if (length(all_genes) == 0) {
    stop("No gene annotations retrieved from Ensembl")
  }

  gene_annotations <- do.call(rbind, all_genes)
  message("Retrieved ", nrow(gene_annotations), " gene records from Ensembl")

  # Filter by biotype if specified
  if (!is.null(gene_biotypes)) {
    message("Filtering for biotypes: ", paste(gene_biotypes, collapse = ", "))
    gene_annotations <- gene_annotations[
      gene_annotations$gene_biotype %in% gene_biotypes,
    ]
    message("Retained ", nrow(gene_annotations), " genes after biotype filtering")
  }

  # Convert to GRanges
  gene_gr <- GenomicRanges::GRanges(
    seqnames = gene_annotations$chromosome_name,
    ranges = IRanges::IRanges(
      start = gene_annotations$start_position,
      end = gene_annotations$end_position
    ),
    strand = ifelse(gene_annotations$strand == 1, "+", "-"),
    gene_id = gene_annotations$ensembl_gene_id,
    gene_name = gene_annotations$external_gene_name,
    biotype = gene_annotations$gene_biotype,
    description = gene_annotations$description,
    entrez_id = gene_annotations$entrezgene_id
  )

  # Find overlaps between probes and genes
  message("Mapping probes to genes...")
  overlaps <- GenomicRanges::findOverlaps(probe_ranges, gene_gr)

  if (length(overlaps) == 0) {
    warning("No overlaps found between probes and genes!")
    return(data.frame())
  }

  # Create annotation table
  result <- data.frame(
    probe_id = probe_ranges$probe_id[queryHits(overlaps)],
    chr = as.character(seqnames(probe_ranges))[queryHits(overlaps)],
    probe_start = start(probe_ranges)[queryHits(overlaps)],
    probe_end = end(probe_ranges)[queryHits(overlaps)],
    probe_strand = as.character(strand(probe_ranges))[queryHits(overlaps)],
    ensembl_gene_id = gene_gr$gene_id[subjectHits(overlaps)],
    gene_name = gene_gr$gene_name[subjectHits(overlaps)],
    gene_biotype = gene_gr$biotype[subjectHits(overlaps)],
    gene_start = start(gene_gr)[subjectHits(overlaps)],
    gene_end = end(gene_gr)[subjectHits(overlaps)],
    gene_strand = as.character(strand(gene_gr))[subjectHits(overlaps)],
    gene_description = gene_gr$description[subjectHits(overlaps)],
    entrez_id = gene_gr$entrez_id[subjectHits(overlaps)],
    stringsAsFactors = FALSE
  )

  message("Annotated ", length(unique(result$probe_id)), " probes to ",
          length(unique(result$ensembl_gene_id)), " genes")

  return(result)
}

#' Annotate Probes Using Bioconductor OrgDb Packages
#'
#' Alternative to biomaRt using local annotation packages.
#' Faster but requires package installation.
#'
#' @param probe_ranges GRanges object with probe locations
#' @param txdb TxDb annotation package (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene)
#' @param orgdb OrgDb annotation package (e.g., org.Hs.eg.db)
#' @return Data frame with gene annotations
#' @export
#' @importFrom GenomicFeatures genes
#' @importFrom AnnotationDbi select mapIds
#' @examples
#' \dontrun{
#'   library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'   library(org.Hs.eg.db)
#'
#'   annotations <- annotate_probes_orgdb(
#'     probe_ranges,
#'     TxDb.Hsapiens.UCSC.hg38.knownGene,
#'     org.Hs.eg.db
#'   )
#' }
annotate_probes_orgdb <- function(probe_ranges, txdb, orgdb) {

  if (!inherits(probe_ranges, "GRanges")) {
    stop("probe_ranges must be a GRanges object")
  }

  message("Extracting gene annotations from TxDb...")

  # Get all genes from TxDb
  all_genes <- GenomicFeatures::genes(txdb)

  # Find overlaps
  message("Finding overlaps between probes and genes...")
  overlaps <- GenomicRanges::findOverlaps(probe_ranges, all_genes)

  if (length(overlaps) == 0) {
    warning("No overlaps found!")
    return(data.frame())
  }

  # Get gene IDs
  gene_ids <- all_genes$gene_id[subjectHits(overlaps)]

  # Map to gene symbols and other info
  message("Retrieving gene symbols and descriptions...")

  gene_symbols <- AnnotationDbi::mapIds(
    orgdb,
    keys = gene_ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  gene_names <- AnnotationDbi::mapIds(
    orgdb,
    keys = gene_ids,
    column = "GENENAME",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  ensembl_ids <- AnnotationDbi::mapIds(
    orgdb,
    keys = gene_ids,
    column = "ENSEMBL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  # Create result data frame
  result <- data.frame(
    probe_id = probe_ranges$probe_id[queryHits(overlaps)],
    chr = as.character(seqnames(probe_ranges))[queryHits(overlaps)],
    probe_start = start(probe_ranges)[queryHits(overlaps)],
    probe_end = end(probe_ranges)[queryHits(overlaps)],
    probe_strand = as.character(strand(probe_ranges))[queryHits(overlaps)],
    entrez_id = gene_ids,
    gene_name = unname(gene_symbols[gene_ids]),
    ensembl_gene_id = unname(ensembl_ids[gene_ids]),
    gene_description = unname(gene_names[gene_ids]),
    gene_start = start(all_genes)[subjectHits(overlaps)],
    gene_end = end(all_genes)[subjectHits(overlaps)],
    gene_strand = as.character(strand(all_genes))[subjectHits(overlaps)],
    stringsAsFactors = FALSE
  )

  message("Annotated ", length(unique(result$probe_id)), " probes to ",
          length(unique(result$entrez_id)), " genes")

  return(result)
}

#' Assign Probes to Genes Based on Overlap Rules
#'
#' Handles complex cases where probes map to multiple genes
#'
#' @param annotations Data frame from update_gene_annotations
#' @param priority Priority rule: "longest_overlap", "protein_coding", "single" (default: "longest_overlap")
#' @param keep_all Keep all mappings or just best match (default: FALSE)
#' @return Data frame with probe-gene assignments
#' @export
assign_probes_to_genes <- function(annotations,
                                  priority = "longest_overlap",
                                  keep_all = FALSE) {

  if (nrow(annotations) == 0) {
    return(annotations)
  }

  message("Assigning probes to genes using priority: ", priority)

  if (keep_all) {
    message("Keeping all probe-gene mappings")
    return(annotations)
  }

  # Calculate overlap length
  annotations$overlap_length <- pmin(annotations$probe_end, annotations$gene_end) -
                                pmax(annotations$probe_start, annotations$gene_start) + 1

  # Group by probe and apply priority rules
  result <- annotations %>%
    dplyr::group_by(probe_id) %>%
    dplyr::arrange(dplyr::desc(overlap_length)) %>%
    dplyr::filter(
      if (priority == "longest_overlap") {
        overlap_length == max(overlap_length)
      } else if (priority == "protein_coding") {
        gene_biotype == "protein_coding" | dplyr::row_number() == 1
      } else if (priority == "single") {
        dplyr::n() == 1
      } else {
        TRUE
      }
    ) %>%
    dplyr::slice(1) %>%  # Take first if still tied
    dplyr::ungroup()

  result <- as.data.frame(result)

  message("Final mapping: ", nrow(result), " probe-gene pairs")

  return(result)
}

#' Retrieve Additional Gene Metadata
#'
#' Fetches extra information like GO terms, pathways, etc.
#'
#' @param gene_ids Vector of Ensembl gene IDs
#' @param species Species name (default: "hsapiens")
#' @param attributes Additional biomaRt attributes to retrieve
#' @return Data frame with gene metadata
#' @export
#' @importFrom biomaRt useMart getBM
#' @examples
#' \dontrun{
#'   metadata <- get_gene_metadata(
#'     c("ENSG00000139618", "ENSG00000141510"),
#'     species = "hsapiens"
#'   )
#' }
get_gene_metadata <- function(gene_ids,
                              species = "hsapiens",
                              attributes = c("go_id", "name_1006",
                                             "namespace_1003")) {

  if (length(gene_ids) == 0) {
    stop("gene_ids is empty")
  }

  # Remove NAs
  gene_ids <- gene_ids[!is.na(gene_ids)]

  if (length(gene_ids) == 0) {
    warning("No valid gene IDs provided")
    return(data.frame())
  }

  message("Retrieving metadata for ", length(gene_ids), " genes...")

  # Connect to Ensembl
  dataset_name <- paste0(species, "_gene_ensembl")
  mart <- biomaRt::useMart("ensembl", dataset = dataset_name)

  # Combine with basic gene info
  all_attributes <- unique(c("ensembl_gene_id", "external_gene_name", attributes))

  # Query in batches to avoid timeout
  batch_size <- 500
  n_batches <- ceiling(length(gene_ids) / batch_size)

  metadata_list <- list()

  for (i in seq_len(n_batches)) {
    batch_start <- (i - 1) * batch_size + 1
    batch_end <- min(i * batch_size, length(gene_ids))
    batch_ids <- gene_ids[batch_start:batch_end]

    message("  Batch ", i, "/", n_batches)

    batch_data <- tryCatch({
      biomaRt::getBM(
        attributes = all_attributes,
        filters = "ensembl_gene_id",
        values = batch_ids,
        mart = mart
      )
    }, error = function(e) {
      warning("Failed to retrieve batch ", i, ": ", e$message)
      return(NULL)
    })

    if (!is.null(batch_data) && nrow(batch_data) > 0) {
      metadata_list[[i]] <- batch_data
    }
  }

  if (length(metadata_list) == 0) {
    warning("No metadata retrieved")
    return(data.frame())
  }

  metadata <- do.call(rbind, metadata_list)

  message("Retrieved metadata for ", length(unique(metadata$ensembl_gene_id)), " genes")

  return(metadata)
}

#' Compare Original and Updated Annotations
#'
#' Compares probe annotations before and after re-annotation to identify changes
#'
#' @param original_annotations Data frame with original probe annotations
#' @param updated_annotations Data frame with updated annotations from update_gene_annotations
#' @param probe_id_col Column name for probe IDs in original data (default: "ID")
#' @param gene_col Column name for gene symbols in original data (default: "Gene Symbol")
#' @return Data frame showing annotation changes
#' @export
#' @examples
#' \dontrun{
#'   comparison <- compare_annotations(
#'     old_annotations,
#'     new_annotations,
#'     probe_id_col = "Probe_ID",
#'     gene_col = "GENE_SYMBOL"
#'   )
#' }
compare_annotations <- function(original_annotations,
                                updated_annotations,
                                probe_id_col = "ID",
                                gene_col = "Gene Symbol") {

  message("Comparing original and updated annotations...")

  # Standardize column names
  original <- data.frame(
    probe_id = original_annotations[[probe_id_col]],
    original_gene = original_annotations[[gene_col]],
    stringsAsFactors = FALSE
  )

  # Aggregate updated annotations (in case of multiple genes per probe)
  updated <- updated_annotations %>%
    dplyr::group_by(probe_id) %>%
    dplyr::summarise(
      updated_gene = paste(unique(gene_name), collapse = ";"),
      n_genes = dplyr::n_distinct(gene_name),.groups = "drop"
    ) %>%
    as.data.frame()

  # Merge
  comparison <- dplyr::full_join(original, updated, by = "probe_id")

  # Clean up gene symbols for comparison
  comparison$original_gene <- trimws(as.character(comparison$original_gene))
  comparison$updated_gene <- trimws(as.character(comparison$updated_gene))

  # Replace empty strings and "---" with NA
  comparison$original_gene[comparison$original_gene %in% c("", "---", "NA")] <- NA
  comparison$updated_gene[comparison$updated_gene %in% c("", "---", "NA")] <- NA

  # Classify changes
  comparison$annotation_status <- dplyr::case_when(
    is.na(comparison$original_gene) & is.na(comparison$updated_gene) ~ "no_annotation",
    is.na(comparison$original_gene) & !is.na(comparison$updated_gene) ~ "newly_annotated",
    !is.na(comparison$original_gene) & is.na(comparison$updated_gene) ~ "annotation_lost",
    comparison$original_gene == comparison$updated_gene ~ "unchanged",
    comparison$original_gene != comparison$updated_gene ~ "changed",
    TRUE ~ "unknown"
  )

  # Add flag for multiple gene mappings
  comparison$multiple_genes <- !is.na(comparison$n_genes) & comparison$n_genes > 1

  # Summary
  message("\nAnnotation Comparison Summary:")
  message("  Total probes: ", nrow(comparison))
  message("  Unchanged: ", sum(comparison$annotation_status == "unchanged"))
  message("  Changed: ", sum(comparison$annotation_status == "changed"))
  message("  Newly annotated: ", sum(comparison$annotation_status == "newly_annotated"))
  message("  Annotation lost: ", sum(comparison$annotation_status == "annotation_lost"))
  message("  No annotation: ", sum(comparison$annotation_status == "no_annotation"))
  message("  Multiple genes: ", sum(comparison$multiple_genes, na.rm = TRUE))

  return(comparison)
}

#' Summarize Probe Annotation Results
#'
#' Creates summary statistics for annotation results
#'
#' @param annotations Data frame from update_gene_annotations
#' @return List with summary statistics
#' @export
summarize_probe_annotations <- function(annotations) {

  if (nrow(annotations) == 0) {
    return(list(
      total_probes = 0,
      annotated_probes = 0,
      unique_genes = 0,
      probes_per_gene = data.frame(),
      genes_per_probe = data.frame(),
      biotype_distribution = data.frame()
    ))
  }

  # Count probes and genes
  n_probes <- length(unique(annotations$probe_id))
  n_genes <- length(unique(annotations$ensembl_gene_id))

  # Probes per gene
  probes_per_gene <- annotations %>%
    dplyr::group_by(ensembl_gene_id, gene_name) %>%
    dplyr::summarise(n_probes = dplyr::n_distinct(probe_id),.groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_probes)) %>%
    as.data.frame()

  # Genes per probe
  genes_per_probe <- annotations %>%
    dplyr::group_by(probe_id) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(ensembl_gene_id),.groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_genes)) %>%
    as.data.frame()

  # Biotype distribution
  if ("gene_biotype" %in% colnames(annotations)) {
    biotype_dist <- annotations %>%
      dplyr::group_by(gene_biotype) %>%
      dplyr::summarise(
        n_genes = dplyr::n_distinct(ensembl_gene_id),
        n_probes = dplyr::n_distinct(probe_id),.groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n_genes)) %>%
      as.data.frame()
  } else {
    biotype_dist <- data.frame()
  }

  summary_list <- list(
    total_probes = n_probes,
    annotated_probes = n_probes,
    unique_genes = n_genes,
    avg_probes_per_gene = round(nrow(annotations) / n_genes, 2),
    avg_genes_per_probe = round(nrow(annotations) / n_probes, 2),
    probes_per_gene = probes_per_gene,
    genes_per_probe = genes_per_probe,
    biotype_distribution = biotype_dist
  )

  class(summary_list) <- c("prism_annotation_summary", "list")

  return(summary_list)
}

#' @export
print.prism_annotation_summary <- function(x,...) {
  cat("PRISM Annotation Summary\n")
  cat("========================\n\n")

  cat("Probes annotated:", x$annotated_probes, "\n")
  cat("Unique genes:", x$unique_genes, "\n")
  cat("Avg probes per gene:", x$avg_probes_per_gene, "\n")
  cat("Avg genes per probe:", x$avg_genes_per_probe, "\n")

  if (nrow(x$biotype_distribution) > 0) {
    cat("\nGene Biotype Distribution:\n")
    print(head(x$biotype_distribution, 10), row.names = FALSE)
    if (nrow(x$biotype_distribution) > 10) {
      cat("... and", nrow(x$biotype_distribution) - 10, "more\n")
    }
  }

  cat("\nTop genes by probe count:\n")
  print(head(x$probes_per_gene, 10), row.names = FALSE)

  multi_gene_probes <- sum(x$genes_per_probe$n_genes > 1)
  if (multi_gene_probes > 0) {
    cat("\nProbes mapping to multiple genes:", multi_gene_probes, "\n")
    cat("Top multi-mapping probes:\n")
    print(head(x$genes_per_probe[x$genes_per_probe$n_genes > 1, ], 5),
          row.names = FALSE)
  }

  invisible(x)
}
