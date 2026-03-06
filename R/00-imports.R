#' @keywords internal
"_PACKAGE"

#' @import methods
#' @importFrom Biostrings DNAStringSet matchPDict PDict
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom dplyr %>% filter mutate group_by summarise
#' @importFrom stats median
NULL

# Suppress R CMD check notes for common variable names
utils::globalVariables(c(
  "probe_id", "gene_id", "gene_name", "chromosome",
  "start", "end", "strand", ".", "n_genes"
))
