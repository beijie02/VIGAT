#' VIGAT: VLM-guided gene activity quantification from single-cell ATAC-seq
#'
#' VIGAT renders TSS-centered single-cell ATAC coverage plots for marker
#' genes, asks a vision-capable large language model to select one compact
#' genomic window per gene, and counts fragments in those windows to build a
#' gene-by-cell ATAC count matrix. The main function is \code{\link{VIGAT}}.
#'
#' @keywords internal
#' @import ggplot2
#' @importFrom SeuratObject DefaultAssay DefaultAssay<- Idents Idents<-
#' @importFrom Signac granges StringToGRanges Annotation CoveragePlot
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom dplyr %>%
#' @importFrom grDevices pdf dev.off
"_PACKAGE"

# Column names used inside dplyr verbs (silences R CMD check notes)
utils::globalVariables(c("gene", "chr", "gene_order", "marker_text"))
