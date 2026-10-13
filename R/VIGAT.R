# ============================================================
# VIGAT.R
#
# VIGAT: VLM-guided gene activity quantification from single-cell ATAC-seq.
# Main entry point: VIGAT()
#
# Usage:
#   library(VIGAT)
#   res <- VIGAT(
#     obj     = obj,
#     genes   = c("INS", "GCG", "SST"),
#     tissue  = "Islet",
#     species = "human"
#   )
#   res$vigat_windows      # GRanges, one window per gene
#   res$count_matrix       # gene x cell ATAC counts in those windows
#
# VIGAT() writes nothing to disk. Intermediate results and highlighted
# coverage plots are returned as R objects when requested.
# ============================================================


# Dependencies are declared in DESCRIPTION / R/VIGAT-package.R


# ============================================================
# HELPER FUNCTIONS  (internal; not intended for direct use)
# ============================================================

.vigat_get_tss_info <- function(gene, annot, flank = 5000) {

  ga_all <- annot[!is.na(annot$gene_name) & annot$gene_name == gene]
  if (length(ga_all) == 0) stop("Gene not found in annotation: ", gene)

  chr_use <- names(sort(table(as.character(seqnames(ga_all))), decreasing = TRUE))[1]
  ga_all  <- ga_all[as.character(seqnames(ga_all)) == chr_use]

  strand_cand <- names(sort(table(as.character(strand(ga_all))), decreasing = TRUE))
  strand_cand <- strand_cand[strand_cand != "*"]
  st <- if (length(strand_cand) > 0) strand_cand[1] else
    names(sort(table(as.character(strand(ga_all))), decreasing = TRUE))[1]
  if (!st %in% c("+", "-")) stop("Cannot determine strand for gene: ", gene)

  gene_start <- min(start(ga_all), na.rm = TRUE)
  gene_end   <- max(end(ga_all),   na.rm = TRUE)
  tss        <- if (st == "+") gene_start else gene_end
  tes        <- if (st == "+") gene_end   else gene_start

  region_start <- max(1L, as.integer(tss - flank))
  region_end   <- as.integer(tss + flank)

  list(
    region       = paste0(chr_use, "-", region_start, "-", region_end),
    tss          = as.integer(tss),
    tes          = as.integer(tes),
    strand       = st,
    chr          = chr_use,
    gene_start   = as.integer(gene_start),
    gene_end     = as.integer(gene_end),
    region_start = region_start,
    region_end   = region_end
  )
}


.vigat_extract_peaks <- function(obj, region_string) {

  DefaultAssay(obj) <- "ATAC"
  peaks_gr  <- granges(obj[["ATAC"]])
  region_gr <- StringToGRanges(region_string)
  hits      <- findOverlaps(peaks_gr, region_gr)

  if (length(hits) == 0) return(data.frame())

  out_gr        <- peaks_gr[queryHits(hits)]
  clipped_start <- pmax(start(out_gr), start(region_gr)[1])
  clipped_end   <- pmin(end(out_gr),   end(region_gr)[1])

  data.frame(
    peak  = paste0(as.character(seqnames(out_gr)), "-", clipped_start, "-", clipped_end),
    chr   = as.character(seqnames(out_gr)),
    start = clipped_start,
    end   = clipped_end,
    width = clipped_end - clipped_start + 1L,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(start <= end) %>%
    dplyr::distinct()
}


.vigat_coverage_theme <- function(p, tss) {
  p &
    scale_x_continuous(
      limits       = c(tss - 5000, tss + 5000),
      breaks       = seq(tss - 4000, tss + 4000, by = 2000),
      minor_breaks = seq(tss - 5000, tss + 5000, by = 200),
      labels       = scales::label_comma(),
      guide        = guide_axis(minor.ticks = TRUE)
    ) &
    theme(
      axis.ticks.x.minor        = element_line(color = "grey40", linewidth = 0.3),
      axis.ticks.length.x.minor = unit(0.18, "cm"),
      axis.ticks.x              = element_line(color = "black",  linewidth = 0.4),
      axis.ticks.length.x       = unit(0.28, "cm")
    )
}


# Coverage plot sent to the model (TSS +/- 5 kb, grouped by ATAC cluster).
# info_row needs: region, tss, chr, strand. title = FALSE drops title/subtitle.
.vigat_plot_coverage <- function(obj, info_row, group_by, tissue, title = FALSE) {
  p <- CoveragePlot(
    object     = obj,
    region     = info_row$region[1],
    group.by   = group_by,
    annotation = TRUE,
    peaks      = TRUE,
    links      = FALSE
  )
  p <- .vigat_coverage_theme(p, tss = info_row$tss[1])
  if (!title) return(p)
  p + patchwork::plot_annotation(
    title = sprintf(
      "%s | TSS \u00b15 kb | %s%s | grouped by ATAC clusters",
      info_row$gene[1], info_row$chr[1],
      ifelse(info_row$strand[1] == "+", "\u2192", "\u2190")
    ),
    subtitle = sprintf(
      "Tissue: %s | Region: %s | TSS: %s",
      tissue, info_row$region[1], scales::comma(info_row$tss[1])
    )
  )
}


# Same coverage plot with the selected window highlighted in blue.
# window_row needs: chr, start, end. title = FALSE drops title/subtitle.
.vigat_plot_highlight <- function(obj, info_row, window_row, group_by, tissue,
                                  title = FALSE) {
  p <- CoveragePlot(
    object     = obj,
    region     = info_row$region[1],
    group.by   = group_by,
    annotation = TRUE,
    peaks      = TRUE,
    links      = FALSE
  )

  p <- p &
    annotate("rect",
             xmin = window_row$start[1], xmax = window_row$end[1],
             ymin = -Inf,                ymax = Inf,
             fill = "#2196F3",           alpha = 0.15) &
    annotate("rect",
             xmin = window_row$start[1], xmax = window_row$end[1],
             ymin = -Inf,                ymax = Inf,
             fill = NA, color = "#0D47A1", linewidth = 0.7)

  p <- .vigat_coverage_theme(p, tss = info_row$tss[1])
  if (!title) return(p)

  p + patchwork::plot_annotation(
    title = sprintf(
      "%s | VIGAT window highlighted | TSS \u00b15 kb | %s%s",
      info_row$gene[1], info_row$chr[1],
      ifelse(info_row$strand[1] == "+", "\u2192", "\u2190")
    ),
    subtitle = sprintf(
      "Tissue: %s | VIGAT window: %s:%s\u2013%s",
      tissue,
      window_row$chr[1], scales::comma(window_row$start[1]),
      scales::comma(window_row$end[1])
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  )
}


.vigat_build_prompt <- function(species, tissue, markers_text, called_peaks_text) {

  context_origin <- if (!is.null(species) && nchar(trimws(species)) > 0) {
    paste0(species, " ", tissue)
  } else {
    tissue
  }

  paste0(
    "The uploaded PDF file(s) contain TSS-centered ATAC CoveragePlots for marker genes.",
    " For each marker gene plot, select one compact absolute genomic window.",
    " The selected window should be suitable for calculating the ATAC signal",
    " of that marker gene later.\n\n",
    "Context: ATAC data from ", context_origin, ".",
    " Cluster identities are unknown.",
    " Use the ATAC coverage signal, gene model, genomic coordinates,",
    " marker gene name, TSS direction provided after each gene name,",
    " and called Peaks if available.\n\n",
    "Markers:\n", markers_text, "\n\n",
    "Called Peaks:\n", called_peaks_text, "\n\n",
    "Instructions:\n\n",
    "1. First, identify the target marker gene for each coverage plot.\n",
    "2. Locate the marker gene's TSS/promoter region based on its strand orientation.",
    ' For a "+" strand gene, the promoter region lies toward decreasing genomic coordinates.',
    ' For a "-" strand gene, the promoter region lies toward increasing genomic coordinates.\n',
    "3. Choose one compact window near the marker gene's TSS/promoter region.\n",
    "    - If a called Peak overlaps with the promoter signal, prioritize a window that",
    " overlaps this Peak and captures the clearest cluster-specific ATAC signal.",
    " The window does not need to span the entire called Peak if only a subset of the",
    " region is informative. You may include a small, contiguous neighboring region if",
    " it helps resolve cluster differences.\n",
    "    - If no called Peak overlaps the promoter signal but a clear, cluster-specific",
    " coverage signal is visible near the TSS, define a compact window based on the",
    " coverage profile alone.\n",
    "4. Do not select a window solely based on the highest signal intensity.",
    " Avoid distal peaks far from the promoter, broad non-specific signals shared across",
    " many clusters, and signals originating from adjacent genes.\n",
    "5. Use a minimum window size of 1 kb by default.",
    " A window smaller than 1 kb should only be used if the promoter-side Peak or coverage",
    " signal is exceptionally sharp and well-isolated.",
    " Otherwise, expand the window modestly around the informative signal while avoiding",
    " unrelated neighboring noise.\n",
    "6. Before final output, verify that each selected start and end coordinate is within",
    " the displayed x-axis range for that gene.\n\n",
    "Output Format:\n",
    "First, provide one row per gene:\n",
    "Gene; TSS direction; chr:start-end; Size;",
    " Anchor evidence (called Peak(s) or coverage-only);",
    " Strongest useful clusters; Brief rationale.\n\n",
    "Then, provide the selected windows as R code for only the genes with a valid selected",
    " range. Keep the same gene order as in the row-format output, and omit genes whose",
    " range is NA or cannot be determined.\n\n",
    'vigat_windows <- GRanges(\n',
    '  seqnames = c(\n',
    '    "chr1", "chr2"\n',
    '  ),\n',
    '  ranges = IRanges(\n',
    '    start = c(\n',
    '      100000, 200000\n',
    '    ),\n',
    '    end = c(\n',
    '      101000, 201000\n',
    '    )\n',
    '  ),\n',
    '  gene = c(\n',
    '    "GENE1", "GENE2"\n',
    '  )\n',
    ')\n\n',
    "============================================================"
  )
}


.vigat_extract_response_text <- function(result) {

  msg_idx <- which(vapply(result$output,
                          function(x) identical(x$type, "message"), logical(1)))
  if (length(msg_idx) == 0) stop("No message output found in API response.")

  all_text <- character(0)
  for (idx in msg_idx) {
    content_i <- result$output[[idx]]$content
    if (is.null(content_i) || length(content_i) == 0) next
    txt_idx <- which(vapply(content_i,
                            function(x) identical(x$type, "output_text"), logical(1)))
    if (length(txt_idx) > 0)
      all_text <- c(all_text,
                    vapply(content_i[txt_idx], function(x) x$text, character(1)))
  }

  if (length(all_text) == 0) stop("No output_text found in API response.")
  paste(all_text, collapse = "\n")
}


.vigat_parse_windows <- function(raw_response) {

  extract_vec <- function(text, pattern) {
    m <- stringr::str_match(text, pattern)
    if (is.na(m[1, 2])) stop("Cannot parse VIGAT window output with pattern: ", pattern)
    m[1, 2]
  }

  seqnames_text <- extract_vec(raw_response,
    'seqnames\\s*=\\s*c\\(([\\s\\S]*?)\\)\\s*,\\s*ranges')
  start_text    <- extract_vec(raw_response,
    'start\\s*=\\s*c\\(([\\s\\S]*?)\\)\\s*,\\s*end')
  end_text      <- extract_vec(raw_response,
    'end\\s*=\\s*c\\(([\\s\\S]*?)\\)\\s*\\)\\s*,\\s*gene')
  gene_text     <- extract_vec(raw_response,
    'gene\\s*=\\s*c\\(([\\s\\S]*?)\\)\\s*\\)\\s*\\)?\\s*```')

  parse_char <- function(x) gsub('"', '', stringr::str_extract_all(x, '"[^"]+"')[[1]])
  parse_num  <- function(x) as.integer(unlist(strsplit(gsub("\\s+", "", x), ",")))

  chr   <- parse_char(seqnames_text)
  start <- parse_num(start_text)
  end   <- parse_num(end_text)
  gene  <- parse_char(gene_text)

  lens <- c(chr = length(chr), start = length(start), end = length(end), gene = length(gene))
  if (length(unique(lens)) != 1)
    stop("Parsed vectors have different lengths: ",
         paste(names(lens), lens, sep = "=", collapse = ", "))

  data.frame(gene  = gene,
             chr   = chr,
             start = start,
             end   = end,
             width = end - start + 1L,
             stringsAsFactors = FALSE)
}

.vigat_count_matrix <- function(obj, window_df, genes_use) {

  DefaultAssay(obj) <- "ATAC"

  window_df <- window_df %>%
    dplyr::filter(gene %in% genes_use) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  if (nrow(window_df) == 0) {
    stop("No VIGAT windows available for selected genes.")
  }

  gr <- GenomicRanges::GRanges(
    seqnames = window_df$chr,
    ranges = IRanges::IRanges(
      start = window_df$start,
      end = window_df$end
    ),
    gene = window_df$gene
  )

  names(gr) <- window_df$gene

  count_mat <- Signac::FeatureMatrix(
    fragments = Signac::Fragments(obj),
    features = gr,
    cells = colnames(obj)
  )

  rownames(count_mat) <- names(gr)

  list(
    windows   = gr,
    window_df = window_df,
    counts   = count_mat
  )
}

.vigat_save_plots_pdf <- function(plots_list, outdir, file_prefix,
                             genes_per_page, width = 22, height = 22) {

  groups    <- split(names(plots_list), ceiling(seq_along(plots_list) / genes_per_page))
  out_files <- character(0)

  for (pp in seq_along(groups)) {
    out_pdf <- file.path(outdir, paste0(file_prefix, "_page", pp, ".pdf"))
    pdf(out_pdf, width = width, height = height)
    print(patchwork::wrap_plots(plots_list[groups[[pp]]], ncol = 2))
    dev.off()
    out_files <- c(out_files, out_pdf)
    message("Saved: ", out_pdf)
  }

  invisible(out_files)
}


# ============================================================
# MAIN FUNCTION
# ============================================================

#' VLM-guided gene activity quantification for one sample
#'
#' Draws TSS-centered ATAC coverage plots for each marker gene, asks a
#' vision-language model to select one compact genomic window per gene, and
#' counts fragments in those windows to build a gene-by-cell count matrix.
#' Nothing is written to disk; everything is returned as R objects.
#'
#' @param obj          A processed Signac Seurat object with an ATAC assay
#'                     and gene annotation. Set Idents(obj) to the desired
#'                     ATAC clustering before calling if needed.
#' @param genes        Character vector of marker gene names.
#' @param tissue       Tissue label injected into the API prompt (e.g. "PBMC").
#' @param species      Species label injected into the API prompt (e.g. "human").
#'                     Set to NULL to omit from prompt.
#' @param cluster_col  Optional. Name of a meta.data column to use as the ATAC
#'                     cluster assignment. If NULL (default), Idents(obj) is used.
#' @param min_atac_cluster_cells  Clusters with fewer cells than this are excluded.
#'                     Default: 50.
#' @param return_intermediates  Logical. If TRUE, also return the intermediate
#'                     results: TSS info, overlapping peaks, the coverage plots
#'                     sent to the model, the prompt, the full API response and
#'                     run information. Default: FALSE.
#' @param return_highlighted_plots  Logical. If TRUE, also return one coverage
#'                     plot per gene with the selected window highlighted.
#'                     Default: FALSE.
#' @param openai_model  OpenAI model name. Default: "gpt-5.5".
#' @section API key:
#' This function reads the OpenAI API key from the \code{OPENAI_API_KEY}
#' environment variable. Set it before calling the function:
#' \preformatted{
#' Sys.setenv(OPENAI_API_KEY = "your_key_here")
#' }
#' or add it to your \code{~/.Renviron} file:
#' \preformatted{
#' OPENAI_API_KEY=your_key_here
#' }
#'
#' @return A named list:
#'   \item{vigat_windows}{GRanges containing VIGAT-selected genomic windows.}
#'   \item{vigat_window_df}{Data frame containing the selected genomic windows.}
#'   \item{count_matrix}{Sparse gene-by-cell raw ATAC count matrix generated
#'     by counting fragments overlapping the VIGAT-selected windows using
#'     \code{Signac::FeatureMatrix()}.}
#'   \item{diagnostics}{List containing \code{raw_response},
#'     \code{genes_found}, and \code{genes_skipped}.}
#'   \item{intermediates}{Only if \code{return_intermediates = TRUE}: list with
#'     \code{tss_info}, \code{overlapping_peaks}, \code{coverage_plots}
#'     (named list of ggplot objects), \code{prompt}, \code{api_response}
#'     and \code{api_run_info}.}
#'   \item{highlighted_plots}{Only if \code{return_highlighted_plots = TRUE}:
#'     named list of ggplot objects, one per gene.}
#'
#' @examples
#' \dontrun{
#' Sys.setenv(OPENAI_API_KEY = "your_key_here")
#' Idents(obj) <- "ATAC_snn_res.0.8"
#' res <- VIGAT(
#'   obj     = obj,
#'   genes   = c("INS", "GCG", "SST"),
#'   tissue  = "Islet",
#'   species = "human"
#' )
#' res$vigat_window_df
#'
#' # Keep the highlighted plots and save them yourself
#' res <- VIGAT(obj, genes = c("INS", "GCG"), tissue = "Islet",
#'              return_highlighted_plots = TRUE)
#' for (g in names(res$highlighted_plots)) {
#'   ggplot2::ggsave(paste0(g, "_vigat.png"), res$highlighted_plots[[g]],
#'                   width = 7, height = 7, dpi = 150)
#' }
#' }
#'
#' @export
VIGAT <- function(
  obj,
  genes,
  tissue,
  species                  = "human",
  cluster_col              = NULL,
  min_atac_cluster_cells   = 50,
  openai_model             = "gpt-5.5",
  return_intermediates     = FALSE,
  return_highlighted_plots = FALSE
) {

  genes_per_page <- 4L   # fixed: max 4 plots per PDF page, 2-column grid

  # ---- Read API key from environment ----
  openai_api_key <- Sys.getenv("OPENAI_API_KEY")
  if (nchar(openai_api_key) == 0)
    stop(
      "OPENAI_API_KEY is not set. Please set your own API key before running:\n",
      "  Sys.setenv(OPENAI_API_KEY = \"your_key_here\")\n",
      "or add it to your ~/.Renviron file:\n",
      "  OPENAI_API_KEY=your_key_here"
    )

  # ---- Validate inputs ----
  if (!inherits(obj, "Seurat"))
    stop("'obj' must be a Seurat object.")
  if (!"ATAC" %in% names(obj@assays))
    stop("'obj' must contain an ATAC assay.")
  if (length(genes) == 0)
    stop("'genes' must be a non-empty character vector.")
  if (nchar(trimws(tissue)) == 0)
    stop("'tissue' must be a non-empty string.")

  cat("=== VIGAT ===\n")
  cat("Tissue:  ", tissue, "\n")
  cat("Species: ", if (!is.null(species)) species else "(not specified)", "\n")
  cat("Genes:   ", paste(genes, collapse = ", "), "\n\n")

  DefaultAssay(obj) <- "ATAC"


  # ----------------------------------------------------------
  # STEP 1: Resolve cluster assignment
  # ----------------------------------------------------------

  if (is.null(cluster_col)) {
    cluster_use <- Idents(obj)
    cat("Using active identity Idents(obj) as ATAC cluster assignment.\n")
  } else {
    if (!cluster_col %in% colnames(obj@meta.data))
      stop("cluster_col '", cluster_col, "' not found in object metadata.")
    cluster_use <- obj[[cluster_col]][, 1]
    cluster_use <- as.character(cluster_use)
    names(cluster_use) <- colnames(obj)

    Idents(obj) <- cluster_use

    cat("Using cluster_col '", cluster_col, "'.\n", sep = "")
  }

  if (length(cluster_use) == 0 || all(is.na(cluster_use)))
    stop(
      "No valid cluster assignment found. Set Idents(obj) to the ATAC cluster ",
      "assignment before calling, or provide cluster_col."
    )


  # ----------------------------------------------------------
  # STEP 2: Filter small clusters
  # ----------------------------------------------------------

  # Remove NA cells
  valid_cells <- names(cluster_use)[!is.na(cluster_use)]
  cluster_use <- cluster_use[valid_cells]
  obj         <- subset(obj, cells = valid_cells)

  # Identify clusters to keep
  cluster_counts <- table(cluster_use)
  clusters_keep  <- names(cluster_counts)[cluster_counts >= min_atac_cluster_cells]
  small_removed  <- setdiff(names(cluster_counts), clusters_keep)

  if (length(small_removed) > 0)
    cat("Filtering out", length(small_removed),
        "small cluster(s) (<", min_atac_cluster_cells, "cells):",
        paste(small_removed, collapse = ", "), "\n")

  if (length(clusters_keep) < 2)
    stop(
      "After filtering clusters with fewer than ", min_atac_cluster_cells,
      " cells, fewer than 2 clusters remain. ",
      "Consider lowering min_atac_cluster_cells."
    )

  cells_keep   <- names(cluster_use)[cluster_use %in% clusters_keep]
  obj_filt     <- subset(obj, cells = cells_keep)

  # Store as a clean unified column used by all downstream steps
  obj_filt$ATAC_cluster_use <- as.character(cluster_use[cells_keep])
  Idents(obj_filt)          <- "ATAC_cluster_use"
  atac_col                  <- "ATAC_cluster_use"

  DefaultAssay(obj_filt) <- "ATAC"

  cat("Cells retained:", ncol(obj_filt), "of", ncol(obj), "\n")
  cat("ATAC cluster sizes:\n")
  print(table(obj_filt[[atac_col]][, 1]))


  # ----------------------------------------------------------
  # STEP 3: TSS info and overlapping peaks
  # ----------------------------------------------------------

  cat("\n-- Computing TSS info and overlapping peaks...\n")

  if (is.null(Annotation(obj_filt[["ATAC"]])))
    stop(
      "ATAC assay annotation is missing. ",
      "Add gene annotation (e.g. via Signac::Annotation<-) before calling."
    )

  annot <- Annotation(obj_filt[["ATAC"]])

  marker_tss_list   <- list()
  marker_peak_list  <- list()
  genes_skipped     <- character(0)

  for (gene_use in genes) {

    cat("  TSS:", gene_use, "\n")

    info <- tryCatch(
      .vigat_get_tss_info(gene = gene_use, annot = annot, flank = 5000),
      error = function(e) {
        cat("    Skipped \u2014", conditionMessage(e), "\n")
        genes_skipped <<- c(genes_skipped, gene_use)
        NULL
      }
    )
    if (is.null(info)) next

    marker_tss_list[[gene_use]] <- data.frame(
      gene         = gene_use,
      chr          = info$chr,
      strand       = info$strand,
      tss          = info$tss,
      region_start = info$region_start,
      region_end   = info$region_end,
      region       = info$region,
      stringsAsFactors = FALSE
    )

    peak_df <- .vigat_extract_peaks(obj = obj_filt, region_string = info$region)
    if (nrow(peak_df) > 0) {
      peak_df$gene   <- gene_use
      peak_df$region <- info$region
    }
    marker_peak_list[[gene_use]] <- peak_df
  }

  marker_tss_df   <- dplyr::bind_rows(marker_tss_list)
  marker_peaks_df <- dplyr::bind_rows(marker_peak_list)
  atac_peak_gr <- Signac::granges(obj_filt[["ATAC"]])
atac_peak_chr <- unique(as.character(GenomicRanges::seqnames(atac_peak_gr)))

unsupported_chr_df <- marker_tss_df %>%
  dplyr::filter(!chr %in% atac_peak_chr)

if (nrow(unsupported_chr_df) > 0) {
  warning(
    "Skipping genes because chromosome is absent from ATAC peak set: ",
    paste0(unsupported_chr_df$gene, "(", unsupported_chr_df$chr, ")", collapse = ", ")
  )
  genes_skipped <- unique(c(genes_skipped, unsupported_chr_df$gene))
}

marker_tss_df <- marker_tss_df %>%
  dplyr::filter(chr %in% atac_peak_chr)

if (is.null(marker_peaks_df) || nrow(marker_peaks_df) == 0 || !"gene" %in% colnames(marker_peaks_df)) {
  marker_peaks_df <- data.frame(
    peak = character(),
    chr = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    gene = character(),
    region = character(),
    stringsAsFactors = FALSE
  )
} else {
  marker_peaks_df <- marker_peaks_df %>%
    dplyr::filter(gene %in% marker_tss_df$gene)
}

genes_found <- marker_tss_df$gene

  cat("Genes found:", paste(genes_found, collapse = ", "), "\n")
  if (length(genes_skipped) > 0)
    cat("Genes skipped:", paste(genes_skipped, collapse = ", "), "\n")

  if (length(genes_found) == 0)
    stop("No marker genes found in the annotation. Cannot proceed.")



  # ----------------------------------------------------------
  # STEP 4: Coverage plots
  #   Always generated (needed for API upload). PDFs go to a
  #   temporary folder that is deleted when VIGAT() returns.
  # ----------------------------------------------------------

  cat("\n-- Generating coverage plots...\n")

  plots <- lapply(genes_found, function(gene_use) {

  info_row <- marker_tss_df %>%
    dplyr::filter(gene == gene_use) %>%
    dplyr::slice(1)

  if (nrow(info_row) == 0) {
    warning("Skipping ", gene_use, ": no TSS info found.")
    genes_skipped <<- unique(c(genes_skipped, gene_use))
    return(NULL)
  }

  p <- tryCatch(
    .vigat_plot_coverage(obj_filt, info_row, group_by = atac_col, tissue = tissue, title = FALSE),
    error = function(e) {
      warning(
        "Skipping ", gene_use,
        ": CoveragePlot failed for region ",
        info_row$region[1],
        ". Error: ",
        conditionMessage(e)
      )
      genes_skipped <<- unique(c(genes_skipped, gene_use))
      return(NULL)
    }
  )

  p
})

names(plots) <- genes_found

plots <- plots[!vapply(plots, is.null, logical(1))]
genes_found <- names(plots)

marker_tss_df <- marker_tss_df %>%
  dplyr::filter(gene %in% genes_found)

marker_peaks_df <- marker_peaks_df %>%
  dplyr::filter(gene %in% genes_found)

if (length(genes_found) == 0) {
  stop("No coverage plots could be generated for any marker gene. Cannot proceed.")
}
  names(plots) <- genes_found

  coverage_subdir <- tempfile("vigat_coverage_")
  dir.create(coverage_subdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(coverage_subdir, recursive = TRUE), add = TRUE)

  coverage_files <- .vigat_save_plots_pdf(
    plots_list     = plots,
    outdir         = coverage_subdir,
    file_prefix    = "coverage_marker_genes_TSS5kb_ATAC_clusters",
    genes_per_page = genes_per_page
  )


  # ----------------------------------------------------------
  # STEP 5: Build prompt
  # ----------------------------------------------------------

  markers_text <- marker_tss_df %>%
    dplyr::mutate(gene_order = match(gene, genes_found)) %>%
    dplyr::arrange(gene_order) %>%
    dplyr::mutate(marker_text = paste0(gene, "(", strand, ")")) %>%
    dplyr::pull(marker_text) %>%
    paste(collapse = "; ")

  called_peaks_lines <- lapply(genes_found, function(gene_use) {
  
  if (is.null(marker_peaks_df) || !"gene" %in% colnames(marker_peaks_df)) {
    peaks_i <- data.frame()
  } else {
    peaks_i <- marker_peaks_df %>%
      dplyr::filter(gene == gene_use)
  }
  
  if (nrow(peaks_i) == 0) return(NULL)
  peak_strings <- paste0(peaks_i$chr, ":", peaks_i$start, "-", peaks_i$end)
  paste0(gene_use, ": ", paste(unique(peak_strings), collapse = "; "))
})

called_peaks_text <- paste(
  unlist(called_peaks_lines[!vapply(called_peaks_lines, is.null, logical(1))]),
  collapse = "\n"
)

if (is.null(called_peaks_text) || length(called_peaks_text) == 0 || called_peaks_text == "") {
  called_peaks_text <- "No called peaks overlap the TSS-centered regions for these marker genes."
}

  prompt <- .vigat_build_prompt(species, tissue, markers_text, called_peaks_text)



  # ----------------------------------------------------------
  # STEP 6: OpenAI API call
  # ----------------------------------------------------------

  pdf_sizes_mb <- file.info(coverage_files)$size / 1024^2
  if (!all(file.exists(coverage_files)))
    stop("Coverage PDF files missing before API call.")
  if (any(pdf_sizes_mb >= 50))
    stop("At least one PDF is >= 50 MB (OpenAI limit).")
  if (sum(pdf_sizes_mb) >= 50)
    stop("Combined PDF size is >= 50 MB (OpenAI limit).")

  cat("\n-- Calling OpenAI API (", openai_model, ")...\n", sep = "")
  cat("   PDFs:", length(coverage_files),
      "| Size:", round(sum(pdf_sizes_mb), 3), "MB\n")

  file_content <- lapply(coverage_files, function(f) {
    list(
      type      = "input_file",
      filename  = basename(f),
      file_data = paste0("data:application/pdf;base64,",
                         base64enc::base64encode(f)),
      detail    = "high"
    )
  })

  request_content <- c(file_content, list(list(type = "input_text", text = prompt)))

  req <- httr2::request("https://api.openai.com/v1/responses") |>
    httr2::req_headers(
      Authorization  = paste("Bearer", openai_api_key),
      `Content-Type` = "application/json"
    ) |>
    httr2::req_body_json(
      list(model = openai_model,
           input = list(list(role = "user", content = request_content))),
      auto_unbox = TRUE
    )

  start_time      <- Sys.time()
  resp            <- req |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_perform()
  elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  status_code     <- httr2::resp_status(resp)

  cat("HTTP status:", status_code,
      "| Elapsed:", round(elapsed_seconds, 1), "seconds\n")

  if (status_code != 200) {
    err_body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "(unreadable)")
    stop("API call failed (HTTP ", status_code, "):\n", err_body)
  }

  result       <- httr2::resp_body_json(resp)
  raw_response <- .vigat_extract_response_text(result)

  cat("\n=== API Response ===\n", raw_response, "\n")

  api_run_info <- data.frame(
    species           = if (!is.null(species)) species else NA_character_,
    tissue            = tissue,
    model             = openai_model,
    n_input_genes     = length(genes),
    n_genes_found     = length(genes_found),
    n_genes_skipped   = length(genes_skipped),
    n_pdf_files       = length(coverage_files),
    total_pdf_size_mb = sum(pdf_sizes_mb),
    http_status       = status_code,
    response_id       = result$id,
    elapsed_seconds   = elapsed_seconds,
    input_tokens      = result$usage$input_tokens ,
    output_tokens     = result$usage$output_tokens,
    total_tokens      = result$usage$total_tokens ,
    stringsAsFactors  = FALSE
  )

# ----------------------------------------------------------
# STEP 7: Parse VIGAT-selected windows
# ----------------------------------------------------------

cat("\n-- Parsing VIGAT-selected windows...\n")

vigat_window_df <- .vigat_parse_windows(raw_response)

cat("VIGAT windows parsed:\n")
print(vigat_window_df)

# ----------------------------------------------------------
# STEP 8: Generate VIGAT gene-by-cell count matrix
# ----------------------------------------------------------

cat("\n-- Generating VIGAT gene-by-cell count matrix...\n")

vigat_quant <- .vigat_count_matrix(
  obj       = obj_filt,
  window_df = vigat_window_df,
  genes_use = genes_found
)

vigat_windows   <- vigat_quant$windows
vigat_window_df <- vigat_quant$window_df
vigat_counts    <- vigat_quant$counts

cat(
  "VIGAT count matrix:",
  nrow(vigat_counts), "genes x",
  ncol(vigat_counts), "cells\n"
)


  # ----------------------------------------------------------
  # STEP 9: Highlighted coverage plots  (optional)
  # ----------------------------------------------------------

  plots_hl <- NULL

  if (return_highlighted_plots) {

    cat("\n-- Generating highlighted coverage plots...\n")

    genes_with_window <- genes_found[genes_found %in% vigat_window_df$gene]

    plots_hl <- lapply(genes_with_window, function(gene_use) {

      cat("  ", gene_use, "\n")

      tryCatch({

        info_row <- marker_tss_df   %>% dplyr::filter(gene == gene_use) %>% dplyr::slice(1)
        gr_row   <- vigat_window_df %>% dplyr::filter(gene == gene_use) %>% dplyr::slice(1)

        .vigat_plot_highlight(obj_filt, info_row, gr_row,
                              group_by = atac_col, tissue = tissue)

      }, error = function(e) {
        cat("  Failed:", gene_use, "|", conditionMessage(e), "\n")
        NULL
      })
    })
    names(plots_hl) <- genes_with_window
    plots_hl        <- plots_hl[!vapply(plots_hl, is.null, logical(1))]
  }


  # ----------------------------------------------------------
  # DONE
  # ----------------------------------------------------------

  cat("\n=== Done ===\n")

  out <- list(
    vigat_windows   = vigat_windows,
    vigat_window_df = vigat_window_df,
    count_matrix    = vigat_counts,
    diagnostics     = list(
      raw_response  = raw_response,
      genes_found   = genes_found,
      genes_skipped = genes_skipped
    )
  )

  if (return_intermediates) {
    out$intermediates <- list(
      tss_info          = marker_tss_df,
      overlapping_peaks = marker_peaks_df,
      coverage_plots    = plots,
      prompt            = prompt,
      api_response      = result,
      api_run_info      = api_run_info
    )
  }

  if (return_highlighted_plots) {
    out$highlighted_plots <- plots_hl
  }

  invisible(out)
}