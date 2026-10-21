# ============================================================================
# DEG LIMMA STUDIO | RStudio / Local Shiny Desktop-Style Application
# Developed by FR | Department of Zoology, University of Rajshahi
# Source: user-provided "DEGs identification R code (1).R"
# Run in RStudio: Source this entire file, or source(file.choose())
# Input: processed, normalized, LOG2-scale microarray expression CSV
#        rows = genes, columns = samples, first column = gene identifier
# Optional: Sample metadata CSV (Sample, Group)
# ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

.install_if_missing <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = NA)
  }
  still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) stop("Could not install: ", paste(still_missing, collapse = ", "))
}

.install_if_missing(c("shiny", "ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "data.table"))
if (!requireNamespace("limma", quietly = TRUE)) {
  .install_if_missing("BiocManager")
  message("Installing limma from Bioconductor; this is a one-time operation.")
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("limma is unavailable. Check your R/Bioconductor compatibility and internet connection.")
}

options(shiny.maxRequestSize = 1024 * 1024^2)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1L])) y else x

get_gse <- function(x) {
  x <- as.character(x)
  vapply(x, function(label) {
    if (is.na(label)) return(NA_character_)
    match <- regexpr("GSE[0-9]+", label, ignore.case = TRUE, perl = TRUE)
    if (match[1L] < 0L) return(NA_character_)
    toupper(substr(label, match[1L], match[1L] + attr(match, "match.length")[1L] - 1L))
  }, character(1), USE.NAMES = FALSE)
}

# The R upload keeps a temporary datapath but retains the original file name.
read_expression <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  tab <- tryCatch({
    if (ext %in% c("tsv", "txt")) {
      data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE,
                        na.strings = c("", "NA", "NaN"), showProgress = FALSE)
    } else {
      data.table::fread(path, data.table = FALSE, check.names = FALSE,
                        na.strings = c("", "NA", "NaN"), showProgress = FALSE)
    }
  }, error = function(e) stop("Could not open expression file: ", conditionMessage(e)))
  if (ncol(tab) < 3L || nrow(tab) < 5L) {
    stop("Expected a matrix with genes in rows and at least two sample columns.")
  }

  # Original uploaded code used the first column as the gene identifier.
  raw_gene <- trimws(as.character(tab[[1L]]))
  valid_gene <- !is.na(raw_gene) & nzchar(raw_gene) & raw_gene != "NA"
  tab <- tab[valid_gene, , drop = FALSE]
  raw_gene <- raw_gene[valid_gene]
  if (!length(raw_gene)) stop("No usable gene identifiers in first column.")

  raw_names <- names(tab)[-1L]
  if (any(!nzchar(raw_names))) stop("One or more sample columns have empty names.")
  # Accept numeric strings while rejecting phenotype/annotation character columns.
  is_sample <- vapply(tab[-1L], function(z) {
    if (is.numeric(z) || is.integer(z)) return(TRUE)
    values <- trimws(as.character(z))
    valid <- !(is.na(values) | values == "" | values == "NA")
    if (!any(valid)) return(FALSE)
    parsed <- suppressWarnings(as.numeric(values[valid]))
    all(is.finite(parsed))
  }, logical(1))
  if (sum(is_sample) < 4L) {
    stop("At least four numeric sample columns are needed (two per biological group).")
  }
  ignored <- raw_names[!is_sample]
  tab <- tab[, c(TRUE, is_sample), drop = FALSE]
  original_sample_names <- names(tab)[-1L]
  if (anyDuplicated(original_sample_names)) {
    stop("Duplicate sample column names found. Supply unique sample identifiers first.")
  }
  values <- lapply(tab[-1L], function(z) suppressWarnings(as.numeric(z)))
  mat <- do.call(cbind, values)
  if (!is.matrix(mat)) mat <- matrix(mat, ncol = length(values))
  storage.mode(mat) <- "double"
  rownames(mat) <- raw_gene
  colnames(mat) <- original_sample_names
  ids <- get_gse(original_sample_names)
  filename_id <- get_gse(name)[1L]
  if (all(is.na(ids)) && !is.na(filename_id)) ids[] <- filename_id
  list(mat = mat, sample_ids = ids, ignored = ignored, file_name = name,
       duplicate_genes = sum(duplicated(raw_gene)))
}

read_groups <- function(path, filename) {
  ext <- tolower(tools::file_ext(filename))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  df <- utils::read.table(path, sep = sep, header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE, quote = '"', comment.char = "",
                          fileEncoding = "UTF-8-BOM")
  normalized <- tolower(gsub("[^a-z0-9]", "", names(df)))
  id_col <- which(normalized %in% c("sample", "sampleid", "gsm", "sampleaccession"))[1]
  group_col <- which(normalized %in% c("group", "condition", "phenotype", "status", "class"))[1]
  if (is.na(id_col) || is.na(group_col)) {
    stop("Metadata CSV needs columns named Sample and Group (or Sample_ID/Condition).")
  }
  ans <- data.frame(Sample = trimws(as.character(df[[id_col]])),
                    Group = trimws(as.character(df[[group_col]])), stringsAsFactors = FALSE)
  if (anyDuplicated(ans$Sample)) stop("Metadata contains duplicate sample identifiers.")
  ans
}

normalize_group <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("control", "normal", "healthy", "non-disease", "nondisease", "0"),
         "Control", ifelse(x %in% c("disease", "case", "patient", "affected", "1"),
                            "Disease", NA_character_))
}

determine_groups <- function(samples, metadata, control_pattern, disease_pattern) {
  if (!is.null(metadata)) {
    idx <- match(samples, metadata$Sample)
    raw_group <- metadata$Group[idx]
    return(normalize_group(raw_group))
  }
  ctl <- grepl(control_pattern, samples, ignore.case = TRUE, perl = TRUE)
  dis <- grepl(disease_pattern, samples, ignore.case = TRUE, perl = TRUE)
  ifelse(ctl & !dis, "Control", ifelse(dis & !ctl, "Disease", NA_character_))
}

preprocess_for_limma <- function(mat, group) {
  if (!is.numeric(mat)) stop("Expression matrix is not numeric.")
  n_input <- nrow(mat)
  clean <- rowSums(!is.finite(mat)) == 0L
  mat <- mat[clean, , drop = FALSE]
  removed_na <- n_input - nrow(mat)
  if (nrow(mat) < 5) stop("Too few complete genes after removing missing/nonfinite rows.")

  # Duplicate gene symbols are combined by mean, with exact provenance in report.
  ngenes <- length(unique(rownames(mat)))
  duplicates <- nrow(mat) - ngenes
  if (duplicates > 0L) {
    group_ids <- rownames(mat)
    summed <- rowsum(mat, group = group_ids, reorder = FALSE)
    denom <- tabulate(match(group_ids, rownames(summed)), nbins = nrow(summed))
    mat <- sweep(summed, 1L, denom, "/")
  }
  variances <- apply(mat, 1L, stats::var)
  keep <- is.finite(variances) & variances > 1e-12
  removed_constant <- sum(!keep)
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 5L) stop("Too few nonconstant genes remain for DEG analysis.")
  if (anyNA(mat) || any(!is.finite(mat))) stop("Nonfinite expression entries remain.")
  if (length(unique(group)) != 2L || min(table(group)) < 2L) {
    stop("At least two Control and two Disease samples are required.")
  }
  list(mat = mat, removed_na = removed_na, duplicates = duplicates,
       removed_constant = removed_constant, n_before = n_input)
}

pca_figure <- function(mat, group, max_genes) {
  variances <- apply(mat, 1L, stats::var)
  nuse <- min(length(variances), as.integer(max_genes))
  selected <- order(variances, decreasing = TRUE)[seq_len(nuse)]
  pca <- stats::prcomp(t(mat[selected, , drop = FALSE]), center = TRUE, scale. = TRUE)
  percentages <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  df <- data.frame(PC1 = pca$x[, 1L], PC2 = pca$x[, 2L], Group = group,
                   Sample = colnames(mat))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 3.0, alpha = 0.88) +
    ggplot2::scale_color_manual(values = c(Control = "#3173AE", Disease = "#D65348")) +
    ggplot2::labs(title = "PCA of Gene Expression Data",
      x = sprintf("PC1 (%.1f%%)", percentages[1L]),
      y = sprintf("PC2 (%.1f%%)", percentages[2L]))
  if (min(table(group)) >= 3L) {
    p <- p + ggplot2::stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.7,
                                    show.legend = FALSE)
  }
  p + publication_theme()
}

publication_theme <- function() {
  ggplot2::theme_bw(base_size = 13) + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 17),
    axis.title = ggplot2::element_text(face = "bold", size = 14),
    axis.text = ggplot2::element_text(face = "bold", color = "#18324B", size = 11),
    legend.text = ggplot2::element_text(face = "bold", size = 11),
    legend.title = ggplot2::element_text(face = "bold", size = 12),
    panel.grid.minor = ggplot2::element_blank())
}

volcano_figure <- function(deg, fc, fdr, nlabels) {
  deg <- as.data.frame(deg)
  deg$negLogFDR <- -log10(pmax(deg$adj.P.Val, .Machine$double.xmin))
  sig <- deg[deg$Regulation != "NotSig" & is.finite(deg$logFC), , drop = FALSE]
  if (nrow(sig)) {
    sig <- sig[order(-abs(sig$logFC), sig$adj.P.Val), , drop = FALSE]
    sig <- head(sig, nlabels)
  }
  p <- ggplot2::ggplot(deg,
          ggplot2::aes(x = logFC, y = negLogFDR, color = Regulation)) +
    ggplot2::geom_point(alpha = 0.67, size = 1.7) +
    ggplot2::scale_color_manual(
       values = c(Up = "#D73027", Down = "#4575B4", NotSig = "grey70"),
       labels = c(Up = "Upregulated", Down = "Downregulated", NotSig = "Not significant"),
       breaks = c("Up", "Down", "NotSig"), drop = FALSE) +
    ggplot2::geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "#263647") +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "#263647") +
    ggplot2::labs(title = "Volcano Plot of Differential Expression (limma)",
      x = "log2 Fold Change (Disease vs Control)", y = "-log10 Adjusted p-value", color = "Status") +
    publication_theme()
  if (nrow(sig)) {
    p <- p + ggrepel::geom_text_repel(
      data = sig, ggplot2::aes(label = Gene), size = 3.6, color = "#1B2C40",
      max.overlaps = Inf, box.padding = 0.32, point.padding = 0.24, segment.size = 0.3,
      show.legend = FALSE)
  }
  p
}

heatmap_figure <- function(mat, group, deg, top_n) {
  sig <- deg[deg$Regulation != "NotSig", , drop = FALSE]
  if (!nrow(sig)) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No significant genes at the selected cutoffs",
                        size = 5, fontface = "bold", color = "#263647") +
      ggplot2::theme_void() + ggplot2::labs(title = "Top DEG Heatmap") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold")))
  }
  sig <- sig[order(sig$adj.P.Val, -abs(sig$logFC)), , drop = FALSE]
  genes <- head(intersect(sig$Gene, rownames(mat)), as.integer(top_n))
  data <- mat[genes, , drop = FALSE]
  z <- t(scale(t(data)))
  z[!is.finite(z)] <- 0
  annotation <- data.frame(Group = group, row.names = colnames(mat))
  anncols <- list(Group = c(Control = "#3173AE", Disease = "#D65348"))
  palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255)
  pheatmap::pheatmap(z, color = palette, cluster_rows = nrow(z) > 1L,
                     cluster_cols = ncol(z) > 1L, annotation_col = annotation,
                     annotation_colors = anncols, show_colnames = FALSE, show_rownames = TRUE,
                     fontsize = 9, fontsize_row = 8.5, fontsize_col = 8.5,
                     main = sprintf("Top %d Differentially Expressed Genes", nrow(z)),
                     border_color = NA, silent = TRUE)
}

run_deg_analysis <- function(mat, group, logfc, fdr, heatmap_n, pca_genes, labels) {
  prepared <- preprocess_for_limma(mat, group)
  expr <- prepared$mat
  design <- stats::model.matrix(~ group)
  colnames(design) <- c("Intercept", "Disease_vs_Control")
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = "Disease_vs_Control", number = Inf,
                        sort.by = "P", adjust.method = "BH")
  tt <- data.frame(Gene = rownames(tt), tt, row.names = NULL, check.names = FALSE)
  need <- c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  for (nm in setdiff(need, names(tt))) tt[[nm]] <- NA_real_
  deg <- tt[, need, drop = FALSE]
  deg$Regulation <- ifelse(is.finite(deg$adj.P.Val) & deg$adj.P.Val < fdr & deg$logFC >= logfc,
                            "Up", ifelse(is.finite(deg$adj.P.Val) & deg$adj.P.Val < fdr &
                            deg$logFC <= -logfc, "Down", "NotSig"))
  deg <- deg[order(deg$P.Value, na.last = TRUE), , drop = FALSE]
  up <- deg[deg$Regulation == "Up", , drop = FALSE]
  down <- deg[deg$Regulation == "Down", , drop = FALSE]
  list(deg = deg, up = up, down = down, expr = expr, group = group,
       pca = pca_figure(expr, group, pca_genes),
       volcano = volcano_figure(deg, logfc, fdr, labels),
       heatmap = heatmap_figure(expr, group, deg, heatmap_n),
       qc = prepared)
}

save_grid_figure <- function(obj, file, fmt, width, height, dpi) {
  if (inherits(obj, "ggplot")) {
    ggplot2::ggsave(filename = file, plot = obj, width = width, height = height,
                    dpi = dpi, units = "in", limitsize = FALSE, bg = "white")
    return(invisible(NULL))
  }
  if (fmt == "png") grDevices::png(file, width = width, height = height,
                                     units = "in", res = dpi, bg = "white")
  else if (fmt %in% c("tif", "tiff")) grDevices::tiff(file, width = width, height = height,
                                     units = "in", res = dpi, compression = "lzw", bg = "white")
  else if (fmt == "svg") grDevices::svg(file, width = width, height = height,
                                        bg = "white", onefile = FALSE)
  else stop("Unsupported figure format: ", fmt)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(obj$gtable)
  invisible(NULL)
}

save_analysis <- function(r, base_dir, identifier, selected_file, dataset_filter, fc, fdr,
                          h_n, p_n, nlabels, formats, dpi, ignored_samples) {
  if (!length(formats)) stop("Select at least one figure format.")
  clean_id <- gsub("[^A-Za-z0-9_-]", "_", identifier)
  parent <- file.path(base_dir, paste0("DEG_", clean_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (!dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create analysis directory: ", parent)
  }
  data_dir <- file.path(parent, "Tables")
  figures_dir <- file.path(parent, "Figures")
  dir.create(data_dir)
  dir.create(figures_dir)

  data.table::fwrite(r$deg[, setdiff(names(r$deg), "Regulation"), drop = FALSE],
                     file.path(data_dir, "deg_all.csv"))
  data.table::fwrite(r$deg, file.path(data_dir, "deg_filtered.csv"))
  data.table::fwrite(r$up, file.path(data_dir, "deg_upregulated.csv"))
  data.table::fwrite(r$down, file.path(data_dir, "deg_downregulated.csv"))
  data.table::fwrite(data.frame(Sample = colnames(r$expr), Group = as.character(r$group)),
                     file.path(data_dir, "sample_metadata_used.csv"))
  settings <- data.frame(
    Parameter = c("Input file", "Selected GEO subset", "log2FC cutoff", "adjusted p-value cutoff",
                  "Heatmap top N", "PCA variable genes", "Volcano labels", "Genes input",
                  "Genes retained", "Removed nonfinite genes", "Collapsed duplicate gene rows",
                  "Removed constant genes", "Control count", "Disease count", "Figure DPI"),
    Value = as.character(c(selected_file, dataset_filter, fc, fdr, h_n, p_n, nlabels,
      r$qc$n_before, nrow(r$expr), r$qc$removed_na, r$qc$duplicates,
      r$qc$removed_constant, sum(r$group == "Control"), sum(r$group == "Disease"), dpi)))
  data.table::fwrite(settings, file.path(parent, "analysis_settings.csv"))

  for (fmt in formats) {
    save_grid_figure(r$pca, file.path(figures_dir, paste0("PCA_before_DEG.", fmt)),
                     fmt, 8.5, 6.5, dpi)
    save_grid_figure(r$volcano, file.path(figures_dir, paste0("Volcano_plot.", fmt)),
                     fmt, 9.2, 6.7, dpi)
    save_grid_figure(r$heatmap, file.path(figures_dir, paste0("Heatmap_top_DEGs.", fmt)),
                     fmt, 8.0, 9.3, dpi)
  }
  summary <- c(
    "Differential Expression Analysis — limma",
    paste("Dataset:", identifier),
    paste("Input:", selected_file),
    paste("Selected within-file GSE:", dataset_filter),
    paste("Contrast: Disease - Control"),
    paste("Controls:", sum(r$group == "Control"), "| Disease:", sum(r$group == "Disease")),
    paste("Retained genes:", nrow(r$expr)),
    paste("Criteria: absolute log2FC >=", fc, "and BH adjusted p-value <", fdr),
    paste("Upregulated genes:", nrow(r$up)),
    paste("Downregulated genes:", nrow(r$down)),
    paste("Total significant genes:", nrow(r$up) + nrow(r$down)),
    paste("Ignored annotation columns:", paste(ignored_samples, collapse = ", ")),
    "Input must already be normalized and appropriate for limma (log2-scale microarray expression).",
    "This tool does not perform batch correction or infer patient independence.",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  writeLines(summary, file.path(parent, "DEG_results_summary.txt"), useBytes = TRUE)
  invisible(parent)
}

studio_css <- "
body {background:#edf3f8 !important; color:#172f4b !important; font-family:'Segoe UI',Arial,sans-serif;}
.container-fluid {padding-left:22px;padding-right:22px;max-width:1650px;}
.hero {position:relative;display:flex;align-items:center;justify-content:center;text-align:center;
  border-radius:13px;background:linear-gradient(110deg,#132d4b,#2665a8);padding:24px 120px;
  min-height:144px;color:#fff;box-shadow:0 5px 16px rgba(19,45,75,.14);
  margin-top:14px;margin-bottom:18px;}
.hero-content {width:100%;max-width:1100px;margin:0 auto;}
.hero h2 {margin:0 0 8px;font-weight:800;color:#fff;font-size:29px;line-height:1.2;text-align:center;}
.hero p {color:#e0edf9;margin:3px 0;font-size:14px;text-align:center;}
.hero-logo {position:absolute;left:28px;top:50%;transform:translateY(-50%);width:84px;height:84px;
  background:#fff;border:3px solid #93c9ed;border-radius:50%;display:flex;align-items:center;
  justify-content:center;box-shadow:0 4px 16px rgba(0,0,0,.16);}
.hero-logo svg {width:72px;height:72px;display:block;}
@media(max-width:900px){.hero{padding:118px 22px 22px;min-height:210px;}
 .hero-logo{top:19px;left:50%;transform:translateX(-50%);width:77px;height:77px;}
 .hero h2{font-size:24px;}.hero p{font-size:12px;}}
.panel, .well {background:white;border:1px solid #d3dfeb;box-shadow:none;border-radius:12px;}
.well {padding:17px;}
.nav-tabs>li>a {font-weight:700;color:#183d64 !important;}
.nav-tabs>li.active>a {color:#0f4d8f !important;}
.form-control,.selectize-input,input[type=number],select {background-color:white!important;
 color:#142e48!important;border:1px solid #9fb7ca!important;}
.btn-primary {background:#206bb0;border-color:#206bb0;font-weight:700;}
.btn-success {background:#0a8067;border-color:#0a8067;font-weight:700;}
.btn {font-weight:700;border-radius:7px;}
label {color:#173c64;font-weight:700;}
.card {background:white;border:1px solid #d1dfec;border-radius:11px;padding:16px;margin:12px 0;}
.card h4 {margin:0 0 9px;font-weight:800;color:#193c60;}
.subtitle {color:#5d7690;font-size:13px;}
.metric {display:inline-block;background:#eaf2fb;border:1px solid #bfd6ee;border-radius:7px;
 padding:11px 15px;margin:4px;color:#173c63;font-weight:800;}
#status {font-weight:700;color:#173c63;}
.table {font-size:12px;}
.modal-content {background:#fff;color:#182e44;border-radius:12px;}
.modal-title {font-weight:800;color:#173d62;}
.guide-heading {font-weight:800;color:#173c63;margin:14px 0 7px;}
.guide-intro {color:#445d77;font-size:14px;line-height:1.7;}
.guide-steps {counter-reset:guide-step;list-style:none;padding-left:0;margin:12px 0;}
.guide-steps li {counter-increment:guide-step;position:relative;padding:13px 15px 13px 57px;
 margin-bottom:9px;border:1px solid #d5e2ef;background:#f7fbff;border-radius:10px;
 font-size:14px;line-height:1.6;color:#233d58;}
.guide-steps li:before {content:counter(guide-step);position:absolute;left:13px;top:13px;
 width:31px;height:31px;border-radius:50%;display:flex;align-items:center;
 justify-content:center;background:#2566a3;color:white;font-weight:800;}
.guide-steps li strong {color:#123d65;font-weight:800;}
.about-banner {padding:17px 19px;border:1px solid #d0e1f1;background:#f2f8ff;border-radius:11px;}
.about-name {font-size:20px;font-weight:800;color:#173d65;margin:0 0 7px;}
.about-table {width:100%;margin-top:15px;border-collapse:collapse;font-size:14px;}
.about-table th {width:118px;color:#173c63;text-align:left;vertical-align:top;padding:11px 8px;border-bottom:1px solid #e1eaf3;}
.about-table td {padding:11px 8px;border-bottom:1px solid #e1eaf3;color:#173248;}
.about-table a {color:#15588d;font-weight:700;text-decoration:underline;}
.about-table a:hover {color:#0b416f;}
.about-note {font-size:12px;color:#61788e;margin-top:12px;}
"

ui <- shiny::fluidPage(
  shiny::tags$head(shiny::tags$style(shiny::HTML(studio_css))),
  shiny::div(class = "hero",
      shiny::div(class = "hero-logo", title = "FR — Fazlay Rabbi",
        shiny::HTML('<svg viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg" aria-label="FR research logo" role="img"><g fill="none" stroke="#3e9fc2" stroke-width="2.5" stroke-linecap="round"><path d="M22 32L34 19L49 27L61 15L77 27M22 32L18 52L31 65L51 58L72 73M49 27L51 58L77 27"/></g><g fill="#2686b8"><circle cx="22" cy="32" r="4"/><circle cx="34" cy="19" r="4"/><circle cx="49" cy="27" r="3.4"/><circle cx="61" cy="15" r="3.7"/><circle cx="77" cy="27" r="4"/><circle cx="18" cy="52" r="3.4"/><circle cx="31" cy="65" r="3.4"/><circle cx="51" cy="58" r="4"/><circle cx="72" cy="73" r="4"/></g><rect x="16" y="41" rx="9" width="68" height="43" fill="white" opacity=".94"/><text x="50" y="75" font-size="43" text-anchor="middle" font-family="Arial,sans-serif" font-weight="900" fill="#153b66">FR</text></svg>')),
      shiny::div(class = "hero-content",
        shiny::tags$h2("DEG Identification & limma Analysis Studio"),
        shiny::tags$p("Dataset selection  |  Disease vs Control  |  Adjustable log2FC/FDR  |  PCA • Volcano • Heatmap"),
        shiny::tags$p("Developed by Md Fazlay Rabbi  |  Department of Zoology, University of Rajshahi"))),
  shiny::sidebarLayout(
    shiny::sidebarPanel(width = 4,
      shiny::h4("1. Expression dataset", style="font-weight:800;color:#163e67;"),
      shiny::fileInput("expression_files", "Upload one or more GSE expression matrices",
          multiple = TRUE, accept = c(".csv", ".tsv", ".txt")),
      shiny::selectInput("which_file", "Select dataset file", choices = character()),
      shiny::selectInput("gse_subset", "Select GSE within file (if present)",
                         choices = c("All samples in selected file" = "ALL")),
      shiny::fileInput("sample_metadata", "Optional sample metadata (Sample, Group)",
                       accept = c(".csv", ".tsv", ".txt")),
      shiny::tags$hr(),
      shiny::h4("2. Biological-group detection", style="font-weight:800;color:#163e67;"),
      shiny::textInput("control_regex", "Control sample-name pattern", "control|healthy|normal"),
      shiny::textInput("disease_regex", "Disease sample-name pattern", "disease|patient|case"),
      shiny::checkboxInput("exclude_unknown", "Explicitly exclude unassigned samples", FALSE),
      shiny::checkboxInput("confirmed_scale", "I confirm these are processed, normalized, log2-scale values", FALSE),
      shiny::tags$hr(),
      shiny::h4("3. DEG threshold settings", style="font-weight:800;color:#163e67;"),
      shiny::numericInput("fc_cutoff", "Absolute log2 Fold Change cutoff", value=1, min=0, max=25, step=0.1),
      shiny::numericInput("fdr_cutoff", "Adjusted p-value cutoff (BH FDR)", value=0.05,
                          min=0.000001, max=0.999999, step=0.01),
      shiny::numericInput("top_heatmap", "Top significant genes for heatmap", value=50, min=1, max=500, step=1),
      shiny::numericInput("pca_genes", "Most-variable genes for PCA", value=2000, min=10, max=40000, step=500),
      shiny::numericInput("volcano_labels", "Number of labelled significant genes", value=15, min=0,max=50),
      shiny::tags$hr(),
      shiny::h4("4. Figure exports", style="font-weight:800;color:#163e67;"),
      shiny::checkboxGroupInput("formats", "Image formats", choices=c("PNG"="png", "SVG"="svg", "TIFF"="tiff"),
                                selected=c("png","svg","tiff"), inline=TRUE),
      shiny::selectInput("dpi", "Raster resolution", choices=c("300 DPI"=300,"600 DPI"=600), selected=600),
      shiny::actionButton("run", "Run DEG Analysis", class = "btn-primary btn-lg", width="100%"),
      shiny::br(), shiny::br(),
      shiny::actionButton("save", "Choose Folder & Save Results", class="btn-success", width="100%"),
      shiny::p("Analysis runs first. After completion, the application asks where to save.", class="subtitle"),
      shiny::tags$div(id="status",shiny::textOutput("status_text"))
    ),
    shiny::mainPanel(width = 8,
      shiny::tabsetPanel(id="tabs",
        shiny::tabPanel("Sample Review",
          shiny::div(class="card",shiny::h4("Dataset and sample summary"),shiny::tableOutput("info_table")),
          shiny::div(class="card",shiny::h4("Control / Disease sample classification"),
              shiny::tableOutput("counts_table"),
              shiny::tableOutput("sample_table"),
              shiny::p("Unassigned samples block analysis unless 'exclude' is enabled.",class="subtitle"))),
        shiny::tabPanel("Results",
          shiny::div(class="card",shiny::h4("Differential expression overview"),
             shiny::uiOutput("results_summary")),
          shiny::div(class="card",shiny::h4("Top genes by nominal p-value"),
             shiny::tableOutput("deg_preview"))),
        shiny::tabPanel("PCA",shiny::div(class="card",shiny::plotOutput("pca_view", height="470px"))),
        shiny::tabPanel("Volcano",shiny::div(class="card",shiny::plotOutput("volcano_view", height="490px"))),
        shiny::tabPanel("Heatmap",shiny::div(class="card",shiny::plotOutput("heatmap_view", height="650px"))),
        shiny::tabPanel("User Guide", shiny::div(class="card",
          shiny::h4("User Guide — Step-by-Step DEG Analysis"),
          shiny::tags$p(class="guide-intro",
             "Follow these steps in order. No edits to R source code are needed between analyses."),
          shiny::tags$ol(class="guide-steps",
            shiny::tags$li(shiny::tags$strong("Upload an expression dataset. "),
              "Select one or more CSV, TSV or TXT files. Rows must be genes, columns must be samples and the first column must contain gene identifiers."),
            shiny::tags$li(shiny::tags$strong("Select the GEO dataset. "),
              "Choose the input file in the dropdown. If sample names include GSE accessions, select the target GSE or choose all samples in that file."),
            shiny::tags$li(shiny::tags$strong("Review Control and Disease groups. "),
              "Open Sample Review to confirm the counts and classifications. You may upload a separate Sample/Group metadata file or adjust the group-name patterns. Do not proceed with unknown or incorrectly assigned samples."),
            shiny::tags$li(shiny::tags$strong("Verify the expression scale. "),
              "Use processed, normalized, log2-scale microarray expression values. Confirm the scale using the checkbox. Raw count data are not suitable for this workflow."),
            shiny::tags$li(shiny::tags$strong("Set your DEG criteria. "),
              "Enter an absolute log2 fold-change cutoff (default 1.0) and BH-adjusted p-value cutoff (default 0.05). Optionally choose the numbers of heatmap genes, PCA genes and volcano labels."),
            shiny::tags$li(shiny::tags$strong("Choose figure outputs. "),
              "Select PNG, SVG and/or TIFF and the raster DPI. Charts have bold titles and axis labels for publication readability."),
            shiny::tags$li(shiny::tags$strong("Run differential-expression analysis. "),
              "Click Run DEG Analysis. The contrast is Disease minus Control, fitted using limma linear models and empirical-Bayes statistics; p-values are BH-adjusted."),
            shiny::tags$li(shiny::tags$strong("Review all results. "),
              "Inspect the Results, PCA, Volcano and Heatmap tabs. Positive log2FC means higher expression in Disease; negative log2FC means lower expression."),
            shiny::tags$li(shiny::tags$strong("Choose the output directory and save. "),
              "After the completion prompt, choose a local Windows folder. The app creates a run-specific subfolder containing DEG CSV tables, selected figures and the analysis summary.")
          ),
          shiny::tags$h4(class="guide-heading", "Important scientific notes"),
          shiny::tags$p(class="guide-intro",
            "This studio does not perform normalization or batch correction. Check biological comparability, confounding and the study design before interpreting the results. No significant genes may be found at stringent cutoffs."))),
        shiny::tabPanel("About", shiny::div(class="card",
          shiny::div(class="about-banner",
            shiny::tags$p(class="about-name", "Md Fazlay Rabbi"),
            shiny::tags$p("Research Fellow | Crop Protection & Toxicology Research Laboratory"),
            shiny::tags$p("Department of Zoology, University of Rajshahi, Rajshahi-6205, Bangladesh")),
          shiny::tags$table(class="about-table",
            shiny::tags$tr(shiny::tags$th("Mobile"),
              shiny::tags$td(shiny::tags$a(href="tel:+8801700823135", "+880 1700-823135"))),
            shiny::tags$tr(shiny::tags$th("Email"),
              shiny::tags$td(shiny::tags$a(href="mailto:rabbi.ms244@gmail.com", "rabbi.ms244@gmail.com"))),
            shiny::tags$tr(shiny::tags$th("Google Scholar"),
              shiny::tags$td(shiny::tags$a(href="https://scholar.google.com/citations?user=pfOKJG0AAAAJ&hl=en&oi=sra",
                target="_blank", rel="noopener noreferrer", "View Google Scholar profile"))),
            shiny::tags$tr(shiny::tags$th("ResearchGate"),
              shiny::tags$td(shiny::tags$a(href="https://www.researchgate.net/profile/Md-Rabbi-51",
                target="_blank", rel="noopener noreferrer", "View ResearchGate profile")))
          ),
          shiny::tags$p(class="about-note",
            "DEG Identification & limma Analysis Studio — developed for research use. Source: the user-provided DEG analysis workflow.")))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- shiny::reactiveValues(result=NULL, message="Upload an expression dataset to begin.", saved=NULL,
                              used_file=NULL, used_subset=NULL, ignored=NULL, used_name=NULL,
                              used_fc=NULL,used_fdr=NULL,used_heatmap=NULL,used_pca=NULL,used_labels=NULL)

  shiny::observeEvent(input$expression_files, {
    files <- input$expression_files
    if (is.null(files) || nrow(files) == 0L) return()
    shiny::updateSelectInput(session, "which_file", choices=stats::setNames(files$name,files$name),
                            selected=files$name[1L])
    rv$result <- NULL
    rv$message <- "Choose a dataset and review sample groups."
  })

  chosen_file <- shiny::reactive({
    shiny::req(input$expression_files, input$which_file)
    f <- input$expression_files
    index <- match(input$which_file, f$name)
    shiny::req(!is.na(index))
    read_expression(f$datapath[index], f$name[index])
  })

  shiny::observeEvent(chosen_file(), {
    obj <- chosen_file()
    ids <- unique(obj$sample_ids[!is.na(obj$sample_ids)])
    choices <- c("All samples in selected file" = "ALL")
    if (length(ids)) choices <- c(choices, stats::setNames(ids, ids))
    shiny::updateSelectInput(session,"gse_subset",choices=choices,
                             selected=if(length(ids)==1L) ids[1L] else "ALL")
    rv$result <- NULL
  })

  optional_metadata <- shiny::reactive({
    if (is.null(input$sample_metadata)) return(NULL)
    read_groups(input$sample_metadata$datapath, input$sample_metadata$name)
  })

  reviewed <- shiny::reactive({
    obj <- chosen_file()
    shiny::req(input$gse_subset)
    included <- if (identical(input$gse_subset,"ALL")) rep(TRUE,ncol(obj$mat))
                else obj$sample_ids == input$gse_subset & !is.na(obj$sample_ids)
    included[is.na(included)] <- FALSE
    if (!any(included)) stop("Selected GEO accession has no matching sample columns.")
    mat <- obj$mat[, included, drop=FALSE]
    sample_ids <- colnames(mat)
    groups <- determine_groups(sample_ids, optional_metadata(), input$control_regex, input$disease_regex)
    group_df <- data.frame(Sample=sample_ids, Group=ifelse(is.na(groups),"Unassigned",groups),
                           stringsAsFactors=FALSE)
    list(mat=mat, groups=groups, group_df=group_df, input=obj,
         all_n=ncol(obj$mat), fileset=input$which_file,
         gse=input$gse_subset, name=if(identical(input$gse_subset,"ALL")) {
           get_gse(input$which_file)[1L] %||% "Combined_Expression"
         } else input$gse_subset)
  })

  output$info_table <- shiny::renderTable({
    dat <- reviewed()
    data.frame(Field=c("Selected CSV", "Selected GEO subset", "Genes (before QC)",
                       "Total samples in file", "Samples included", "Duplicate gene IDs", "Ignored text columns"),
       Value=c(dat$fileset, dat$gse, nrow(dat$mat), dat$all_n,
               ncol(dat$mat), dat$input$duplicate_genes,
               if(length(dat$input$ignored)) paste(dat$input$ignored,collapse=", ") else "None"))
  }, striped=TRUE)

  output$counts_table <- shiny::renderTable({
    d <- reviewed()
    tab <- table(factor(d$group_df$Group,levels=c("Control","Disease","Unassigned")))
    data.frame(Group=names(tab),Samples=as.integer(tab))
  }, striped=TRUE)
  output$sample_table <- shiny::renderTable({
    head(reviewed()$group_df,40)
  }, striped=TRUE)

  output$status_text <- shiny::renderText(rv$message)

  shiny::observeEvent(input$run, {
    tryCatch({
      if (!is.finite(input$fc_cutoff) || input$fc_cutoff < 0) stop("log2FC cutoff must be >= 0.")
      if (!is.finite(input$fdr_cutoff) || input$fdr_cutoff <= 0 || input$fdr_cutoff >= 1) {
        stop("Adjusted p-value cutoff must be strictly between 0 and 1.")
      }
      if (!isTRUE(input$confirmed_scale)) {
        stop("Please confirm that your matrix is processed, normalized and log2-scale. This app does not transform raw counts.")
      }
      dat <- reviewed()
      mat <- dat$mat
      groups <- dat$groups
      unassigned <- is.na(groups)
      if (any(unassigned) && !isTRUE(input$exclude_unknown)) {
        stop(sum(unassigned), " sample(s) are unassigned. Add Sample/Group metadata or adjust patterns.")
      }
      if (any(unassigned)) {
        mat <- mat[, !unassigned, drop=FALSE]
        groups <- groups[!unassigned]
      }
      group <- factor(groups, levels=c("Control","Disease"))
      if (any(table(group) < 2L)) {
        stop("Need at least 2 Control and 2 Disease samples after group review.")
      }
      if (!length(input$formats)) stop("Choose at least one output figure format.")
      rv$message <- "Running limma differential-expression analysis..."
      shiny::withProgress(message="DEG analysis", value=0, {
        shiny::incProgress(.15, detail="Preparing and validating expression data")
        result <- run_deg_analysis(mat, group, input$fc_cutoff, input$fdr_cutoff,
                     as.integer(input$top_heatmap), as.integer(input$pca_genes),
                     as.integer(input$volcano_labels))
        shiny::incProgress(.65, detail="limma, results and visualizations completed")
        rv$result <- result
        rv$used_file <- dat$fileset
        rv$used_subset <- dat$gse
        rv$used_name <- dat$name
        rv$ignored <- dat$input$ignored
        rv$used_fc <- input$fc_cutoff
        rv$used_fdr <- input$fdr_cutoff
        rv$used_heatmap <- input$top_heatmap
        rv$used_pca <- input$pca_genes
        rv$used_labels <- input$volcano_labels
        shiny::incProgress(.20, detail="Ready to save")
      })
      rv$message <- sprintf("Analysis completed: %d Up + %d Down DEGs. Choose folder to save.",
                            nrow(rv$result$up), nrow(rv$result$down))
      shiny::updateTabsetPanel(session,"tabs",selected="Results")
      shiny::showModal(shiny::modalDialog(
         title="DEG analysis complete",
         shiny::tags$p(sprintf("Upregulated: %d  |  Downregulated: %d",
             nrow(rv$result$up), nrow(rv$result$down))),
         shiny::tags$p("Would you like to choose the output directory now?"),
         footer=shiny::tagList(shiny::modalButton("Later"),
                   shiny::actionButton("save_after", "Choose folder & save",class="btn-success")),
         easyClose=FALSE))
    }, error=function(e) {
      rv$message <- paste("Error:", conditionMessage(e))
      shiny::showModal(shiny::modalDialog(title="Unable to complete analysis",
        shiny::tags$pre(style="white-space:pre-wrap;color:#233e57;", conditionMessage(e)),
        easyClose=TRUE))
    })
  })

  output$results_summary <- shiny::renderUI({
    shiny::req(rv$result)
    r <- rv$result
    shiny::tagList(shiny::span(class="metric",paste("Genes:",nrow(r$expr))),
      shiny::span(class="metric",paste("Upregulated:",nrow(r$up))),
      shiny::span(class="metric",paste("Downregulated:",nrow(r$down))),
      shiny::span(class="metric",paste("Control:",sum(r$group=="Control"))),
      shiny::span(class="metric",paste("Disease:",sum(r$group=="Disease"))))
  })
  output$deg_preview <- shiny::renderTable({
    shiny::req(rv$result)
    utils::head(rv$result$deg, 15L)
  },digits=4,striped=TRUE)
  output$pca_view <- shiny::renderPlot({shiny::req(rv$result);print(rv$result$pca)},res=110)
  output$volcano_view <- shiny::renderPlot({shiny::req(rv$result);print(rv$result$volcano)},res=110)
  output$heatmap_view <- shiny::renderPlot({
    shiny::req(rv$result)
    if (inherits(rv$result$heatmap,"ggplot")) print(rv$result$heatmap)
    else {grid::grid.newpage();grid::grid.draw(rv$result$heatmap$gtable)}
  },res=110)

  save_to_chosen_folder <- function() {
    if (is.null(rv$result)) {
      shiny::showNotification("Run the DEG analysis before saving.",type="warning"); return()
    }
    if (.Platform$OS.type == "windows") {
      folder <- utils::choose.dir(caption="Choose folder for DEG results")
    } else {
      folder <- NA_character_
    }
    if (is.na(folder) || !nzchar(folder)) {
      # Cross-platform fallback: an explicit directory path in a modal.
      if (.Platform$OS.type != "windows") {
        shiny::showModal(shiny::modalDialog(title="Enter output directory",
             shiny::textInput("manual_folder","Existing local folder path",getwd()),
             footer=shiny::tagList(shiny::modalButton("Cancel"),
                                   shiny::actionButton("save_manual","Save results",class="btn-success"))))
      } else {
        shiny::showNotification("No output folder selected; results remain in the app.",type="message")
      }
      return()
    }
    tryCatch({
      path <- save_analysis(rv$result, folder, rv$used_name %||% "GSE_Selected",
           rv$used_file, rv$used_subset, rv$used_fc, rv$used_fdr,
           rv$used_heatmap, rv$used_pca, rv$used_labels,
           input$formats,as.integer(input$dpi),rv$ignored)
      rv$saved <- path
      rv$message <- paste("Results saved to:",path)
      shiny::showModal(shiny::modalDialog(title="Results saved successfully",
          shiny::tags$p("Output directory:"),shiny::tags$code(path),easyClose=TRUE))
    },error=function(e){
      rv$message <- paste("Saving failed:",conditionMessage(e))
      shiny::showModal(shiny::modalDialog(title="Saving failed",conditionMessage(e)))
    })
  }
  shiny::observeEvent(input$save, {save_to_chosen_folder()})
  shiny::observeEvent(input$save_after, {shiny::removeModal();save_to_chosen_folder()})
  shiny::observeEvent(input$save_manual, {
    shiny::req(rv$result, input$manual_folder)
    folder <- path.expand(input$manual_folder)
    if (!dir.exists(folder)) {
      shiny::showNotification("That folder does not exist.",type="error"); return()
    }
    tryCatch({
      path <- save_analysis(rv$result,folder,rv$used_name %||% "GSE_Selected",
         rv$used_file, rv$used_subset, rv$used_fc,rv$used_fdr,rv$used_heatmap,
         rv$used_pca,rv$used_labels,input$formats,as.integer(input$dpi),rv$ignored)
      rv$saved <- path
      rv$message <- paste("Results saved to:",path)
      shiny::removeModal()
      shiny::showNotification(paste("Saved:",path),type="message",duration=12)
    },error=function(e) shiny::showNotification(conditionMessage(e),type="error",duration=NULL))
  })
}

# Source this complete .R file in RStudio; a LOCAL browser window opens.
# The user does not need to modify this code between analyses.
shiny::runApp(shiny::shinyApp(ui = ui, server = server), launch.browser = TRUE,
              host = "127.0.0.1")
