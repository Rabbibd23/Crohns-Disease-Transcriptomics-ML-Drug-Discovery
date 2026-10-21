# ============================================================================
# GEO Independent Gene Validation Studio | RStudio / local Shiny application
# Developed by Md Fazlay Rabbi | Department of Zoology, University of Rajshahi
# Based on: STAT3_6_Dataset_Independent_Validation_R_code.R
# Version 1.3 | Exact GSE labels, gene-first processing, automatic groups and downloads
# ============================================================================

options(shiny.maxRequestSize = 700 * 1024^2)

# Check first, install only missing packages; no reinstall on later runs.
required_packages <- c("shiny", "shinyFiles", "ggplot2", "dplyr", "tidyr", "data.table", "zip")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
still_missing <- required_packages[!vapply(required_packages, requireNamespace,
                                             quietly = TRUE, FUN.VALUE = logical(1))]
if (length(still_missing)) stop("Packages unavailable: ", paste(still_missing, collapse = ", "))

suppressPackageStartupMessages({
  library(shiny)
  library(shinyFiles)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, alternative) if (is.null(x) || length(x) == 0L) alternative else x

# Keep the original six-panel order when those datasets are selected.
PREFERRED_GSE_ORDER <- c("GSE66407", "GSE179285", "GSE95095",
                         "GSE75214", "GSE102133", "GSE57945")
COLORS <- c(Control = "#4DBBD5", Disease = "#E64B35")

clean_gene <- function(gene) toupper(trimws(as.character(gene)))


# User-friendly fixed group detection. No regular expressions to edit in the GUI.
# read.csv(check.names=TRUE) may convert names such as "control 1" to "control.1".
# Accepted examples: control, Control.1, control_2, control-3, disease, disease.1.
identify_sample_groups <- function(column_names) {
  control_cols <- grep("^control($|[._-]|[0-9])", column_names,
                       value = TRUE, ignore.case = TRUE, perl = TRUE)
  disease_cols <- grep("^disease($|[._-]|[0-9])", column_names,
                       value = TRUE, ignore.case = TRUE, perl = TRUE)
  list(control = control_cols, disease = disease_cols)
}

# Inspect only the CSV header/first row for responsive dataset previews.
# Full-expression data are read only when the user starts the analysis.
preview_sample_groups <- function(path) {
  tryCatch({
    sample_row <- utils::read.csv(path, nrows = 1L, check.names = TRUE,
                                  stringsAsFactors = FALSE)
    labels <- identify_sample_groups(names(sample_row))
    nc <- length(labels$control)
    nd <- length(labels$disease)
    list(n_control = nc, n_disease = nd,
         status = if (nc > 0L && nd > 0L) "Ready" else "Missing Control or Disease columns",
         control_columns = labels$control, disease_columns = labels$disease)
  }, error = function(e) list(n_control = 0L, n_disease = 0L,
                               status = paste("Could not read CSV:", conditionMessage(e)),
                               control_columns = character(), disease_columns = character()))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

# Dataset identity is taken only from the selected filename.
# Never relabel GSE102134 as GSE102133 without independently verified metadata.
get_dataset_label <- function(file_name) {
  basename_ <- basename(file_name)
  hit <- regmatches(basename_, regexpr("GSE[0-9]+", basename_, ignore.case = TRUE))
  if (!length(hit) || !nzchar(hit)) hit <- tools::file_path_sans_ext(basename_)
  toupper(hit)
}

# FAST gene-first workflow. Only requested target-gene rows undergo numeric
# conversion and duplicate averaging. This reduces the major bottleneck in
# large processed GEO expression matrices while retaining the SAME values
# for selected genes as whole-matrix duplicate-gene averaging.
# Users may separately opt in to full gene-averaged matrix export.
read_and_prepare_dataset <- function(path, original_name, dataset_id, target_genes,
                                     save_averaged = FALSE, averaged_dir = NULL) {
  df <- data.table::fread(path, data.table = FALSE, check.names = TRUE,
                          showProgress = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) < 1L || ncol(df) < 3L) {
    stop("Expression matrix is empty or lacks sample columns: ", original_name)
  }
  possible_gene_cols <- c("gene_symbol", "Gene_symbol", "GENE_SYMBOL",
                          "official_gene_symbol", "Gene", "gene", "Symbol", "symbol")
  matched <- possible_gene_cols[possible_gene_cols %in% names(df)]
  gene_col <- if (length(matched)) matched[[1L]] else names(df)[[1L]]
  groups <- identify_sample_groups(names(df))
  control_cols <- groups$control
  disease_cols <- groups$disease
  if (!length(control_cols) || !length(disease_cols)) {
    stop("Control and Disease columns were not both found in ", original_name,
         ". Expected columns such as control, control.1, disease and disease.1.")
  }
  sample_cols <- c(control_cols, disease_cols)
  df[[gene_col]] <- clean_gene(df[[gene_col]])
  df <- df[!is.na(df[[gene_col]]) & nzchar(df[[gene_col]]) &
             df[[gene_col]] != "NA", , drop = FALSE]
  if (!nrow(df)) stop("No valid gene symbols found in ", original_name)
  original_rows <- nrow(df)
  unique_genes <- length(unique(df[[gene_col]]))

  convert_sample_columns <- function(dat) {
    for (col in sample_cols) {
      raw <- as.character(dat[[col]])
      converted <- suppressWarnings(as.numeric(raw))
      bad <- sum(!is.na(raw) & nzchar(trimws(raw)) & is.na(converted))
      if (bad) warning(dataset_id, ": ", bad, " nonnumeric entries in ", col,
                       " were converted to NA")
      dat[[col]] <- converted
    }
    dat
  }
  average_duplicates <- function(dat) {
    if (!nrow(dat)) {
      result <- data.frame(gene_symbol = character(), stringsAsFactors = FALSE)
      for (col in sample_cols) result[[col]] <- numeric()
      return(result)
    }
    dat <- convert_sample_columns(dat)
    dat %>%
      dplyr::group_by(gene_symbol = .data[[gene_col]]) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols),
        ~ if (all(is.na(.x))) NA_real_ else base::mean(.x, na.rm = TRUE)),
        .groups = "drop")
  }
  target_genes <- unique(clean_gene(target_genes))
  requested <- df[df[[gene_col]] %in% target_genes,
                  c(gene_col, sample_cols), drop = FALSE]
  averaged <- average_duplicates(requested)

  # The optional full export performs all-gene averaging, as in the
  # original script. Off by default for efficient validation runs.
  if (isTRUE(save_averaged)) {
    if (is.null(averaged_dir)) stop("Specify a folder for full averaged matrix export.")
    dir.create(averaged_dir, recursive = TRUE, showWarnings = FALSE)
    full <- average_duplicates(df[, c(gene_col, sample_cols), drop = FALSE])
    utils::write.csv(full,
      file.path(averaged_dir,
                paste0(safe_name(dataset_id), "_unique_gene_averaged_expression.csv")),
      row.names = FALSE)
    rm(full)
  }
  rm(df)
  list(averaged = averaged, sample_cols = sample_cols,
       control_cols = control_cols, disease_cols = disease_cols,
       genes_before = original_rows, genes_after = unique_genes,
       dataset_id = dataset_id, input_name = original_name)
}

# Build long-format expression for the requested gene; no gene is imputed.
extract_gene <- function(prepared, gene) {
  gene_row <- prepared$averaged[prepared$averaged$gene_symbol == gene, , drop = FALSE]
  if (nrow(gene_row) == 0L) return(NULL)
  if (nrow(gene_row) != 1L) stop("Gene row is not unique after averaging: ", gene)
  cols <- prepared$sample_cols
  values <- as.numeric(unlist(gene_row[1L, cols, drop = FALSE], use.names = FALSE))
  tbl <- data.frame(Gene = gene, Dataset = prepared$dataset_id,
                    Sample = cols,
                    Group = ifelse(cols %in% prepared$control_cols, "Control", "Disease"),
                    Expression = values, stringsAsFactors = FALSE)
  tbl <- tbl[is.finite(tbl$Expression), , drop = FALSE]
  if (!nrow(tbl) || !all(c("Control", "Disease") %in% tbl$Group)) return(NULL)
  tbl
}

# The same two-sided Wilcoxon rank-sum test as the uploaded R script.
wilcox_p <- function(values, groups) {
  if (!all(c("Control", "Disease") %in% as.character(groups))) return(NA_real_)
  tryCatch(suppressWarnings(stats::wilcox.test(values ~ groups,
     alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value),
     error = function(e) NA_real_)
}

make_statistics <- function(long_df, fdr_scope = "per_gene") {
  stat <- long_df %>%
    dplyr::group_by(Gene, Dataset) %>%
    dplyr::summarise(
      n_control = sum(Group == "Control"),
      n_disease = sum(Group == "Disease"),
      median_control = stats::median(Expression[Group == "Control"], na.rm = TRUE),
      median_disease = stats::median(Expression[Group == "Disease"], na.rm = TRUE),
      mean_control = base::mean(Expression[Group == "Control"], na.rm = TRUE),
      mean_disease = base::mean(Expression[Group == "Disease"], na.rm = TRUE),
      median_difference = median_disease - median_control,
      p_value = wilcox_p(Expression, Group), .groups = "drop")
  if (identical(fdr_scope, "global")) {
    stat$p_adj_BH <- stats::p.adjust(stat$p_value, method = "BH")
    correction_description <- "BH correction across all gene-by-dataset tests"
  } else {
    stat <- stat %>% dplyr::group_by(Gene) %>%
      dplyr::mutate(p_adj_BH = stats::p.adjust(p_value, method = "BH")) %>%
      dplyr::ungroup()
    correction_description <- "BH correction separately per gene across selected datasets"
  }
  stat <- stat %>% dplyr::mutate(
    direction = dplyr::case_when(
      median_difference > 0 ~ "Upregulated in disease",
      median_difference < 0 ~ "Downregulated in disease", TRUE ~ "No median change"),
    significance = dplyr::case_when(is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***", p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*", TRUE ~ "ns"),
    fdr_label = dplyr::case_when(is.na(p_adj_BH) ~ "FDR = NA",
      p_adj_BH < 0.001 ~ "FDR < 0.001",
      TRUE ~ paste0("FDR = ", signif(p_adj_BH, 3))),
    final_label = paste0(significance, "\n", fdr_label))
  attr(stat, "correction_description") <- correction_description
  stat
}

# Original colors, original point/box layout, annotation placement and text theme.
make_original_validation_plot <- function(df, stat, gene, dataset_order) {
  df <- df[df$Gene == gene, , drop = FALSE]
  stat <- stat[stat$Gene == gene, , drop = FALSE]
  present <- intersect(dataset_order, unique(as.character(df$Dataset)))
  df$Dataset <- factor(as.character(df$Dataset), levels = present)
  df$Group <- factor(df$Group, levels = c("Control", "Disease"))
  stat$Dataset <- factor(as.character(stat$Dataset), levels = present)
  position <- df %>% dplyr::group_by(Dataset) %>%
    dplyr::summarise(y_max = max(Expression, na.rm = TRUE),
      y_min = min(Expression, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(y_range = y_max - y_min,
      y_range = ifelse(y_range == 0, abs(y_max) + 1, y_range),
      y_pos = y_max + 0.18 * y_range, x_pos = 1.5) %>%
    dplyr::left_join(stat[, c("Dataset", "final_label")], by = "Dataset")
  ncol_plot <- min(3L, length(present))
  ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Expression, fill = Group)) +
    ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.92,
                          linewidth = 0.65, color = "black") +
    ggplot2::geom_jitter(ggplot2::aes(color = Group), width = 0.12,
                         size = 1.65, alpha = 0.65, show.legend = FALSE) +
    ggplot2::facet_wrap(~ Dataset, ncol = ncol_plot, scales = "free_y") +
    ggplot2::geom_label(data = position,
      ggplot2::aes(x = x_pos, y = y_pos, label = final_label),
      inherit.aes = FALSE, fontface = "bold", size = 4.7, lineheight = 0.88,
      label.size = 0, label.padding = grid::unit(0.18, "lines"),
      fill = "white", color = "black") +
    ggplot2::scale_fill_manual(values = COLORS) +
    ggplot2::scale_color_manual(values = COLORS) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.35))) +
    ggplot2::labs(title = paste0(gene, " expression validation across independent GEO datasets"),
                  x = "", y = paste0(gene, " expression")) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 20,
                                        hjust = 0.5, color = "black"),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      axis.text.x = ggplot2::element_text(face = "bold", size = 16, color = "black"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 13, color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black", linewidth = 0.8),
      strip.text = ggplot2::element_text(face = "bold", size = 15, color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(face = "bold", size = 16, color = "black"),
      legend.key.width = grid::unit(1.20, "cm"),
      legend.key.height = grid::unit(0.82, "cm"),
      legend.spacing.x = grid::unit(0.32, "cm"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(alpha = 1, linewidth = 0.8),
      nrow = 1,
      label.theme = ggplot2::element_text(face = "bold", size = 16, color = "black")
    ))
}

# Optional complementary plot: only significant direction, not pooled effects.
make_multigene_overview <- function(stat, gene_order, dataset_order) {
  grid <- expand.grid(Gene = gene_order, Dataset = dataset_order,
                      stringsAsFactors = FALSE)
  dat <- dplyr::left_join(grid, stat[, c("Gene", "Dataset", "p_adj_BH", "median_difference")],
                          by = c("Gene", "Dataset"))
  dat$signal <- ifelse(is.na(dat$p_adj_BH), "Missing",
                      ifelse(dat$p_adj_BH >= 0.05, "Not significant",
                             ifelse(dat$median_difference > 0, "Up in disease",
                                    ifelse(dat$median_difference < 0, "Down in disease", "Not significant"))))
  dat$label <- ifelse(is.na(dat$p_adj_BH), "NA",
                      ifelse(dat$p_adj_BH < 0.001, "***",
                             ifelse(dat$p_adj_BH < 0.01, "**",
                                    ifelse(dat$p_adj_BH < 0.05, "*", "ns"))))
  dat$Gene <- factor(dat$Gene, levels = rev(gene_order))
  dat$Dataset <- factor(dat$Dataset, levels = dataset_order)
  ggplot2::ggplot(dat, ggplot2::aes(x = Dataset, y = Gene, fill = signal)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), fontface = "bold", size = 4) +
    ggplot2::scale_fill_manual(values = c("Up in disease" = "#E64B35",
      "Down in disease" = "#4DBBD5", "Not significant" = "#DDE5EC",
      "Missing" = "#FFFFFF"), drop = FALSE) +
    ggplot2::labs(title = "Independent validation overview",
      subtitle = "Direction is shown only for BH-significant comparisons; * FDR < 0.05, ** < 0.01, *** < 0.001",
      x = "GEO dataset", y = "Validated genes", fill = "Interpretation") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", color = "#152C42"),
      legend.text = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank())
}

save_validation_plot <- function(plot, file_prefix, width, height, dpi,
                                 svg = TRUE) {
  ggplot2::ggsave(paste0(file_prefix, ".png"), plot = plot,
    width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
  ggplot2::ggsave(paste0(file_prefix, ".tiff"), plot = plot,
    width = width, height = height, units = "in", dpi = dpi,
    compression = "lzw", limitsize = FALSE)
  ggplot2::ggsave(paste0(file_prefix, ".pdf"), plot = plot,
    width = width, height = height, units = "in", limitsize = FALSE)
  if (isTRUE(svg)) {
    # Optional SVG device may be unsupported in minimal R installations.
    tryCatch(ggplot2::ggsave(paste0(file_prefix, ".svg"), plot = plot,
      device = grDevices::svg, width = width, height = height,
      units = "in", limitsize = FALSE),
      error = function(e) warning("SVG could not be saved: ", conditionMessage(e)))
  }
}

# Single entry point used by the GUI; calculates first into a temporary folder.
perform_validation <- function(info, selected_rows, genes, settings, progress = NULL) {
  stopifnot(length(genes) >= 1L, length(selected_rows) >= 1L)
  files <- info[selected_rows, , drop = FALSE]
  if (anyDuplicated(files$dataset_id)) {
    stop("Two files resolve to the same GSE accession. Select only one per dataset.")
  }
  run_dir <- tempfile("GEO_Validation_Workspace_")
  dir.create(run_dir, recursive = TRUE)
  figures_dir <- file.path(run_dir, "Figures")
  data_dir <- file.path(run_dir, "Expression_Data")
  dir.create(figures_dir)
  dir.create(data_dir)

  dataset_order <- c(PREFERRED_GSE_ORDER[PREFERRED_GSE_ORDER %in% files$dataset_id],
                     setdiff(files$dataset_id, PREFERRED_GSE_ORDER))
  missing <- list()
  combined <- list()
  dataset_qc <- list()
  for (j in seq_len(nrow(files))) {
    if (is.function(progress)) progress(j, nrow(files), files$dataset_id[[j]], "start")
    x <- read_and_prepare_dataset(files$datapath[[j]], files$filename[[j]],
      files$dataset_id[[j]], target_genes = genes,
      save_averaged = settings$save_averaged, averaged_dir = data_dir)
    if (is.function(progress)) progress(j, nrow(files), files$dataset_id[[j]], "done")
    dataset_qc[[j]] <- data.frame(Dataset = x$dataset_id,
      File = x$input_name, Control_samples = length(x$control_cols),
      Disease_samples = length(x$disease_cols), Gene_rows_before = x$genes_before,
      Unique_genes_in_input = x$genes_after)
    for (gene in genes) {
      val <- extract_gene(x, gene)
      if (is.null(val)) {
        missing[[length(missing) + 1L]] <- data.frame(Gene = gene,
          Dataset = x$dataset_id, Reason = "Absent gene or no usable control/disease values")
      } else {
        combined[[length(combined) + 1L]] <- val
      }
    }
    rm(x)
  }
  if (!length(combined)) stop("None of the requested genes could be validated in the selected files.")
  long_df <- dplyr::bind_rows(combined)
  long_df$Group <- factor(long_df$Group, levels = c("Control", "Disease"))
  long_df$Gene <- factor(long_df$Gene, levels = genes)
  long_df$Dataset <- factor(long_df$Dataset, levels = dataset_order)
  stats <- make_statistics(long_df, settings$fdr_scope)
  stats$Gene <- as.character(stats$Gene)
  stats$Dataset <- as.character(stats$Dataset)
  qc <- dplyr::bind_rows(dataset_qc)
  missing_df <- if (length(missing)) dplyr::bind_rows(missing) else
    data.frame(Gene = character(), Dataset = character(), Reason = character())
  utils::write.csv(long_df, file.path(run_dir, "Combined_Gene_Expression_Long.csv"),
                   row.names = FALSE)
  utils::write.csv(stats, file.path(run_dir, "Wilcoxon_BH_Statistics.csv"),
                   row.names = FALSE)
  utils::write.csv(qc, file.path(run_dir, "Dataset_QC_Summary.csv"), row.names = FALSE)
  utils::write.csv(missing_df, file.path(run_dir, "Missing_Gene_Dataset_Comparisons.csv"),
                   row.names = FALSE)
  plots <- list()
  genes_found <- intersect(genes, unique(as.character(long_df$Gene)))
  for (gene in genes_found) {
    p <- make_original_validation_plot(long_df, stats, gene, dataset_order)
    plots[[gene]] <- p
    present_count <- length(unique(as.character(long_df$Dataset[long_df$Gene == gene])))
    columns <- min(3L, present_count)
    rows <- ceiling(present_count / columns)
    width <- if (isTRUE(settings$auto_dimensions)) max(6.2, 4.15 * columns) else settings$width
    height <- if (isTRUE(settings$auto_dimensions)) max(4.8, 4.5 * rows) else settings$height
    save_validation_plot(p, file.path(figures_dir,
      paste0(safe_name(gene), "_Independent_Validation")), width, height,
      settings$dpi, svg = settings$save_svg)
  }
  overview <- NULL
  if (length(genes_found) > 1L && isTRUE(settings$overview)) {
    overview <- make_multigene_overview(stats, genes, dataset_order)
    ggplot2::ggsave(file.path(figures_dir, "Multi_Gene_Validation_Overview.png"),
      overview, width = max(8, 1.5 * length(dataset_order) + 2),
      height = max(5, 0.65 * length(genes) + 3), dpi = settings$dpi,
      units = "in", limitsize = FALSE)
    ggplot2::ggsave(file.path(figures_dir, "Multi_Gene_Validation_Overview.pdf"),
      overview, width = max(8, 1.5 * length(dataset_order) + 2),
      height = max(5, 0.65 * length(genes) + 3), units = "in", limitsize = FALSE)
  }
  method_lines <- c(
    "Independent GEO Expression Validation — Methods Summary",
    paste("Genes requested:", paste(genes, collapse = ", ")),
    paste("Dataset accessions:", paste(files$dataset_id, collapse = ", ")),
    "Expression input: processed GEO expression matrices with Control/Disease columns.",
    "Only requested gene-symbol rows were selected prior to sample-wise duplicate averaging.",
    if (isTRUE(settings$save_averaged)) "Full duplicate-averaged expression matrices were additionally exported." else "Full-matrix averaging was not performed.",
    "Dataset identifiers were taken from input filenames; no GSE accession was automatically remapped.",
    "Two-sided Wilcoxon rank-sum tests with continuity correction were performed within each dataset.",
    paste("Multiple-testing correction:",
          if (identical(settings$fdr_scope, "global"))
            "BH correction across all gene-by-dataset tests" else
            "BH correction separately per gene across selected datasets"),
    "For each dataset, means, medians, median differences and sample counts were reported.",
    "Expression values were not log-transformed, renormalized or batch-corrected by this tool.",
    "Samples were treated as independent; paired/repeated-measure designs require a different analysis.",
    "Missing gene-dataset combinations were reported, not imputed.",
    paste("Figure export:", settings$dpi, "DPI raster; PNG, TIFF, PDF",
          if (isTRUE(settings$save_svg)) "and SVG (where supported)" else ""),
    paste("BH correction scope:", settings$fdr_scope))
  writeLines(method_lines, file.path(run_dir, "Methods.txt"), useBytes = TRUE)
  result_lines <- c("Independent GEO Expression Validation — Results Summary",
    paste("Validated genes:", paste(genes_found, collapse = ", ")),
    paste("Datasets processed:", length(dataset_order)),
    paste("Gene-dataset tests available:", nrow(stats)),
    paste("Comparisons FDR < 0.05:", sum(stats$p_adj_BH < 0.05, na.rm = TRUE)),
    paste("Missing comparisons:", nrow(missing_df)), "", "Per gene and dataset:")
  for (i in seq_len(nrow(stats))) {
    s <- stats[i, ]
    result_lines <- c(result_lines, sprintf(
      "%s | %s | n Control=%d, Disease=%d | median difference=%.4g | p=%.4g | BH-FDR=%.4g | %s",
      s$Gene, s$Dataset, s$n_control, s$n_disease,
      s$median_difference, s$p_value, s$p_adj_BH, s$direction))
  }
  writeLines(result_lines, file.path(run_dir, "Results.txt"), useBytes = TRUE)
  list(run_dir = run_dir, stats = stats, long_df = long_df,
       qc = qc, missing = missing_df, plots = plots, overview = overview,
       genes_found = genes_found, dataset_order = dataset_order, settings = settings)
}

# Safe save after analysis: never overwrite a previous run folder.
save_run_to_directory <- function(run_dir, target_parent) {
  if (!dir.exists(target_parent)) stop("Selected output directory does not exist.")
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  target <- file.path(target_parent, paste0("GEO_Gene_Validation_", stamp))
  index <- 2L
  while (file.exists(target)) {
    target <- file.path(target_parent, paste0("GEO_Gene_Validation_", stamp, "_", index))
    index <- index + 1L
  }
  if (!dir.create(target, recursive = TRUE)) stop("Cannot create output folder: ", target)
  entries <- list.files(run_dir, all.files = FALSE, full.names = TRUE)
  copied <- file.copy(entries, target, recursive = TRUE, overwrite = FALSE)
  if (!all(copied)) stop("Some results could not be copied to: ", target)
  target
}

# User-visible local file system roots for choosing output folder AFTER analysis.
create_roots <- function() {
  roots <- c(Home = normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE))
  working <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (dir.exists(working)) roots <- c(roots, Working = working)
  if (.Platform$OS.type == "windows") {
    drive_letters <- LETTERS
    paths <- paste0(drive_letters, ":/")
    valid <- dir.exists(paths)
    if (any(valid)) roots <- c(roots, stats::setNames(paths[valid],
      paste0("Drive_", drive_letters[valid])))
  }
  roots[!duplicated(roots)]
}

ui <- shiny::fluidPage(
  tags$head(tags$meta(charset = "utf-8"),
    tags$style(HTML('
      body {background:#ECF2F4; color:#173247; font-family:"Segoe UI",Arial,sans-serif;}
      .hero {position:relative; background:linear-gradient(112deg,#123F51,#176F78 63%,#268A82);
        border-radius:13px; padding:24px 96px; color:white; text-align:center;
        box-shadow:0 7px 18px rgba(22,64,75,.12); margin:9px 0 20px; min-height:143px;}
      .frmark {position:absolute; left:20px; top:26px; width:77px; height:77px; border-radius:50%;
        background:white; border:4px solid #77E0D0; color:#155768; display:flex;
        align-items:center; justify-content:center; font-size:31px; font-weight:900;}
      .hero h1 {font-size:30px; font-weight:800; margin:5px 0 8px; color:white;}
      .hero p {margin:5px 0; color:#E7FEFA; font-size:13px;}
      .panelbox {background:white; border:1px solid #CFDFE5; border-radius:12px;
        padding:18px 18px 11px; margin-bottom:15px; box-shadow:0 3px 12px rgba(31,74,83,.04);}
      .section-title {color:#115568; font-weight:800; font-size:17px; margin:0 0 13px;}
      label, .control-label {font-weight:700; color:#1D4155;}
      input, .form-control, .selectize-input, .selectize-dropdown, select {
        background:#FFF !important; color:#162C39 !important; border-color:#A6C6D0 !important;}
      .btn-run {background:#07806C !important; border-color:#07806C !important;
        color:#FFFFFF !important; font-weight:800; padding:12px 17px; width:100%;}
      .btn-run:hover {background:#056F60 !important; color:#FFFFFF !important;}
      .auto-group {border:1px solid #A5DCCF; background:#E9F8F3;
        color:#124D47; border-radius:8px; padding:12px; font-size:13px; line-height:1.65;}
      .downloads-row {display:flex; flex-wrap:wrap; gap:9px; margin:12px 0;}
      .downloads-row a, .downloads-row button {font-weight:750;}
      .btn-download {background:#197B88 !important; border:1px solid #197B88 !important;
        color:white !important; font-weight:750; margin:4px 6px 4px 0;}
      .btn-download:hover {background:#0C626E !important; color:white !important;}
      .nav-tabs > li > a {color:#145567; font-weight:700;}
      .nav-tabs > li.active > a {color:#096759 !important; font-weight:800;}
      .table th {background:#E5F2F2;color:#183F4D;}
      .small-note {font-size:12px; color:#4A6D7A; line-height:1.55;}
      .modal-content {background:#FFF;color:#182F3D;border-radius:12px;}
      .modal-header,.modal-footer {border-color:#E0EBED;}
      .modal-body,.modal-body p {color:#193A4B;}
      #status_text {font-weight:700; color:#076C60;}
      .helpsteps li {margin-bottom:12px; line-height:1.55;}
      .helpsteps {font-size:14px;}
      .navbar {margin-bottom:0;}
      @media(max-width:700px){.hero{padding:19px 15px 19px 93px; text-align:left;}
       .hero h1{font-size:21px;} .frmark{width:60px;height:60px;font-size:25px;left:16px;}}
    '))),
  div(class = "hero",
    div(class = "frmark", "FR"),
    tags$h1("Independent GEO Gene Validation Studio"),
    tags$p("Gene-wise validation  |  One or more GEO datasets  |  Wilcoxon + BH-FDR  |  Publication figures"),
    tags$p("Developed by Md Fazlay Rabbi  |  Department of Zoology, University of Rajshahi, Bangladesh")),
  shiny::fluidRow(
    shiny::column(width = 4,
      div(class = "panelbox",
        div(class = "section-title", "1. Select genes"),
        shiny::numericInput("gene_count", "How many genes do you want to validate?",
          value = 1, min = 1, max = 50, step = 1),
        shiny::uiOutput("gene_inputs"),
        tags$p(class = "small-note", "Enter official gene symbols, e.g. STAT3, IL6, JAK2. A separate original-style plot will be generated for each gene.")),
      div(class = "panelbox",
        div(class = "section-title", "2. Select GEO expression files"),
        shiny::fileInput("expression_files", "Upload one or multiple CSV datasets",
          accept = ".csv", multiple = TRUE),
        shiny::checkboxGroupInput("datasets_to_use", "Datasets to analyze", choices = NULL),
        tags$p(class = "small-note", "Expression CSVs must contain a gene-symbol column and columns beginning with Control or Disease. No regex or manual grouping is needed.")),
      div(class = "panelbox",
        div(class = "section-title", "3. Statistical and plot settings"),
        div(class = "auto-group",
          tags$b("Control / Disease: automatically detected"),
          tags$div("Your CSV sample columns can be named control, control.1, control_2, disease, disease.1, etc."),
          tags$div("Check the Dataset Review tab to confirm how many samples were detected in each group.")),
        shiny::selectInput("fdr_scope", "BH multiple-testing scope",
          choices = c("Per gene across GEO datasets (original rule)" = "per_gene",
                      "All gene × dataset tests together" = "global"),
          selected = "per_gene"),
        shiny::numericInput("dpi", "PNG/TIFF DPI", value = 600, min = 150, max = 1200, step = 50),
        shiny::checkboxInput("auto_dimensions", "Automatic figure size", TRUE),
        shiny::conditionalPanel("input.auto_dimensions == false",
          shiny::fluidRow(shiny::column(6, shiny::numericInput("plot_width", "Width (in)", 12.5, min = 4, max = 30)),
                         shiny::column(6, shiny::numericInput("plot_height", "Height (in)", 9, min = 3, max = 30)))),
        shiny::checkboxInput("overview_plot", "Include multi-gene validation overview", TRUE),
        shiny::checkboxInput("save_svg", "Also export SVG where supported", TRUE),
        shiny::checkboxInput("save_averaged",
          "Optional: export full duplicate-averaged expression matrices (slower)", FALSE),
        tags$p(class = "small-note",
          "Fast default: selected genes are filtered before duplicate averaging. GSE labels come from filenames without automatic rewriting.")),
      div(class = "panelbox",
        div(class = "section-title", "4. Run and save"),
        shiny::actionButton("run_validation", "Run Independent Validation", class = "btn-run"),
        tags$hr(),
        shiny::textOutput("status_text"),
        shiny::uiOutput("save_again_ui"),
        tags$p(class = "small-note", "After analysis, download a ZIP or individual tables directly. You may also choose a Windows folder for saving everything."))
    ),
    shiny::column(width = 8,
      div(class = "panelbox",
        shiny::tabsetPanel(id = "main_tabs",
          shiny::tabPanel("Dataset Review", tags$br(), shiny::tableOutput("dataset_review"),
            tags$hr(), tags$p(class = "auto-group", "Groups are recognized automatically from Control/Disease column names. Rows marked Ready have samples in both groups. Please verify counts before running.")),
          shiny::tabPanel("Statistics", tags$br(), shiny::tableOutput("stats_view")),
          shiny::tabPanel("Download Results", tags$br(),
            tags$h3("Download completed analysis"),
            tags$p("After an analysis completes, download the full results as one ZIP file. You can also download individual tables, reports and the selected gene figure."),
            shiny::uiOutput("download_panel"),
            tags$hr(), tags$p("To save all files directly into a selected folder on this computer, use Choose Save Folder / Save Again on the left.")),
          shiny::tabPanel("Validation Figure", tags$br(),
            shiny::selectInput("selected_gene_plot", "Choose gene", choices = NULL),
            shiny::plotOutput("validation_figure", height = "640px")),
          shiny::tabPanel("Multi-Gene Overview", tags$br(),
            shiny::plotOutput("overview_figure", height = "560px"),
            tags$p(class = "small-note", "Overview shows the direction of statistically significant changes, not a pooled effect size. Dataset expression scales may differ.")),
          shiny::tabPanel("User Guide", tags$br(),
            tags$h3("User Guide — Step-by-Step"),
            tags$ol(class = "helpsteps",
              tags$li(tags$b("Choose gene count."), " Start with 1 for STAT3 or choose a larger number for independent validation of several targets."),
              tags$li(tags$b("Provide gene names."), " Enter one official gene symbol in each text box."),
              tags$li(tags$b("Upload GEO CSV matrices."), " Select one or several processed expression datasets. The file names should contain their correct GSE accessions."),
              tags$li(tags$b("Select datasets."), " Check only the datasets you want to process. The software uses the GSE number found in each filename without silently changing it."),
              tags$li(tags$b("Verify automatically detected groups."), " Sample columns named control, control.1, control_2, disease, disease.1 and similar names are recognized automatically. In Dataset Review, check the Control and Disease sample counts and Ready status. No manual pattern editing is needed."),
              tags$li(tags$b("Review statistics."), " The original two-sided Wilcoxon rank-sum test and BH-FDR adjustment are retained. With multiple genes, the default BH correction is applied separately per gene across datasets."),
              tags$li(tags$b("Run analysis."), " Click Run Independent Validation. Fast processing averages duplicate entries only for the selected genes. Full-matrix average export is optional and slower."),
              tags$li(tags$b("Download your results."), " After analysis, choose Download Complete ZIP to receive all output, or download individual tables/figures. A folder-selection option is also provided."),
              tags$li(tags$b("Save to your computer folder."), " Optionally, after analysis click Choose output directory and Save Results to write the complete results into a selected local folder."),
              tags$li(tags$b("Interpret correctly."), " A missing gene-dataset comparison is reported, not imputed. Distinct tissues or sample designs should not automatically be pooled.")),
            tags$p("Saved outputs: combined expression table, Wilcoxon + BH statistics, dataset QC, missing-comparison list, gene-specific plots (PNG, TIFF, PDF and optional SVG), optional overview, Methods.txt and Results.txt.")),
          shiny::tabPanel("About", tags$br(),
            tags$h3("Independent GEO Gene Validation Studio"),
            tags$p(tags$b("Version: "), "1.2"),
            tags$p(tags$b("Developed by: "), "Md Fazlay Rabbi"),
            tags$p(tags$b("Affiliation: "), "Department of Zoology, University of Rajshahi, Bangladesh"),
            tags$p(tags$b("Email: "), tags$a(href = "mailto:rabbi.ms244@gmail.com", "rabbi.ms244@gmail.com")),
            tags$p(tags$b("Research objective: "), "Validation of selected gene expression differences across independent GEO cohorts."),
            tags$p("This software extends the original six-dataset STAT3 analysis while preserving its dataset-wise Wilcoxon/BH workflow and cyan/red publication boxplots."))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  roots <- create_roots()
  shinyFiles::shinyDirChoose(input, "save_dir", roots = roots, session = session)
  chosen_save_dir <- shiny::reactiveVal(NULL)
  result <- shiny::reactiveVal(NULL)
  status <- shiny::reactiveVal("Ready — select genes and GEO expression files.")

  output$status_text <- shiny::renderText(status())
  output$save_again_ui <- shiny::renderUI({
    if (is.null(result())) return(NULL)
    shiny::actionButton("open_save_dialog", "Choose Save Folder / Save Again",
                         class = "btn btn-default")
  })

  output$download_panel <- shiny::renderUI({
    res <- result()
    if (is.null(res)) return(tags$p(class = "small-note", "Run the analysis first to enable downloads."))
    shiny::tagList(
      tags$div(class = "downloads-row",
        shiny::downloadButton("download_complete_zip", "Download Complete ZIP", class = "btn-download"),
        shiny::downloadButton("download_statistics", "Wilcoxon + BH CSV", class = "btn-download"),
        shiny::downloadButton("download_expression", "Expression Data CSV", class = "btn-download")),
      tags$div(class = "downloads-row",
        shiny::downloadButton("download_methods", "Methods.txt", class = "btn-download"),
        shiny::downloadButton("download_results_text", "Results.txt", class = "btn-download"),
        shiny::downloadButton("download_selected_figure", "Selected Gene Figure (PNG)", class = "btn-download")),
      tags$p(class = "small-note", "ZIP includes all completed analysis outputs: figures, statistical tables, QC, missing-gene report and text summaries.")
    )
  })

  make_complete_zip <- function(destination) {
    res <- result()
    if (is.null(res)) stop("No finished analysis is available.")
    contents <- list.files(res$run_dir, recursive = TRUE, full.names = FALSE,
                           all.files = FALSE, include.dirs = FALSE)
    if (!length(contents)) stop("Analysis output folder is empty.")
    zip::zipr(zipfile = destination, files = contents,
              root = res$run_dir, mode = "mirror", include_directories = FALSE)
  }

  output$download_complete_zip <- shiny::downloadHandler(
    filename = function() paste0("GEO_Gene_Validation_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = make_complete_zip, contentType = "application/zip")
  output$download_complete_zip_modal <- shiny::downloadHandler(
    filename = function() paste0("GEO_Gene_Validation_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = make_complete_zip, contentType = "application/zip")

  make_file_download <- function(input_name, output_name, mime_type = "text/plain") {
    shiny::downloadHandler(
      filename = function() output_name,
      content = function(file) {
        res <- result()
        if (is.null(res)) stop("Run an analysis before downloading.")
        source_file <- file.path(res$run_dir, input_name)
        if (!file.exists(source_file)) stop("Output not found: ", input_name)
        if (!file.copy(source_file, file, overwrite = TRUE)) stop("Could not copy download file.")
      },
      contentType = mime_type)
  }
  output$download_statistics <- make_file_download("Wilcoxon_BH_Statistics.csv",
    "Wilcoxon_BH_Statistics.csv", "text/csv")
  output$download_expression <- make_file_download("Combined_Gene_Expression_Long.csv",
    "Combined_Gene_Expression_Long.csv", "text/csv")
  output$download_methods <- make_file_download("Methods.txt", "Methods.txt")
  output$download_results_text <- make_file_download("Results.txt", "Results.txt")
  output$download_selected_figure <- shiny::downloadHandler(
    filename = function() {
      gene <- input$selected_gene_plot %||% "Gene"
      paste0(safe_name(gene), "_Independent_Validation.png")
    },
    content = function(file) {
      res <- result()
      if (is.null(res)) stop("Run an analysis before downloading.")
      gene <- input$selected_gene_plot
      if (is.null(gene) || !(gene %in% names(res$plots)))
        stop("Select an available gene in the Validation Figure tab.")
      source_file <- file.path(res$run_dir, "Figures",
        paste0(safe_name(gene), "_Independent_Validation.png"))
      if (!file.exists(source_file) || !file.copy(source_file, file, overwrite = TRUE))
        stop("Could not download the requested figure.")
    }, contentType = "image/png")

  show_save_prompt <- function() {
    shiny::showModal(shiny::modalDialog(
      title = "Analysis finished — select output folder",
      tags$p("Calculations and figures are ready. Download the complete ZIP directly, or choose a directory to save everything."),
      shiny::downloadButton("download_complete_zip_modal", "Download Complete ZIP", class = "btn-download"),
      tags$hr(),
      shinyFiles::shinyDirButton("save_dir", "Choose output directory", "Select a folder",
        class = "btn btn-primary"),
      tags$br(), tags$br(), shiny::textOutput("chosen_dir_text"),
      footer = shiny::tagList(shiny::modalButton("Save later"),
        shiny::actionButton("confirm_save", "Save Results", class = "btn-success")),
      easyClose = FALSE, size = "m"))
  }

  shiny::observeEvent(input$open_save_dialog, {
    if (!is.null(result())) show_save_prompt()
  })

  output$gene_inputs <- shiny::renderUI({
    count <- min(50L, max(1L, as.integer(input$gene_count %||% 1L)))
    lapply(seq_len(count), function(i) {
      id <- paste0("gene_", i)
      old_value <- isolate(input[[id]])
      if (is.null(old_value)) old_value <- if (i == 1L) "STAT3" else ""
      shiny::textInput(id, paste("Gene", i), value = old_value,
                       placeholder = "e.g. STAT3")
    })
  })

  uploaded_index <- shiny::reactive({
    upload <- input$expression_files
    if (is.null(upload) || !nrow(upload)) return(NULL)
    data.frame(row = seq_len(nrow(upload)), filename = upload$name,
      datapath = upload$datapath,
      dataset_id = vapply(upload$name, get_dataset_label, character(1)),
      stringsAsFactors = FALSE)
  })

  dataset_preview <- shiny::reactive({
    info <- uploaded_index()
    if (is.null(info)) return(NULL)
    checks <- lapply(info$datapath, preview_sample_groups)
    data.frame(Dataset = info$dataset_id,
      File = info$filename,
      Control_samples = vapply(checks, `[[`, integer(1), "n_control"),
      Disease_samples = vapply(checks, `[[`, integer(1), "n_disease"),
      Status = vapply(checks, `[[`, character(1), "status"),
      stringsAsFactors = FALSE)
  })

  shiny::observeEvent(uploaded_index(), {
    info <- uploaded_index()
    if (is.null(info)) return()
    labels <- paste0(info$dataset_id, " — ", info$filename)
    choices <- stats::setNames(as.character(info$row), labels)
    shiny::updateCheckboxGroupInput(session, "datasets_to_use", choices = choices,
      selected = as.character(info$row))
  })

  output$dataset_review <- shiny::renderTable({
    checks <- dataset_preview()
    if (is.null(checks)) return(data.frame(Status = "Upload one or more GEO expression CSV files."))
    checks
  }, striped = TRUE, bordered = TRUE)

  output$invalid_dataset_review <- shiny::renderTable({
    checks <- dataset_preview()
    selected <- suppressWarnings(as.integer(input$datasets_to_use))
    if (is.null(checks) || !length(selected)) return(NULL)
    bad <- checks[selected, , drop = FALSE]
    bad[bad$Status != "Ready", , drop = FALSE]
  }, striped = TRUE, bordered = TRUE)

  output$stats_view <- shiny::renderTable({
    res <- result()
    if (is.null(res)) return(data.frame(Status = "Run analysis to view Wilcoxon + BH results."))
    stat <- res$stats
    keep <- c("Gene", "Dataset", "n_control", "n_disease", "median_control",
              "median_disease", "median_difference", "p_value", "p_adj_BH", "significance")
    stat[, keep, drop = FALSE]
  }, digits = 5, striped = TRUE, bordered = TRUE)

  output$validation_figure <- shiny::renderPlot({
    res <- result()
    if (is.null(res)) return(NULL)
    gene <- input$selected_gene_plot
    shiny::req(gene, gene %in% names(res$plots))
    print(res$plots[[gene]])
  }, res = 100)

  output$overview_figure <- shiny::renderPlot({
    res <- result()
    if (is.null(res) || is.null(res$overview)) return(NULL)
    print(res$overview)
  }, res = 100)

  shiny::observeEvent(input$run_validation, {
    info <- uploaded_index()
    if (is.null(info)) {
      shiny::showNotification("First upload at least one expression CSV.", type = "error")
      return()
    }
    selected <- suppressWarnings(as.integer(input$datasets_to_use))
    if (!length(selected) || anyNA(selected)) {
      shiny::showNotification("Select at least one GEO dataset.", type = "error")
      return()
    }
    count <- as.integer(input$gene_count)
    if (!is.finite(count) || count < 1L || count > 50L) {
      shiny::showNotification("Gene number must be between 1 and 50.", type = "error")
      return()
    }
    genes <- vapply(seq_len(count), function(i) clean_gene(input[[paste0("gene_", i)]] %||% ""),
                    character(1))
    if (any(!nzchar(genes)) || anyNA(genes) || anyDuplicated(genes)) {
      shiny::showNotification("Enter unique, non-empty gene symbols for every field.", type = "error")
      return()
    }
    checks <- dataset_preview()
    selected_checks <- checks[selected, , drop = FALSE]
    bad <- selected_checks[selected_checks$Status != "Ready", , drop = FALSE]
    if (nrow(bad)) {
      shiny::showModal(shiny::modalDialog(
        title = "Please check Control / Disease columns",
        tags$p("These datasets do not have recognizable Control and Disease columns:"),
        shiny::tableOutput("invalid_dataset_review"),
        tags$p("Use names such as control, control.1, disease, disease.1 and upload again."),
        easyClose = TRUE))
      return()
    }
    settings <- list(fdr_scope = input$fdr_scope,
      dpi = as.integer(input$dpi),
      auto_dimensions = isTRUE(input$auto_dimensions),
      width = input$plot_width %||% 12.5,
      height = input$plot_height %||% 9,
      save_svg = isTRUE(input$save_svg),
      save_averaged = isTRUE(input$save_averaged),
      overview = isTRUE(input$overview_plot))
    result(NULL)
    status("Running independent validation...")
    tryCatch({
      res <- shiny::withProgress(message = "Analyzing selected GEO datasets...", value = 0, {
        perform_validation(info, selected, genes, settings,
          progress = function(j, n, id, stage) {
            shiny::setProgress(value = if (identical(stage, "start")) (j - 1) / n else j / n,
              detail = paste(if (identical(stage, "start")) "Processing" else "Completed",
                             id, sprintf("(%d/%d)", j, n)))
          })
      })
      result(res)
      shiny::updateSelectInput(session, "selected_gene_plot",
        choices = names(res$plots), selected = names(res$plots)[1L])
      shiny::updateTabsetPanel(session, "main_tabs", selected = "Download Results")
      status(paste("Finished:", length(res$genes_found), "genes across",
                   length(res$dataset_order), "datasets. Download results or save to a folder."))
      chosen_save_dir(NULL)
      show_save_prompt()
    }, error = function(e) {
      status(paste("Error:", conditionMessage(e)))
      shiny::showModal(shiny::modalDialog(title = "Validation could not finish",
        tags$p(style = "color:#9D2323;font-weight:bold;", conditionMessage(e)),
        footer = shiny::modalButton("Close"), easyClose = TRUE))
    })
  })

  shiny::observeEvent(input$save_dir, {
    parsed <- tryCatch(shinyFiles::parseDirPath(roots, input$save_dir),
                       error = function(e) character(0))
    if (length(parsed) && nzchar(parsed[[1L]]) && dir.exists(parsed[[1L]])) {
      chosen_save_dir(parsed[[1L]])
    }
  })

  output$chosen_dir_text <- shiny::renderText({
    dir <- chosen_save_dir()
    if (is.null(dir)) "No output directory chosen yet." else paste("Selected:", dir)
  })

  shiny::observeEvent(input$confirm_save, {
    res <- result()
    if (is.null(res)) return()
    target <- chosen_save_dir()
    if (is.null(target)) {
      shiny::showNotification("Choose an output directory first.", type = "warning")
      return()
    }
    tryCatch({
      saved <- save_run_to_directory(res$run_dir, target)
      shiny::removeModal()
      status(paste("Results saved to:", saved))
      shiny::showModal(shiny::modalDialog(title = "Results saved successfully",
        tags$p("All statistical tables and publication figures were saved in:"),
        tags$code(saved), footer = shiny::modalButton("Close"), easyClose = TRUE))
    }, error = function(e) {
      shiny::showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })
}

# Run in RStudio by clicking Source or: source(file.choose())
shiny::runApp(shiny::shinyApp(ui = ui, server = server), launch.browser = TRUE)
