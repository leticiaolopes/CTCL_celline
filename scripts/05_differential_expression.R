# differential expression - ctcl cell lines
# deg cutoff: padj < 0.01 && abs(log2fc) >= 2

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(apeglm)
library(rtracklayer)
library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Analysis scope. The default preserves the original all-gene workflow.
gene_scope <- getOption("ctcl.gene_scope", "all")
if (!gene_scope %in% c("all", "protein_coding")) {
  stop("Unsupported ctcl.gene_scope: ", gene_scope)
}
analysis_suffix <- if (gene_scope == "protein_coding") "_protein_coding" else ""

# paths
deseq2_dir <- file.path(paths$results, paste0("deseq2", analysis_suffix))

de_dir <- file.path(paths$results, paste0("differential_expression", analysis_suffix))

figures_de_dir <- file.path(
  paths$figures,
  "drafts",
  paste0("deseq2", analysis_suffix),
  "differential_expression")

dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(figures_de_dir, recursive = TRUE, showWarnings = FALSE)

# load deseq2 object
dds <- readRDS(file.path(deseq2_dir, "dds.rds"))

if (gene_scope == "protein_coding") {
  protein_coding_ids <- readRDS(
    file.path(deseq2_dir, "protein_coding_gene_ids.rds"))
  if (!all(rownames(dds) %in% protein_coding_ids)) {
    stop("The protein-coding DESeq2 object contains genes outside the selected biotype.")
  }
}

resultsNames(dds)

# raw pairwise deseq2 results
get_de_result <- function(dds, numerator, denominator) {
  res <- results(dds, contrast = c("cell_line_model", numerator, denominator))

  as.data.frame(res) |>
    rownames_to_column("Geneid") |>
    mutate(comparison = paste0(numerator, "_vs_", denominator)) |>
    arrange(padj)}

res_myla_vs_hh <- get_de_result(dds, numerator = "MyLa", denominator = "HH")

res_hut78_vs_hh <- get_de_result(dds, numerator = "HuT78", denominator = "HH")

res_seax_vs_hh <- get_de_result(dds, numerator = "SeAx", denominator = "HH")

res_hut78_vs_myla <- get_de_result(dds, numerator = "HuT78", denominator = "MyLa")

res_seax_vs_myla <- get_de_result(dds, numerator = "SeAx", denominator = "MyLa")

res_seax_vs_hut78 <- get_de_result(dds, numerator = "SeAx", denominator = "HuT78")

# save raw deseq2 results
write_csv(res_myla_vs_hh, file.path(de_dir, "MyLa_vs_HH.csv"))

write_csv(res_hut78_vs_hh, file.path(de_dir, "HuT78_vs_HH.csv"))

write_csv(res_seax_vs_hh, file.path(de_dir, "SeAx_vs_HH.csv"))

write_csv(res_hut78_vs_myla, file.path(de_dir, "HuT78_vs_MyLa.csv"))

write_csv(res_seax_vs_myla, file.path(de_dir, "SeAx_vs_MyLa.csv"))

write_csv(res_seax_vs_hut78, file.path(de_dir, "SeAx_vs_HuT78.csv"))

# apeglm shrinkage: hh reference
res_myla_vs_hh_shrunk <- lfcShrink(dds, coef = "cell_line_model_MyLa_vs_HH", type = "apeglm")

res_hut78_vs_hh_shrunk <- lfcShrink(dds, coef = "cell_line_model_HuT78_vs_HH", type = "apeglm")

res_seax_vs_hh_shrunk <- lfcShrink(dds, coef = "cell_line_model_SeAx_vs_HH", type = "apeglm")

res_myla_vs_hh_shrunk_tbl <- as.data.frame(res_myla_vs_hh_shrunk) |>
  rownames_to_column("Geneid") |>
  mutate(comparison = "MyLa_vs_HH") |>
  arrange(padj)

res_hut78_vs_hh_shrunk_tbl <- as.data.frame(res_hut78_vs_hh_shrunk) |>
  rownames_to_column("Geneid") |>
  mutate(comparison = "HuT78_vs_HH") |>
  arrange(padj)

res_seax_vs_hh_shrunk_tbl <- as.data.frame(res_seax_vs_hh_shrunk) |>
  rownames_to_column("Geneid") |>
  mutate(comparison = "SeAx_vs_HH") |>
  arrange(padj)

# apeglm shrinkage: myla reference
dds_myla_ref <- dds

dds_myla_ref$cell_line_model <- relevel(dds_myla_ref$cell_line_model, ref = "MyLa")

design(dds_myla_ref) <- ~cell_line_model

dds_myla_ref <- DESeq(dds_myla_ref, fitType = "local")

resultsNames(dds_myla_ref)

res_hut78_vs_myla_shrunk <- lfcShrink(dds_myla_ref, coef = "cell_line_model_HuT78_vs_MyLa", type = "apeglm")

res_seax_vs_myla_shrunk <- lfcShrink(dds_myla_ref, coef = "cell_line_model_SeAx_vs_MyLa", type = "apeglm")

res_hut78_vs_myla_shrunk_tbl <- as.data.frame(res_hut78_vs_myla_shrunk) |>
  rownames_to_column("Geneid") |>
  mutate(comparison = "HuT78_vs_MyLa") |>
  arrange(padj)

res_seax_vs_myla_shrunk_tbl <- as.data.frame(res_seax_vs_myla_shrunk) |>
  rownames_to_column("Geneid") |>
  mutate(comparison = "SeAx_vs_MyLa") |>
  arrange(padj)

# apeglm shrinkage: hut78 reference
dds_hut78_ref <- dds

dds_hut78_ref$cell_line_model <- relevel(dds_hut78_ref$cell_line_model, ref = "HuT78")

design(dds_hut78_ref) <- ~cell_line_model

dds_hut78_ref <- DESeq(dds_hut78_ref, fitType = "local")

resultsNames(dds_hut78_ref)

res_seax_vs_hut78_shrunk <- lfcShrink(dds_hut78_ref, coef = "cell_line_model_SeAx_vs_HuT78", type = "apeglm")

res_seax_vs_hut78_shrunk_tbl <- as.data.frame(res_seax_vs_hut78_shrunk) |>
  rownames_to_column("Geneid") |>
  mutate(comparison = "SeAx_vs_HuT78") |>
  arrange(padj)

# combine shrunken results
all_shrunk_results <- bind_rows(res_myla_vs_hh_shrunk_tbl, res_hut78_vs_hh_shrunk_tbl, res_seax_vs_hh_shrunk_tbl,
  res_hut78_vs_myla_shrunk_tbl, res_seax_vs_myla_shrunk_tbl, res_seax_vs_hut78_shrunk_tbl)

saveRDS(all_shrunk_results, file.path(de_dir, "all_shrunk_results.rds"))

write_csv(all_shrunk_results, file.path(de_dir, "all_shrunk_results.csv"))

# gene annotation
annotation_rds <- file.path(deseq2_dir, "gencode_v50_gene_annotation.rds")

if (file.exists(annotation_rds)) {
  gene_annotation <- readRDS(annotation_rds)
} else {
  gtf <- import(file.path(paths$data, "reference", "gencode.v50.annotation.gtf.gz"))

  gene_annotation <- as.data.frame(gtf) |>
    filter(type == "gene") |>
    transmute(Geneid = gene_id, Ensembl = str_remove(gene_id, "\\.\\d+$"), gene_symbol = gene_name, gene_type = gene_type) |>
    distinct(Geneid, .keep_all = TRUE)
}

all_shrunk_results_annotated <- all_shrunk_results |>
  left_join(gene_annotation, by = "Geneid")

# annotation qc
missing_gene_symbols <- sum(is.na(all_shrunk_results_annotated$gene_symbol))

cat("\nGenes without gene symbol:", missing_gene_symbols, "\n")

# deg classification
padj_cutoff <- 0.01
lfc_cutoff <- 2

all_shrunk_results_annotated <- all_shrunk_results_annotated |>
  mutate(DEG_class = case_when(!is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "Up",
    !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Down", TRUE ~ "NS"))

# DEG summary
de_summary <- all_shrunk_results_annotated |>
  dplyr::count(comparison, DEG_class) |>
  tidyr::pivot_wider(names_from = DEG_class, values_from = n, values_fill = 0)

de_summary

write_csv(de_summary, file.path(de_dir, "DEG_summary_padj001_log2FC2.csv"))

# save annotated results
saveRDS(all_shrunk_results_annotated, file.path(de_dir, "all_shrunk_results_annotated.rds"))

write_csv(all_shrunk_results_annotated, file.path(de_dir, "all_shrunk_results_annotated.csv"))

# volcano plot function
plot_volcano <- function(data, comparison_label, lfc_cutoff = 2, padj_cutoff = 0.01, n_labels = 5, y_max = 100) {
  plot_data <- data |>
    filter(comparison == comparison_label, !is.na(padj), !is.na(log2FoldChange)) |>
    mutate(minus_log10_padj = -log10(padj))

  max_finite_y <- max(plot_data$minus_log10_padj[is.finite(plot_data$minus_log10_padj)], na.rm = TRUE)

  plot_data <- plot_data |>
    mutate(minus_log10_padj = case_when(is.infinite(minus_log10_padj) ~ max_finite_y + 1, TRUE ~
      minus_log10_padj))

  label_genes <- bind_rows(plot_data |>
    filter(DEG_class == "Up") |>
    arrange(desc(log2FoldChange)) |>
    slice_head(n = n_labels), plot_data |>
    filter(DEG_class == "Down") |>
    arrange(log2FoldChange) |>
    slice_head(n = n_labels))

  ggplot() + geom_point(data = plot_data |>
    filter(DEG_class == "NS"), aes(x = log2FoldChange, y = minus_log10_padj, color = DEG_class),
    alpha = 0.3, size = 0.8) + geom_point(data = plot_data |>
    filter(DEG_class != "NS"), aes(x = log2FoldChange, y = minus_log10_padj, color = DEG_class),
    alpha = 0.7, size = 1) + ggrepel::geom_text_repel(data = label_genes, aes(x = log2FoldChange,
    y = minus_log10_padj, label = gene_symbol), size = 2.3, color = black, box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0, segment.size = 0.25, max.overlaps = Inf) + scale_color_up_down() + geom_vline(xintercept = c(-lfc_cutoff,
    lfc_cutoff), linetype = "dashed", linewidth = 0.3, color = mid_grey) + geom_hline(yintercept = -log10(padj_cutoff),
    linetype = "dashed", linewidth = 0.3, color = mid_grey) + labs(x = "Shrunken log2 fold change",
    y = expression(-log[10]("adjusted p-value")), color = NULL, title = gsub("_", " ", comparison_label)) +
    coord_cartesian(ylim = c(0, y_max)) + theme_clean() + ggplot2::theme(legend.position = "top",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))}

# volcano plots
p_volcano_myla_hh <- plot_volcano(all_shrunk_results_annotated, "MyLa_vs_HH")

p_volcano_hut78_hh <- plot_volcano(all_shrunk_results_annotated, "HuT78_vs_HH")

p_volcano_seax_hh <- plot_volcano(all_shrunk_results_annotated, "SeAx_vs_HH")

p_volcano_hut78_myla <- plot_volcano(all_shrunk_results_annotated, "HuT78_vs_MyLa")

p_volcano_seax_myla <- plot_volcano(all_shrunk_results_annotated, "SeAx_vs_MyLa")

p_volcano_seax_hut78 <- plot_volcano(all_shrunk_results_annotated, "SeAx_vs_HuT78")

# save volcano plots
save_single(file.path(figures_de_dir, "volcano_MyLa_vs_HH.pdf"), p_volcano_myla_hh, height = height_medium)

save_single(file.path(figures_de_dir, "volcano_HuT78_vs_HH.pdf"), p_volcano_hut78_hh, height = height_medium)

save_single(file.path(figures_de_dir, "volcano_SeAx_vs_HH.pdf"), p_volcano_seax_hh, height = height_medium)

save_single(file.path(figures_de_dir, "volcano_HuT78_vs_MyLa.pdf"), p_volcano_hut78_myla, height = height_medium)

save_single(file.path(figures_de_dir, "volcano_SeAx_vs_MyLa.pdf"), p_volcano_seax_myla, height = height_medium)

save_single(file.path(figures_de_dir, "volcano_SeAx_vs_HuT78.pdf"), p_volcano_seax_hut78, height = height_medium)

# top deg tables
top_degs <- all_shrunk_results_annotated |>
  filter(DEG_class %in% c("Up", "Down")) |>
  group_by(comparison, DEG_class) |>
  arrange(desc(abs(log2FoldChange)), .by_group = TRUE) |>
  slice_head(n = 100) |>
  ungroup()

write_csv(top_degs, file.path(de_dir, "top100_DEGs_per_direction_per_comparison.csv"))

# heatmap setup

# load vst matrix
vst_matrix <- readRDS(file.path(deseq2_dir, "vst_matrix.rds"))

# sample metadata
coldata_heatmap <- SummarizedExperiment::colData(dds) |>
  as.data.frame() |>
  tibble::rownames_to_column("original_sample_name")

sample_name_map <- coldata_heatmap |>
  dplyr::select(original_sample_name, analysis_sample_name, analysis_cell_line)

# sample order for final heatmaps
sample_order <- sample_name_map |>
  dplyr::mutate(analysis_cell_line = factor(analysis_cell_line, levels = c("HH", "MyLa", "HuT 78",
    "SeAx"))) |>
  dplyr::arrange(analysis_cell_line, analysis_sample_name) |>
  dplyr::pull(analysis_sample_name)

print(sample_order)

# top-deg heatmap
heatmap_genes_tbl <- all_shrunk_results_annotated |>
  dplyr::filter(DEG_class %in% c("Up", "Down")) |>
  dplyr::group_by(comparison, DEG_class) |>
  dplyr::arrange(dplyr::desc(abs(log2FoldChange)), .by_group = TRUE) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup()

heatmap_genes <- heatmap_genes_tbl |>
  dplyr::pull(Geneid) |>
  unique()

cat("\nUnique genes in top DEG heatmap:", length(heatmap_genes), "\n")

# gene labels
gene_labels <- all_shrunk_results_annotated |>
  dplyr::select(Geneid, Ensembl, gene_symbol) |>
  dplyr::distinct(Geneid, .keep_all = TRUE)

gene_labels$orgdb_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = gene_labels$Ensembl, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first")

gene_labels <- gene_labels |>
  dplyr::mutate(manual_symbol = dplyr::case_when(Ensembl == "ENSG00000243276" ~ "LOC105374059", Ensembl ==
    "ENSG00000299902" ~ "HSALNG0033462", TRUE ~ NA_character_), heatmap_label = dplyr::case_when(!is.na(orgdb_symbol) &
    orgdb_symbol != "" & !stringr::str_detect(orgdb_symbol, "^ENSG") ~ orgdb_symbol, !is.na(manual_symbol) ~
    manual_symbol, !is.na(gene_symbol) & gene_symbol != "" & !stringr::str_detect(gene_symbol, "^ENSG") ~
    gene_symbol, !is.na(Ensembl) ~ Ensembl, TRUE ~ stringr::str_remove(Geneid, "\\.\\d+$")))

gene_labels_vector <- gene_labels$heatmap_label

names(gene_labels_vector) <- gene_labels$Geneid

heatmap_matrix <- vst_matrix[heatmap_genes, , drop = FALSE]

rownames(heatmap_matrix) <- gene_labels_vector[rownames(heatmap_matrix)]

heatmap_matrix_z <- t(scale(t(heatmap_matrix)))

# color scale
heatmap_z_limit <- 1.5

heatmap_breaks <- seq(-heatmap_z_limit, heatmap_z_limit, length.out = length(heatmap_deg) + 1)

# sample identity mapping
original_matrix_names <- colnames(heatmap_matrix_z)

new_sample_names <- sample_name_map$analysis_sample_name[match(original_matrix_names, sample_name_map$original_sample_name)]

if (any(is.na(new_sample_names))) {
  stop("Could not map all original sample names to analysis sample names.")}

colnames(heatmap_matrix_z) <- unname(new_sample_names)

# sample annotation
annotation_col <- sample_name_map |>
  dplyr::select(analysis_sample_name, analysis_cell_line) |>
  dplyr::rename(`Cell line` = analysis_cell_line) |>
  tibble::column_to_rownames("analysis_sample_name")

annotation_colors <- list(`Cell line` = cell_lines)

# apply final sample order
heatmap_matrix_z <- heatmap_matrix_z[, sample_order, drop = FALSE]

annotation_col <- annotation_col[sample_order, , drop = FALSE]

heatmap_sample_check <- data.frame(heatmap_sample = colnames(heatmap_matrix_z), cell_line = annotation_col[colnames(heatmap_matrix_z),
  "Cell line"])

print(heatmap_sample_check)

# heatmap
top_deg_heatmap <- ComplexHeatmap::pheatmap(heatmap_matrix_z, color = heatmap_deg, breaks = heatmap_breaks,
  annotation_col = annotation_col, annotation_colors = annotation_colors, cluster_rows = TRUE, cluster_cols = FALSE,
  clustering_distance_rows = "correlation", clustering_method = "complete", show_rownames = TRUE, border_color = NA,
  fontsize = 6.5, fontsize_row = 5.5, fontsize_col = 6.5, angle_col = "45", treeheight_row = 20, main = "Top differentially expressed genes",
  name = "Z-score", heatmap_legend_param = list(at = c(-1.5, 0, 1.5), labels = c("-1.5", "0", "1.5"),
    legend_height = grid::unit(22, "mm"), legend_width = grid::unit(5, "mm"), border = black, title_position = "topleft"))

grDevices::pdf(file.path(figures_de_dir, "top_DEGs_global_heatmap.pdf"), width = width_heatmap, height = height_large,
  useDingbats = FALSE)
ComplexHeatmap::draw(top_deg_heatmap, heatmap_legend_side = "right")
grDevices::dev.off()

# all-deg heatmap
# gene selection
all_deg_genes <- all_shrunk_results_annotated |>
  dplyr::filter(DEG_class %in% c("Up", "Down")) |>
  dplyr::pull(Geneid) |>
  unique()

n_top_heatmap_genes <- length(heatmap_genes)

n_all_deg_heatmap_genes <- length(all_deg_genes)

cat("\nUnique DEGs in all-DEG heatmap:", n_all_deg_heatmap_genes, "\n")

# expression matrix
large_heatmap_matrix <- vst_matrix[all_deg_genes, , drop = FALSE]

# gene-wise z scores
large_heatmap_matrix_z <- t(scale(t(large_heatmap_matrix)))

large_heatmap_matrix_z <- large_heatmap_matrix_z[apply(large_heatmap_matrix_z, 1, function(x) {
  all(is.finite(x))
}), , drop = FALSE]

# sample identity mapping
large_original_names <- colnames(large_heatmap_matrix_z)

large_new_names <- sample_name_map$analysis_sample_name[match(large_original_names, sample_name_map$original_sample_name)]

if (any(is.na(large_new_names))) {
  stop("could not map all original sample names to analysis sample names...")}

colnames(large_heatmap_matrix_z) <- unname(large_new_names)

# sample order
large_heatmap_matrix_z <- large_heatmap_matrix_z[, sample_order, drop = FALSE]

# sample annotation
sample_cell_line <- sample_name_map$analysis_cell_line[match(sample_order, sample_name_map$analysis_sample_name)]

names(sample_cell_line) <- sample_order

top_annotation <- ComplexHeatmap::HeatmapAnnotation(`Cell line` = sample_cell_line, col = list(`Cell line` = cell_lines),
  show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

# color scale
large_heatmap_colors <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("#3F648C", "#F8F8F8", "#B55350"))

# heatmap
large_deg_heatmap <- ComplexHeatmap::Heatmap(large_heatmap_matrix_z, name = "Z-score", col = large_heatmap_colors,
  top_annotation = top_annotation, cluster_rows = TRUE, cluster_columns = FALSE, clustering_distance_rows = "pearson",
  clustering_method_rows = "complete", show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
  show_row_dend = TRUE, row_dend_width = grid::unit(7, "mm"), border = FALSE, use_raster = TRUE, raster_quality = 3,
  heatmap_legend_param = list(at = c(-1.5, 0, 1.5), labels = c("-1.5", "0", "1.5"), legend_height = grid::unit(22,
    "mm"), legend_width = grid::unit(5, "mm"), border = black, title_position = "topleft"), column_title = "All differentially expressed genes",
  column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"))

# export
graphics.off()

pdf(file.path(figures_de_dir, "all_DEGs_large_heatmap.pdf"), width = 5, height = 6.5, useDingbats = FALSE)

ComplexHeatmap::draw(large_deg_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")

dev.off()

# final summary
cat("\nDifferential expression analysis completed\n", "------------------------------------------\n",
  "Comparisons: 6\n", "Adjusted p-value cutoff: ", padj_cutoff, "\n", "Absolute shrunken log2FC cutoff: ",
  lfc_cutoff, "\n", "Top DEG heatmap genes: ", n_top_heatmap_genes, "\n", "All-DEG heatmap genes: ",
  n_all_deg_heatmap_genes, "\n", sep = "")
