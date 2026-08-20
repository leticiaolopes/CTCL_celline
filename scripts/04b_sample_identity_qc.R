# evaluate suspected sample swaps using internal expression signatures

# IMPORTANT: this script intentionally uses original sample identities 
# to investigate suspected sample swaps
# do not replace original labels with inferred labels in this QC 

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

library(DESeq2)
library(dplyr)
library(tibble)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

deseq2_dir <- file.path(paths$results, "deseq2")

qc_dir <- file.path(paths$qc, "sample_identity")

figures_dir <- file.path(paths$figures, "drafts", "deseq2", "sample_identity")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# load vst matrix && sample metadata
vst_matrix <- readRDS(file.path(deseq2_dir, "vst_matrix.rds"))

dds <- readRDS(file.path(deseq2_dir, "dds.rds"))

coldata <- SummarizedExperiment::colData(dds) |>
  as.data.frame()

if (!identical(colnames(vst_matrix),
  rownames(coldata))) {stop("vst matrix columns & coldata row names are not aligned!!!")}

# define trusted reference samples
hh_ref <- c("HH-1", "HH-2")

myla_ref <- c("M-1", "M-2")

hut78_ref <- c("78-1", "78-2")

seax_ref <- c("S-1", "S-2")

# build internal identity signatures
hh_mean <- rowMeans(vst_matrix[, hh_ref, drop = FALSE])

myla_mean <- rowMeans(vst_matrix[, myla_ref, drop = FALSE])

hh_myla_signature <- tibble(Geneid = rownames(vst_matrix), HH_mean = hh_mean, MyLa_mean = myla_mean,
  difference = HH_mean - myla_mean, abs_difference = abs(difference)) |>
  dplyr::arrange(dplyr::desc(abs_difference))

top_hh_myla <- hh_myla_signature |>
  dplyr::slice_head(n = 50) |>
  dplyr::pull(Geneid)

hut78_mean <- rowMeans(vst_matrix[, hut78_ref, drop = FALSE])

seax_mean <- rowMeans(vst_matrix[, seax_ref, drop = FALSE])

hut78_seax_signature <- tibble(Geneid = rownames(vst_matrix), HuT78_mean = hut78_mean, SeAx_mean = seax_mean,
  difference = HuT78_mean - SeAx_mean, abs_difference = abs(difference)) |>
  dplyr::arrange(dplyr::desc(abs_difference))

top_hut78_seax <- hut78_seax_signature |>
  dplyr::slice_head(n = 50) |>
  dplyr::pull(Geneid)

identity_genes <- unique(c(top_hh_myla, top_hut78_seax))

identity_matrix <- vst_matrix[identity_genes, , drop = FALSE]

identity_matrix_z <- t(scale(t(identity_matrix)))

identity_matrix_z <- identity_matrix_z[apply(identity_matrix_z, 1, function(x) {
  all(is.finite(x))
}), , drop = FALSE]

# sample annotation
sample_annotation <- coldata |>
  dplyr::transmute(cell_line = dplyr::recode(as.character(cell_line_model), HuT78 = "HuT 78"))

rownames(sample_annotation) <- rownames(coldata)

sample_annotation <- sample_annotation[colnames(identity_matrix_z), , drop = FALSE]

annotation_colors <- list(cell_line = cell_lines)

# internal identity signature heatmap
pheatmap::pheatmap(identity_matrix_z, color = heatmap_deg, breaks = seq(-1.5, 1.5, length.out = length(heatmap_deg) +
  1), legend_breaks = c(-1.5, 0, 1.5), annotation_col = sample_annotation, annotation_colors = annotation_colors,
  cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, border_color = NA, fontsize = 7,
  fontsize_col = 7, angle_col = 45, filename = file.path(figures_dir, "internal_identity_signature_heatmap.pdf"),
  width = width_double, height = height_large)

# paired heatmaps: before and after identity correction
export_identity_heatmap <- function(matrix, display_names, groups, filename) {
  display_matrix <- matrix
  colnames(display_matrix) <- display_names

  annotation <- data.frame(`Cell line` = groups, row.names = display_names, check.names = FALSE)

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    `Cell line` = groups, col = list(`Cell line` = cell_lines),
    show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

  heatmap <- ComplexHeatmap::Heatmap(display_matrix, name = "Z-score", 
                                     col = circlize::colorRamp2(c(-1.5, 0, 1.5), 
                                     c("#3F648C", "#F7F7F7", "#B55350")), 
                                     top_annotation = top_annotation, cluster_rows = TRUE,
    cluster_columns = TRUE, show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 7), border = FALSE, heatmap_legend_param = list(at = c(-1.5,
      0, 1.5), labels = c("-1.5", "0", "1.5"), legend_height = grid::unit(22, "mm"), 
      legend_width = grid::unit(5, "mm"), border = black, title_position = "topleft"))

  grDevices::pdf(file.path(figures_dir, filename), width = 5.6, height = 5.2, useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()}

original_identity_groups <- dplyr::recode(as.character(coldata$original_cell_line), HuT78 = "HuT 78")

export_identity_heatmap(identity_matrix_z, colnames(identity_matrix_z), original_identity_groups, "internal_identity_signature_original.pdf")

export_identity_heatmap(identity_matrix_z, coldata$analysis_sample_name, as.character(coldata$analysis_cell_line), "internal_identity_signature_inferred.pdf")

# reference centroid correlations
reference_centroids <- cbind(HH = rowMeans(vst_matrix[, hh_ref, drop = FALSE]), MyLa = rowMeans(vst_matrix[,
  myla_ref, drop = FALSE]), HuT78 = rowMeans(vst_matrix[, hut78_ref, drop = FALSE]), SeAx = rowMeans(vst_matrix[,
  seax_ref, drop = FALSE]))

identity_cor <- cor(vst_matrix, reference_centroids, method = "pearson")

identity_cor <- identity_cor[colnames(vst_matrix), c("HH", "MyLa", "HuT78", "SeAx"), drop = FALSE]

print(round(identity_cor, 3))

# predicted identity from centroid correlation
identity_prediction <- tibble(sample_name = rownames(identity_cor), HH = identity_cor[, "HH"], MyLa = identity_cor[,
  "MyLa"], HuT78 = identity_cor[, "HuT78"], SeAx = identity_cor[, "SeAx"], predicted_cell_line = colnames(identity_cor)[max.col(identity_cor,
  ties.method = "first")]) |>
  dplyr::left_join(coldata |>
    tibble::rownames_to_column("sample_name") |>
    dplyr::select(sample_name, annotated_cell_line = cell_line_model), by = "sample_name")

print(identity_prediction)

write.csv(identity_prediction, file.path(qc_dir, "sample_identity_prediction.csv"), row.names = FALSE)

# centroid correlation heatmap
pheatmap::pheatmap(identity_cor, color = heatmap_blue, cluster_rows = TRUE, cluster_cols = FALSE, border_color = NA,
  fontsize = 7, fontsize_row = 7, fontsize_col = 7, angle_col = 0, filename = file.path(figures_dir,
    "sample_to_reference_correlation.pdf"), width = width_single, height = height_medium)

# reference-only deg validation
reference_samples <- c("HH-1", "HH-2", "M-1", "M-2", "78-1", "78-2", "S-1", "S-2")

dds_ref <- dds[, reference_samples]

# hh vs myla
dds_hh_myla <- dds_ref[, SummarizedExperiment::colData(dds_ref)$cell_line_model %in% c("HH", "MyLa")]

dds_hh_myla$cell_line_model <- droplevels(dds_hh_myla$cell_line_model)

dds_hh_myla <- DESeq2::DESeq(dds_hh_myla, fitType = "local")

res_hh_myla <- DESeq2::results(dds_hh_myla, contrast = c("cell_line_model", "MyLa", "HH"))

res_hh_myla_tbl <- as.data.frame(res_hh_myla) |>
  tibble::rownames_to_column("Geneid") |>
  dplyr::filter(!is.na(padj)) |>
  dplyr::arrange(padj)

myla_markers <- res_hh_myla_tbl |>
  dplyr::filter(padj < 0.01, log2FoldChange > 2) |>
  dplyr::arrange(dplyr::desc(log2FoldChange)) |>
  dplyr::slice_head(n = 15) |>
  dplyr::pull(Geneid)

hh_markers <- res_hh_myla_tbl |>
  dplyr::filter(padj < 0.01, log2FoldChange < -2) |>
  dplyr::arrange(log2FoldChange) |>
  dplyr::slice_head(n = 15) |>
  dplyr::pull(Geneid)

# hut78 vs seax
dds_hut_seax <- dds_ref[, SummarizedExperiment::colData(dds_ref)$cell_line_model %in% c("HuT78", "SeAx")]

dds_hut_seax$cell_line_model <- droplevels(dds_hut_seax$cell_line_model)

dds_hut_seax <- DESeq2::DESeq(dds_hut_seax, fitType = "local")

res_hut_seax <- DESeq2::results(dds_hut_seax, contrast = c("cell_line_model", "HuT78", "SeAx"))

res_hut_seax_tbl <- as.data.frame(res_hut_seax) |>
  tibble::rownames_to_column("Geneid") |>
  dplyr::filter(!is.na(padj)) |>
  dplyr::arrange(padj)

hut78_markers <- res_hut_seax_tbl |>
  dplyr::filter(padj < 0.01, log2FoldChange > 2) |>
  dplyr::arrange(dplyr::desc(log2FoldChange)) |>
  dplyr::slice_head(n = 15) |>
  dplyr::pull(Geneid)

seax_markers <- res_hut_seax_tbl |>
  dplyr::filter(padj < 0.01, log2FoldChange < -2) |>
  dplyr::arrange(log2FoldChange) |>
  dplyr::slice_head(n = 15) |>
  dplyr::pull(Geneid)

# combined validation signature
validation_genes <- unique(c(hh_markers, myla_markers, hut78_markers, seax_markers))

print(length(validation_genes))

validation_matrix <- vst_matrix[validation_genes, , drop = FALSE]

validation_matrix_z <- t(scale(t(validation_matrix)))

validation_matrix_z <- validation_matrix_z[apply(validation_matrix_z, 1, function(x) {
  all(is.finite(x))
}), , drop = FALSE]

# gene and sample annotations
gene_class <- tibble(Geneid = validation_genes, identity_marker = dplyr::case_when(Geneid %in% hh_markers ~
  "HH", Geneid %in% myla_markers ~ "MyLa", Geneid %in% hut78_markers ~ "HuT 78", Geneid %in% seax_markers ~
  "SeAx", TRUE ~ NA_character_))

gene_annotation <- gene_class |>
  tibble::column_to_rownames("Geneid")

gene_annotation <- gene_annotation[rownames(validation_matrix_z), , drop = FALSE]

validation_sample_annotation <- coldata |>
  dplyr::select(cell_line_model)

rownames(validation_sample_annotation) <- rownames(coldata)

validation_sample_annotation <- validation_sample_annotation[colnames(validation_matrix_z), , drop = FALSE]

validation_annotation_colors <- list(cell_line = cell_lines, identity_marker = cell_lines)

# deg validation heatmap
pheatmap::pheatmap(validation_matrix_z, color = heatmap_deg, breaks = seq(-1.5, 1.5, length.out = length(heatmap_deg) +
  1), legend_breaks = c(-1.5, 0, 1.5), annotation_col = validation_sample_annotation, annotation_row = gene_annotation,
  annotation_colors = validation_annotation_colors, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames = FALSE,
  border_color = NA, fontsize = 7, fontsize_col = 7, angle_col = 45, filename = file.path(figures_dir,
    "sample_identity_DEG_validation.pdf"), width = width_double, height = height_large)

# marker scores per sample
marker_score <- function(genes) {
  colMeans(validation_matrix_z[genes, , drop = FALSE])}

identity_scores <- tibble(sample_name = colnames(vst_matrix), HH_score = marker_score(hh_markers), MyLa_score = marker_score(myla_markers),
  HuT78_score = marker_score(hut78_markers), SeAx_score = marker_score(seax_markers))

score_matrix <- identity_scores |>
  dplyr::select(HH_score, MyLa_score, HuT78_score, SeAx_score) |>
  as.matrix()

identity_scores$predicted_identity <- c("HH", "MyLa", "HuT 78", "SeAx")[max.col(score_matrix, ties.method = "first")]

print(identity_scores, width = Inf)

write.csv(identity_scores, file.path(qc_dir, "sample_identity_marker_scores.csv"), row.names = FALSE)

write.csv(res_hh_myla_tbl, file.path(qc_dir, "HH_vs_MyLa_reference_DEGs.csv"), row.names = FALSE)

write.csv(res_hut_seax_tbl, file.path(qc_dir, "HuT78_vs_SeAx_reference_DEGs.csv"), row.names = FALSE)

