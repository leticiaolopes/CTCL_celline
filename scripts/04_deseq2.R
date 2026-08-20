# build the deseq2 dataset, calculate normalized and variance counts

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# directories
deseq2_dir <- file.path(paths$results, "deseq2")
qc_dir <- file.path(paths$qc, "deseq2")
figures_dir <- file.path(paths$figures, "drafts", "deseq2")

dir.create(deseq2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# load counts & metadata
count_matrix <- readRDS(file.path(paths$featurecounts, "count_matrix_raw.rds"))
coldata <- readRDS(file.path(paths$featurecounts, "coldata.rds"))

dim(count_matrix)
coldata

# model metadata
coldata <- coldata |>
  rownames_to_column("sample_name") |>
  mutate(original_cell_line = factor(original_cell_line, levels = c("HH", "MyLa", "HuT 78", "SeAx")),
    inferred_cell_line = factor(inferred_cell_line, levels = c("HH", "MyLa", "HuT 78", "SeAx")),
    analysis_cell_line = factor(analysis_cell_line, levels = c("HH", "MyLa", "HuT 78", "SeAx")),
    cell_line_model = case_when(analysis_cell_line == "HH" ~ "HH", analysis_cell_line == "MyLa" ~
      "MyLa", analysis_cell_line == "HuT 78" ~ "HuT78", analysis_cell_line == "SeAx" ~ "SeAx",
      TRUE ~ NA_character_), cell_line_model = factor(cell_line_model, levels = c("HH", "MyLa",
      "HuT78", "SeAx"))) |>
  column_to_rownames("sample_name")

coldata

# raw-count qc
library_size <- colSums(count_matrix)
detected_genes <- colSums(count_matrix > 0)

qc_samples <- tibble(sample_name = colnames(count_matrix), library_size = library_size, detected_genes = detected_genes) |>
  left_join(coldata |>
    as.data.frame() |>
    rownames_to_column("sample_name") |>
    dplyr::select(sample_name, analysis_sample_name, analysis_cell_line, original_cell_line), by = "sample_name")

qc_samples

write.csv(qc_samples, file.path(qc_dir, "sample_count_qc.csv"), row.names = FALSE)

ymax_library <- max(qc_samples$library_size)
ymax_library_plot <- ceiling(ymax_library/5e+06) * 5e+06

p_library <- ggplot(qc_samples, aes(x = analysis_sample_name, y = library_size, fill = analysis_cell_line)) +
  geom_col(width = 0.75) + scale_fill_cell_line() + scale_y_continuous(labels = scales::comma, breaks = seq(0,
  ymax_library_plot, by = 5e+06), limits = c(0, ymax_library_plot), expand = expansion(mult = c(0.015,
  0))) + labs(x = NULL, y = "Assigned fragments", fill = "Cell line") + theme_clean() + ggplot2::theme(axis.text.x = element_text(angle = 45,
  hjust = 1), legend.position = "top")

p_library

save_double(file.path(figures_dir, "library_size.pdf"), p_library, height = height_small)

ymax_genes <- max(qc_samples$detected_genes)
ymax_genes_plot <- ceiling(ymax_genes/5000) * 5000

p_genes <- ggplot(qc_samples, aes(x = analysis_sample_name, y = detected_genes, fill = analysis_cell_line)) +
  geom_col(width = 0.75) + scale_fill_cell_line() + scale_y_continuous(labels = scales::comma, breaks = seq(0,
  ymax_genes_plot, by = 5000), limits = c(0, ymax_genes_plot), expand = expansion(mult = c(0.015, 0))) +
  labs(x = NULL, y = "Genes detected", fill = "Cell line") + theme_clean() + ggplot2::theme(axis.text.x = element_text(angle = 45,
  hjust = 1), legend.position = "top")

p_genes

save_double(file.path(figures_dir, "detected_genes.pdf"), p_genes, height = height_small)

# filter low count
keep <- rowSums(count_matrix >= 10) >= 3
table(keep)

count_matrix_filtered <- count_matrix[keep, ]
dim(count_matrix)
dim(count_matrix_filtered)

filter_summary <- tibble(metric = c("Genes before filtering", "Genes after filtering", "Genes removed"),
  value = c(nrow(count_matrix), nrow(count_matrix_filtered), nrow(count_matrix) - nrow(count_matrix_filtered)))
filter_summary

write.csv(filter_summary, file.path(qc_dir, "gene_filter_summary.csv"), row.names = FALSE)

# deseq2 
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, colData = coldata, design = ~cell_line_model)
dds

# normalization
dds <- DESeq(dds)
sizeFactors(dds)
resultsNames(dds)

# export normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(data.frame(Geneid = rownames(normalized_counts), normalized_counts, check.names = FALSE), file.path(deseq2_dir,
  "normalized_counts.csv"), row.names = FALSE)

saveRDS(normalized_counts, file.path(deseq2_dir, "normalized_counts.rds"))

# vst
vsd <- vst(dds, blind = TRUE)

vst_matrix <- assay(vsd)

saveRDS(vsd, file.path(deseq2_dir, "vsd.rds"))

saveRDS(vst_matrix, file.path(deseq2_dir, "vst_matrix.rds"))

# pca
pca_data <- plotPCA(vsd, intgroup = "analysis_cell_line", returnData = TRUE)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_data <- pca_data |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_name") |>
  dplyr::select(sample_name, PC1, PC2) |>
  dplyr::left_join(coldata |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_name") |>
    dplyr::select(sample_name, analysis_sample_name, analysis_cell_line, original_cell_line), by = "sample_name")

pca_data

p_pca_original <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = original_cell_line, label = sample_name)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.45, color = black) + geom_text(vjust = -0.9, size = 2.2,
  show.legend = FALSE) + scale_fill_manual(values = cell_lines, breaks = names(cell_lines)) + labs(x = paste0("PC1 (",
  percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), fill = "Cell line") + theme_clean() +
  ggplot2::theme(legend.position = "top")

p_pca_original

save_single(file.path(figures_dir, "PCA_vst_original_identity.pdf"), p_pca_original, height = height_medium)

p_pca_inferred <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = analysis_cell_line, label = analysis_sample_name)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.45, color = black) + geom_text(vjust = -0.9, size = 2.2,
  show.legend = FALSE) + scale_fill_cell_line() + labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (",
  percent_var[2], "%)"), fill = "Cell line") + theme_clean() + ggplot2::theme(legend.position = "top")

p_pca_inferred

save_single(file.path(figures_dir, "PCA_vst_inferred_identity.pdf"), p_pca_inferred, height = height_medium)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = analysis_cell_line, label = analysis_sample_name)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.45, color = black) + geom_text(vjust = -0.9, size = 2.2,
  show.legend = FALSE) + scale_fill_cell_line() + labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (",
  percent_var[2], "%)"), fill = "Cell line") + theme_clean() + ggplot2::theme(legend.position = "top")

p_pca

save_single(file.path(figures_dir, "PCA_vst.pdf"), p_pca, height = height_medium)

write.csv(pca_data, file.path(qc_dir, "PCA_coordinates.csv"), row.names = FALSE)

# sample 2 sample correlation
sample_cor <- cor(vst_matrix, method = "pearson")

analysis_names <- coldata$analysis_sample_name
names(analysis_names) <- rownames(coldata)

rownames(sample_cor) <- analysis_names[rownames(sample_cor)]
colnames(sample_cor) <- analysis_names[colnames(sample_cor)]

annotation_col <- coldata |>
  as.data.frame() |>
  dplyr::select(analysis_cell_line)

rownames(annotation_col) <- coldata$analysis_sample_name

annotation_colors <- list(analysis_cell_line = cell_lines)

pheatmap::pheatmap(sample_cor, color = heatmap_blue, annotation_col = annotation_col, annotation_row = annotation_col,
  annotation_colors = annotation_colors, border_color = NA, fontsize = 7, fontsize_row = 6.5, fontsize_col = 6.5,
  angle_col = 45, filename = file.path(figures_dir, "sample_correlation_vst.pdf"), width = width_double,
  height = width_double)

# sample distance heatmap
sample_dist <- dist(t(vst_matrix))
sample_dist_matrix <- as.matrix(sample_dist)

rownames(sample_dist_matrix) <- analysis_names[rownames(sample_dist_matrix)]
colnames(sample_dist_matrix) <- analysis_names[colnames(sample_dist_matrix)]

pheatmap::pheatmap(sample_dist_matrix, color = heatmap_purple, annotation_col = annotation_col, annotation_row = annotation_col,
  annotation_colors = annotation_colors, border_color = NA, fontsize = 7, fontsize_row = 6.5, fontsize_col = 6.5,
  angle_col = 45, filename = file.path(figures_dir, "sample_distance_vst.pdf"), width = width_double,
  height = width_double)

# exports -> original vs identity
correlation_colors <- (grDevices::colorRampPalette(c("#F7F7F7", "#C6DED9", "#2A7F7F")))(101)

distance_colors <- (grDevices::colorRampPalette(c("#FFF7EC", "#F4A582", "#B55350")))(101)

sample_cor_raw <- cor(vst_matrix, method = "pearson")
sample_dist_raw <- as.matrix(dist(t(vst_matrix)))

original_names <- rownames(coldata)
inferred_names <- coldata$analysis_sample_name

original_groups <- dplyr::recode(as.character(coldata$original_cell_line), HuT78 = "HuT 78")
inferred_groups <- as.character(coldata$analysis_cell_line)

export_sample_heatmap <- function(matrix, display_names, groups, colors, filename, legend_breaks) {
  display_matrix <- matrix
  rownames(display_matrix) <- display_names
  colnames(display_matrix) <- display_names

  annotation <- data.frame(`Cell line` = groups, row.names = display_names, check.names = FALSE)

  pheatmap::pheatmap(display_matrix, color = colors, breaks = seq(min(legend_breaks), max(legend_breaks),
    length.out = length(colors) + 1), legend_breaks = legend_breaks, annotation_col = annotation,
    annotation_row = annotation, annotation_colors = list(`Cell line` = cell_lines), border_color = NA,
    fontsize = 7, fontsize_row = 6.5, fontsize_col = 6.5, angle_col = 45, filename = file.path(figures_dir,
      filename), width = 5.6, height = 5.2)}

cor_limits <- c(floor(min(sample_cor_raw) * 20)/20, 1)
dist_limits <- c(0, ceiling(max(sample_dist_raw)/25) * 25)

export_sample_heatmap(sample_cor_raw, original_names, original_groups, correlation_colors, "sample_correlation_original_identity.pdf",
  pretty(cor_limits, n = 4))

export_sample_heatmap(sample_cor_raw, inferred_names, inferred_groups, correlation_colors, "sample_correlation_inferred_identity.pdf",
  pretty(cor_limits, n = 4))

export_sample_heatmap(sample_dist_raw, original_names, original_groups, distance_colors, "sample_distance_original_identity.pdf",
  pretty(dist_limits, n = 4))

export_sample_heatmap(sample_dist_raw, inferred_names, inferred_groups, distance_colors, "sample_distance_inferred_identity.pdf",
  pretty(dist_limits, n = 4))

write.csv(sample_dist_matrix, file.path(qc_dir, "sample_distance_vst.csv"))

# deseq2 objects
saveRDS(dds, file.path(deseq2_dir, "dds.rds"))

saveRDS(coldata, file.path(deseq2_dir, "coldata_deseq2.rds"))

saveRDS(count_matrix_filtered, file.path(deseq2_dir, "count_matrix_filtered.rds"))

# summary
cat("\nDESeq2 preparation completed\n")
cat("----------------------------\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(dds), "\n")
cat("Design: ~ cell_line_model\n\n")

resultsNames(dds)
