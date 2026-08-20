# compare ctcl cell lines with healthy t cell references

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(stringr)
library(ExpressionAtlas)
library(DESeq2)
library(edgeR)

external_dir <- file.path(paths$data, "external_reference")
gse197067_dir <- file.path(external_dir, "GSE197067")
blueprint_dir <- file.path(external_dir, "BLUEPRINT")
results_external_dir <- file.path(paths$results, "external_reference")
figures_external_dir <- file.path(paths$figures, "drafts", "external_reference")
deseq2_dir <- file.path(paths$results, "deseq2")

dir.create(gse197067_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(blueprint_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(results_external_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(figures_external_dir, recursive = TRUE, showWarnings = FALSE)

# gse197067
gse197067_counts_file <- file.path(gse197067_dir, "GSE197067_HTSeq_counts.csv.gz")

if (!file.exists(gse197067_counts_file)) {
  download.file(url = paste0("https://www.ncbi.nlm.nih.gov/geo/download/", "?acc=GSE197067&file=GSE197067_HTSeq_counts.csv.gz&format=file"),
    destfile = gse197067_counts_file, mode = "wb")}

gse197067_counts_raw <- readr::read_csv(gse197067_counts_file, show_col_types = FALSE)

dim(gse197067_counts_raw)

colnames(gse197067_counts_raw)[1:10]

head(gse197067_counts_raw)

# select healthy pan t 0h samples
gse197067_0h_cols <- grep("_Pan_T_cells_0h_D[1-4]$", colnames(gse197067_counts_raw), value = TRUE)

print(gse197067_0h_cols)

gse197067_0h_counts <- gse197067_counts_raw |>
  dplyr::select(EnsemblID, dplyr::all_of(gse197067_0h_cols)) |>
  dplyr::mutate(Ensembl = stringr::str_remove(EnsemblID, "\\.\\d+$")) |>
  dplyr::select(Ensembl, dplyr::everything(), -EnsemblID)

gse197067_0h_counts <- gse197067_0h_counts |>
  dplyr::rename_with(~stringr::str_extract(.x, "D[1-4]$"), -Ensembl)

dim(gse197067_0h_counts)

colnames(gse197067_0h_counts)

head(gse197067_0h_counts)

sum(duplicated(gse197067_0h_counts$Ensembl))

saveRDS(gse197067_0h_counts, file.path(gse197067_dir, "GSE197067_PanT_0h_counts.rds"))

readr::write_csv(gse197067_0h_counts, file.path(gse197067_dir, "GSE197067_PanT_0h_counts.csv"))

# metadata
gse197067_0h_metadata <- tibble::tibble(sample = c("D1", "D2", "D3", "D4"), donor = c("D1", "D2", "D3",
  "D4"), cell_type = "Pan T", condition = "Healthy", activation = "Resting", time = "0h", dataset = "GSE197067")

readr::write_csv(gse197067_0h_metadata, file.path(gse197067_dir, "GSE197067_PanT_0h_metadata.csv"))

# qc
library_size <- gse197067_0h_counts |>
  dplyr::summarise(dplyr::across(-Ensembl, sum)) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "sample", values_to = "library_size")

print(library_size)

detected_genes <- gse197067_0h_counts |>
  dplyr::summarise(dplyr::across(-Ensembl, ~sum(.x > 0))) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "sample", values_to = "detected_genes")

print(detected_genes)

gse197067_0h_matrix <- gse197067_0h_counts |>
  tibble::column_to_rownames("Ensembl") |>
  as.matrix()

gse197067_0h_dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(gse197067_0h_matrix), colData = data.frame(row.names = colnames(gse197067_0h_matrix)),
  design = ~1)

gse197067_0h_vst <- DESeq2::vst(gse197067_0h_dds, blind = TRUE)

gse197067_0h_vst_matrix <- SummarizedExperiment::assay(gse197067_0h_vst)

# pca
gse197067_pca <- prcomp(t(gse197067_0h_vst_matrix))

gse197067_pca_data <- as.data.frame(gse197067_pca$x) |>
  tibble::rownames_to_column("sample")

percent_var <- (gse197067_pca$sdev^2/sum(gse197067_pca$sdev^2)) * 100

p_gse197067_pca <- ggplot2::ggplot(gse197067_pca_data, ggplot2::aes(x = PC1, y = PC2, label = sample)) +
  ggplot2::geom_point(shape = 21, size = 3.4, stroke = 0.5, color = black, fill = external_group_colors["Healthy Pan T"]) +
  ggrepel::geom_text_repel(size = 3) + ggplot2::labs(x = paste0("PC1 (", round(percent_var[1], 1),
  "%)"), y = paste0("PC2 (", round(percent_var[2], 1), "%)")) + theme_clean() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
  face = "bold"))

p_gse197067_pca

saveRDS(gse197067_0h_vst_matrix, file.path(gse197067_dir, "GSE197067_PanT_0h_vst_matrix.rds"))

save_single(file.path(figures_external_dir, "GSE197067_PanT_0h_PCA.pdf"), p_gse197067_pca, height = height_medium)

gse197067_cor <- cor(gse197067_0h_vst_matrix, method = "pearson")

correlation_heatmap_colors <- (grDevices::colorRampPalette(c("#F7F7F7", "#C6DED9", "#2A7F7F")))(101)

round(gse197067_cor, 3)

pheatmap::pheatmap(gse197067_cor, color = correlation_heatmap_colors, breaks = seq(min(gse197067_cor),
  1, length.out = 102), legend_breaks = pretty(c(min(gse197067_cor), 1), n = 4), cluster_rows = TRUE,
  cluster_cols = TRUE, border_color = NA, display_numbers = TRUE, number_format = "%.3f", fontsize = 8,
  fontsize_number = 7, main = "GSE197067 Pan T cells 0h")

gse197067_cor_no_diag <- gse197067_cor

diag(gse197067_cor_no_diag) <- NA

min(gse197067_cor_no_diag, na.rm = TRUE)

which(gse197067_cor_no_diag == min(gse197067_cor_no_diag, na.rm = TRUE), arr.ind = TRUE)

readr::write_csv(tibble::rownames_to_column(as.data.frame(gse197067_cor), "sample"), file.path(gse197067_dir,
  "GSE197067_PanT_0h_sample_correlations.csv"))

saveRDS(gse197067_cor, file.path(gse197067_dir, "GSE197067_PanT_0h_sample_correlations.rds"))

pdf(file.path(figures_external_dir, "GSE197067_PanT_0h_correlation_heatmap.pdf"), width = 4, height = 4,
  useDingbats = FALSE)

pheatmap::pheatmap(gse197067_cor, color = correlation_heatmap_colors, breaks = seq(min(gse197067_cor),
  1, length.out = 102), legend_breaks = pretty(c(min(gse197067_cor), 1), n = 4), cluster_rows = TRUE,
  cluster_cols = TRUE, border_color = NA, display_numbers = TRUE, number_format = "%.3f", fontsize = 8,
  fontsize_number = 7, main = "GSE197067 Pan T cells 0h")

dev.off()

gse197067_qc_summary <- dplyr::left_join(library_size, detected_genes, by = "sample")

readr::write_csv(gse197067_qc_summary, file.path(gse197067_dir, "GSE197067_PanT_0h_QC_summary.csv"))

# blueprint
blueprint_data <- ExpressionAtlas::getAtlasData("E-MTAB-3827")

blueprint_data

names(blueprint_data)
class(blueprint_data[[1]])
blueprint_data[[1]]

blueprint_rnaseq <- blueprint_data[["E-MTAB-3827"]][["rnaseq"]]

class(blueprint_rnaseq)

blueprint_rnaseq

colnames(SummarizedExperiment::colData(blueprint_rnaseq))

blueprint_metadata <- SummarizedExperiment::colData(blueprint_rnaseq) |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_id")

dim(blueprint_metadata)

head(blueprint_metadata)

SummarizedExperiment::assayNames(blueprint_rnaseq)

dim(blueprint_rnaseq)

blueprint_counts <- SummarizedExperiment::assay(blueprint_rnaseq, "counts")

blueprint_cell_types <- blueprint_metadata |>
  dplyr::count(cell_type, organism_part, disease, name = "n") |>
  dplyr::arrange(dplyr::desc(n), cell_type)

print(blueprint_cell_types)

blueprint_t_cells <- blueprint_metadata |>
  dplyr::filter(stringr::str_detect(stringr::str_to_lower(cell_type), "t cell")) |>
  dplyr::select(sample_id, AtlasAssayGroup, cell_type, organism_part, disease, age, sex, ethnic_group)

print(blueprint_t_cells)

blueprint_t_cell_summary <- blueprint_t_cells |>
  dplyr::count(cell_type, organism_part, disease, name = "n") |>
  dplyr::arrange(dplyr::desc(n), cell_type)

print(blueprint_t_cell_summary)

blueprint_cd4_metadata <- blueprint_metadata |>
  dplyr::filter(cell_type == "CD4-positive, alpha-beta T cell", organism_part == "venous blood", disease ==
    "normal") |>
  dplyr::select(sample_id, AtlasAssayGroup, cell_type, organism_part, disease, age, sex, ethnic_group)

print(blueprint_cd4_metadata)

blueprint_cd4_samples <- blueprint_cd4_metadata$sample_id

blueprint_cd4_counts <- blueprint_counts[, blueprint_cd4_samples, drop = FALSE]

dim(blueprint_cd4_counts)

blueprint_cd4_counts_tbl <- blueprint_cd4_counts |>
  as.data.frame() |>
  tibble::rownames_to_column("Ensembl")

head(blueprint_cd4_counts_tbl)

dim(blueprint_cd4_counts_tbl)

blueprint_cd4_library_size <- blueprint_cd4_counts_tbl |>
  dplyr::summarise(dplyr::across(-Ensembl, sum)) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "sample", values_to = "library_size")

print(blueprint_cd4_library_size)

blueprint_cd4_detected_genes <- blueprint_cd4_counts_tbl |>
  dplyr::summarise(dplyr::across(-Ensembl, ~sum(.x > 0))) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "sample", values_to = "detected_genes")

print(blueprint_cd4_detected_genes)

saveRDS(blueprint_cd4_counts_tbl, file.path(blueprint_dir, "BLUEPRINT_CD4_venous_blood_counts.rds"))

readr::write_csv(blueprint_cd4_counts_tbl, file.path(blueprint_dir, "BLUEPRINT_CD4_venous_blood_counts.csv"))

readr::write_csv(blueprint_cd4_metadata, file.path(blueprint_dir, "BLUEPRINT_CD4_venous_blood_metadata.csv"))

blueprint_cd4_dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(blueprint_cd4_counts), colData = data.frame(row.names = colnames(blueprint_cd4_counts)),
  design = ~1)

blueprint_cd4_vst <- DESeq2::vst(blueprint_cd4_dds, blind = TRUE)

blueprint_cd4_vst_matrix <- SummarizedExperiment::assay(blueprint_cd4_vst)

blueprint_cd4_pca <- prcomp(t(blueprint_cd4_vst_matrix))

blueprint_cd4_pca_data <- as.data.frame(blueprint_cd4_pca$x) |>
  tibble::rownames_to_column("sample")

blueprint_cd4_percent_var <- (blueprint_cd4_pca$sdev^2/sum(blueprint_cd4_pca$sdev^2)) * 100

p_blueprint_cd4_pca <- ggplot2::ggplot(blueprint_cd4_pca_data, ggplot2::aes(x = PC1, y = PC2, label = sample)) +
  ggplot2::geom_point(shape = 21, size = 3.4, stroke = 0.5, color = black, fill = external_group_colors["Healthy CD4"]) +
  ggrepel::geom_text_repel(size = 3) + ggplot2::labs(x = paste0("PC1 (", round(blueprint_cd4_percent_var[1],
  1), "%)"), y = paste0("PC2 (", round(blueprint_cd4_percent_var[2], 1), "%)")) + theme_clean() + ggplot2::labs(title = "BLUEPRINT healthy CD4 T cells") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

p_blueprint_cd4_pca

# qc
blueprint_cd4_cor <- cor(blueprint_cd4_vst_matrix, method = "pearson")

round(blueprint_cd4_cor, 3)

pheatmap::pheatmap(blueprint_cd4_cor, color = correlation_heatmap_colors, breaks = seq(min(blueprint_cd4_cor),
  1, length.out = 102), legend_breaks = pretty(c(min(blueprint_cd4_cor), 1), n = 4), cluster_rows = TRUE,
  cluster_cols = TRUE, border_color = NA, display_numbers = TRUE, number_format = "%.3f", fontsize = 8,
  fontsize_number = 7, main = "BLUEPRINT healthy CD4 T cells")

blueprint_cd4_cor_no_diag <- blueprint_cd4_cor

diag(blueprint_cd4_cor_no_diag) <- NA

min(blueprint_cd4_cor_no_diag, na.rm = TRUE)

which(blueprint_cd4_cor_no_diag == min(blueprint_cd4_cor_no_diag, na.rm = TRUE), arr.ind = TRUE)

blueprint_cd4_metadata |>
  dplyr::arrange(sample_id) |>
  dplyr::select(sample_id, AtlasAssayGroup, age, sex, ethnic_group, organism_part, cell_type) |>
  print()

blueprint_cd4_gene_var <- apply(blueprint_cd4_vst_matrix, 1, var)

blueprint_cd4_top_variable <- names(sort(blueprint_cd4_gene_var, decreasing = TRUE))[1:500]

blueprint_cd4_pca_top500 <- prcomp(t(blueprint_cd4_vst_matrix[blueprint_cd4_top_variable, , drop = FALSE]))

blueprint_cd4_pca_top500_data <- as.data.frame(blueprint_cd4_pca_top500$x) |>
  tibble::rownames_to_column("sample")

blueprint_cd4_top500_percent_var <- (blueprint_cd4_pca_top500$sdev^2/sum(blueprint_cd4_pca_top500$sdev^2)) *
  100

p_blueprint_cd4_pca_top500 <- ggplot2::ggplot(blueprint_cd4_pca_top500_data, ggplot2::aes(x = PC1, y = PC2,
  label = sample)) + ggplot2::geom_point(shape = 21, size = 3.4, stroke = 0.5, color = black, fill = external_group_colors["Healthy CD4"]) +
  ggrepel::geom_text_repel(size = 3) + ggplot2::labs(x = paste0("PC1 (", round(blueprint_cd4_top500_percent_var[1],
  1), "%)"), y = paste0("PC2 (", round(blueprint_cd4_top500_percent_var[2], 1), "%)")) + theme_clean() +
  ggplot2::labs(title = "BLUEPRINT healthy CD4 T cells - top 500 variable genes") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
  face = "bold"))

p_blueprint_cd4_pca_top500

saveRDS(blueprint_cd4_vst_matrix, file.path(blueprint_dir, "BLUEPRINT_CD4_venous_blood_vst_matrix.rds"))

saveRDS(blueprint_cd4_cor, file.path(blueprint_dir, "BLUEPRINT_CD4_venous_blood_sample_correlations.rds"))

readr::write_csv(tibble::rownames_to_column(as.data.frame(blueprint_cd4_cor), "sample"), file.path(blueprint_dir,
  "BLUEPRINT_CD4_venous_blood_sample_correlations.csv"))

blueprint_cd4_qc_summary <- dplyr::left_join(blueprint_cd4_library_size, blueprint_cd4_detected_genes,
  by = "sample")

readr::write_csv(blueprint_cd4_qc_summary, file.path(blueprint_dir, "BLUEPRINT_CD4_venous_blood_QC_summary.csv"))

save_single(file.path(figures_external_dir, "BLUEPRINT_CD4_PCA.pdf"), p_blueprint_cd4_pca, height = height_medium)

save_single(file.path(figures_external_dir, "BLUEPRINT_CD4_PCA_top500_variable_genes.pdf"), p_blueprint_cd4_pca_top500,
  height = height_medium)

pdf(file.path(figures_external_dir, "BLUEPRINT_CD4_correlation_heatmap.pdf"), width = 5, height = 5,
  useDingbats = FALSE)

pheatmap::pheatmap(blueprint_cd4_cor, color = correlation_heatmap_colors, breaks = seq(min(blueprint_cd4_cor),
  1, length.out = 102), legend_breaks = pretty(c(min(blueprint_cd4_cor), 1), n = 4), cluster_rows = TRUE,
  cluster_cols = TRUE, border_color = NA, display_numbers = TRUE, number_format = "%.3f", fontsize = 8,
  fontsize_number = 7, main = "BLUEPRINT healthy CD4 T cells")

dev.off()

# cross-dataset expression comparison
dds <- readRDS(file.path(deseq2_dir, "dds.rds"))
dds

# cell-line raw counts from fitted deseq2 object
cell_line_counts <- DESeq2::counts(dds, normalized = FALSE)

rownames(cell_line_counts) <- stringr::str_remove(rownames(cell_line_counts), "\\.\\d+$")

cell_line_metadata <- SummarizedExperiment::colData(dds) |>
  as.data.frame() |>
  tibble::rownames_to_column("original_sample_name")

cell_line_names <- cell_line_metadata$analysis_sample_name[match(colnames(cell_line_counts), cell_line_metadata$original_sample_name)]

colnames(cell_line_counts) <- cell_line_names
dim(cell_line_counts)
colnames(cell_line_counts)

# gse197067 counts
gse197067_counts_matrix <- gse197067_0h_counts |>
  tibble::column_to_rownames("Ensembl") |>
  as.matrix()

# blueprint counts
blueprint_counts_matrix <- blueprint_cd4_counts_tbl |>
  tibble::column_to_rownames("Ensembl") |>
  as.matrix()

common_genes <- Reduce(intersect, list(rownames(cell_line_counts), rownames(gse197067_counts_matrix),
  rownames(blueprint_counts_matrix)))

length(common_genes)
c(cell_lines = nrow(cell_line_counts), GSE197067 = nrow(gse197067_counts_matrix), BLUEPRINT = nrow(blueprint_counts_matrix),
  shared = length(common_genes))

cell_line_cpm <- edgeR::cpm(cell_line_counts)

gse197067_cpm <- edgeR::cpm(gse197067_counts_matrix)

blueprint_cpm <- edgeR::cpm(blueprint_counts_matrix)

cell_line_cpm <- cell_line_cpm[common_genes, , drop = FALSE]

gse197067_cpm <- gse197067_cpm[common_genes, , drop = FALSE]

blueprint_cpm <- blueprint_cpm[common_genes, , drop = FALSE]

expressed_cell_lines <- rowSums(cell_line_cpm >= 1) >= 2

expressed_gse197067 <- rowSums(gse197067_cpm >= 1) >= 2

expressed_blueprint <- rowSums(blueprint_cpm >= 1) >= 2

comparison_genes <- common_genes[expressed_cell_lines & expressed_gse197067 & expressed_blueprint]

length(comparison_genes)

comparison_gene_summary <- tibble::tibble(metric = c("Cell-line genes", "GSE197067 genes", "BLUEPRINT genes",
  "Shared genes", "Shared expressed genes"), n = c(nrow(cell_line_counts), nrow(gse197067_counts_matrix),
  nrow(blueprint_counts_matrix), length(common_genes), length(comparison_genes)))

print(comparison_gene_summary)

cell_line_dge <- edgeR::DGEList(counts = cell_line_counts)
cell_line_dge <- edgeR::calcNormFactors(cell_line_dge, method = "TMM")
cell_line_logcpm <- edgeR::cpm(cell_line_dge, log = TRUE, prior.count = 1)
gse197067_dge <- edgeR::DGEList(counts = gse197067_counts_matrix)
gse197067_dge <- edgeR::calcNormFactors(gse197067_dge, method = "TMM")
gse197067_logcpm <- edgeR::cpm(gse197067_dge, log = TRUE, prior.count = 1)
blueprint_dge <- edgeR::DGEList(counts = blueprint_counts_matrix)
blueprint_dge <- edgeR::calcNormFactors(blueprint_dge, method = "TMM")
blueprint_logcpm <- edgeR::cpm(blueprint_dge, log = TRUE, prior.count = 1)

cell_line_logcpm <- cell_line_logcpm[comparison_genes, , drop = FALSE]
gse197067_logcpm <- gse197067_logcpm[comparison_genes, , drop = FALSE]
blueprint_logcpm <- blueprint_logcpm[comparison_genes, , drop = FALSE]

colnames(gse197067_logcpm) <- paste0("GSE_", colnames(gse197067_logcpm))
colnames(blueprint_logcpm) <- paste0("BP_", colnames(blueprint_logcpm))

external_comparison_matrix <- cbind(cell_line_logcpm, gse197067_logcpm, blueprint_logcpm)

dim(external_comparison_matrix)

external_cor <- cor(external_comparison_matrix, method = "spearman")

round(external_cor[13:16, 17:24], 3)

cell_line_annotation <- cell_line_metadata |>
  dplyr::select(analysis_sample_name, analysis_cell_line) |>
  dplyr::distinct()

cell_line_group <- cell_line_annotation$analysis_cell_line[match(colnames(cell_line_logcpm), cell_line_annotation$analysis_sample_name)]

external_annotation <- data.frame(row.names = colnames(external_comparison_matrix), Dataset = c(rep("Cell lines",
  12), rep("GSE197067", 4), rep("BLUEPRINT", 8)), Group = c(cell_line_group, rep("Healthy Pan T", 4),
  rep("Healthy CD4", 8)))

pheatmap::pheatmap(external_cor, color = correlation_heatmap_colors, breaks = seq(min(external_cor),
  1, length.out = 102), legend_breaks = pretty(c(min(external_cor), 1), n = 4), annotation_col = external_annotation,
  annotation_row = external_annotation, annotation_colors = list(Dataset = c(`Cell lines` = "#777777",
    GSE197067 = unname(external_group_colors["Healthy Pan T"]), BLUEPRINT = unname(external_group_colors["Healthy CD4"])),
    Group = external_group_colors), cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA,
  show_colnames = TRUE, show_rownames = TRUE, fontsize = 6, main = "Cross-dataset transcriptomic similarity",
  filename = file.path(figures_external_dir, "cross_dataset_correlation_heatmap.pdf"), width = 6.2,
  height = 6.2)

cell_line_samples <- colnames(cell_line_logcpm)
gse_samples <- colnames(gse197067_logcpm)
blueprint_samples <- colnames(blueprint_logcpm)

cell_line_similarity <- tibble::tibble(sample = cell_line_samples, mean_cor_GSE197067 = rowMeans(external_cor[cell_line_samples,
  gse_samples, drop = FALSE]), mean_cor_BLUEPRINT = rowMeans(external_cor[cell_line_samples, blueprint_samples,
  drop = FALSE]))

cell_line_similarity <- cell_line_similarity |>
  dplyr::left_join(cell_line_annotation, by = c(sample = "analysis_sample_name"))

cell_line_similarity |>
  dplyr::arrange(dplyr::desc(mean_cor_BLUEPRINT))

cell_line_similarity <- cell_line_similarity |>
  dplyr::mutate(mean_cor_healthy = rowMeans(cbind(mean_cor_GSE197067, mean_cor_BLUEPRINT)))

cell_line_similarity |>
  dplyr::arrange(dplyr::desc(mean_cor_healthy))

cell_line_similarity_summary <- cell_line_similarity |>
  dplyr::group_by(analysis_cell_line) |>
  dplyr::summarise(mean_GSE197067 = mean(mean_cor_GSE197067), sd_GSE197067 = sd(mean_cor_GSE197067),
    mean_BLUEPRINT = mean(mean_cor_BLUEPRINT), sd_BLUEPRINT = sd(mean_cor_BLUEPRINT), mean_healthy = mean(mean_cor_healthy),
    sd_healthy = sd(mean_cor_healthy), .groups = "drop")

cell_line_similarity_summary |>
  dplyr::arrange(dplyr::desc(mean_healthy))

similarity_plot_data <- cell_line_similarity |>
  tidyr::pivot_longer(cols = c(mean_cor_GSE197067, mean_cor_BLUEPRINT), names_to = "reference", values_to = "correlation")

similarity_plot_data <- similarity_plot_data |>
  dplyr::mutate(reference = dplyr::recode(reference, mean_cor_GSE197067 = "Healthy Pan T", mean_cor_BLUEPRINT = "Healthy CD4"))

p_external_similarity <- ggplot2::ggplot(similarity_plot_data, ggplot2::aes(x = analysis_cell_line, y = correlation,
  color = analysis_cell_line)) + ggplot2::geom_jitter(width = 0.08, size = 2, alpha = 0.8) + ggplot2::stat_summary(fun = mean,
  geom = "crossbar", width = 0.35, linewidth = 0.4, color = black) + ggplot2::facet_wrap(~reference) +
  ggplot2::scale_color_manual(values = cell_lines) + ggplot2::labs(x = NULL, y = "Mean Spearman correlation") +
  theme_clean() + ggplot2::theme(legend.position = "none")

save_double(file.path(figures_external_dir, "cell_line_similarity_to_healthy_references.pdf"), p_external_similarity,
  height = height_medium)

external_rank_matrix <- apply(external_comparison_matrix, 2, rank, ties.method = "average")

saveRDS(list(external_rank_matrix = external_rank_matrix, gse_samples = gse_samples, blueprint_samples = blueprint_samples),
  file.path(results_external_dir, "external_rank_reference.rds"))

dim(external_rank_matrix)
external_rank_var <- apply(external_rank_matrix, 1, var)
external_rank_top2000 <- names(sort(external_rank_var, decreasing = TRUE))[1:2000]

external_pca <- prcomp(t(external_rank_matrix[external_rank_top2000, , drop = FALSE]), center = TRUE,
  scale. = TRUE)

external_pca_data <- as.data.frame(external_pca$x) |>
  tibble::rownames_to_column("sample")

external_pca_percent_var <- (external_pca$sdev^2/sum(external_pca$sdev^2)) * 100

external_pca_data <- external_pca_data |>
  dplyr::left_join(tibble::rownames_to_column(external_annotation, "sample"), by = "sample")

external_pca_data |>
  dplyr::select(sample, Dataset, Group, PC1, PC2)

cell_line_group <- as.character(cell_line_annotation$analysis_cell_line[match(colnames(cell_line_logcpm),
  cell_line_annotation$analysis_sample_name)])

external_annotation <- data.frame(row.names = colnames(external_comparison_matrix), Dataset = c(rep("Cell lines",
  12), rep("GSE197067", 4), rep("BLUEPRINT", 8)), Group = c(cell_line_group, rep("Healthy Pan T", 4),
  rep("Healthy CD4", 8)))

external_pca_data <- as.data.frame(external_pca$x) |>
  tibble::rownames_to_column("sample") |>
  dplyr::left_join(tibble::rownames_to_column(external_annotation, "sample"), by = "sample")

p_external_pca <- ggplot2::ggplot(external_pca_data, ggplot2::aes(x = PC1, y = PC2, color = Group, shape = Dataset)) +
  ggplot2::geom_point(size = 3, alpha = 0.9) + ggplot2::labs(x = paste0("PC1 (", round(external_pca_percent_var[1],
  1), "%)"), y = paste0("PC2 (", round(external_pca_percent_var[2], 1), "%)"), color = NULL, shape = NULL) +
  theme_clean() + ggplot2::theme(legend.position = "right")
p_external_pca

external_pca_all <- prcomp(t(external_rank_matrix), center = TRUE, scale. = TRUE)
external_pca_all_percent_var <- (external_pca_all$sdev^2/sum(external_pca_all$sdev^2)) * 100
external_pca_all_data <- as.data.frame(external_pca_all$x) |>
  tibble::rownames_to_column("sample") |>
  dplyr::left_join(tibble::rownames_to_column(external_annotation, "sample"), by = "sample")
p_external_pca_all <- ggplot2::ggplot(external_pca_all_data, ggplot2::aes(x = PC1, y = PC2, color = Group,
  shape = Dataset)) + ggplot2::geom_point(size = 3, alpha = 0.9) + ggplot2::labs(x = paste0("PC1 (",
  round(external_pca_all_percent_var[1], 1), "%)"), y = paste0("PC2 (", round(external_pca_all_percent_var[2],
  1), "%)"), color = NULL, shape = NULL) + theme_clean()
p_external_pca_all

external_spearman_dist <- as.dist(1 - external_cor)
external_mds <- cmdscale(external_spearman_dist, k = 2, eig = TRUE)
external_mds_data <- as.data.frame(external_mds$points) |>
  tibble::rownames_to_column("sample")
colnames(external_mds_data)[2:3] <- c("MDS1", "MDS2")

external_mds_data <- external_mds_data |>
  dplyr::left_join(tibble::rownames_to_column(external_annotation, "sample"), by = "sample")

p_external_mds <- ggplot2::ggplot(external_mds_data, ggplot2::aes(x = MDS1, y = MDS2, color = Group,
  shape = Dataset)) + ggplot2::geom_point(size = 3, alpha = 0.9) + ggplot2::labs(x = "MDS1", y = "MDS2",
  color = NULL, shape = NULL) + theme_clean()
p_external_mds

gene_modules <- readr::read_csv(file.path(paths$results, "functional_enrichment", "gene_module_assignments.csv"),
  show_col_types = FALSE)
colnames(gene_modules)
head(gene_modules)

gene_modules |>
  dplyr::count(module)

module_genes_shared <- intersect(gene_modules$Ensembl, rownames(external_comparison_matrix))
length(module_genes_shared)

gene_modules |>
  dplyr::filter(Ensembl %in% module_genes_shared) |>
  dplyr::count(module)

module_map <- gene_modules |>
  dplyr::filter(Ensembl %in% module_genes_shared) |>
  dplyr::select(Ensembl, module) |>
  dplyr::distinct()

exists("external_rank_matrix")

module_rank_matrix <- external_rank_matrix[module_map$Ensembl, , drop = FALSE]
module_scores <- module_map |>
  dplyr::group_by(module) |>
  dplyr::group_modify(~tibble::tibble(sample = colnames(module_rank_matrix), score = colMeans(module_rank_matrix[.x$Ensembl,
    , drop = FALSE]))) |>
  dplyr::ungroup()

module_scores <- module_scores |>
  dplyr::left_join(tibble::rownames_to_column(external_annotation, "sample"), by = "sample")

head(module_scores)

module_scores <- module_scores |>
  dplyr::group_by(module) |>
  dplyr::mutate(score_z = as.numeric(scale(score))) |>
  dplyr::ungroup()

module_scores_summary <- module_scores |>
  dplyr::group_by(module, Group) |>
  dplyr::summarise(mean_score_z = mean(score_z), sd_score_z = sd(score_z), .groups = "drop")
module_scores_summary

module_heatmap_df <- module_scores_summary |>
  dplyr::select(module, Group, mean_score_z) |>
  tidyr::pivot_wider(names_from = Group, values_from = mean_score_z)

module_heatmap_matrix <- module_heatmap_df |>
  tibble::column_to_rownames("module") |>
  as.matrix()

module_heatmap_matrix <- module_heatmap_matrix[, c("HH", "MyLa", "HuT 78", "SeAx", "Healthy Pan T", "Healthy CD4"),
  drop = FALSE]

pheatmap::pheatmap(module_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap_deg,
  border_color = NA, fontsize = 8, main = "Gene-module CTCL cell lines vs healthy T-cell references")

module_coverage <- gene_modules |>
  dplyr::count(module, name = "total_genes") |>
  dplyr::left_join(gene_modules |>
    dplyr::filter(Ensembl %in% module_genes_shared) |>
    dplyr::count(module, name = "shared_genes"), by = "module") |>
  dplyr::mutate(coverage = 100 * shared_genes/total_genes)
module_coverage

module_heatmap_matrix <- module_heatmap_matrix[c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5",
  "Module 6"), c("HH", "MyLa", "HuT 78", "SeAx", "Healthy Pan T", "Healthy CD4"), drop = FALSE]

pheatmap::pheatmap(module_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap_deg,
  breaks = seq(-2, 2, length.out = length(heatmap_deg) + 1), legend_breaks = c(-2, 0, 2), border_color = NA,
  fontsize = 8, fontsize_row = 8, fontsize_col = 8, angle_col = 45, main = "Gene-module activity: CTCL vs healthy T-cell references",
  filename = file.path(figures_external_dir, "gene_module_activity_healthy_references.pdf"), width = 6.2,
  height = 3.8)

module_go_reduced <- readr::read_csv(file.path(paths$results, "functional_enrichment", "GO_BP_by_gene_module_reduced.csv"),
  show_col_types = FALSE)
colnames(module_go_reduced)
head(module_go_reduced)

module_go_summary <- module_go_reduced |>
  dplyr::select(module, representative_ID, representative_term, p.adjust, Count) |>
  dplyr::distinct(module, representative_ID, .keep_all = TRUE) |>
  dplyr::arrange(module, p.adjust)

module_go_summary |>
  dplyr::group_by(module) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup() |>
  dplyr::select(module, representative_term, p.adjust, Count) |>
  print(n = 37)

module_coverage <- gene_modules |>
  dplyr::count(module, name = "total_genes") |>
  dplyr::left_join(gene_modules |>
    dplyr::filter(Ensembl %in% module_genes_shared) |>
    dplyr::count(module, name = "shared_genes"), by = "module") |>
  dplyr::mutate(coverage_percent = 100 * shared_genes/total_genes)
module_coverage

module_scores <- module_scores |>
  dplyr::mutate(Group = factor(Group, levels = c("HH", "MyLa", "HuT 78", "SeAx", "Healthy Pan T", "Healthy CD4")))

p_module_scores <- ggplot2::ggplot(module_scores, ggplot2::aes(x = Group, y = score_z, color = Group)) +
  ggplot2::geom_jitter(width = 0.08, size = 1.8, alpha = 0.75) + ggplot2::stat_summary(fun = mean,
  geom = "point", size = 3) + ggplot2::stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
  geom = "errorbar", width = 0.18, linewidth = 0.4) + ggplot2::facet_wrap(~module, ncol = 3, scales = "free_y") +
  ggplot2::scale_color_manual(values = external_group_colors) + ggplot2::geom_hline(yintercept = 0,
  linewidth = 0.3, color = light_grey) + ggplot2::labs(x = NULL, y = "Module score (Z-score)", color = NULL) +
  theme_clean() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45,
  hjust = 1), strip.text = ggplot2::element_text(face = "bold"))
p_module_scores

save_double(file.path(figures_external_dir, "gene_module_scores_CTCL_vs_healthy.pdf"), p_module_scores,
  height = height_large)

module_labels <- c(`Module 1` = "Module 1\nAdhesion / morphogenesis", `Module 2` = "Module 2\nT-cell activation / cytokine signaling",
  `Module 3` = "Module 3\nDevelopmental program", `Module 4` = "Module 4\nIon transport", `Module 5` = "Module 5\nAdaptive immunity / T-cell migration",
  `Module 6` = "Module 6\nImmune proliferation")

p_module_scores <- ggplot2::ggplot(module_scores, ggplot2::aes(x = Group, y = score_z, color = Group)) +
  ggplot2::geom_jitter(width = 0.08, size = 1.8, alpha = 0.75) + ggplot2::stat_summary(fun = mean,
  geom = "point", size = 3) + ggplot2::stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
  geom = "errorbar", width = 0.18, linewidth = 0.4) + ggplot2::facet_wrap(~module, ncol = 3, scales = "free_y",
  labeller = ggplot2::as_labeller(module_labels)) + ggplot2::scale_color_manual(values = external_group_colors) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = light_grey) + ggplot2::labs(x = NULL,
  y = "Module score (Z-score)", color = NULL) + theme_clean() + ggplot2::theme(legend.position = "none",
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), strip.text = ggplot2::element_text(face = "bold"))

save_double(file.path(figures_external_dir, "gene_module_scores_CTCL_vs_healthy.pdf"), p_module_scores,
  height = height_large)

module_scores_export <- module_scores_summary |>
  dplyr::left_join(module_coverage, by = "module")
readr::write_csv(module_scores_export, file.path(results_external_dir, "gene_module_scores_CTCL_vs_healthy_summary.csv"))

external_group_levels <- c("HH", "MyLa", "HuT 78", "SeAx", "Healthy Pan T", "Healthy CD4")
external_pca_data <- external_pca_data |>
  dplyr::mutate(Group = factor(Group, levels = external_group_levels))
external_mds_data <- external_mds_data |>
  dplyr::mutate(Group = factor(Group, levels = external_group_levels))
p_external_pca <- ggplot2::ggplot(external_pca_data, ggplot2::aes(x = PC1, y = PC2, color = Group, shape = Dataset)) +
  ggplot2::geom_point(size = 3, alpha = 0.9) + ggplot2::scale_color_manual(values = external_group_colors) +
  ggplot2::scale_shape_manual(values = c(`Cell lines` = 17, GSE197067 = 15, BLUEPRINT = 16)) + ggplot2::labs(x = paste0("PC1 (",
  round(external_pca_percent_var[1], 1), "%)"), y = paste0("PC2 (", round(external_pca_percent_var[2],
  1), "%)"), color = NULL, shape = NULL) + theme_clean() + ggplot2::theme(legend.position = "right")
save_double(file.path(figures_external_dir, "cross_dataset_rank_PCA.pdf"), p_external_pca, height = height_medium)

p_external_mds <- ggplot2::ggplot(external_mds_data, ggplot2::aes(x = MDS1, y = MDS2, color = Group,
  shape = Dataset)) + ggplot2::geom_point(size = 3, alpha = 0.9) + ggplot2::scale_color_manual(values = external_group_colors) +
  ggplot2::scale_shape_manual(values = c(`Cell lines` = 17, GSE197067 = 15, BLUEPRINT = 16)) + ggplot2::labs(x = "MDS1",
  y = "MDS2", color = NULL, shape = NULL) + theme_clean() + ggplot2::theme(legend.position = "right")
save_double(file.path(figures_external_dir, "cross_dataset_spearman_MDS.pdf"), p_external_mds, height = height_medium)

module_scores <- module_scores |>
  dplyr::mutate(Group = factor(Group, levels = external_group_levels))
module_labels <- c(`Module 1` = "Module 1\nAdhesion / morphogenesis", `Module 2` = "Module 2\nT-cell activation / cytokine signaling",
  `Module 3` = "Module 3\nDevelopmental program", `Module 4` = "Module 4\nIon transport", `Module 5` = "Module 5\nAdaptive immunity / T-cell migration",
  `Module 6` = "Module 6\nImmune proliferation")

p_module_scores <- ggplot2::ggplot(module_scores, ggplot2::aes(x = Group, y = score_z, color = Group)) +
  ggplot2::geom_jitter(width = 0.08, size = 1.8, alpha = 0.75) + ggplot2::stat_summary(fun = mean,
  geom = "point", size = 3) + ggplot2::stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
  geom = "errorbar", width = 0.18, linewidth = 0.4) + ggplot2::facet_wrap(~module, ncol = 3, scales = "free_y",
  labeller = ggplot2::as_labeller(module_labels)) + ggplot2::scale_color_manual(values = external_group_colors) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = light_grey) + ggplot2::labs(x = NULL,
  y = "Module score (Z-score)", color = NULL) + theme_clean() + ggplot2::theme(legend.position = "none",
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), strip.text = ggplot2::element_text(face = "bold"))

save_double(file.path(figures_external_dir, "gene_module_scores_CTCL_vs_healthy.pdf"), p_module_scores,
  height = height_large)
