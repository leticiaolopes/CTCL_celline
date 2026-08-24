# rebuild supplementary panels from saved rds & csv files

source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(edgeR)})

gene_scope <- getOption("ctcl.gene_scope", "all")
if (!gene_scope %in% c("all", "protein_coding")) {
  stop("Unsupported ctcl.gene_scope: ", gene_scope)
}
analysis_suffix <- if (gene_scope == "protein_coding") "_protein_coding" else ""

ensure_file <- function(path, fallback_script) {
  if (!file.exists(path)) {
    message("missing object: ", path, "; running ", fallback_script)
    source(fallback_script, local = new.env(parent = globalenv()))}
  if (!file.exists(path))
    stop("required object is still missing: ", path)
  path}

deseq2_dir <- file.path(paths$results, paste0("deseq2", analysis_suffix))
de_dir <- file.path(paths$results, paste0("differential_expression", analysis_suffix))
functional_dir <- file.path(paths$results, paste0("functional_enrichment", analysis_suffix))
external_dir <- file.path(paths$results, paste0("external_reference", analysis_suffix))
external_data_dir <- file.path(paths$data, "external_reference")

figures_deseq2_dir <- file.path(paths$figures, "drafts", paste0("deseq2", analysis_suffix))
figures_de_dir <- file.path(figures_deseq2_dir, "differential_expression")
figures_external_dir <- file.path(paths$figures, "drafts", paste0("external_reference", analysis_suffix))
figures_functional_dir <- file.path(paths$figures, "drafts", paste0("functional_enrichment", analysis_suffix))

dir.create(figures_de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_external_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_functional_dir, recursive = TRUE, showWarnings = FALSE)

dds_file <- ensure_file(file.path(deseq2_dir, "dds.rds"), "scripts/04_deseq2.R")
vst_file <- ensure_file(file.path(deseq2_dir, "vst_matrix.rds"), "scripts/04_deseq2.R")
de_file <- ensure_file(file.path(de_dir, "all_shrunk_results_annotated.rds"), "scripts/05_differential_expression.R")

dds <- readRDS(dds_file)
vst_matrix <- readRDS(vst_file)
all_de <- readRDS(de_file)

coldata <- SummarizedExperiment::colData(dds) |>
  as.data.frame() |>
  tibble::rownames_to_column("original_sample_name")

# s1g-h: identity signature heatmaps 
identity_figures_dir <- file.path(figures_deseq2_dir, "sample_identity")
dir.create(identity_figures_dir, recursive = TRUE, showWarnings = FALSE)

signature_difference <- function(group_a, group_b, n = 50) {
  difference <- rowMeans(vst_matrix[, group_a, drop = FALSE]) - rowMeans(vst_matrix[, group_b, drop = FALSE])
  names(sort(abs(difference), decreasing = TRUE))[seq_len(n)]}

identity_genes <- unique(c(signature_difference(c("HH-1", "HH-2"), c("M-1", "M-2")), signature_difference(c("78-1",
  "78-2"), c("S-1", "S-2"))))

identity_z <- t(scale(t(vst_matrix[identity_genes, , drop = FALSE])))
identity_z <- identity_z[apply(identity_z, 1, function(x) all(is.finite(x))), , drop = FALSE]

export_identity_heatmap <- function(display_names, groups, filename) {
  matrix <- identity_z
  colnames(matrix) <- display_names

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(`Cell line` = groups, col = list(`Cell line` = cell_lines),
    show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

  heatmap <- ComplexHeatmap::Heatmap(matrix, name = "Z-score", col = circlize::colorRamp2(c(-1.5, 0,
    1.5), c("#3F648C", "#F7F7F7", "#B55350")), top_annotation = top_annotation, cluster_rows = TRUE,
    cluster_columns = TRUE, show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 7), border = FALSE, heatmap_legend_param = list(at = c(-1.5,
      0, 1.5), labels = c("-1.5", "0", "1.5"), legend_height = grid::unit(22, "mm"), legend_width = grid::unit(5,
      "mm"), border = black, title_position = "topleft"))

  grDevices::pdf(filename, width = 5.6, height = 5.2, useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "right")
  grDevices::dev.off()}

original_groups <- dplyr::recode(as.character(coldata$original_cell_line), HuT78 = "HuT 78")

export_identity_heatmap(colnames(identity_z), original_groups, file.path(identity_figures_dir, "internal_identity_signature_original.pdf"))

export_identity_heatmap(coldata$analysis_sample_name, as.character(coldata$analysis_cell_line), file.path(identity_figures_dir,
  "internal_identity_signature_inferred.pdf"))

# s1c-f
sample_cor <- cor(vst_matrix, method = "pearson")
sample_dist <- as.matrix(dist(t(vst_matrix)))

export_sample_matrix <- function(matrix, display_names, groups, filename, legend_title, palette, legend_values,
  legend_labels) {
  display_matrix <- matrix
  rownames(display_matrix) <- display_names
  colnames(display_matrix) <- display_names

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(`Cell line` = groups, col = list(`Cell line` = cell_lines),
    show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

  left_annotation <- ComplexHeatmap::rowAnnotation(`Cell line` = groups, col = list(`Cell line` = cell_lines),
    show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

  heatmap <- ComplexHeatmap::Heatmap(display_matrix, name = legend_title, col = circlize::colorRamp2(legend_values,
    palette), top_annotation = top_annotation, left_annotation = left_annotation, cluster_rows = TRUE,
    cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 6.5),
    column_names_gp = grid::gpar(fontsize = 6.5), column_names_rot = 45, border = FALSE, heatmap_legend_param = list(at = legend_values,
      labels = legend_labels, legend_height = grid::unit(22, "mm"), legend_width = grid::unit(5,
        "mm"), border = black, title_position = "topleft"))

  grDevices::pdf(filename, width = 5.6, height = 5.2, useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()}

cor_min <- floor(min(sample_cor) * 100)/100
cor_mid <- (cor_min + 1)/2
dist_max <- ceiling(max(sample_dist)/25) * 25

for (identity_version in c("original", "inferred")) {
  if (identity_version == "original") {
    display_names <- colnames(vst_matrix)
    groups <- original_groups
  } else {
    display_names <- coldata$analysis_sample_name
    groups <- as.character(coldata$analysis_cell_line)}

  export_sample_matrix(sample_cor, display_names, groups, file.path(figures_deseq2_dir,
    paste0("sample_correlation_", identity_version, "_identity.pdf")), "Pearson r", c("#F7F7F7",
    "#C6DED9", "#2A7F7F"), c(cor_min, cor_mid, 1), format(c(cor_min, cor_mid, 1), digits = 2))

  export_sample_matrix(sample_dist, display_names, groups, file.path(figures_deseq2_dir,
    paste0("sample_distance_", identity_version, "_identity.pdf")), "Distance", c("#FFF7EC", "#F4A582",
    "#B55350"), c(0, dist_max/2, dist_max), format(c(0, dist_max/2, dist_max), trim = TRUE))}

# s2a
top_gene_ids <- all_de |>
  dplyr::filter(DEG_class %in% c("Up", "Down")) |>
  dplyr::group_by(comparison, DEG_class) |>
  dplyr::arrange(dplyr::desc(abs(log2FoldChange)), .by_group = TRUE) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup() |>
  dplyr::pull(Geneid) |>
  unique()

top_gene_ids <- intersect(top_gene_ids, rownames(vst_matrix))
top_z <- t(scale(t(vst_matrix[top_gene_ids, , drop = FALSE])))
top_z <- top_z[apply(top_z, 1, function(x) all(is.finite(x))), , drop = FALSE]

gene_label_table <- all_de |>
  dplyr::select(Geneid, Ensembl, gene_symbol) |>
  dplyr::distinct(Geneid, .keep_all = TRUE)
gene_labels <- gene_label_table$gene_symbol[match(rownames(top_z), gene_label_table$Geneid)]
fallback_labels <- gene_label_table$Ensembl[match(rownames(top_z), gene_label_table$Geneid)]
missing_labels <- is.na(gene_labels) | gene_labels == ""
gene_labels[missing_labels] <- fallback_labels[missing_labels]
missing_labels <- is.na(gene_labels) | gene_labels == ""
gene_labels[missing_labels] <- rownames(top_z)[missing_labels]
rownames(top_z) <- make.unique(gene_labels)

top_names <- coldata$analysis_sample_name[match(colnames(top_z), coldata$original_sample_name)]
colnames(top_z) <- top_names
top_order <- c("HH-1", "HH-2", "HH-3", "M-1", "M-2", "M-3", "78-1", "78-2", "78-3", "S-1", "S-2", "S-3")
top_z <- top_z[, top_order, drop = FALSE]
top_groups <- coldata$analysis_cell_line[match(top_order, coldata$analysis_sample_name)]

top_annotation <- ComplexHeatmap::HeatmapAnnotation(`Cell line` = top_groups, col = list(`Cell line` = cell_lines),
  show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

top_heatmap <- ComplexHeatmap::Heatmap(top_z, name = "Z-score", col = circlize::colorRamp2(c(-1.5, 0,
  1.5), c("#3F648C", "#F7F7F7", "#B55350")), top_annotation = top_annotation, cluster_rows = TRUE,
  cluster_columns = FALSE, clustering_distance_rows = "pearson", clustering_method_rows = "complete",
  show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 5.5), show_column_names = TRUE, column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 6.5), column_title = "Top differentially expressed genes",
  column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"), border = FALSE, heatmap_legend_param = list(at = c(-1.5,
    0, 1.5), labels = c("-1.5", "0", "1.5"), legend_height = grid::unit(22, "mm"), legend_width = grid::unit(5,
    "mm"), border = black, title_position = "topleft"))

grDevices::pdf(file.path(figures_de_dir, "top_DEGs_global_heatmap.pdf"), width = width_heatmap, height = height_large,
  useDingbats = FALSE)
ComplexHeatmap::draw(top_heatmap, heatmap_legend_side = "right")
grDevices::dev.off()

all_deg_heatmap_file <- file.path(figures_de_dir, "all_DEGs_large_heatmap.pdf")

if (file.exists(all_deg_heatmap_file) && file.info(all_deg_heatmap_file)$size > 1000) {
  message("Reusing valid all-DEG heatmap: ", all_deg_heatmap_file)
} else {
  all_deg_genes <- all_de |>
    dplyr::filter(DEG_class %in% c("Up", "Down")) |>
    dplyr::pull(Geneid) |>
    unique()

  all_deg_genes <- intersect(all_deg_genes, rownames(vst_matrix))
  large_z <- t(scale(t(vst_matrix[all_deg_genes, , drop = FALSE])))
  large_z <- large_z[apply(large_z, 1, function(x) all(is.finite(x))), , drop = FALSE]

  coldata <- SummarizedExperiment::colData(dds) |>
    as.data.frame() |>
    tibble::rownames_to_column("original_sample_name")

  analysis_names <- coldata$analysis_sample_name[match(colnames(large_z), coldata$original_sample_name)]
  colnames(large_z) <- analysis_names

  sample_order <- c("HH-1", "HH-2", "HH-3", "M-1", "M-2", "M-3", "78-1", "78-2", "78-3", "S-1", "S-2",
    "S-3")
  large_z <- large_z[, sample_order, drop = FALSE]

  sample_groups <- coldata$analysis_cell_line[match(sample_order, coldata$analysis_sample_name)]

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(`Cell line` = sample_groups, col = list(`Cell line` = cell_lines),
    show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

  large_heatmap <- ComplexHeatmap::Heatmap(large_z, name = "Row Z-score", col = circlize::colorRamp2(c(-1.5,
    0, 1.5), c("#3F648C", "#F7F7F7", "#B55350")), top_annotation = top_annotation, cluster_rows = TRUE,
    cluster_columns = FALSE, clustering_distance_rows = "pearson", clustering_method_rows = "complete",
    show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45, row_dend_width = grid::unit(6,
      "mm"), border = FALSE, use_raster = TRUE, raster_quality = 3, column_title = "All differentially expressed genes",
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"), heatmap_legend_param = list(at = c(-1.5,
      0, 1.5), labels = c("-1.5", "0", "1.5"), legend_height = grid::unit(22, "mm"), legend_width = grid::unit(5,
      "mm"), border = black, title_position = "topleft"))

  grDevices::pdf(all_deg_heatmap_file, width = 5, height = 6.5, useDingbats = FALSE)
  ComplexHeatmap::draw(large_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()}

# s3b-d
gse_vst_file <- ensure_file(file.path(external_data_dir, "GSE197067", "GSE197067_PanT_0h_vst_matrix.rds"),
  "scripts/07_external_healthy_reference.R")
bp_vst_file <- ensure_file(file.path(external_data_dir, "BLUEPRINT", "BLUEPRINT_CD4_venous_blood_vst_matrix.rds"),
  "scripts/07_external_healthy_reference.R")

gse_vst <- readRDS(gse_vst_file)
bp_vst <- readRDS(bp_vst_file)

plot_reference_pca <- function(matrix, group, title, filename, top_variable = NULL) {
  if (!is.null(top_variable)) {
    gene_var <- apply(matrix, 1, var)
    keep <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(top_variable, length(gene_var)))]
    matrix <- matrix[keep, , drop = FALSE]}

  pca <- stats::prcomp(t(matrix))
  pct <- 100 * pca$sdev^2/sum(pca$sdev^2)
  data <- as.data.frame(pca$x) |>
    tibble::rownames_to_column("sample")

  p <- ggplot2::ggplot(data, ggplot2::aes(PC1, PC2, label = sample)) + ggplot2::geom_hline(yintercept = 0,
    color = light_grey, linewidth = 0.3) + ggplot2::geom_vline(xintercept = 0, color = light_grey,
    linewidth = 0.3) + ggplot2::geom_point(shape = 21, size = 3.4, stroke = 0.5, color = black, fill = external_group_colors[group]) +
    ggrepel::geom_text_repel(size = 2.3, color = black) + ggplot2::labs(title = title, x = paste0("PC1 (",
    round(pct[1], 1), "%)"), y = paste0("PC2 (", round(pct[2], 1), "%)")) + theme_clean() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
    face = "bold"))

  save_single(filename, p, height = height_medium)}

plot_reference_pca(bp_vst, "Healthy CD4", "BLUEPRINT healthy CD4 T cells", file.path(figures_external_dir,
  "BLUEPRINT_CD4_PCA.pdf"))

plot_reference_pca(bp_vst, "Healthy CD4", "BLUEPRINT healthy CD4 T cells - top 500 variable genes", file.path(figures_external_dir,
  "BLUEPRINT_CD4_PCA_top500_variable_genes.pdf"), top_variable = 500)

plot_reference_pca(gse_vst, "Healthy Pan T", "GSE197067 healthy Pan T cells at 0 h", file.path(figures_external_dir,
  "GSE197067_PanT_0h_PCA.pdf"))

# s3a
gse_counts_file <- ensure_file(file.path(external_data_dir, "GSE197067", "GSE197067_PanT_0h_counts.rds"),
  "scripts/07_external_healthy_reference.R")
bp_counts_file <- ensure_file(file.path(external_data_dir, "BLUEPRINT", "BLUEPRINT_CD4_venous_blood_counts.rds"),
  "scripts/07_external_healthy_reference.R")

cell_counts <- DESeq2::counts(dds, normalized = FALSE)
rownames(cell_counts) <- stringr::str_remove(rownames(cell_counts), "\\.\\d+$")
colnames(cell_counts) <- coldata$analysis_sample_name[match(colnames(cell_counts), coldata$original_sample_name)]

to_ensembl_matrix <- function(object) {
  if (is.matrix(object)) {
    return(object)}
  if ("Ensembl" %in% colnames(object)) {
    return(object |>
      tibble::column_to_rownames("Ensembl") |>
      as.matrix())}
  as.matrix(object)}

gse_counts <- to_ensembl_matrix(readRDS(gse_counts_file))
bp_counts <- to_ensembl_matrix(readRDS(bp_counts_file))

common <- Reduce(intersect, list(rownames(cell_counts), rownames(gse_counts), rownames(bp_counts)))

cell_dge <- edgeR::calcNormFactors(edgeR::DGEList(cell_counts), method = "TMM")
gse_dge <- edgeR::calcNormFactors(edgeR::DGEList(gse_counts), method = "TMM")
bp_dge <- edgeR::calcNormFactors(edgeR::DGEList(bp_counts), method = "TMM")

cell_cpm <- edgeR::cpm(cell_dge)[common, , drop = FALSE]
gse_cpm <- edgeR::cpm(gse_dge)[common, , drop = FALSE]
bp_cpm <- edgeR::cpm(bp_dge)[common, , drop = FALSE]

keep <- rowSums(cell_cpm >= 1) >= 2 & rowSums(gse_cpm >= 1) >= 2 & rowSums(bp_cpm >= 1) >= 2
comparison_genes <- common[keep]

cell_log <- edgeR::cpm(cell_dge, log = TRUE, prior.count = 1)[comparison_genes, , drop = FALSE]
gse_log <- edgeR::cpm(gse_dge, log = TRUE, prior.count = 1)[comparison_genes, , drop = FALSE]
bp_log <- edgeR::cpm(bp_dge, log = TRUE, prior.count = 1)[comparison_genes, , drop = FALSE]

colnames(gse_log) <- paste0("GSE_", colnames(gse_log))
colnames(bp_log) <- paste0("BP_", colnames(bp_log))
external_matrix <- cbind(cell_log, gse_log, bp_log)
external_cor <- cor(external_matrix, method = "spearman")

cell_groups <- coldata$analysis_cell_line[match(colnames(cell_log), coldata$analysis_sample_name)]
external_annotation <- data.frame(Dataset = c(rep("Cell lines", 12), rep("GSE197067", ncol(gse_log)),
  rep("BLUEPRINT", ncol(bp_log))), Group = c(as.character(cell_groups), rep("Healthy Pan T", ncol(gse_log)),
  rep("Healthy CD4", ncol(bp_log))), row.names = colnames(external_matrix), check.names = FALSE)

cor_colors <- (grDevices::colorRampPalette(c("#F7F7F7", "#C6DED9", "#2A7F7F")))(101)
external_dataset_colors <- c(`Cell lines` = "#777777", GSE197067 = unname(external_group_colors["Healthy Pan T"]),
  BLUEPRINT = unname(external_group_colors["Healthy CD4"]))
external_annotation_colors <- list(Dataset = external_dataset_colors, Group = external_group_colors)

external_top_annotation <- ComplexHeatmap::HeatmapAnnotation(df = external_annotation, col = external_annotation_colors,
  show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))
external_left_annotation <- ComplexHeatmap::rowAnnotation(df = external_annotation, col = external_annotation_colors,
  show_annotation_name = FALSE, simple_anno_size = grid::unit(3, "mm"))

external_cor_min <- floor(min(external_cor) * 100)/100
external_cor_mid <- (external_cor_min + 1)/2
external_cor_marks <- c(external_cor_min, external_cor_mid, 1)

external_cor_heatmap <- ComplexHeatmap::Heatmap(external_cor, name = "Spearman rho", col = circlize::colorRamp2(external_cor_marks,
  c("#F7F7F7", "#C6DED9", "#2A7F7F")), top_annotation = external_top_annotation, left_annotation = external_left_annotation,
  cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 5.5),
  column_names_gp = grid::gpar(fontsize = 5.5), column_names_rot = 45, column_title = "Cross-dataset transcriptomic similarity",
  column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"), border = FALSE, heatmap_legend_param = list(at = external_cor_marks,
    labels = format(external_cor_marks, digits = 2), legend_height = grid::unit(22, "mm"), legend_width = grid::unit(5,
      "mm"), border = black, title_position = "topleft"))

grDevices::pdf(file.path(figures_external_dir, "cross_dataset_correlation_heatmap.pdf"), width = 6.2,
  height = 6.2, useDingbats = FALSE)
ComplexHeatmap::draw(external_cor_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
grDevices::dev.off()

# recreate venn panels
external_rank_matrix <- apply(external_matrix, 2, rank, ties.method = "average")
gse_samples <- colnames(gse_log)
blueprint_samples <- colnames(bp_log)
results_external_dir <- external_dir
source("scripts/08_venn_diagrams.R")

# s4a-b
module_expression_file <- ensure_file(file.path(functional_dir, "gene_module_expression_by_cell_line.csv"),
  "scripts/06_functional_enrichment.R")
module_external_file <- ensure_file(file.path(external_dir, "gene_module_scores_CTCL_vs_healthy_summary.csv"),
  "scripts/07_external_healthy_reference.R")

module_labels <- c(`Module 1` = "Adhesion / morphogenesis", `Module 2` = "T-cell activation / cytokine signaling",
  `Module 3` = "Developmental program", `Module 4` = "Ion transport", `Module 5` = "Adaptive immunity / T-cell migration",
  `Module 6` = "Immune proliferation")

plot_module_heatmap <- function(data, group_column, value_column, title, filename, groups) {
  data <- data |>
    dplyr::mutate(module = factor(module, levels = names(module_labels)), module_label = factor(unname(module_labels[as.character(module)]),
      levels = unname(module_labels)), plot_group = factor(.data[[group_column]], levels = groups))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = plot_group, y = module_label, fill = .data[[value_column]])) +
    ggplot2::geom_tile(color = white, linewidth = 0.35) + ggplot2::scale_fill_gradient2(low = gradient_diverging["low"],
    mid = gradient_diverging["mid"], high = gradient_diverging["high"], midpoint = 0, limits = c(-2,
      2), oob = scales::squish, name = "Z-score", breaks = c(-2, 0, 2), labels = c("-2", "0", "2"),
    guide = ggplot2::guide_colourbar(barheight = grid::unit(22, "mm"), barwidth = grid::unit(5, "mm"),
      frame.colour = black, frame.linewidth = 0.6, ticks = FALSE, title.position = "top")) + ggplot2::labs(title = title,
    x = NULL, y = NULL) + theme_clean() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40,
    hjust = 1), axis.text.y = ggplot2::element_text(size = 6.2), axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  ggplot2::ggsave(filename, p, width = 5.6, height = 3.5, units = "in", dpi = dpi)}

plot_module_heatmap(readr::read_csv(module_expression_file, show_col_types = FALSE), "analysis_cell_line",
  "mean_z_score", "Cell-line gene-module expression", file.path(figures_functional_dir, "gene_module_expression_patterns.pdf"),
  c("HH", "MyLa", "HuT 78", "SeAx"))

plot_module_heatmap(readr::read_csv(module_external_file, show_col_types = FALSE), "Group", "mean_score_z",
  "Gene-module activity: CTCL vs healthy T-cell references", file.path(figures_external_dir, "gene_module_activity_healthy_references.pdf"),
  c("HH", "MyLa", "HuT 78", "SeAx", "Healthy Pan T", "Healthy CD4"))
