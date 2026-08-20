# top 50 DEG heatmap

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

figures_external_dir <- file.path(paths$figures, "drafts", "external_reference")

# load upstream results
one_vs_rest <- readr::read_csv(file.path(paths$results, "functional_enrichment", "one_vs_rest_DESeq2_results.csv"),
  show_col_types = FALSE)

dds <- readRDS(file.path(paths$results, "deseq2", "dds.rds"))

cell_line_metadata <- SummarizedExperiment::colData(dds) |>
  as.data.frame() |>
  tibble::rownames_to_column("original_sample_name")

# select and scale the strongest degs
top50_genes <- one_vs_rest |>
  dplyr::filter(!is.na(padj), padj < 0.01, abs(log2FoldChange) >= 2) |>
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) |>
  dplyr::mutate(Ensembl = stringr::str_remove(Geneid, "\\.\\d+$")) |>
  dplyr::distinct(Ensembl, .keep_all = TRUE) |>
  dplyr::slice_head(n = 50)

vst_matrix <- readRDS(file.path(paths$results, "deseq2", "vst_matrix.rds"))

rownames(vst_matrix) <- stringr::str_remove(rownames(vst_matrix), "\\.\\d+$")

top50_matrix <- vst_matrix[top50_genes$Ensembl, , drop = FALSE]

new_sample_names <- cell_line_metadata$analysis_sample_name[match(colnames(top50_matrix), cell_line_metadata$original_sample_name)]
colnames(top50_matrix) <- new_sample_names

sample_order <- c("HH-1", "HH-2", "HH-3", "M-1", "M-2", "M-3", "78-1", "78-2", "78-3", "S-1", "S-2",
  "S-3")

top50_matrix <- top50_matrix[, sample_order, drop = FALSE]

top50_matrix_z <- t(scale(t(top50_matrix)))

complete_genes <- apply(top50_matrix_z, 1, function(x) all(is.finite(x)))

top50_matrix_z <- top50_matrix_z[complete_genes, , drop = FALSE]

all_shrunk_results_annotated <- readRDS(file.path(paths$results, "differential_expression", "all_shrunk_results_annotated.rds"))

gene_label_map <- all_shrunk_results_annotated |>
  dplyr::select(Ensembl, gene_symbol) |>
  dplyr::distinct(Ensembl, .keep_all = TRUE)

gene_labels <- gene_label_map$gene_symbol[match(rownames(top50_matrix_z), gene_label_map$Ensembl)]

missing_gene_labels <- is.na(gene_labels) | gene_labels == ""
gene_labels[missing_gene_labels] <- rownames(top50_matrix_z)[missing_gene_labels]

rownames(top50_matrix_z) <- gene_labels

annotation_col <- data.frame(`Cell line` = rep(c("HH", "MyLa", "HuT 78", "SeAx"), each = 3), row.names = sample_order,
  check.names = FALSE)

annotation_colors <- list(`Cell line` = cell_lines)

rownames(top50_matrix_z)[rownames(top50_matrix_z) == "ENSG00000287597"] <- "HSALNG0026834"

# order genes && samples within cell line groups
cell_line_means <- cbind(HH = rowMeans(top50_matrix_z[, c("HH-1", "HH-2", "HH-3"), drop = FALSE]), MyLa = rowMeans(top50_matrix_z[,
  c("M-1", "M-2", "M-3"), drop = FALSE]), `HuT 78` = rowMeans(top50_matrix_z[, c("78-1", "78-2", "78-3"),
  drop = FALSE]), SeAx = rowMeans(top50_matrix_z[, c("S-1", "S-2", "S-3"), drop = FALSE]))

dominant_cell_line <- colnames(cell_line_means)[max.col(cell_line_means, ties.method = "first")]

dominant_cell_line <- factor(dominant_cell_line, levels = c("HH", "MyLa", "HuT 78", "SeAx"))
row_order <- order(dominant_cell_line, -apply(cell_line_means, 1, max))
top50_matrix_z <- top50_matrix_z[row_order, , drop = FALSE]
dominant_cell_line <- dominant_cell_line[row_order]
gaps_row <- cumsum(table(dominant_cell_line))
gaps_row <- gaps_row[-length(gaps_row)]

column_groups <- list(HH = c("HH-1", "HH-2", "HH-3"), MyLa = c("M-1", "M-2", "M-3"), `HuT 78` = c("78-1",
  "78-2", "78-3"), SeAx = c("S-1", "S-2", "S-3"))

clustered_column_order <- unlist(lapply(column_groups, function(samples) {
  sample_correlations <- cor(top50_matrix_z[, samples, drop = FALSE], method = "pearson")

  samples[hclust(as.dist(1 - sample_correlations), method = "complete")$order]
}), use.names = FALSE)

top50_matrix_z <- top50_matrix_z[, clustered_column_order, drop = FALSE]
annotation_col <- annotation_col[clustered_column_order, , drop = FALSE]

row_group_levels <- c("HH", "MyLa", "HuT 78", "SeAx")
clustered_row_order <- unlist(lapply(row_group_levels, function(group) {
  idx <- which(dominant_cell_line == group)
  if (length(idx) <= 1) {
    return(idx)}
  gene_correlations <- cor(t(top50_matrix_z[idx, , drop = FALSE]), method = "pearson")

  idx[hclust(as.dist(1 - gene_correlations), method = "complete")$order]
}), use.names = FALSE)
top50_matrix_z <- top50_matrix_z[clustered_row_order, , drop = FALSE]
row_group_sizes <- table(dominant_cell_line)
gaps_row <- cumsum(row_group_sizes)
gaps_row <- gaps_row[-length(gaps_row)]

p_top50_heatmap <- pheatmap::pheatmap(top50_matrix_z, color = heatmap_deg, annotation_col = annotation_col,
  annotation_colors = annotation_colors, cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps_row,
  gaps_col = c(3, 6, 9), show_rownames = TRUE, border_color = NA, fontsize = 7, fontsize_row = 6, fontsize_col = 7,
  angle_col = 45, annotation_names_col = FALSE, main = "Top 50 differentially expressed genes")

# export heatmap
pdf(file.path(figures_external_dir, "top50_differentially_expressed_genes_heatmap.pdf"), width = width_heatmap,
  height = height_large, useDingbats = FALSE)
grid::grid.newpage()
grid::grid.draw(p_top50_heatmap$gtable)
dev.off()
