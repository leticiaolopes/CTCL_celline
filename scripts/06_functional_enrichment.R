# gene modules & functional enrichment

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(rrvgo)
library(BiocParallel)
library(GOSemSim)
library(tidytext)

BiocParallel::register(BiocParallel::SerialParam())

# Analysis scope. The default preserves the original all-gene workflow.
gene_scope <- getOption("ctcl.gene_scope", "all")
if (!gene_scope %in% c("all", "protein_coding")) {
  stop("Unsupported ctcl.gene_scope: ", gene_scope)
}
analysis_suffix <- if (gene_scope == "protein_coding") "_protein_coding" else ""

de_dir <- file.path(paths$results, paste0("differential_expression", analysis_suffix))

deseq2_dir <- file.path(paths$results, paste0("deseq2", analysis_suffix))

enrichment_dir <- file.path(paths$results, paste0("functional_enrichment", analysis_suffix))

figures_enrichment_dir <- file.path(
  paths$figures,
  "drafts",
  paste0("functional_enrichment", analysis_suffix))

dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(figures_enrichment_dir, recursive = TRUE, showWarnings = FALSE)

module_k <- 6
go_padj_cutoff <- 0.05
gsea_padj_cutoff <- 0.05
rrvgo_similarity_threshold <- 0.7
go_terms_per_module <- 6

go_semdata_bp <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP")

# load de results and vst matrix
dds <- readRDS(file.path(deseq2_dir, "dds.rds"))

all_shrunk_results_annotated <- readRDS(file.path(de_dir, "all_shrunk_results_annotated.rds"))

vst_matrix <- readRDS(file.path(deseq2_dir, "vst_matrix.rds"))

if (gene_scope == "protein_coding") {
  protein_coding_ids <- readRDS(
    file.path(deseq2_dir, "protein_coding_gene_ids.rds"))
  if (!all(rownames(dds) %in% protein_coding_ids) ||
      !all(rownames(vst_matrix) %in% protein_coding_ids)) {
    stop("Protein-coding enrichment inputs contain genes outside the selected biotype.")
  }
}

coldata_enrichment <- colData(dds) |>
  as.data.frame() |>
  rownames_to_column("original_sample_name")

sample_name_map <- coldata_enrichment |>
  dplyr::select(original_sample_name, analysis_sample_name, analysis_cell_line)

sample_order <- sample_name_map |>
  dplyr::mutate(analysis_cell_line = factor(analysis_cell_line, levels = c("HH", "MyLa", "HuT 78",
    "SeAx"))) |>
  dplyr::arrange(analysis_cell_line, analysis_sample_name) |>
  dplyr::pull(analysis_sample_name)

# reconstruct all DEG heatmap clustering
all_deg_genes <- all_shrunk_results_annotated |>
  dplyr::filter(DEG_class %in% c("Up", "Down")) |>
  dplyr::pull(Geneid) |>
  unique()

large_heatmap_matrix <- vst_matrix[all_deg_genes, , drop = FALSE]

large_heatmap_matrix_z <- t(scale(t(large_heatmap_matrix)))

large_heatmap_matrix_z <- large_heatmap_matrix_z[apply(large_heatmap_matrix_z, 1, function(x) {
  all(is.finite(x))
}), , drop = FALSE]

new_sample_names <- sample_name_map$analysis_sample_name[match(colnames(large_heatmap_matrix_z), sample_name_map$original_sample_name)]

if (any(is.na(new_sample_names))) {
  stop("Could not map all original sample names to analysis sample names.")}

colnames(large_heatmap_matrix_z) <- unname(new_sample_names)

large_heatmap_matrix_z <- large_heatmap_matrix_z[, sample_order, drop = FALSE]

gene_distance <- as.dist(1 - cor(t(large_heatmap_matrix_z), method = "pearson"))

gene_hclust <- hclust(gene_distance, method = "complete")

# define gene modules
gene_modules <- cutree(gene_hclust, k = module_k)

module_assignments <- tibble(Geneid = names(gene_modules), module = paste0("Module ", gene_modules))

gene_annotation <- all_shrunk_results_annotated |>
  dplyr::select(Geneid, Ensembl, gene_symbol, gene_type) |>
  dplyr::distinct(Geneid, .keep_all = TRUE)

module_assignments <- module_assignments |>
  dplyr::left_join(gene_annotation, by = "Geneid")

module_sizes <- module_assignments |>
  dplyr::count(module, name = "n_genes") |>
  dplyr::arrange(module)

print(module_sizes)

# summarize module expression by cell line
expression_long <- large_heatmap_matrix_z |>
  as.data.frame() |>
  rownames_to_column("Geneid") |>
  pivot_longer(cols = -Geneid, names_to = "analysis_sample_name", values_to = "z_score") |>
  left_join(module_assignments |>
    dplyr::select(Geneid, module), by = "Geneid") |>
  left_join(sample_name_map |>
    dplyr::select(analysis_sample_name, analysis_cell_line), by = "analysis_sample_name")

module_expression <- expression_long |>
  dplyr::group_by(module, analysis_cell_line) |>
  dplyr::summarise(mean_z_score = mean(z_score), .groups = "drop")

module_expression_wide <- module_expression |>
  tidyr::pivot_wider(names_from = analysis_cell_line, values_from = mean_z_score)

print(module_expression_wide)

module_expression_matrix <- module_expression |>
  dplyr::mutate(analysis_cell_line = factor(analysis_cell_line, levels = c("HH", "MyLa", "HuT 78",
    "SeAx"))) |>
  dplyr::arrange(module, analysis_cell_line) |>
  tidyr::pivot_wider(names_from = analysis_cell_line, values_from = mean_z_score) |>
  tibble::column_to_rownames("module") |>
  as.matrix()

p_module_expression <- pheatmap::pheatmap(module_expression_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
  color = heatmap_deg, border_color = NA, fontsize = 7, fontsize_row = 7, fontsize_col = 7, angle_col = 45,
  main = "Cell-line gene-module expression", silent = TRUE)

matrix_index <- which(p_module_expression$gtable$layout$name == "matrix")

matrix_layout <- p_module_expression$gtable$layout[matrix_index, ]

p_module_expression$gtable <- gtable::gtable_add_grob(p_module_expression$gtable, grobs = grid::rectGrob(gp = grid::gpar(fill = NA,
  col = black, lwd = 0.5)), t = matrix_layout$t, l = matrix_layout$l, b = matrix_layout$b, r = matrix_layout$r,
  z = Inf, name = "matrix_border")

pdf(file.path(figures_enrichment_dir, "gene_module_expression_patterns.pdf"), width = 4.4, height = 3.6,
  useDingbats = FALSE)

grid::grid.newpage()

grid::grid.draw(p_module_expression$gtable)

dev.off()

# GO biological process enrichment by module
tested_genes <- all_shrunk_results_annotated |>
  dplyr::select(Geneid, Ensembl) |>
  dplyr::distinct(Geneid, .keep_all = TRUE)

tested_genes$ENTREZID <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = tested_genes$Ensembl, column = "ENTREZID",
  keytype = "ENSEMBL", multiVals = "first")

go_universe <- tested_genes |>
  dplyr::filter(!is.na(ENTREZID)) |>
  dplyr::pull(ENTREZID) |>
  unique()

module_gene_mapping <- module_assignments |>
  dplyr::mutate(ENTREZID = AnnotationDbi::mapIds(org.Hs.eg.db, keys = Ensembl, column = "ENTREZID",
    keytype = "ENSEMBL", multiVals = "first")) |>
  dplyr::filter(!is.na(ENTREZID))

module_go_results <- lapply(sort(unique(module_gene_mapping$module)), function(module_name) {
  module_genes <- module_gene_mapping |>
    dplyr::filter(module == module_name) |>
    dplyr::pull(ENTREZID) |>
    unique()

  go_result <- clusterProfiler::enrichGO(gene = module_genes, universe = go_universe, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

  go_table <- as.data.frame(go_result)

  if (nrow(go_table) == 0) {
    return(tibble::tibble())}

  go_table |>
    dplyr::mutate(module = module_name)})

module_go_table <- dplyr::bind_rows(module_go_results) |>
  dplyr::filter(p.adjust < go_padj_cutoff)

module_mapping_summary <- module_assignments |>
  dplyr::left_join(module_gene_mapping |>
    dplyr::select(Geneid, ENTREZID), by = "Geneid") |>
  dplyr::group_by(module) |>
  dplyr::summarise(total_genes = dplyr::n(), mapped_to_entrez = sum(!is.na(ENTREZID)), mapping_percent = 100 *
    mapped_to_entrez/total_genes, .groups = "drop")

print(module_mapping_summary)

go_terms <- module_go_table |>
  dplyr::pull(ID) |>
  unique()

go_scores <- module_go_table |>
  dplyr::group_by(ID) |>
  dplyr::summarise(score = max(-log10(p.adjust), na.rm = TRUE), .groups = "drop")

go_score_vector <- go_scores$score

names(go_score_vector) <- go_scores$ID

if (any(is.infinite(go_score_vector))) {
  max_finite_go_score <- max(go_score_vector[is.finite(go_score_vector)], na.rm = TRUE)

  go_score_vector[is.infinite(go_score_vector)] <- max_finite_go_score + 1}

go_similarity_matrix <- rrvgo::calculateSimMatrix(go_terms, orgdb = "org.Hs.eg.db", semdata = go_semdata_bp,
  ont = "BP", method = "Rel")

saveRDS(go_similarity_matrix, file.path(enrichment_dir, "GO_BP_similarity_matrix.rds"))

go_rrvgo <- rrvgo::reduceSimMatrix(go_similarity_matrix, scores = go_score_vector, threshold = rrvgo_similarity_threshold,
  orgdb = "org.Hs.eg.db")

go_rrvgo_mapping <- go_rrvgo |>
  dplyr::select(go, term, parent, parentTerm, score, size)

module_go_reduced <- module_go_table |>
  dplyr::left_join(go_rrvgo_mapping, by = c(ID = "go")) |>
  dplyr::group_by(module, parent, parentTerm) |>
  dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::rename(representative_ID = parent, representative_term = parentTerm)

selected_representative_terms <- module_go_reduced |>
  dplyr::group_by(module) |>
  dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) |>
  dplyr::slice_head(n = go_terms_per_module) |>
  dplyr::ungroup() |>
  dplyr::pull(representative_ID) |>
  unique()

representative_term_table <- module_go_reduced |>
  dplyr::filter(representative_ID %in% selected_representative_terms) |>
  dplyr::select(representative_ID, representative_term) |>
  dplyr::distinct() |>
  dplyr::left_join(go_rrvgo_mapping |>
    dplyr::select(go, representative_score = score), by = c(representative_ID = "go")) |>
  dplyr::arrange(dplyr::desc(representative_score))

module_levels <- paste0("Module ", seq_len(module_k))

module_go_plot_grid <- tidyr::expand_grid(representative_ID = representative_term_table$representative_ID,
  module = module_levels) |>
  dplyr::left_join(representative_term_table, by = "representative_ID") |>
  dplyr::left_join(module_go_reduced |>
    dplyr::select(representative_ID, module, p.adjust, Count, GeneRatio), by = c("representative_ID",
    "module")) |>
  dplyr::mutate(status = dplyr::if_else(!is.na(p.adjust) & p.adjust < go_padj_cutoff, "Significant",
    "NS"), minus_log10_padj = dplyr::if_else(status == "Significant", -log10(p.adjust), NA_real_),
    module = factor(module, levels = module_levels), representative_term = stringr::str_wrap(representative_term,
      width = 45))

term_levels <- representative_term_table |>
  dplyr::mutate(representative_term = stringr::str_wrap(representative_term, width = 45)) |>
  dplyr::pull(representative_term)

module_go_plot_grid <- module_go_plot_grid |>
  dplyr::mutate(representative_term = factor(representative_term, levels = rev(term_levels)))

# one-vs-rest cell-line contrasts
coef_names <- resultsNames(dds)

print(coef_names)

expected_coef_names <- c("Intercept", "cell_line_model_MyLa_vs_HH", "cell_line_model_HuT78_vs_HH", "cell_line_model_SeAx_vs_HH")

if (!all(expected_coef_names %in% coef_names)) {
  stop("Unexpected DESeq2 coefficient names. Check resultsNames(dds).")}

contrast_hh_vs_rest <- c(0, -1/3, -1/3, -1/3)

contrast_myla_vs_rest <- c(0, 1, -1/3, -1/3)

contrast_hut78_vs_rest <- c(0, -1/3, 1, -1/3)

contrast_seax_vs_rest <- c(0, -1/3, -1/3, 1)

names(contrast_hh_vs_rest) <- coef_names

names(contrast_myla_vs_rest) <- coef_names

names(contrast_hut78_vs_rest) <- coef_names

names(contrast_seax_vs_rest) <- coef_names

res_hh_vs_rest <- results(dds, contrast = contrast_hh_vs_rest) |>
  as.data.frame() |>
  rownames_to_column("Geneid") |>
  mutate(cell_line = "HH")

res_myla_vs_rest <- results(dds, contrast = contrast_myla_vs_rest) |>
  as.data.frame() |>
  rownames_to_column("Geneid") |>
  mutate(cell_line = "MyLa")

res_hut78_vs_rest <- results(dds, contrast = contrast_hut78_vs_rest) |>
  as.data.frame() |>
  rownames_to_column("Geneid") |>
  mutate(cell_line = "HuT 78")

res_seax_vs_rest <- results(dds, contrast = contrast_seax_vs_rest) |>
  as.data.frame() |>
  rownames_to_column("Geneid") |>
  mutate(cell_line = "SeAx")

one_vs_rest_results <- dplyr::bind_rows(res_hh_vs_rest, res_myla_vs_rest, res_hut78_vs_rest, res_seax_vs_rest) |>
  dplyr::left_join(gene_annotation, by = "Geneid")

# go-bp gsea by cell line
run_go_gsea <- function(result_table, cell_line_name) {
  ranked_genes <- result_table |>
    dplyr::filter(cell_line == cell_line_name, !is.na(stat), !is.na(Ensembl)) |>
    dplyr::mutate(ENTREZID = AnnotationDbi::mapIds(org.Hs.eg.db, keys = Ensembl, column = "ENTREZID",
      keytype = "ENSEMBL", multiVals = "first")) |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::group_by(ENTREZID) |>
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(stat))

  gene_list <- ranked_genes$stat

  names(gene_list) <- ranked_genes$ENTREZID

  gene_list <- sort(gene_list, decreasing = TRUE)

  gsea_result <- clusterProfiler::gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
    ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", by = "fgsea",
    seed = TRUE, verbose = FALSE)

  gsea_table <- as.data.frame(gsea_result)

  if (nrow(gsea_table) == 0) {
    return(tibble::tibble())}

  gsea_table |>
    dplyr::mutate(cell_line = cell_line_name, direction = dplyr::case_when(NES > 0 ~ "Up", NES <
      0 ~ "Down", TRUE ~ "Neutral"))}

gsea_hh <- run_go_gsea(one_vs_rest_results, "HH")

gsea_myla <- run_go_gsea(one_vs_rest_results, "MyLa")

gsea_hut78 <- run_go_gsea(one_vs_rest_results, "HuT 78")

gsea_seax <- run_go_gsea(one_vs_rest_results, "SeAx")

gsea_go_table <- dplyr::bind_rows(gsea_hh, gsea_myla, gsea_hut78, gsea_seax) |>
  dplyr::filter(p.adjust < gsea_padj_cutoff, direction %in% c("Up", "Down"))

gsea_groups <- split(gsea_go_table, interaction(gsea_go_table$cell_line, gsea_go_table$direction, drop = TRUE))

gsea_rrvgo_results <- lapply(gsea_groups, function(group_data) {
  if (nrow(group_data) == 0) {
    return(tibble::tibble())}

  group_data <- group_data |>
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) |>
    dplyr::distinct(ID, .keep_all = TRUE)

  go_terms_group <- group_data$ID

  go_scores_group <- -log10(group_data$p.adjust)

  names(go_scores_group) <- group_data$ID

  if (any(is.infinite(go_scores_group))) {
    max_finite_score <- max(go_scores_group[is.finite(go_scores_group)], na.rm = TRUE)

    go_scores_group[is.infinite(go_scores_group)] <- max_finite_score + 1}

  if (length(go_terms_group) == 1) {
    return(group_data |>
      dplyr::mutate(representative_ID = ID, representative_term = Description))}

  similarity_matrix <- rrvgo::calculateSimMatrix(go_terms_group, orgdb = "org.Hs.eg.db", semdata = go_semdata_bp,
    ont = "BP", method = "Rel")

  reduced_terms <- rrvgo::reduceSimMatrix(similarity_matrix, scores = go_scores_group, threshold = rrvgo_similarity_threshold,
    orgdb = "org.Hs.eg.db")

  reduced_mapping <- reduced_terms |>
    dplyr::select(go, parent, parentTerm)

  group_data |>
    dplyr::left_join(reduced_mapping, by = c(ID = "go")) |>
    dplyr::group_by(parent, parentTerm) |>
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES)), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::rename(representative_ID = parent, representative_term = parentTerm)})

gsea_go_reduced <- dplyr::bind_rows(gsea_rrvgo_results)

gsea_plot_data <- gsea_go_reduced |>
  dplyr::filter(!is.na(representative_term), !is.na(NES), !is.na(p.adjust)) |>
  dplyr::mutate(cell_line = factor(cell_line, levels = c("HH", "MyLa", "HuT 78", "SeAx")), direction = dplyr::case_when(NES >
    0 ~ "Up", NES < 0 ~ "Down", TRUE ~ NA_character_)) |>
  dplyr::filter(!is.na(direction)) |>
  dplyr::group_by(cell_line, direction) |>
  dplyr::arrange(dplyr::desc(abs(NES)), p.adjust, .by_group = TRUE) |>
  dplyr::slice_head(n = 5) |>
  dplyr::ungroup() |>
  dplyr::mutate(minus_log10_padj = -log10(p.adjust), representative_term = stringr::str_wrap(representative_term,
    width = 42))

gsea_plot_data <- gsea_plot_data |>
  dplyr::group_by(cell_line) |>
  dplyr::arrange(dplyr::desc(direction == "Up"), dplyr::desc(abs(NES)), .by_group = TRUE) |>
  dplyr::mutate(plot_order = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(term_plot = tidytext::reorder_within(representative_term, -plot_order, cell_line))

gsea_plot_summary <- gsea_plot_data |>
  dplyr::count(cell_line, direction)

print(gsea_plot_summary)

# export tables
write_csv(module_assignments, file.path(enrichment_dir, "gene_module_assignments.csv"))

write_csv(module_sizes, file.path(enrichment_dir, "gene_module_sizes.csv"))

write_csv(module_expression, file.path(enrichment_dir, "gene_module_expression_by_cell_line.csv"))

write_csv(module_mapping_summary, file.path(enrichment_dir, "gene_module_GO_mapping_summary.csv"))

write_csv(module_go_table, file.path(enrichment_dir, "GO_BP_by_gene_module.csv"))

write_csv(go_rrvgo_mapping, file.path(enrichment_dir, "GO_BP_rrvgo_reduction.csv"))

write_csv(module_go_reduced, file.path(enrichment_dir, "GO_BP_by_gene_module_reduced.csv"))

write_csv(module_go_plot_grid, file.path(enrichment_dir, "GO_BP_gene_module_plot_grid.csv"))

write_csv(one_vs_rest_results, file.path(enrichment_dir, "one_vs_rest_DESeq2_results.csv"))

write_csv(gsea_go_table, file.path(enrichment_dir, "GO_BP_GSEA_one_vs_rest.csv"))

write_csv(gsea_go_reduced, file.path(enrichment_dir, "GO_BP_GSEA_one_vs_rest_rrvgo.csv"))

write_csv(gsea_plot_data, file.path(enrichment_dir, "GO_BP_GSEA_top5_up_top5_down.csv"))

# generate enrichment figures
raw_selected_terms <- module_go_table |>
  dplyr::group_by(module) |>
  dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) |>
  dplyr::slice_head(n = go_terms_per_module) |>
  dplyr::ungroup() |>
  dplyr::pull(ID) |>
  unique()

raw_term_table <- module_go_table |>
  dplyr::filter(ID %in% raw_selected_terms) |>
  dplyr::select(ID, Description) |>
  dplyr::distinct()

raw_module_go_grid <- tidyr::expand_grid(ID = raw_selected_terms, module = module_levels) |>
  dplyr::left_join(raw_term_table, by = "ID") |>
  dplyr::left_join(module_go_table |>
    dplyr::select(ID, module, p.adjust, Count), by = c("ID", "module")) |>
  dplyr::mutate(status = dplyr::if_else(!is.na(p.adjust) & p.adjust < go_padj_cutoff, "Significant",
    "NS"), minus_log10_padj = dplyr::if_else(status == "Significant", -log10(p.adjust), NA_real_),
    Description = stringr::str_wrap(Description, width = 45), module = factor(module, levels = module_levels))

p_module_go_raw <- ggplot2::ggplot(raw_module_go_grid, ggplot2::aes(x = module, y = Description)) + ggplot2::geom_point(data = raw_module_go_grid |>
  dplyr::filter(status == "Significant"), ggplot2::aes(size = Count, color = minus_log10_padj), shape = 16,
  alpha = 0.9) + ggplot2::geom_point(data = raw_module_go_grid |>
  dplyr::filter(status == "NS"), ggplot2::aes(shape = status), size = 1.7, color = mid_grey) + ggplot2::scale_color_gradient(low = gradient_diverging["low"],
  high = gradient_diverging["high"], name = expression(-log[10](padj))) + ggplot2::scale_size_continuous(name = "Gene count",
  range = c(1.8, 4.5)) + ggplot2::scale_shape_manual(values = c(NS = 18), name = NULL) + ggplot2::labs(x = NULL,
  y = NULL) + theme_clean(base_size = 8) + ggplot2::theme(panel.border = ggplot2::element_rect(colour = black,
  fill = NA, linewidth = 0.4), panel.grid.major = ggplot2::element_line(colour = light_grey, linewidth = 0.2),
  panel.grid.minor = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
  axis.text.y = ggplot2::element_text(size = 6.5), legend.position = "right", plot.margin = ggplot2::margin(2,
    2, 2, 2))

ggplot2::ggsave(file.path(figures_enrichment_dir, "GO_BP_gene_modules_dotplot.pdf"), p_module_go_raw,
  width = 7, height = 7, units = "in")

p_module_go <- ggplot2::ggplot(module_go_plot_grid, ggplot2::aes(x = module, y = representative_term)) +
  ggplot2::geom_point(data = module_go_plot_grid |>
    dplyr::filter(status == "Significant"), ggplot2::aes(size = Count, color = minus_log10_padj),
    shape = 16, alpha = 0.9) + ggplot2::geom_point(data = module_go_plot_grid |>
  dplyr::filter(status == "NS"), ggplot2::aes(shape = status), size = 1.7, color = mid_grey) + ggplot2::scale_color_gradient(low = gradient_diverging["low"],
  high = gradient_diverging["high"], name = expression(-log[10](padj))) + ggplot2::scale_size_continuous(name = "Gene count",
  range = c(1.8, 4.5)) + ggplot2::scale_shape_manual(values = c(NS = 18), name = NULL) + ggplot2::labs(x = NULL,
  y = NULL) + theme_clean(base_size = 8) + ggplot2::theme(panel.border = ggplot2::element_rect(colour = black,
  fill = NA, linewidth = 0.4), panel.grid.major = ggplot2::element_line(colour = light_grey, linewidth = 0.2),
  panel.grid.minor = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
  axis.text.y = ggplot2::element_text(size = 6.5), legend.position = "right", plot.margin = ggplot2::margin(2,
    2, 2, 2))

ggplot2::ggsave(file.path(figures_enrichment_dir, "GO_BP_gene_modules_rrvgo_dotplot.pdf"), p_module_go,
  width = 7, height = 7, units = "in")

p_gsea <- ggplot2::ggplot(gsea_plot_data, ggplot2::aes(x = NES, y = term_plot)) + ggplot2::geom_point(ggplot2::aes(size = minus_log10_padj,
  color = NES), shape = 16, alpha = 0.95) + ggplot2::facet_wrap(~cell_line, scales = "free_y", ncol = 2) +
  tidytext::scale_y_reordered() + ggplot2::scale_color_gradient2(low = gradient_diverging["low"], mid = gradient_diverging["mid"],
  high = gradient_diverging["high"], midpoint = 0, name = "NES") + ggplot2::scale_size_continuous(name = expression(-log[10](padj)),
  range = c(2.5, 6), guide = ggplot2::guide_legend(reverse = TRUE)) + ggplot2::labs(x = "Normalized enrichment score",
  y = NULL) + theme_clean(base_size = 8) + ggplot2::theme(panel.border = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_line(colour = light_grey, linewidth = 0.25), panel.grid.minor = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_text(size = 6), strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "plain"),
  legend.position = "right")

ggplot2::ggsave(file.path(figures_enrichment_dir, "GO_BP_GSEA_one_vs_rest.pdf"), p_gsea, width = 7, height = 4.5,
  units = "in")
