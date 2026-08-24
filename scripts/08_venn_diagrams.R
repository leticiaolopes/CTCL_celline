# compare cell line DEG overlaps with healthy T cell refs

# setup
source("scripts/00_config.R")
source("scripts/00_aesthetics.R")

library(ggVennDiagram)
library(patchwork)

# Analysis scope. The default preserves the original all-gene workflow.
gene_scope <- getOption("ctcl.gene_scope", "all")
if (!gene_scope %in% c("all", "protein_coding")) {
  stop("Unsupported ctcl.gene_scope: ", gene_scope)
}
analysis_suffix <- if (gene_scope == "protein_coding") "_protein_coding" else ""

# input and output paths
results_external_dir <- file.path(paths$results, paste0("external_reference", analysis_suffix))
figures_external_dir <- file.path(
  paths$figures,
  "drafts",
  paste0("external_reference", analysis_suffix))
functional_enrichment_dir <- file.path(
  paths$results,
  paste0("functional_enrichment", analysis_suffix))

dir.create(results_external_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_external_dir, recursive = TRUE, showWarnings = FALSE)

reference_cache_file <- file.path(results_external_dir, "external_rank_reference.rds")

if (!file.exists(reference_cache_file)) {
  prerequisite <- if (gene_scope == "protein_coding") {
    "scripts/07b_external_healthy_reference_protein_coding.R"
  } else {
    "scripts/07_external_healthy_reference.R"
  }
  stop("Missing external-reference cache. Run ", prerequisite, " first.")}

reference_cache <- readRDS(reference_cache_file)
external_rank_matrix <- reference_cache$external_rank_matrix
gse_samples <- reference_cache$gse_samples
blueprint_samples <- reference_cache$blueprint_samples

# cell-line overlaps
one_vs_rest <- readr::read_csv(file.path(functional_enrichment_dir, "one_vs_rest_DESeq2_results.csv"),
  show_col_types = FALSE)

if (gene_scope == "protein_coding" &&
    "gene_type" %in% colnames(one_vs_rest) &&
    any(one_vs_rest$gene_type != "protein_coding", na.rm = TRUE)) {
  stop("The one-vs-rest input contains non-protein-coding genes.")
}

colnames(one_vs_rest)

venn_up <- one_vs_rest |>
  dplyr::filter(!is.na(padj), padj < 0.01, log2FoldChange >= 2) |>
  dplyr::group_by(cell_line) |>
  dplyr::summarise(genes = list(unique(stringr::str_remove(Geneid, "\\.\\d+$"))), .groups = "drop")

venn_down <- one_vs_rest |>
  dplyr::filter(!is.na(padj), padj < 0.01, log2FoldChange <= -2) |>
  dplyr::group_by(cell_line) |>
  dplyr::summarise(genes = list(unique(stringr::str_remove(Geneid, "\\.\\d+$"))), .groups = "drop")

venn_up_list <- stats::setNames(venn_up$genes, venn_up$cell_line)
venn_down_list <- stats::setNames(venn_down$genes, venn_down$cell_line)

p_venn_up <- ggVennDiagram::ggVennDiagram(venn_up_list, label_alpha = 0, edge_size = 0.2, set_size = 3.5,
  label_size = 3.5) + ggplot2::scale_fill_gradient(low = very_light_grey, high = gradient_red["high"]) +
  ggplot2::scale_x_continuous(breaks = NULL) + ggplot2::scale_y_continuous(breaks = NULL) + ggplot2::labs(title = "Upregulated genes") +
  theme_clean() + ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(hjust = 0.5))

p_venn_down <- ggVennDiagram::ggVennDiagram(venn_down_list, label_alpha = 0, edge_size = 0.2, set_size = 3.5,
  label_size = 3.5) + ggplot2::scale_fill_gradient(low = very_light_grey, high = gradient_blue["high"]) +
  ggplot2::scale_x_continuous(breaks = NULL) + ggplot2::scale_y_continuous(breaks = NULL) + ggplot2::labs(title = "Downregulated genes") +
  theme_clean() + ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5))

p_venn_cell_lines <- p_venn_up + p_venn_down

ggplot2::ggsave(file.path(figures_external_dir, "Venn_cell_line_up_down_genes.pdf"), plot = p_venn_cell_lines,
  width = 12, height = height_large)

# healthy-reference overlap
healthy_panT_rank_mean <- rowMeans(external_rank_matrix[, gse_samples, drop = FALSE])
healthy_cd4_rank_mean <- rowMeans(external_rank_matrix[, blueprint_samples, drop = FALSE])

healthy_panT_high <- names(healthy_panT_rank_mean)[healthy_panT_rank_mean >= stats::quantile(healthy_panT_rank_mean,
  0.8)]

healthy_cd4_high <- names(healthy_cd4_rank_mean)[healthy_cd4_rank_mean >= stats::quantile(healthy_cd4_rank_mean,
  0.8)]

c(Healthy_PanT = length(healthy_panT_high), Healthy_CD4 = length(healthy_cd4_high), Healthy_consensus = length(intersect(healthy_panT_high,
  healthy_cd4_high)))

make_healthy_venn <- function(cell_line, direction = c("Up", "Down")) {
  direction <- match.arg(direction)
  gene_list <- if (direction == "Up")
    venn_up_list else venn_down_list
  fill_high <- if (direction == "Up")
    gradient_red["high"] else gradient_blue["high"]

  venn_sets <- list(`  CTCL` = gene_list[[cell_line]], `Healthy Pan T               ` = healthy_panT_high,
    `Healthy CD4` = healthy_cd4_high)

  ggVennDiagram::ggVennDiagram(venn_sets, label_alpha = 0, edge_size = 0.2, set_size = 3.2, label_size = 3.2) +
    ggplot2::scale_fill_gradient(low = very_light_grey, high = fill_high) + ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = NULL) + ggplot2::labs(title = paste(cell_line, direction)) +
    ggplot2::guides(fill = "none") + theme_clean() + ggplot2::theme(axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(), plot.margin = ggplot2::margin(8,
      12, 8, 28), plot.title = ggplot2::element_text(hjust = 0.5))}

p_healthy_hh_up <- make_healthy_venn("HH", "Up")
p_healthy_myla_up <- make_healthy_venn("MyLa", "Up")
p_healthy_hut78_up <- make_healthy_venn("HuT 78", "Up")
p_healthy_seax_up <- make_healthy_venn("SeAx", "Up")

p_healthy_hh_down <- make_healthy_venn("HH", "Down")
p_healthy_myla_down <- make_healthy_venn("MyLa", "Down")
p_healthy_hut78_down <- make_healthy_venn("HuT 78", "Down")
p_healthy_seax_down <- make_healthy_venn("SeAx", "Down")

p_venn_healthy_refs <- (p_healthy_hh_up + p_healthy_myla_up + p_healthy_hut78_up + p_healthy_seax_up)/(p_healthy_hh_down +
  p_healthy_myla_down + p_healthy_hut78_down + p_healthy_seax_down)

ggplot2::ggsave(file.path(figures_external_dir, "Venn_cell_lines_healthy_reference_overlap.pdf"), plot = p_venn_healthy_refs,
  width = 11, height = 8, units = "in")

summarize_healthy_overlap <- function(cell_line, direction) {
  gene_list <- if (direction == "Up")
    venn_up_list else venn_down_list
  cell_line_genes <- gene_list[[cell_line]]

  tibble::tibble(cell_line = cell_line, direction = direction, overlap_PanT = length(intersect(cell_line_genes,
    healthy_panT_high)), overlap_CD4 = length(intersect(cell_line_genes, healthy_cd4_high)), overlap_both = length(Reduce(intersect,
    list(cell_line_genes, healthy_panT_high, healthy_cd4_high))))}

overlap_comparisons <- tidyr::expand_grid(direction = c("Up", "Down"), cell_line = c("HH", "MyLa", "HuT 78",
  "SeAx"))

healthy_overlap_summary <- purrr::map2_dfr(overlap_comparisons$cell_line, overlap_comparisons$direction,
  summarize_healthy_overlap)

healthy_overlap_summary

readr::write_csv(healthy_overlap_summary, file.path(results_external_dir, "cell_line_healthy_reference_gene_overlap_summary.csv"))
