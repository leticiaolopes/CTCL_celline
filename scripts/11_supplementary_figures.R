# assemble supplementary figures 

# setup
source("scripts/00_config.R")

if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required. Install it with install.packages('magick').")}

gene_scope <- getOption("ctcl.gene_scope", "all")
if (!gene_scope %in% c("all", "protein_coding")) {
  stop("Unsupported ctcl.gene_scope: ", gene_scope)
}
analysis_suffix <- if (gene_scope == "protein_coding") "_protein_coding" else ""

supplementary_dir <- file.path(paths$figures, paste0("supplementary", analysis_suffix))
dir.create(supplementary_dir, recursive = TRUE, showWarnings = FALSE)

deseq2_figures_dir <- file.path(paths$figures, "drafts", paste0("deseq2", analysis_suffix))
external_figures_dir <- file.path(paths$figures, "drafts", paste0("external_reference", analysis_suffix))
functional_figures_dir <- file.path(paths$figures, "drafts", paste0("functional_enrichment", analysis_suffix))

sample_identity_dir <- file.path(deseq2_figures_dir, "sample_identity")

differential_expression_dir <- file.path(deseq2_figures_dir, "differential_expression")

render_density <- 300

a4_portrait <- c(width = 2480, height = 3508)
a4_landscape <- c(width = 3508, height = 2480)

# helpers
read_pdf_panel <- function(filename, density = render_density) {
  if (!file.exists(filename)) {
    stop("Source figure not found: ", filename)}

  pdftoppm_bin <- Sys.which("pdftoppm")

  if (pdftoppm_bin == "") {
    stop("Poppler 'pdftoppm' was not found on PATH.")}

  temporary_prefix <- tempfile(pattern = "supplementary_panel_")
  temporary_png <- paste0(temporary_prefix, ".png")

  render_status <- system2(command = pdftoppm_bin, args = c("-f", "1", "-singlefile", "-png", "-r",
    as.character(density), shQuote(normalizePath(filename, winslash = "/")), shQuote(temporary_prefix)),
    stdout = FALSE, stderr = FALSE)
  
  if (!file.exists(temporary_png)) {
    stop("Could not render source PDF: ", filename)}

  panel <- magick::image_read(temporary_png)[1]
  unlink(temporary_png)

  panel <- magick::image_trim(panel, fuzz = 2)

  if (basename(filename) == "Venn_cell_lines_healthy_reference_overlap.pdf") {
    panel <- magick::image_border(panel, color = "white", geometry = "120x55")}

  panel}

place_panel <- function(canvas, panel, column, row, ncol, nrow, label, title = "", row_weights = rep(1,
  nrow), outer_margin = 45, inner_margin = 28) {
  canvas_info <- magick::image_info(canvas)
  cell_width <- floor((canvas_info$width - 2 * outer_margin)/ncol)
  usable_height <- canvas_info$height - 2 * outer_margin
  row_heights <- floor(usable_height * row_weights/sum(row_weights))
  row_heights[length(row_heights)] <- usable_height - sum(row_heights[-length(row_heights)])
  cell_height <- row_heights[row]

  available_width <- cell_width - 2 * inner_margin
  available_height <- cell_height - 2 * inner_margin

  panel_info <- magick::image_info(panel)
  resize_factor <- min(available_width/panel_info$width, available_height/panel_info$height)

  resized_width <- max(1, floor(panel_info$width * resize_factor))
  resized_height <- max(1, floor(panel_info$height * resize_factor))

  panel <- magick::image_resize(panel, geometry = paste0(resized_width, "x", resized_height, "!"))

  cell_left <- outer_margin + (column - 1) * cell_width
  cell_top <- outer_margin + if (row == 1)
    0 else sum(row_heights[seq_len(row - 1)])

  panel_left <- cell_left + floor((cell_width - resized_width)/2)
  panel_top <- cell_top + floor((cell_height - resized_height)/2)

  canvas <- magick::image_composite(canvas, panel, offset = paste0("+", panel_left, "+", panel_top),
    operator = "over")

  label_image <- magick::image_blank(width = 100, height = 100, color = "transparent") |>
    magick::image_annotate(text = label, gravity = "northwest", location = "+0+0", size = 44, font = "Arial",
      weight = 700, color = "#222222")

  canvas <- magick::image_composite(canvas, label_image, offset = paste0("+", cell_left + 2, "+", cell_top +
    2), operator = "over")

  if (nzchar(title)) {
    title_image <- magick::image_blank(width = max(200, cell_width - 120), height = 64, color = "transparent") |>
      magick::image_annotate(text = title, gravity = "north", location = "+0+0", size = 27, font = "Arial",
        weight = 700, color = "#222222")

    canvas <- magick::image_composite(canvas, title_image, offset = paste0("+", cell_left + 70, "+",
      cell_top + 2), operator = "over")}

  canvas}

assemble_figure <- function(source_files, output_stem, ncol, nrow, page_size = a4_portrait, panel_titles = rep("",
  length(source_files)), row_weights = rep(1, nrow), show_labels = TRUE) {
  expected_panels <- ncol * nrow

  if (length(source_files) != expected_panels) {
    stop(output_stem, " requires ", expected_panels, " panels; received ", length(source_files),
      ".")}

  canvas <- magick::image_blank(width = unname(page_size["width"]), height = unname(page_size["height"]),
    color = "white")

  panel_labels <- LETTERS[seq_along(source_files)]

  for (index in seq_along(source_files)) {
    panel <- read_pdf_panel(source_files[[index]])
    row <- ceiling(index/ncol)
    column <- index - (row - 1) * ncol

    canvas <- place_panel(canvas = canvas, panel = panel, column = column, row = row, ncol = ncol,
      nrow = nrow, label = if (show_labels)
        panel_labels[index] else "", title = panel_titles[index], row_weights = row_weights)}

  png_file <- file.path(supplementary_dir, paste0(output_stem, ".png"))

  pdf_file <- file.path(supplementary_dir, paste0(output_stem, ".pdf"))

  magick::image_write(canvas, path = png_file, format = "png", density = "300x300")

  page_width_in <- unname(page_size["width"])/300
  page_height_in <- unname(page_size["height"])/300

  grDevices::cairo_pdf(filename = pdf_file, width = page_width_in, height = page_height_in, family = "sans",
    onefile = FALSE)

  grid::grid.newpage()
  grid::grid.raster(as.raster(canvas), x = 0.5, y = 0.5, width = 1, height = 1, interpolate = FALSE)

  grDevices::dev.off()

  if (!file.exists(pdf_file) || file.info(pdf_file)$size == 0) {
    stop("Failed to create supplementary PDF: ", pdf_file)}

  message("Created: ", pdf_file)
  invisible(pdf_file)}

# sup fig 1: qc
s1_files <- c(file.path(deseq2_figures_dir, "PCA_vst_original_identity.pdf"), file.path(deseq2_figures_dir,
  "PCA_vst_inferred_identity.pdf"), file.path(deseq2_figures_dir,
  "sample_correlation_original_identity.pdf"), file.path(deseq2_figures_dir, "sample_correlation_inferred_identity.pdf"),
  file.path(deseq2_figures_dir, "sample_distance_original_identity.pdf"), file.path(deseq2_figures_dir,
    "sample_distance_inferred_identity.pdf"), file.path(sample_identity_dir,
    "internal_identity_signature_original.pdf"), file.path(sample_identity_dir, "internal_identity_signature_inferred.pdf"))

assemble_figure(source_files = s1_files, output_stem = "Supplementary_Figure_S1_QC_sample_identity",
  ncol = 2, nrow = 4, page_size = a4_portrait, panel_titles = c("Before identity correction", "After identity correction",
    "Sample correlation - before", "Sample correlation - after", "Sample distance - before", "Sample distance - after",
    "Identity signature - before", "Identity signature - after"))

# sup figure s2: differential expression
s2_files <- c(file.path(differential_expression_dir, "top_DEGs_global_heatmap.pdf"), file.path(differential_expression_dir,
  "all_DEGs_large_heatmap.pdf"), file.path(differential_expression_dir, "volcano_MyLa_vs_HH.pdf"),
  file.path(differential_expression_dir, "volcano_HuT78_vs_HH.pdf"), file.path(differential_expression_dir,
    "volcano_SeAx_vs_HH.pdf"), file.path(differential_expression_dir, "volcano_HuT78_vs_MyLa.pdf"),
  file.path(differential_expression_dir, "volcano_SeAx_vs_MyLa.pdf"), file.path(differential_expression_dir,
    "volcano_SeAx_vs_HuT78.pdf"))

assemble_figure(source_files = s2_files, output_stem = "Supplementary_Figure_S2_differential_expression",
  ncol = 2, nrow = 4, page_size = a4_portrait)

# sup figure s3: external healthy reference validation
s3_files <- c(file.path(external_figures_dir, "cross_dataset_correlation_heatmap.pdf"),
  file.path(external_figures_dir, "BLUEPRINT_CD4_PCA.pdf"), file.path(external_figures_dir,
    "BLUEPRINT_CD4_PCA_top500_variable_genes.pdf"), file.path(external_figures_dir,
    "GSE197067_PanT_0h_PCA.pdf"), file.path(external_figures_dir,
    "cross_dataset_rank_PCA.pdf"), file.path(external_figures_dir,
    "cross_dataset_spearman_MDS.pdf"))

assemble_figure(source_files = s3_files, output_stem = "Supplementary_Figure_S3_external_reference_validation",
  ncol = 3, nrow = 2, page_size = a4_landscape, panel_titles = c("", "", "", "", "Cross-dataset rank PCA",
    "Cross-dataset Spearman MDS"))

# sup figure s4: modules and functional enrichment
s4_files <- c(file.path(functional_figures_dir, "gene_module_expression_patterns.pdf"),
  file.path(external_figures_dir, "gene_module_activity_healthy_references.pdf"),
  file.path(functional_figures_dir, "GO_BP_gene_modules_rrvgo_dotplot.pdf"),
  file.path(functional_figures_dir, "GO_BP_GSEA_one_vs_rest.pdf"))

assemble_figure(source_files = s4_files, output_stem = "Supplementary_Figure_S4_modules_functional_enrichment",
  ncol = 2, nrow = 2, page_size = a4_landscape, panel_titles = rep("", 4), row_weights = c(0.3, 0.7))

# sup figure s5: gene overlaps
s5_files <- c(file.path(external_figures_dir, "Venn_cell_lines_healthy_reference_overlap.pdf"))

assemble_figure(source_files = s5_files, output_stem = "Supplementary_Figure_S5_gene_overlaps", ncol = 1,
  nrow = 1, page_size = a4_landscape, show_labels = FALSE)

message("All supplementary figures were saved in: ", supplementary_dir)
