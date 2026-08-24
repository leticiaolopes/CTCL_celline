# export final ctcl sequencing figure from the approved png
source("scripts/00_config.R")

gene_scope <- getOption("ctcl.gene_scope", "all")
if (!gene_scope %in% c("all", "protein_coding")) {
  stop("Unsupported ctcl.gene_scope: ", gene_scope)
}
analysis_suffix <- if (gene_scope == "protein_coding") "_protein_coding" else ""

figure_png <- file.path(paths$figures, paste0("CTCL_sequencing", analysis_suffix, ".png"))
figure_pdf <- file.path(paths$figures, paste0("CTCL_sequencing", analysis_suffix, ".pdf"))

if (!file.exists(figure_png)) {
  stop("Approved main-figure PNG not found: ", figure_png)
}
figure_image <- magick::image_read(figure_png)
figure_info <- magick::image_info(figure_image)
expected_ratio <- 210/297
observed_ratio <- figure_info$width/figure_info$height

grDevices::cairo_pdf(filename = figure_pdf, width = 210/25.4, height = 297/25.4, family = "sans", onefile = FALSE)

grid::grid.newpage()
grid::grid.raster(as.raster(figure_image), x = 0.5, y = 0.5, width = 1, height = 1, interpolate = FALSE)

grDevices::dev.off()

pdf_size <- file.info(figure_pdf)$size

