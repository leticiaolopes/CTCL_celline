library(ggplot2)
library(scales)
library(paletteer)

# typography
font <- "sans"

# figure dimensions
width_single <- 90 / 25.4
width_onehalf <- 130 / 25.4
width_double <- 180 / 25.4
width_heatmap <- 120 / 25.4

height_small <- 70 / 25.4
height_medium <- 100 / 25.4
height_large <- 150 / 25.4

dpi <- 300

# core colors
black <- "#222222"
dark_grey <- "#4D4D4D"
mid_grey <- "#8A8A8A"
light_grey <- "#D9D9D9"
very_light_grey <- "#F2F2F2"
white <- "#FFFFFF"

# main categorical and continuous palettes
palette <- paletteer_d("wesanderson::Zissou1")

palette_continuous <- paletteer_d(
  "wesanderson::Zissou1Continuous")

# cell-line palette
cell_lines <- c(
  "HH" = "#3B9AB2FF",
  "MyLa" = "#F21A00FF",
  "HuT 78" = "#BDC881FF",
  "SeAx" = "#EBCC2AFF")

# healthy palette 
external_group_colors <- c(cell_lines, "Healthy Pan T" = "#507d69", "Healthy CD4" = "#052162")

# 2 group palettes

binary <- c("Group 1" = "#7F9DB9", "Group 2" = "#D9918F")

up_down <- c("Down" = "#3F648C", "NS" = "#BDBDBD", "Up" = "#B55350")

# continuous palettes
gradient_blue <- c(
  low = "#F3F6F8",
  mid = "#B6CADB",
  high = "#4F789D")

gradient_red <- c(
  low = "#FAF3F3",
  mid = "#E4B5B3",
  high = "#B45F5B")

gradient_diverging <- c(
  low = "#3F648C",
  mid = "#F7F7F7",
  high = "#B55350")

heatmap_blue <- colorRampPalette(
  c("#F3F6F8", "#C4D4E1", "#7F9DB9", "#4F789D"))(100)

heatmap_purple <- colorRampPalette(
  c("#F7F7F7", "#D9CFE3", "#B19BC7", "#78678F"))(100)

heatmap_deg <- colorRampPalette(
  c("#3F648C",
    "#F7F7F7",
    "#B55350"))(101)

annotation_colors <- list(cell_line = cell_lines)

# main ggplot theme
theme_clean <- function(base_size = 7, base_family = font) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family) +
    ggplot2::theme(text = ggplot2::element_text(
        family = base_family, colour = black),
      axis.title = ggplot2::element_text(
        size = base_size, colour = black),
      axis.text = ggplot2::element_text(
        size = base_size - 0.5, colour = black),
      axis.text.x = ggplot2::element_text(colour = black),
      axis.text.y = ggplot2::element_text(colour = black),
      axis.line = ggplot2::element_line(colour = black, linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = black, linewidth = 0.35),
      axis.ticks.length = grid::unit(1.5, "mm"),
      legend.title = ggplot2::element_text(size = base_size, face = "plain"),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = white, colour = NA),
      plot.background = ggplot2::element_rect(fill = white, colour = NA),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, colour = dark_grey),
      plot.caption = ggplot2::element_text(size = base_size - 1, colour = dark_grey, hjust = 0),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "plain"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 4, 4, 4))}

# minimal version for PCA, scatterplots, UMAPs
theme_minimal_clean <- function(base_size = 7, base_family = font) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(text = ggplot2::element_text(family = base_family, colour = black),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = black),
      panel.grid.major = ggplot2::element_line(colour = "#EAEAEA", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "plain"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain"),
      plot.background = ggplot2::element_rect(fill = white, colour = NA))}

# scale helpers
scale_color_leti <- function(...) {scale_color_manual(values = palette, ...)}

scale_fill_leti <- function(...) {scale_fill_manual(values = palette, ...)}

scale_color_cell_line <- function(...) {scale_color_manual(values = cell_lines, ...)}

scale_fill_cell_line <- function(...) {scale_fill_manual(values = cell_lines, ...)}

scale_color_cell_line_model <- function(...) {scale_color_manual(values = cell_lines_model, ...)}

scale_fill_cell_line_model <- function(...) {scale_fill_manual(values = cell_lines_model, ...)}

scale_color_up_down <- function(...) {scale_color_manual(values = up_down, ...)}

scale_fill_up_down <- function(...) {scale_fill_manual(values = up_down, ...)}

scale_color_diverging <- function(midpoint = 0, limits = NULL, ...) {
  scale_color_gradient2(
    low = gradient_diverging["low"],
    mid = gradient_diverging["mid"],
    high = gradient_diverging["high"],
    midpoint = midpoint,
    limits = limits,
    ...)}

scale_fill_diverging <- function(midpoint = 0,  limits = NULL,
    ...) {
  scale_fill_gradient2(
    low = gradient_diverging["low"],
    mid = gradient_diverging["mid"],
    high = gradient_diverging["high"],
    midpoint = midpoint,
    limits = limits,
    ...)}

# export helpers
save_single <- function(filename, plot, height = height_medium) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width_single,
    height = height,
    dpi = dpi,
    units = "in")}

save_double <- function(filename, plot, height = height_medium) {
  ggplot2::ggsave(filename = filename,
    plot = plot,
    width = width_double,
    height = height,
    dpi = dpi,
    units = "in")}

# default ggplot settings
ggplot2::theme_set(theme_clean())
