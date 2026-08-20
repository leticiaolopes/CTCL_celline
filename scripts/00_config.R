project_dir <- "D:/Collabs/eva/CTCL_bulk_RNAseq"

paths <- list(
  data = file.path(project_dir, "data"),
  metadata = file.path(project_dir, "metadata"),
  scripts = file.path(project_dir, "scripts"),
  qc = file.path(project_dir, "qc"),
  counts = file.path(project_dir, "counts"),
  featurecounts = file.path(project_dir, "counts", "featurecounts"),
  featurecounts_raw = file.path(project_dir, "counts", "featurecounts", "Counts"),
  results = file.path(project_dir, "results"),
  figures = file.path(project_dir, "figures"))

sample_metadata <- tibble::tribble(
  ~sample_id,      ~sample_name, ~original_cell_line,
  "25060a001_01",  "HH-1",       "HH",
  "25060a002_01",  "HH-2",       "HH",
  "25060a003_01",  "HH-3",       "HH",
  "25060a004_01",  "M-1",        "MyLa",
  "25060a005_01",  "M-2",        "MyLa",
  "25060a006_01",  "M-3",        "MyLa",
  "25060a007_01",  "78-1",       "HuT 78",
  "25060a008_01",  "78-2",       "HuT 78",
  "25060a009_01",  "78-3",       "HuT 78",
  "25060a010_01",  "S-1",        "SeAx",
  "25060a011_01",  "S-2",        "SeAx",
  "25060a012_01",  "S-3",        "SeAx") |>
  dplyr::mutate(inferred_cell_line = dplyr::case_when(
      sample_name == "HH-3" ~ "MyLa",
      sample_name == "M-3" ~ "HH",
      sample_name == "78-3" ~ "SeAx",
      sample_name == "S-3" ~ "HuT 78",
      TRUE ~ original_cell_line),
    identity_status = dplyr::case_when(
      sample_name %in% c("HH-3", "M-3", "78-3", "S-3") ~ "suspected_swap",
      TRUE ~ "consistent"))

use_inferred_identity <- TRUE

sample_metadata <- sample_metadata |>
  dplyr::mutate(
    analysis_cell_line = if (use_inferred_identity) {
      inferred_cell_line} else {original_cell_line})

sample_metadata <- sample_metadata |>
  dplyr::mutate(analysis_sample_name = dplyr::case_when(
      sample_name == "HH-3" ~ "M-3",
      sample_name == "M-3" ~ "HH-3",
      sample_name == "78-3" ~ "S-3",
      sample_name == "S-3" ~ "78-3",
      TRUE ~ sample_name))
