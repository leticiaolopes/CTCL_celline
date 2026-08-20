# build sample & metadata from the sequencing

# setup
source("scripts/00_config.R")

library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(tidyr)
library(readr)

metadata_dir <- file.path(paths$data, "metadata_json")

metadata_files <- list.files(metadata_dir, pattern = "\\.metadata\\.json$", full.names = TRUE)

metadata_files
length(metadata_files)

# parse sample metadata
read_sample_metadata <- function(file) {
  x <- fromJSON(file)

  qc_table <- x$sample1$qc

  read_count <- qc_table |>
    filter(name == "read count") |>
    pull(value)

  tibble(sample_id = x$sample1$id_genetics, prepared_sample_id = x$sample1$id_qbic, data_type = x$type,
    processing_system = x$sample1$processing_system, genome = x$sample1$genome, tumor_core_annotation = x$sample1$tumor,
    facility_read_count = as.numeric(read_count), metadata_file = basename(file))}

sample_metadata_raw <- map_dfr(metadata_files, read_sample_metadata)

sample_metadata_raw

nrow(sample_metadata_raw)

sample_metadata_raw |>
  select(sample_id, prepared_sample_id, processing_system, genome, facility_read_count)

sample_metadata_raw |>
  arrange(sample_id)

experimental_metadata <- tribble(~sample_id, ~prepared_sample_id_expected, ~sample_name, ~cell_line,
  "25060a001_01", "QNJNX013AD", "HH-1", "HH", "25060a002_01", "QNJNX014AL", "HH-2", "HH", "25060a003_01",
  "QNJNX015AT", "HH-3", "HH", "25060a004_01", "QNJNX016A3", "M-1", "M", "25060a005_01", "QNJNX017AB",
  "M-2", "M", "25060a006_01", "QNJNX018AJ", "M-3", "M", "25060a007_01", "QNJNX019AR", "78-1", "78",
  "25060a008_01", "QNJNX020AU", "78-2", "78", "25060a009_01", "QNJNX021A4", "78-3", "78", "25060a010_01",
  "QNJNX022AC", "S-1", "S", "25060a011_01", "QNJNX023AK", "S-2", "S", "25060a012_01", "QNJNX024AS",
  "S-3", "S")

# reconcile experimental sample identifiers
sample_metadata <- experimental_metadata |>
  left_join(sample_metadata_raw, by = "sample_id")

sample_metadata |>
  select(sample_id, sample_name, cell_line, prepared_sample_id_expected, prepared_sample_id, facility_read_count)

sample_metadata |>
  filter(!is.na(prepared_sample_id), prepared_sample_id != prepared_sample_id_expected)

missing_metadata <- sample_metadata |>
  filter(is.na(metadata_file))

missing_metadata

# export sample metadata
write_csv(sample_metadata, file.path(paths$metadata, "sample_metadata.csv"))

read_file_metadata <- function(file) {
  x <- fromJSON(file)

  tibble(sample_id = x$sample1$id_genetics, filename = x$files)}

file_metadata_raw <- metadata_files |>
  map_dfr(read_file_metadata)

file_metadata_raw

file_metadata <- file_metadata_raw |>
  mutate(sequencing_id = str_extract(filename, "S[0-9]+"), lane = str_extract(filename, "L[0-9]{3}"),
    file_type = case_when(str_detect(filename, "_R1_") ~ "R1", str_detect(filename, "_R2_") ~ "R2",
      str_detect(filename, "_index_") ~ "index", TRUE ~ "unknown")) |>
  arrange(sample_id, lane, file_type)

# check fastq
file_metadata |>
  print(n = Inf, width = Inf)

file_check <- file_metadata |>
  dplyr::count(sample_id, file_type) |>
  tidyr::pivot_wider(names_from = file_type, values_from = n, values_fill = 0)

file_check

lane_check <- file_metadata |>
  dplyr::distinct(sample_id, lane) |>
  dplyr::count(sample_id, name = "n_lanes")

lane_check

# export file metadata & completeness checks
write_csv(file_metadata, file.path(paths$metadata, "file_metadata.csv"))

write_csv(file_check, file.path(paths$metadata, "file_check.csv"))

write_csv(lane_check, file.path(paths$metadata, "lane_check.csv"))

cat("\nSamples expected:", nrow(sample_metadata), "\nSamples with metadata:", sum(!is.na(sample_metadata$metadata_file)),
  "\nSamples missing metadata:", sum(is.na(sample_metadata$metadata_file)), "\nFASTQ records:", nrow(file_metadata),
  "\n")
