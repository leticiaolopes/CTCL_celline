# CTCL cell line bulk RNA-seq

bulk RNA-seq analysis of the CTCL cell lines HH, MyLa, HuT 78 and SeAx.

this repo has the scripts used for metadata checks, FASTQ QC, count processing, sample identity checks, differential expression, functional enrichment and figure generation.

## quick notes

- scripts should be run from the project root;
- the main project paths and sample identities are defined in `scripts/00_config.R`;
- plot colors and themes are in `scripts/00_aesthetics.R`;
- the analysis currently uses the inferred sample identities defined in `00_config.R`;
- raw FASTQ files and the large GENCODE GTF are not stored in GitHub.

## scripts

| script | what it does |
| --- | --- |
| `00_config.R` | project paths, sample names and sample identity settings |
| `00_aesthetics.R` | colors, themes and figure helpers |
| `01_build_metadata.R` | builds sample and file metadata from the facility JSON files |
| `02_fastq_qc.R` | summarizes FASTQ QC from MultiQC |
| `03_build_count_matrix.R` | combines featureCounts files into the count matrix |
| `04_deseq2.R` | normalization, VST, PCA and sample 2 sample QC |
| `04b_sample_identity_qc.R` | checks suspected sample swaps |
| `05_differential_expression.R` | pairwise differential expression and DEG heatmaps |
| `06_functional_enrichment.R` | gene modules, GO enrichment and GSEA |
| `07_external_healthy_reference.R` | comparison with healthy T cell references |
| `08_venn_diagrams.R` | DEG and healthy-reference overlaps |
| `09_heatmap.R` | top 50 DEG heatmap |
| `10_CTCL_sequencing.R` | exports the main sequencing figure |
| `10b_replot_supplementary_panels.R` | rebuilds supplementary panels from saved results |
| `11_supplementary_figures.R` | assembles the supplementary figures |

## folders

- `metadata/`: sample and file metadata;
- `qc/`: sequencing and sample-level QC;
- `counts/`: count matrices and featureCounts outputs;
- `results/`: analysis tables and saved R objects;
- `figures/`: final figures and draft panels;

## running the analysis

open R in the project root, check `project_dir` in `scripts/00_config.R`, install the packages loaded by the scripts and run the scripts in numerical order.

some external-reference steps download public data and can take a while. the raw sequencing files and the GENCODE annotation need to be available locally before running the relevant steps.
