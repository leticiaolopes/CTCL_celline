# Compare CTCL cell lines with healthy references using the protein-coding
# DESeq2 object and protein-coding-derived gene modules. Cross-dataset analyses
# inherit the protein-coding universe from the fitted cell-line object.
# Outputs are written to results/external_reference_protein_coding and
# figures/drafts/external_reference_protein_coding.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/07_external_healthy_reference.R", local = TRUE)
})
