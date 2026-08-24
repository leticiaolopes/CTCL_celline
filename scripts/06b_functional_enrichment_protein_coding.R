# Build gene modules and functional enrichment from protein-coding DE results.
# Outputs are written to results/functional_enrichment_protein_coding and
# figures/drafts/functional_enrichment_protein_coding.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/06_functional_enrichment.R", local = TRUE)
})
