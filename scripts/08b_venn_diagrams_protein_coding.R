# Generate Venn overlaps from protein-coding one-vs-rest and external-reference
# results. Outputs are kept under the protein-coding result and figure folders.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/08_venn_diagrams.R", local = TRUE)
})
