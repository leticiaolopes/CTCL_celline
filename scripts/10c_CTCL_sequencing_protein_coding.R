# Export the approved protein-coding main-figure PNG as the final PDF.
# The PNG must first be assembled and approved as
# figures/CTCL_sequencing_protein_coding.png.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/10_CTCL_sequencing.R", local = TRUE)
})
