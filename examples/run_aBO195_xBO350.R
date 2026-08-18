# Reproduces the xBO350 panel of aBO195's pre/post doublet-removal UMAP,
# reading straight from the raw/filtered PCA exports that
# aBO195/code/prepare_adaptive_cell_qc.R already writes to
# aBO195/ints/adaptive_qc/. Doublet flag: spatial_multiplet (bead-based
# spatial-multiplet call). Keep flag: global_min_keep (read-depth cutoff).

library(scrubadub)

ints_dir <- "/n/data1/hms/scrb/chen/lab/bco/analyses/aBO195/ints/adaptive_qc"

meta <- read.csv(gzfile(file.path(ints_dir, "raw_metadata.csv.gz")), stringsAsFactors = FALSE)
meta$keep <- as.logical(meta$global_min_keep)
meta$doublet <- as.logical(meta$spatial_multiplet)

fig <- sad_prepost_umap(
  raw_pcs      = file.path(ints_dir, "raw_rna_pcs.csv.gz"),
  filtered_pcs = file.path(ints_dir, "filtered_rna_pcs.csv.gz"),
  metadata     = meta,
  dataset      = "xBO350"
)

ggplot2::ggsave("xBO350_prepost_doublet_umap.png", fig, width = 12, height = 5.5, dpi = 240, bg = "white")
