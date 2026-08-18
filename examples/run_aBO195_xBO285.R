# Reproduces the xBO285 panel of aBO195's pre/post doublet-removal UMAP.
# Doublet flag here is stdbl_call, the audited `stdoublet` classifier's
# doublet call (aBO195/ints/adaptive_qc_xBO285/doublet_scores.csv.gz),
# not the older spatial_multiplet heuristic used for xBO350.

library(scrubadub)

ints_dir <- "/n/data1/hms/scrb/chen/lab/bco/analyses/aBO195/ints/adaptive_qc_xBO285"

meta <- read.csv(gzfile(file.path(ints_dir, "raw_metadata.csv.gz")), stringsAsFactors = FALSE)
meta$keep <- as.logical(meta$global_min_keep)

dbl <- read.csv(gzfile(file.path(ints_dir, "doublet_scores.csv.gz")), stringsAsFactors = FALSE)
dbl$doublet <- toupper(as.character(dbl$stdbl_call)) == "TRUE"
meta <- merge(meta, dbl[, c("barcode", "doublet")], by = "barcode", all.x = TRUE)
meta$doublet[is.na(meta$doublet)] <- FALSE

fig <- sad_prepost_umap(
  raw_pcs      = file.path(ints_dir, "raw_rna_pcs.csv.gz"),
  filtered_pcs = file.path(ints_dir, "filtered_rna_pcs.csv.gz"),
  metadata     = meta,
  dataset      = "xBO285"
)

ggplot2::ggsave("xBO285_prepost_doublet_umap.png", fig, width = 12, height = 5.5, dpi = 240, bg = "white")
