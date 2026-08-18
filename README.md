# scrub-a-dub

Pre- and post-filtering UMAPs for a doublet- or QC-filtering step: one panel
of every raw cell (colored by why it was later dropped), next to one panel
of the cells that survive.

scrub-a-dub does not call doublets or decide what "low-quality" means -- it
draws the before/after picture once you already have a `keep` flag and a
`doublet` flag per cell, from whatever classifier or cutoff you used (a
lab example: the `stdoublet` calls and `global_min_keep` depth cutoff used
in `aBO195`). For calling low-quality cells in the first place, see the
companion package [`duddrop`](../duddrop).

## Install

```r
# from an R session
devtools::install("/n/data1/hms/scrb/chen/lab/bco/packages/scrub-a-dub")
```

## Usage

```r
library(scrubadub)

fig <- sad_prepost_umap(
  raw_pcs      = "raw_rna_pcs.csv.gz",       # barcode + PC columns, every raw cell
  filtered_pcs = "filtered_rna_pcs.csv.gz",  # barcode + PC columns, refit post-QC
  metadata     = metadata,                   # data.frame(barcode, keep, doublet)
  dataset      = "xBO350"
)
ggplot2::ggsave("prepost_doublet_umap.png", fig, width = 12, height = 5.5, dpi = 240)
```

`raw_pcs`/`filtered_pcs` accept either a path (read with `read.csv`) or an
already-loaded data.frame -- any non-`barcode` column is treated as a PC.
`metadata` needs `barcode`, a logical `keep` column (passed general QC),
and a logical `doublet` column (flagged by whatever doublet caller you
used). Points where `doublet` is `TRUE` are dropped from the post panel
regardless of `keep`.

If you already have your own UMAP coordinates (e.g. from a WNN graph
instead of RNA PCA), skip the embedding step and call
`sad_plot_prepost(pre_frame, post_frame, dataset)` directly with
`data.frame(barcode, umap_1, umap_2, category)` /
`data.frame(barcode, umap_1, umap_2)`.

See `examples/` for the runs this was built for: `aBO195`'s pre/post
doublet-removal UMAPs for xBO350 and xBO285.

## Tests

```bash
Rscript -e 'devtools::test("/n/data1/hms/scrb/chen/lab/bco/packages/scrub-a-dub")'
```
