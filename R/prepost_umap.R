#' Read a barcode + PCA-columns table
#'
#' Expects a data.frame (or a path readable by [utils::read.csv]) with a
#' `barcode` column and one column per PC, in any name/order -- every
#' non-`barcode` column is treated as a PC.
#'
#' @param x a data.frame, or a path/connection passed to [utils::read.csv]
#' @param n_pcs number of leading PC columns to keep
#' @return list(barcode = character vector, mat = numeric matrix)
#' @export
sad_read_pcs <- function(x, n_pcs = 20L) {
  df <- if (is.data.frame(x)) x else utils::read.csv(x, stringsAsFactors = FALSE)
  if (!"barcode" %in% colnames(df)) stop("expected a 'barcode' column")
  pc_cols <- setdiff(colnames(df), "barcode")
  if (!length(pc_cols)) stop("no PC columns found besides 'barcode'")
  n_pcs <- min(n_pcs, length(pc_cols))
  list(barcode = df$barcode, mat = as.matrix(df[, pc_cols[seq_len(n_pcs)], drop = FALSE]))
}

#' UMAP embed a PC matrix
#'
#' Thin wrapper around [uwot::umap()] with the defaults scrubadub's figures
#' were tuned against; pass `...` to override any of them.
#'
#' @param mat numeric matrix, cells x PCs
#' @param seed random seed for reproducibility
#' @param ... additional arguments forwarded to [uwot::umap()]
#' @return two-column matrix of UMAP coordinates
#' @export
sad_umap <- function(mat, seed = 1L, ...) {
  set.seed(seed)
  args <- utils::modifyList(
    list(n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose = FALSE),
    list(...)
  )
  do.call(uwot::umap, c(list(X = mat), args))
}

#' Pre- and post-filtering UMAP, side by side
#'
#' Embeds `raw` (every cell that entered the QC step) and `filtered` (the
#' cells you actually kept) independently, and returns both panels as one
#' [patchwork] figure. The left panel colors every raw cell by why it was
#' later dropped so you can see where doublets/low-quality cells sit
#' relative to the clusters they blur; the right panel is the resulting
#' embedding once they're gone.
#'
#' scrubadub does not call doublets or decide what "low-quality" means --
#' `metadata` must already carry those calls (e.g. from your own doublet
#' classifier, `duddrop`, or a manual cutoff).
#'
#' @param raw_pcs list(barcode, mat) from [sad_read_pcs()], or a data.frame
#'   with a `barcode` column plus PC columns, covering every raw cell
#' @param filtered_pcs same shape as `raw_pcs`, but only the retained cells,
#'   refit however you refit it (fresh PCA on the subset is typical)
#' @param metadata data.frame with a `barcode` column plus `keep` and
#'   `doublet` logical columns, covering (at least) every cell in `raw_pcs`
#' @param dataset label used in panel titles
#' @param keep_col name of the logical "passed non-doublet QC" column
#' @param doublet_col name of the logical "flagged doublet" column
#' @param n_pcs number of leading PCs to embed on
#' @param seed random seed, reused (and +1'd for the post panel) for
#'   reproducibility
#' @return a [patchwork] object (pre panel | post panel); `print()` or
#'   [ggplot2::ggsave()] it
#' @export
sad_prepost_umap <- function(raw_pcs, filtered_pcs, metadata,
                              dataset = "dataset",
                              keep_col = "keep", doublet_col = "doublet",
                              n_pcs = 20L, seed = 1L) {
  is_parsed_pcs <- function(x) is.list(x) && all(c("barcode", "mat") %in% names(x))
  if (!is_parsed_pcs(raw_pcs)) raw_pcs <- sad_read_pcs(raw_pcs, n_pcs)
  if (!is_parsed_pcs(filtered_pcs)) filtered_pcs <- sad_read_pcs(filtered_pcs, n_pcs)

  required <- c("barcode", keep_col, doublet_col)
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) stop("metadata is missing column(s): ", paste(missing, collapse = ", "))

  meta <- metadata[, required]
  colnames(meta) <- c("barcode", "keep", "doublet")
  meta$keep <- as.logical(meta$keep)
  meta$doublet <- as.logical(meta$doublet)
  meta$doublet[is.na(meta$doublet)] <- FALSE

  meta$category <- ifelse(
    meta$doublet, "doublet",
    ifelse(!meta$keep, "low-quality (excluded)", "retained")
  )
  meta$category <- factor(meta$category,
                           levels = c("retained", "low-quality (excluded)", "doublet"))

  pre_emb <- sad_umap(raw_pcs$mat, seed = seed)
  pre_frame <- data.frame(barcode = raw_pcs$barcode, umap_1 = pre_emb[, 1], umap_2 = pre_emb[, 2])
  pre_frame <- merge(pre_frame, meta[, c("barcode", "category")], by = "barcode")
  if (!nrow(pre_frame)) stop("no raw_pcs barcodes matched metadata$barcode")

  keep_barcodes <- meta$barcode[meta$keep & !meta$doublet]
  keep_idx <- filtered_pcs$barcode %in% keep_barcodes
  post_mat <- filtered_pcs$mat[keep_idx, , drop = FALSE]
  if (!nrow(post_mat)) stop("no filtered_pcs barcodes were both kept and non-doublet")
  post_emb <- sad_umap(post_mat, seed = seed + 1L)
  post_frame <- data.frame(
    barcode = filtered_pcs$barcode[keep_idx],
    umap_1 = post_emb[, 1], umap_2 = post_emb[, 2]
  )

  sad_plot_prepost(pre_frame, post_frame, dataset = dataset)
}

#' Assemble the pre/post panels from already-embedded coordinates
#'
#' Exposed separately from [sad_prepost_umap()] so you can supply your own
#' UMAP (e.g. computed on a WNN graph) instead of the RNA-PCA UMAP that
#' function runs internally.
#'
#' @param pre_frame data.frame(barcode, umap_1, umap_2, category) for every
#'   raw cell; `category` a factor with levels "retained",
#'   "low-quality (excluded)", "doublet"
#' @param post_frame data.frame(barcode, umap_1, umap_2) for the retained cells
#' @param dataset label used in panel titles
#' @return a [patchwork] object
#' @importFrom ggplot2 .data
#' @export
sad_plot_prepost <- function(pre_frame, post_frame, dataset = "dataset") {
  cat_colors <- c(
    "retained" = "#7F7F7F",
    "low-quality (excluded)" = "#EE7733",
    "doublet" = "#CC3311"
  )
  base_theme <- ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  n_doublet <- sum(pre_frame$category == "doublet")
  n_lowq <- sum(pre_frame$category == "low-quality (excluded)")

  p_pre <- ggplot2::ggplot(pre_frame, ggplot2::aes(.data$umap_1, .data$umap_2)) +
    ggplot2::geom_point(
      data = pre_frame[pre_frame$category == "retained", ],
      color = cat_colors[["retained"]], size = 0.35, alpha = 0.4
    ) +
    ggplot2::geom_point(
      data = pre_frame[pre_frame$category != "retained", ],
      ggplot2::aes(color = .data$category), size = 0.55, alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = cat_colors, drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s: pre-filtering (n=%d)", dataset, nrow(pre_frame)),
      subtitle = sprintf("%d doublet, %d low-quality flagged", n_doublet, n_lowq),
      color = NULL
    ) +
    ggplot2::coord_equal() + base_theme

  p_post <- ggplot2::ggplot(post_frame, ggplot2::aes(.data$umap_1, .data$umap_2)) +
    ggplot2::geom_point(color = "#0072B2", size = 0.4, alpha = 0.55) +
    ggplot2::labs(
      title = sprintf("%s: post-filtering (n=%d)", dataset, nrow(post_frame)),
      subtitle = "low-quality + doublet cells removed"
    ) +
    ggplot2::coord_equal() + base_theme

  patchwork::wrap_plots(p_pre, p_post)
}
