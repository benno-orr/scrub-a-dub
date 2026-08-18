make_toy <- function(n_retained = 40, n_lowq = 8, n_doublet = 6, n_pcs = 5, seed = 1) {
  set.seed(seed)
  n <- n_retained + n_lowq + n_doublet
  barcode <- sprintf("cell_%03d", seq_len(n))
  category <- c(
    rep("retained", n_retained),
    rep("low-quality (excluded)", n_lowq),
    rep("doublet", n_doublet)
  )
  mat <- matrix(rnorm(n * n_pcs), nrow = n)
  colnames(mat) <- paste0("PC_", seq_len(n_pcs))
  raw_df <- data.frame(barcode = barcode, mat, check.names = FALSE)

  keep <- category != "low-quality (excluded)"
  doublet <- category == "doublet"
  metadata <- data.frame(barcode = barcode, keep = keep, doublet = doublet)

  filtered_barcodes <- barcode[keep & !doublet]
  filt_df <- raw_df[raw_df$barcode %in% filtered_barcodes, ]

  list(raw_df = raw_df, filt_df = filt_df, metadata = metadata,
       n_retained = n_retained, n_lowq = n_lowq, n_doublet = n_doublet)
}

test_that("sad_read_pcs parses a barcode + PC data.frame", {
  toy <- make_toy()
  pcs <- sad_read_pcs(toy$raw_df, n_pcs = 3)
  expect_equal(length(pcs$barcode), nrow(toy$raw_df))
  expect_equal(ncol(pcs$mat), 3)
})

test_that("sad_read_pcs requires a barcode column", {
  expect_error(sad_read_pcs(data.frame(PC_1 = 1:3)), "barcode")
})

test_that("sad_prepost_umap embeds raw cells and drops doublets + low-quality from post", {
  toy <- make_toy()
  fig <- sad_prepost_umap(
    toy$raw_df, toy$filt_df, toy$metadata,
    dataset = "toy", n_pcs = 5
  )
  expect_s3_class(fig, "patchwork")

  built <- ggplot2::ggplot_build(fig[[1]])
  n_pre_points <- sum(vapply(built$data, nrow, integer(1)))
  expect_equal(n_pre_points, toy$n_retained + toy$n_lowq + toy$n_doublet)

  built_post <- ggplot2::ggplot_build(fig[[2]])
  expect_equal(nrow(built_post$data[[1]]), toy$n_retained)
})

test_that("sad_prepost_umap errors when metadata lacks required columns", {
  toy <- make_toy()
  bad_meta <- toy$metadata[, "barcode", drop = FALSE]
  expect_error(
    sad_prepost_umap(toy$raw_df, toy$filt_df, bad_meta, n_pcs = 5),
    "missing column"
  )
})

test_that("sad_prepost_umap accepts file paths, not just data.frames", {
  toy <- make_toy()
  raw_path <- tempfile(fileext = ".csv")
  filt_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(raw_path, filt_path)), add = TRUE)
  utils::write.csv(toy$raw_df, raw_path, row.names = FALSE)
  utils::write.csv(toy$filt_df, filt_path, row.names = FALSE)

  fig <- sad_prepost_umap(raw_path, filt_path, toy$metadata, dataset = "toy", n_pcs = 5)
  expect_s3_class(fig, "patchwork")
})

test_that("sad_prepost_umap errors when no barcodes overlap", {
  toy <- make_toy()
  bad_meta <- toy$metadata
  bad_meta$barcode <- paste0("nomatch_", bad_meta$barcode)
  expect_error(
    sad_prepost_umap(toy$raw_df, toy$filt_df, bad_meta, n_pcs = 5),
    "no raw_pcs barcodes matched"
  )
})
