
context("Testing DGE analysis in microarray data")


test_that("OMA DGE anlayisis produces similar results to running DGE by hand", {
  expect_equal(
    d_dge$GSE10761$coeff, limma_manual$coefficients[, 1] %>% unname,
  )
  expect_equal(
    d_dge$GSE10761$p_value, limma_manual$p.value[, 1] %>% unname,
  )
})

