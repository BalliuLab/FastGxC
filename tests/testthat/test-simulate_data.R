library(FastGxC)
test_that("simulate_data creates expected output files", {
  data_dir <- file.path(tempfile(), "")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(data_dir, recursive = TRUE, force = TRUE), add = TRUE)

  simulate_data(
    data_dir = data_dir,
    N = 10,
    n_genes = 5,
    n_snps_per_gene = 10,
    n_contexts = 3,
    seed = 42
  )

  expect_true(file.exists(file.path(data_dir, "expression.txt")))
  expect_true(file.exists(file.path(data_dir, "SNPs.txt")))
  expect_true(file.exists(file.path(data_dir, "snpsloc.txt")))
  expect_true(file.exists(file.path(data_dir, "geneloc.txt")))
})

test_that("simulate_data respects seed for reproducibility", {
  data_dir1 <- file.path(tempfile(), "")
  data_dir2 <- file.path(tempfile(), "")
  dir.create(data_dir1, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir2, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(c(data_dir1, data_dir2), recursive = TRUE, force = TRUE), add = TRUE)

  simulate_data(data_dir = data_dir1, N = 10, n_genes = 5,
                n_snps_per_gene = 10, n_contexts = 3, seed = 42)
  simulate_data(data_dir = data_dir2, N = 10, n_genes = 5,
                n_snps_per_gene = 10, n_contexts = 3, seed = 42)

  exp1 <- read.table(file.path(data_dir1, "expression.txt"), header = TRUE, sep = "\t")
  exp2 <- read.table(file.path(data_dir2, "expression.txt"), header = TRUE, sep = "\t")
  expect_equal(exp1, exp2)
})

test_that("simulate_data fails with invalid parameters", {
  expect_error(simulate_data(data_dir = tempdir(), N = 0, n_genes = 5,
                             n_snps_per_gene = 10, n_contexts = 3))
  expect_error(simulate_data(data_dir = tempdir(), N = 10, n_genes = 5,
                             n_snps_per_gene = 10, n_contexts = 1))
})
