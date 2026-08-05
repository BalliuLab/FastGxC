library(FastGxC)
test_that("simulate_data creates expected output files", {
  data_dir <- tempfile()
  dir.create(data_dir)
  data_dir <- paste0(data_dir, "/")

  simulate_data(
    data_dir = data_dir,
    N = 10,
    n_genes = 5,
    n_snps_per_gene = 10,
    n_contexts = 3,
    seed = 42
  )

  expect_true(file.exists(paste0(data_dir, "expression.txt")))
  expect_true(file.exists(paste0(data_dir, "SNPs.txt")))
  expect_true(file.exists(paste0(data_dir, "snpsloc.txt")))
  expect_true(file.exists(paste0(data_dir, "geneloc.txt")))
})

test_that("simulate_data respects seed for reproducibility", {
  data_dir1 <- paste0(tempfile(), "/")
  data_dir2 <- paste0(tempfile(), "/")
  dir.create(data_dir1)
  dir.create(data_dir2)

  simulate_data(data_dir = data_dir1, N = 10, n_genes = 5,
                n_snps_per_gene = 10, n_contexts = 3, seed = 42)
  simulate_data(data_dir = data_dir2, N = 10, n_genes = 5,
                n_snps_per_gene = 10, n_contexts = 3, seed = 42)

  exp1 <- read.table(paste0(data_dir1, "expression.txt"), header = TRUE, sep = "\t")
  exp2 <- read.table(paste0(data_dir2, "expression.txt"), header = TRUE, sep = "\t")
  expect_equal(exp1, exp2)
})

test_that("simulate_data fails with invalid parameters", {
  expect_error(simulate_data(data_dir = tempdir(), N = 0, n_genes = 5,
                             n_snps_per_gene = 10, n_contexts = 3))
  expect_error(simulate_data(data_dir = tempdir(), N = 10, n_genes = 5,
                             n_snps_per_gene = 10, n_contexts = 1))
})
