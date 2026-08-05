library(FastGxC)

test_that("decomposition_step creates expected output files", {
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

  decomposition_step(
    exp_mat_filename = file.path(data_dir, "expression.txt"),
    data_dir = data_dir
  )

  expect_true(file.exists(file.path(data_dir, "context_shared_expression.txt")))
  expect_true(file.exists(file.path(data_dir, "context1_specific_expression.txt")))
  expect_true(file.exists(file.path(data_dir, "context2_specific_expression.txt")))
  expect_true(file.exists(file.path(data_dir, "context3_specific_expression.txt")))
})
