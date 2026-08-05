library(FastGxC)

test_that("full pipeline runs end to end", {
  data_dir <- paste0(tempfile(), "/")
  out_dir <- paste0(tempfile(), "/")
  dir.create(data_dir)
  dir.create(out_dir)

  # Step 1: simulate data
  simulate_data(
    data_dir = data_dir,
    N = 20,
    n_genes = 5,
    n_snps_per_gene = 10,
    n_contexts = 3,
    seed = 42
  )

  # Step 2: decompose expression
  decomposition_step(
    exp_mat_filename = paste0(data_dir, "expression.txt"),
    data_dir = data_dir
  )

  # Step 3: eQTL mapping for each context and shared
  context_names <- c("context1", "context2", "context3")
  all_contexts <- c(context_names, "shared")

  for (ctx in all_contexts) {
    tag <- ifelse(ctx == "shared", "shared", "specific")
    exp_file <- ifelse(ctx == "shared",
                       paste0(data_dir, "context_shared_expression.txt"),
                       paste0(data_dir, ctx, "_specific_expression.txt"))
    eQTL_mapping_step(
      SNP_file_name = paste0(data_dir, "SNPs.txt"),
      snps_location_file_name = paste0(data_dir, "snpsloc.txt"),
      expression_file_name = exp_file,
      gene_location_file_name = paste0(data_dir, "geneloc.txt"),
      context = ctx,
      out_dir = out_dir,
      output_file_name_cis = paste0(out_dir, ctx, "_", tag, ".cis_pairs.txt"),
      output_file_name_tra = paste0(out_dir, ctx, "_", tag, ".trans_pairs.txt")
    )
  }

  # Check eQTL output files exist
  for (ctx in all_contexts) {
    tag <- ifelse(ctx == "shared", "shared", "specific")
    expect_true(file.exists(paste0(out_dir, ctx, "_", tag, ".cis_pairs.txt")))
  }

  # Step 4: treeQTL multiple testing correction
  treeQTL_step(
    data_dir = out_dir,
    snps_location_file_name = paste0(data_dir, "snpsloc.txt"),
    gene_location_file_name = paste0(data_dir, "geneloc.txt"),
    context_names = context_names,
    out_dir = out_dir
  )

  expect_true(file.exists(file.path(out_dir, "specific_eGenes.txt")))
  expect_true(file.exists(file.path(out_dir, "shared_eGenes.txt")))
})
