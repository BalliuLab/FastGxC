#' eQTL Mapping Step
#'
#' Function to map cis eQTLs using either Matrix eQTL or TensorQTL.
#'
#' @param  SNP_file_name Full path to SNP genotype matrix.
#' @param  snps_location_file_name Full path to SNP location file.
#' @param  expression_file_name Full path to expression matrix.
#' @param  gene_location_file_name Full path to gene location file.
#' @param  context Context name for labeling output.
#' @param  out_dir Output directory.
#' @param  output_file_name_cis Path to write cis-eQTL output.
#' @param  output_file_name_tra Path to write trans-eQTL output.
#' @param  method Either "MatrixEQTL" or "tensorQTL".
#' @param  use_model MatrixEQTL model (default: modelLINEAR).
#' @param  cis_dist Distance threshold for cis-eQTLs.
#' @param  pv_threshold_cis P-value threshold for cis.
#' @param  pv_threshold_tra P-value threshold for trans.
#' @param  error_covariance Covariance matrix (or numeric()).
#'
#' @return Writes cis-eQTLs and trans-eQTLs (optional) to disk.
#' @export
eQTL_mapping_step = function(
    SNP_file_name,
    snps_location_file_name,
    expression_file_name,
    gene_location_file_name,
    context,
    out_dir,
    output_file_name_cis = output_file_name_cis,
    output_file_name_tra = output_file_name_tra,
    method = "MatrixEQTL",
    use_model = modelLINEAR,
    cis_dist = 1e6,
    pv_threshold_cis = 1,
    pv_threshold_tra = 0,
    error_covariance = numeric()
) {
  if (method == "MatrixEQTL") {
    setDTthreads(1)
    
    expression_mat <- as.matrix(data.frame(fread(expression_file_name, header = TRUE), row.names = 1, check.names = FALSE))
    genepos <- read.table(gene_location_file_name, header = TRUE)[, c("geneid", "chr", "s1", "s2")]
    
    snps <- SlicedData$new()
    snps$fileDelimiter <- "\t"
    snps$fileOmitCharacters <- "NA"
    snps$fileSkipRows <- 1
    snps$fileSkipColumns <- 1
    snps$fileSliceSize <- 2000
    snps$LoadFile(SNP_file_name)
    
    snpspos <- read.table(snps_location_file_name, header = TRUE)[, c("snpid", "chr", "pos")]
    snps$ColumnSubsample(which(snps$columnNames %in% colnames(expression_mat)))
    expression_mat <- expression_mat[, snps$columnNames]
    
    gene <- SlicedData$new()
    gene$CreateFromMatrix(expression_mat)
    
    Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = SlicedData$new(),
      output_file_name = output_file_name_tra,
      pvOutputThreshold = pv_threshold_tra,
      useModel = use_model,
      errorCovariance = error_covariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pv_threshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cis_dist,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    
  } else if (method == "tensorQTL") {
    SNP_file_name <- path.expand(SNP_file_name)
    snps_location_file_name <- path.expand(snps_location_file_name)
    gene_location_file_name <- path.expand(gene_location_file_name)
    expression_file_name <- path.expand(expression_file_name)
    out_dir <- path.expand(out_dir)
    
    py_run_string("
import os
import torch
import pandas as pd
import numpy as np
from tensorqtl import cis
import builtins

os.environ['CUDA_VISIBLE_DEVICES'] = ''
torch.set_num_threads(1)

def run_tensorqtl(phenotype_df, phenotype_pos_df, genotype_df, variant_df, prefix=''):
    if 's1' in phenotype_pos_df.columns and 's2' in phenotype_pos_df.columns:
        phenotype_pos_df.rename(columns={'s1': 'start', 's2': 'end'}, inplace=True)
    if 'chr' in variant_df.columns:
        variant_df.rename(columns={'chr': 'chrom'}, inplace=True)
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix=prefix)

builtins.run_tensorqtl = run_tensorqtl
")

    phenotype_df <- read.table(expression_file_name, header = TRUE, row.names = 1)
    phenotype_pos_df <- read.table(gene_location_file_name, header = TRUE, row.names = 1)
    genotype_df <- read.table(SNP_file_name, header = TRUE, row.names = 1)
    variant_df <- read.table(snps_location_file_name, header = TRUE, row.names = 1)
    
    py$phenotype_df <- phenotype_df
    py$phenotype_pos_df <- phenotype_pos_df
    py$genotype_df <- genotype_df
    py$variant_df <- variant_df
    py$prefix <- file.path(out_dir, paste0(context, "_tensor"))

    tryCatch({
      py$run_tensorqtl(
        py$phenotype_df,
        py$phenotype_pos_df,
        py$genotype_df,
        py$variant_df,
        prefix = py$prefix
      )
    }, error = function(e) {
      message("TensorQTL failed for context: ", context)
      message(e$message)
    })
    
    parquet_path <- output_file_name_cis
    if (!file.exists(parquet_path)) stop("TensorQTL failed: .parquet output not found")

    pd <- import("pandas", convert = FALSE)
    pq <- pd$read_parquet(parquet_path)
    py$pq <- pq

    pq$rename(columns = dict(
      phenotype_id = "gene",
      variant_id = "SNP",
      slope = "beta",
      pval_nominal = "p.value"
    ), inplace = TRUE)

    pq_df <- py_eval("pq.loc[:, ['SNP', 'gene', 'beta', 'p.value']]", convert = TRUE)
    pq_df$FDR <- p.adjust(pq_df$p.value, method = "BH")
    pq_df <- pq_df[order(pq_df$p.value), ]
    
    write.table(pq_df, file = output_file_name_cis, sep = '\t', row.names = FALSE, quote = FALSE)
  } else {
    stop("Invalid method. Use 'MatrixEQTL' or 'tensorQTL'")
  }
}