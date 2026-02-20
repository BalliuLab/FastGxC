#' eQTL Mapping Step
#'
#' Function to map cis eQTLs - cis window is defined as 1Mb
#'
#' @param  SNP_file_name Full path to SNP genotype matrix.
#' @param  snps_location_file_name Full path to SNP location file.
#' @param  expression_file_name Full path to expression matrix.
#' @param  gene_location_file_name Full path to gene location file.
#' @param  context Context name for labeling output.
#' @param  shared_specific Context specific or context shared
#' @param  out_dir Output directory.
#' @param  output_file_name_cis Path to write cis-eQTL output.
#' @param  output_file_name_tra Path to write trans-eQTL output.
#' @param  method Either "MatrixEQTL" or "tensorqtl".
#' @param  use_model MatrixEQTL model (default: modelLINEAR).
#' @param  cis_dist Distance threshold for cis-eQTLs.
#' @param  pv_threshold_cis P-value threshold for cis.
#' @param  pv_threshold_tra P-value threshold for trans.
#' @param  error_covariance Covariance matrix (or numeric()).
#' 
#' @return Writes cis-eQTLs and trans-eQTLs (optional) to file.
#'
#' @export
eQTL_mapping_step = function(SNP_file_name, 
                             snps_location_file_name, 
                             expression_file_name, 
                             gene_location_file_name, 
                             context, 
                             shared_specific,
                             out_dir,
                             output_file_name_cis = file.path(out_dir, paste0(context, "_", shared_specific, ".cis_pairs.txt")),
                             output_file_name_tra = file.path(out_dir, paste0(context, "_", shared_specific, ".trans_pairs.txt")),
                             method = "MatrixEQTL",
                             use_model = modelLINEAR,
                             cis_dist = 1e6,
                             pv_threshold_cis = 1,
                             pv_threshold_tra = 0,
                             error_covariance = numeric()){

  setDTthreads(1)
  
  string1 = sprintf("Running analysis for %s \n", expression_file_name)
  cat(string1)
  
  if(method == "MatrixEQTL"){
    
    #%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%% Read files
    #%%%%%%%%%%%%%%%%%%%%%%%%
    
    ## Raw gene expression data with gene position
    expression_mat=as.matrix(data.frame(data.table::fread(input = expression_file_name, header = T),row.names = 1, check.names = F))
    genos = data.frame(fread(file = SNP_file_name, sep = '\t'),row.names = 1)
    
    genepos <- data.table::fread(gene_location_file_name, sep = "\t", header = TRUE, data.table = FALSE)
    names(genepos) <- tolower(names(genepos)) 
    
    genepos <- genepos |>
      dplyr::rename(
        geneid = dplyr::coalesce(names(genepos)[grepl("gene", names(genepos))][1], "geneid"),
        s1     = dplyr::coalesce(names(genepos)[grepl("start|s1", names(genepos))][1], "s1"),
        s2     = dplyr::coalesce(names(genepos)[grepl("end|s2", names(genepos))][1], "s2")
      ) |>
      dplyr::select(geneid, chr, s1, s2)
    
    snpspos <- read.table(snps_location_file_name, header = TRUE)[, c("snpid", "chr", "pos")]
    
    # Filter individuals with all NAs
    expression_mat = data.frame(expression_mat) %>% select_if(~ !all(is.na(.)))
    # Filter individuals from genotypes
    genos = genos[,colnames(expression_mat)]
    
    ## Load genotype data
    snps = SlicedData$new();
    snps$CreateFromMatrix(as.matrix(genos))
    
    snps$ColumnSubsample(which(snps$columnNames %in% colnames(expression_mat)))
    expression_mat <- expression_mat[, snps$columnNames]
    
    gene = SlicedData$new();
    gene$CreateFromMatrix(as.matrix(expression_mat))
    
    #%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%% Run the analysis
    #%%%%%%%%%%%%%%%%%%%%%%%%
    
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
  }
  else if(method == "tensorqtl"){
    SNP_file_name <- path.expand(SNP_file_name)
    snps_location_file_name <- path.expand(snps_location_file_name)
    gene_location_file_name <- path.expand(gene_location_file_name)
    expression_file_name <- path.expand(expression_file_name)
    out_dir <- path.expand(out_dir)
    suppressWarnings({
      py_run_string("
import os
import torch
import pandas as pd
import numpy as np
from tensorqtl import cis
import builtins

torch.set_num_threads(1)

def run_tensorqtl(phenotype_df, phenotype_pos_df, genotype_df, variant_df, prefix=''):
        
    if 's1' in phenotype_pos_df.columns and 's2' in phenotype_pos_df.columns:
        phenotype_pos_df.rename(columns={'s1': 'start', 's2': 'end'}, inplace=True)
    if 'chr' in variant_df.columns:
        variant_df.rename(columns={'chr': 'chrom'}, inplace=True)
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix=prefix)

builtins.run_tensorqtl = run_tensorqtl
        ")})
    
    expr <- read.table(expression_file_name, header = TRUE, sep = "\t", check.names = FALSE)
    
    rownames(expr) <- expr[, 1]
    expr <- expr[, -1, drop = FALSE]
    
    expr <- t(apply(expr, 1, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    }))
    
    write.table(
      expr,
      file = expression_file_name,
      sep = "\t",
      row.names = TRUE,
      col.names = NA,   
      quote = FALSE
    )
    phenotype_df <- suppressWarnings(read.table(expression_file_name, header = TRUE, row.names = 1))
    phenotype_pos_df <- suppressWarnings(read.table(gene_location_file_name, header = TRUE, row.names = 1))
    genotype_df <- suppressWarnings(read.table(SNP_file_name, header = TRUE, row.names = 1))
    variant_df <- suppressWarnings(read.table(snps_location_file_name, header = TRUE, row.names = 1))
    
    py$phenotype_df <- phenotype_df
    py$phenotype_pos_df <- phenotype_pos_df
    py$genotype_df <- genotype_df
    py$variant_df <- variant_df
    tensorqtl_prefix <- file.path(out_dir, paste0(context, "_", shared_specific))
    py$prefix <- tensorqtl_prefix
    
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
    
    
    # 1. Identify all matching parquet files in the folder
    all_parquet_files <- list.files(path = dirname(py$prefix), 
                                    pattern = paste0(basename(py$prefix), ".*\\.parquet$"), 
                                    full.names = TRUE)
    
    py$all_files <- all_parquet_files
    
    py_run_string("
import pandas as pd
# Read every file in the list and concatenate them immediately
pq_all = pd.concat([pd.read_parquet(f) for f in all_files]).reset_index(drop=True)

# Rename columns to your desired R format
pq_all = pq_all.rename(columns={
    'variant_id': 'SNP', 
    'phenotype_id': 'gene', 
    'slope': 'beta', 
    'pval_nominal': 'p-value'
})
")
    
    pq_df <- py$pq_all
    
    pq_df <- py_eval("pq.reset_index().loc[:, ['SNP', 'gene', 'beta', 'p-value']]", convert = TRUE)
    pq_df$FDR <- p.adjust(pq_df$`p-value`, method = "BH")
    pq_df <- pq_df[order(pq_df$`p-value`), ]
    
    write.table(pq_df, file = output_file_name_cis, sep = '\t', row.names = FALSE, quote = FALSE)
    
    unlink(all_parquet_files)
    } else {
      stop("Invalid method. Use 'MatrixEQTL' or 'tensorqtl'")
    }
}

