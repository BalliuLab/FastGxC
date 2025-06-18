#' eQTL Mapping Step
#'
#' Function to map cis eQTLs using either Matrix eQTL or TensorQTL. The cis window is defined as 1Mb.
#'
#' @param  SNP_file_name Full file path with genotypes of all individuals.
#' @param  snps_location_file_name Full file path with SNP IDs, start, and end positions.
#' @param  expression_file_name Full file path with expression matrix of individuals per gene.
#' @param  gene_location_file_name Full file path with gene IDs, start, and end positions.
#' @param  context Name of context for output file naming purposes.
#' @param  shared_specific Either "shared" if mapping eQTLs with a shared component, or "specific" if mapping eQTLs with a specific component.
#' @param  out_dir Full path of the output directory where eQTL result files will be written.
#' @param  method Mapping method to use: either "matrix" (Matrix eQTL) or "tensor" (TensorQTL). Default is "matrix".
#'
#' @return Writes one file of mapped cis eQTLs to the output directory.
#' @export
eQTL_mapping_step = function(SNP_file_name, snps_location_file_name, expression_file_name, gene_location_file_name, context, shared_specific, out_dir, method = "matrix"){
  if(method == "matrix"){
    # Output file name
    output_file_name_cis = paste0(out_dir, context, "_", shared_specific, ".all_pairs.txt"); 
    output_file_name_tra = tempfile();
    
    setDTthreads(1)
    print(paste0("data.table getDTthreads(): ",getDTthreads()))
    
    string1 = sprintf("Running analysis for %s \n", expression_file_name)
    cat(string1)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%% MatrixEQTL parameters
    #%%%%%%%%%%%%%%%%%%%%%%%%
    
    # Linear model to use
    useModel = modelLINEAR; 
    
    # Only associations significant at this level will be saved
    pvOutputThreshold_cis = 1; 
    pvOutputThreshold_tra = 0;
    
    # Error covariance matrix
    errorCovariance = numeric();
    
    # Distance for local gene-SNP pairs
    cisDist = 1e6;
    
    #%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%% Read files
    #%%%%%%%%%%%%%%%%%%%%%%%%
    
    ## Raw gene expression data with gene position
    expression_mat=as.matrix(data.frame(data.table::fread(input = expression_file_name, header = T),row.names = 1, check.names = F))
    genepos = read.table(file = gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)[,1:4];
    
    ## Genotype data with snp position
    snps = MatrixEQTL::SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name);
    
    snpspos = read.table(file = snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)[,1:3];
    
    print(dim(snps))
    print(snps$columnNames)
    
    ## Make sure order of individuals is the same in gene expression and genotype matrices  
    snps$ColumnSubsample(which(snps$columnNames %in% colnames(expression_mat))) # Match SNP and expression individuals
    expression_mat=expression_mat[,snps$columnNames]
    gene = SlicedData$new();
    gene$CreateFromMatrix(expression_mat)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%% Run the analysis
    #%%%%%%%%%%%%%%%%%%%%%%%%
    
    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = SlicedData$new(),
      output_file_name  = output_file_name_tra,
      pvOutputThreshold  = pvOutputThreshold_tra,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    
    ## Results:
    cat('Analysis finished in: ', me$time.in.sec, ' seconds', '\n')
  }
  else if (method == "tensor") {
    use_python(Sys.which("python"), required = TRUE)
    
    py_run_string("
import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post

def run_tensorqtl(phenotype_df, phenotype_pos_df, genotype_df, variant_df, prefix=''):
    phenotype_pos_df.rename(columns={'s1': 'start', 's2': 'end'}, inplace=True)
    variant_df.rename(columns={'chr': 'chrom'}, inplace=True)
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix=prefix)
    pairs_df = pd.read_parquet(f'{prefix}.cis_qtl_pairs.chr1.parquet')
    pairs_df.rename(columns={'phenotype_id': 'gene', 'variant_id': 'SNP'}, inplace=True)
    return pairs_df
")
    
    phenotype_df = read.table(expression_file_name, header = TRUE, row.names = 1)
    phenotype_pos_df = read.table(gene_location_file_name, header = TRUE, row.names = 1)
    genotype_df = read.table(SNP_file_name, header = TRUE, row.names = 1)
    variant_df = read.table(snps_location_file_name, header = TRUE, row.names = 1)
    
    py$phenotype_df <- phenotype_df
    py$phenotype_pos_df <- phenotype_pos_df
    py$genotype_df <- genotype_df
    py$variant_df <- variant_df
    py$prefix <- paste0(out_dir, '/', context, '_tensor')
    
    py_run_string("result = run_tensorqtl(phenotype_df, phenotype_pos_df, genotype_df, variant_df, prefix)")
    tensorqtl_result <- py$result
    write.table(tensorqtl_result, file = paste0(out_dir, '/', context, '_tensorqtl_results.txt'),
                sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    stop("Invalid method. Use 'matrix' or 'tensor'")
  }
}