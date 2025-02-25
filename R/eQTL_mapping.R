#' eQTL Mapping Step
#'
#' Function to map cis eQTLs - cis window is defined as 1Mb
#'
#' @param  SNP_file_name - full file path with genotypes of all individuals
#' @param  snps_location_file_name - full file path with snp ids, start, and end positions
#' @param  expression_file_name - full file path with expression matrix of individuals per gene
#' @param  gene_location_file_name - full file path with gene ids, start, and end positions
#' @param  context - name of context for output file naming purposes
#' @param  shared_specific - either "shared" if mapping eQTLs with shared component or "specific" if mapping eQTLs with specific component
#' @param  out_dir - full file path of output directory where eQTL file will be written out
#' @return outputs one file of mapped cis eQTLs
#'
#' @export
eQTL_mapping_step = function(SNP_file_name, snps_location_file_name, expression_file_name, gene_location_file_name, context, shared_specific, out_dir){

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
