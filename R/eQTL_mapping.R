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
eQTL_mapping_step = function(SNP_file_name, 
                             snps_location_file_name, 
                             expression_file_name, 
                             ene_location_file_name, 
                             context, 
                             shared_specific, 
                             out_dir,
                             output_file_name_cis = file.path(out_dir, paste0(context, "_", shared_specific, ".cis_pairs.txt")),
                             output_file_name_tra = file.path(out_dir, paste0(context, "_", shared_specific, ".trans_pairs.txt")),
                             use_model = modelLINEAR,
                             cis_dist = 1e6,
                             pv_threshold_cis = 1,
                             pv_threshold_tra = 0,
                             error_covariance = numeric()){

setDTthreads(1)

string1 = sprintf("Running analysis for %s \n", expression_file_name)
cat(string1)

#%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%% Read files
#%%%%%%%%%%%%%%%%%%%%%%%%

## Raw gene expression data with gene position
expression_mat=as.matrix(data.frame(data.table::fread(input = expression_file_name, header = T),row.names = 1, check.names = F))
genepos <- data.table::fread(gene_location_file_name, sep = NULL, header = TRUE, data.table = FALSE)
names(genepos) <- tolower(names(genepos)) 

genepos <- genepos |>
  dplyr::rename(
    geneid = dplyr::coalesce(names(genepos)[grepl("gene", names(genepos))][1], "geneid"),
    s1     = dplyr::coalesce(names(genepos)[grepl("start|s1", names(genepos))][1], "s1"),
    s2     = dplyr::coalesce(names(genepos)[grepl("end|s2", names(genepos))][1], "s2")
  ) |>
  dplyr::select(geneid, chr, s1, s2)

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

gene = SlicedData$new();
gene$CreateFromMatrix(expression_mat)

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
