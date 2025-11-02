#' eQTL Mapping Step
#'
#' Function to map cis eQTLs - cis window is defined as 1Mb
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
#' @return Writes cis-eQTLs and trans-eQTLs (optional) to file.
#'
#' @export
eQTL_mapping_step = function(SNP_file_name, 
                             snps_location_file_name, 
                             expression_file_name, 
                             gene_location_file_name, 
                             context, 
                             out_dir,
                             output_file_name_cis = file.path(out_dir, paste0(context, "_", shared_specific, ".cis_pairs.txt")),
                             output_file_name_tra = file.path(out_dir, paste0(context, "_", shared_specific, ".trans_pairs.txt")),
                             use_model = modelLINEAR,
                             cis_dist = 1e6,
                             pv_threshold_cis = 1,
                             pv_threshold_tra = 0,
                             error_covariance = numeric()){

  setDTthreads(1)
  
  expression_file_name = file.path(input_dir, paste0(context, "_specific_expression.txt"))
  string1 = sprintf("Running analysis for %s \n", expression_file_name)
  cat(string1)
  
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
