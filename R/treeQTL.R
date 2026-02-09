#' Multiple Testing Correction
#'
#' Function to adjust for hierarchical multiple testing correction using TreeQTL. Runs multiple testing correction for both FastGxC shared and specific eQTLs. 
#'
#' @param data_dir - full filepath of the directory where eQTL output files are stored. This function assumes that files are named as outputted by FastGxC's eQTL mapping function
#' @param snps_location_file_name - full filepath of the snpsloc file used in the eQTL mapping step
#' @param gene_location_file_name - full filepath of the geneloc file used in the eQTL mapping step
#' @param out_dir - full filepath of the output directory that FDR adjusted eQTLs should be written out to.
#' @param cisDist - numeric value of defined cis window. Note: this should match the cisDist set in the eQTL mapping step.
#' @param context_names - vector of all context names in the format c("tissue1", "tissue2", ..., etc.)
#' @param fdr_thresh - value between 0 and 1 that signifies what FDR threshold for multiple testing correction. The same value will be used across all hierarchical levels.
#' @param four_level - boolean value (T or F) that signifies whether to use the 4-level hierarchy (set this parameter to R and test for a global eQTL across shared and specific components) or a 3-level hierarchy (this parameter is default set to F to test for shared vs specific eQTLs)
#' @param qtl_type - string value "cis" or "trans" denoting the type of eQTL mapped. Default is set to "cis"
#' @param treeBH_method - character string specifying which TreeBH implementation to use when four_level = TRUE. Options: "original" (TreeBH package), "datatable" (optimized R), "cpp" (fast C++ implementation, default). Ignored when four_level = FALSE.
#' @param treeBH_test - character string specifying the p-value aggregation method for TreeBH when four_level = TRUE. Options: "simes" (Simes' method, default), "fisher" (Fisher's method). Ignored when four_level = FALSE.
#' @return outputs one file of specific eGenes across all contexts and one file of shared eGenes. Outputs an eAssociation file for each context and one for shared eQTLs with snp-gene pairs and FDR adjusted p-values.
#'
#' @export
treeQTL_step = function(data_dir, snps_location_file_name, gene_location_file_name, context_names, out_dir, cisDist = 1e6, fdr_thresh = 0.05, four_level = F, qtl_type = "cis", treeBH_method = "cpp", treeBH_test = "simes"){

# use a single thread
print(paste0("data.table getDTthreads(): ",getDTthreads()))
setDTthreads(1)
print(paste0("after setting as 1 thread; data.table getDTthreads(): ",getDTthreads()))

# Display all warnings as they occur
options(warn=1)

# FDR thresholds
level1=fdr_thresh
level2=fdr_thresh
level3=fdr_thresh

# Distance for local gene-SNP pairs
nearby = TRUE
if(qtl_type != "cis"){
  nearby = FALSE
  print(paste("Performing multiple testing adjustment for FastGxC trans-eQTLs."))
}else{
  print(paste("Performing multiple testing adjustment for FastGxC cis-eQTLs."))
}

snpspos = fread(file = snps_location_file_name);
genepos = fread(file = gene_location_file_name);

## get context names
if (is.vector(context_names)) {
    print(paste("Input for context names is a valid vector."))
}else{
    stop(print(paste0("No valid input for context names.")))
}

# If four_level, use TreeBH for hierarchical multiple testing and return early
if (four_level) {
  shared_file <- list.files(path = data_dir, pattern = "shared_shared.cis_pairs.txt", full.names = TRUE)
  if (length(shared_file) == 0) stop("No shared_shared.cis_pairs.txt file found.")

  to_TreeBH_input(data_dir = data_dir, shared_file = shared_file, context_names = context_names, out_dir = out_dir)
  df <- read.table(file.path(out_dir, "treeBH_input.txt"), header = TRUE)
  treeBH_step(matrix = df, fdr_thres = fdr_thresh, out_dir = out_dir, method = treeBH_method, test = treeBH_test)

  message("TreeBH multiple testing finished. Output written to ", out_dir)
  return(invisible(NULL))
}

# Use treeQTL to perform hierarchical FDR and get specific_eGenes, i.e. genes with at least one context-specific eQTL, and shared_eGenes, i.e. genes with at least one context-shared eQTL
  

names(genepos) <- tolower(names(genepos)) 
names(snpspos) <- tolower(names(snpspos))

genepos <- genepos |>
  dplyr::rename(
    geneid = dplyr::coalesce(names(genepos)[grepl("gene", names(genepos))][1], "geneid"),
    s1     = dplyr::coalesce(names(genepos)[grepl("start|s1", names(genepos))][1], "s1"),
    s2     = dplyr::coalesce(names(genepos)[grepl("end|s2", names(genepos))][1], "s2")
  ) 
snpspos <- snpspos |>
  dplyr::rename(
    snpid = dplyr::coalesce(names(snpspos)[grepl("SNP|snp", names(snpspos))][1], "snpid"),
    pos     = dplyr::coalesce(names(snpspos)[grepl("position|pos", names(snpspos))][1], "pos")
  ) 


shared_n_tests_per_gene = get_n_tests_per_gene(snp_map = snpspos %>% dplyr::select(snpid, chr, pos), gene_map = genepos %>% dplyr::select(geneid, chr, s1, s2), 
                                            nearby = nearby, dist = cisDist)
shared_n_tests_per_gene = data.frame(shared_n_tests_per_gene)

if (ncol(shared_n_tests_per_gene) != 2) {
  stop("unable to compute number of tests per gene.")
}

specific_eGenes=get_eGenes_multi_tissue_mod(
                              m_eqtl_out_dir = data_dir, 
                              treeQTL_dir = out_dir, 
                              tissue_names = context_names,
                              level1 = level1, level2 = level2, level3 = level3, 
                              exp_suffix = "specific",
                              four_level = four_level,
                              qtl_type = qtl_type,
                              shared_n_tests_per_gene = shared_n_tests_per_gene)
write.table(x = specific_eGenes, file = paste0(out_dir,"specific_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')

pattern=(paste0("shared.", qtl_type, "_pairs.txt"))
shared_eGenes = get_eGenes(n_tests_per_gene = shared_n_tests_per_gene, 
                            m_eqtl_out = list.files(data_dir, pattern = pattern,full.names = T), 
                            method = "BH",
                            level1 = level1, level2 = level2,
                            slice_size = 1e+05,
                            silent = FALSE)

write.table(x = shared_eGenes, file = paste0(out_dir,"shared_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')


eAssociations = get_eAssociations(eDiscoveries = shared_eGenes, n_tests = shared_n_tests_per_gene, 
                m_eqtl_out = list.files(data_dir, pattern = pattern,full.names = T),
                out_file = paste0(out_dir,"eAssoc_by_gene.context_shared.txt"), 
                by_snp = F, slice_size = 1e+05,
                silent = FALSE)


  
}
