#' Multiple Testing Correction
#'
#' Function to adjust for hierarchical multiple testing correction using TreeQTL. Runs multiple testing correction for both FastGxC shared and specific eQTLs. 
#'
#' @param data_dir - full filepath of the directory where eQTL output files are stored. This function assumes that files are named as outputted by FastGxC's eQTL mapping function
#' @param snps_location_file_name - full filepath of the snpsloc file used in the eQTL mapping step
#' @param gene_location_file_name - full filepath of the geneloc file used in the eQTL mapping step
#' @param out_dir - full filepath of the output directory that FDR adjusted eQTLs should be written out to.
#' @param context_names - vector of all context names in the format c("tissue1", "tissue2", ..., etc.)
#' @param fdr_thresh - value between 0 and 1 that signifies what FDR threshold for multiple testing correction. The same value will be used across all hierarchical levels.
#' @param four_level - boolean value (T or F) that signifies whether to use the 4-level hierarchy (set this parameter to R and test for a global eQTL across shared and specific components) or a 3-level hierarchy (this parameter is default set to F to test for shared vs specific eQTLs)
#' @return outputs one file of specific eGenes across all contexts and one file of shared eGenes. Outputs an eAssociation file for each context and one for shared eQTLs with snp-gene pairs and FDR adjusted p-values. 
#'
#' @export
treeQTL_step = function(data_dir, snps_location_file_name, gene_location_file_name, out_dir, fdr_thresh = 0.05, four_level = F){

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
cisDist = 1e6;

snpspos = fread(file = snps_location_file_name);
genepos = fread(file = gene_location_file_name);

## get context names
if (is.vector(context_names)) {
    print(paste("Input for context names is a valid vector."))
}else{
    stop(print(paste0("No valid input for context names.")))
}
# Use treeQTL to perform hierarchical FDR and get specific_eGenes, i.e. genes with at least one context-specific eQTL, and shared_eGenes, i.e. genes with at least one context-shared eQTL
  
shared_n_tests_per_gene = get_n_tests_per_gene(snp_map = snpspos[,1:3], gene_map = genepos[,1:4], 
                                            nearby = TRUE, dist = cisDist)

if(four_level){
    specific_eGenes=get_eGenes_multi_tissue_mod(
                                  m_eqtl_out_dir = data_dir, 
                                  treeQTL_dir = out_dir, 
                                  tissue_names = context_names,
                                  level1 = level1, level2 = level2, level3 = level3, 
                                  exp_suffix = "specific",
                                  four_level = four_level,
                                  shared_n_tests_per_gene = shared_n_tests_per_gene)
    write.table(x = specific_eGenes, file = paste0(work_dir,"specific_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
}else{
    specific_eGenes=get_eGenes_multi_tissue_mod(
                                  m_eqtl_out_dir = data_dir, 
                                  treeQTL_dir = out_dir, 
                                  tissue_names = context_names,
                                  level1 = level1, level2 = level2, level3 = level3, 
                                  exp_suffix = "specific",
                                  four_level = four_level)
    write.table(x = specific_eGenes, file = paste0(work_dir,"specific_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')

    pattern=("shared.all_pairs.txt")
    shared_eGenes = get_eGenes(n_tests_per_gene = n_tests_per_gene, 
                                m_eqtl_out = list.files(data_dir, pattern = pattern,full.names = T), 
                                method = "BH",
                                level1 = level1, level2 = level2,
                                slice_size = 1e+05,
                                silent = FALSE)

    write.table(x = shared_eGenes, file = paste0(work_dir,"shared_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')


    eAssociations = get_eAssociations(eDiscoveries = shared_eGenes, n_tests = n_tests_per_gene, 
                    m_eqtl_out = list.files(data_dir, pattern = pattern,full.names = T),
                    out_file = paste0(work_dir,"eAssoc_by_gene.context_shared.txt"), 
                    by_snp = F, slice_size = 1e+05,
                    silent = FALSE)
}

  
}
