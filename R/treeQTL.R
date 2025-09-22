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
treeQTL_step = function(data_dir, snps_location_file_name, gene_location_file_name, context_names, out_dir, fdr_thresh = 0.05, four_level = FALSE){
  
  print(paste0("data.table getDTthreads(): ", getDTthreads()))
  setDTthreads(1)
  print(paste0("after setting as 1 thread; data.table getDTthreads(): ", getDTthreads()))
  
  options(warn=1)
  
  level1 = fdr_thresh
  level2 = fdr_thresh
  level3 = fdr_thresh
  
  cisDist = 1e6;
  
  snpspos = fread(file = snps_location_file_name);
  genepos = fread(file = gene_location_file_name);
  
  if (is.vector(context_names)) {
    print(paste("Input for context names is a valid vector."))
  } else {
    stop(print(paste0("No valid input for context names.")))
  }
  
  shared_n_tests_per_gene = get_n_tests_per_gene(snp_map = snpspos[,1:3], gene_map = genepos[,1:4], 
                                                 nearby = TRUE, dist = cisDist)
  shared_n_tests_per_gene = data.frame(shared_n_tests_per_gene)
  shared_n_tests_per_gene = data.frame(gene = rownames(shared_n_tests_per_gene), shared_n_tests_per_gene)
  names(shared_n_tests_per_gene) = c("family", "n_tests")
  
  get_eGenes_multi_tissue_mod <- function(m_eqtl_out_dir, treeQTL_dir, tissue_names, level1 = 0.05, level2 = 0.05, level3 = 0.05, exp_suffix, four_level = FALSE, shared_n_tests_per_gene) {
    pattern <- paste0("context", 1:length(tissue_names), "_", exp_suffix, ".cis_pairs.txt")
    m_eqtl_outfiles <- file.path(m_eqtl_out_dir, pattern)
    print("Constructed MatrixEQTL output files:")
    print(m_eqtl_outfiles)
    
    if (!all(file.exists(m_eqtl_outfiles))) {
      stop("Some expected MatrixEQTL files do not exist.")
    }
    
    n_tissue <- length(tissue_names)
    for (i in 1:n_tissue) {
      cur_tissue_name <- tissue_names[i]
      m_eqtl_out <- m_eqtl_outfiles[i]
      print(paste("Computing summary statistics for tissue", cur_tissue_name))
      n_SNPs_per_gene_this_tissue <- fread(m_eqtl_out) %>% select(gene) %>% group_by(gene) %>% mutate(n = n()) %>% distinct()
      colnames(n_SNPs_per_gene_this_tissue) <- c("family", "n_tests")
      gene_simes_cur_tissue <- get_eGenes(n_tests_per_gene = n_SNPs_per_gene_this_tissue, m_eqtl_out = m_eqtl_out, method = "BH", level1 = 1, level2 = 1, silent = TRUE)
      gene_simes_cur_tissue <- merge(gene_simes_cur_tissue, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
      gene_simes_cur_tissue$fam_p[which(is.na(gene_simes_cur_tissue$fam_p))] <- 1
      gene_simes_cur_tissue$eGene <- as.integer(gene_simes_cur_tissue$fam_p <= level1)
      
      # Output eAssociations for specific context
      context_eGenes <- gene_simes_cur_tissue[gene_simes_cur_tissue$eGene == 1, , drop = FALSE]
      if (nrow(context_eGenes) > 0) {
        colnames(context_eGenes)[1] <- "family"
        context_eGenes$pval <- 0
        context_eGenes$n_sel <- 1
        
        eAssoc_file <- paste0(treeQTL_dir, "/eAssoc_by_gene.context_", cur_tissue_name, ".txt")
        get_eAssociations(
          eDiscoveries = context_eGenes,
          n_tests = n_SNPs_per_gene_this_tissue,
          m_eqtl_out = m_eqtl_out,
          out_file = eAssoc_file,
          by_snp = FALSE,
          slice_size = 1e+05,
          silent = FALSE
        )
      }
      
      if (i == 1) {
        eGene_pvals <- gene_simes_cur_tissue[, c("family", "fam_p")]
        n_SNPs_per_gene_xT <- n_SNPs_per_gene_this_tissue
      } else {
        eGene_pvals <- merge(eGene_pvals, gene_simes_cur_tissue[, c("family", "fam_p")], by = "family", all = TRUE)
        n_SNPs_per_gene_xT <- merge(n_SNPs_per_gene_xT, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
      }
      names(eGene_pvals)[i + 1] <- cur_tissue_name
      names(n_SNPs_per_gene_xT)[i + 1] <- cur_tissue_name
    }
    names(eGene_pvals)[1] <- "gene"
    col_ind_pvals <- 2:(n_tissue + 1)
    eGene_pvals$simes_p <- apply(eGene_pvals[, col_ind_pvals], 1, TreeQTL:::get_simes_p)
    eGene_xT_qvals <- qvalue(eGene_pvals$simes_p, lambda = 0)$qvalue
    R_G <- sum(eGene_xT_qvals <= level1)
    print(paste("Number of cross-tissue eGenes = ", R_G))
    
    if (R_G == 0) {
      warning("No eGenes selected at level1; skipping lower-level selections.")
      return(list(
        qvalue_table = data.frame(),
        raw_pval_table = eGene_pvals,
        qvals = eGene_xT_qvals
      ))
    }
    
    q2_adj <- R_G * level2 / nrow(eGene_pvals)
    ind_sel_simes <- which(eGene_xT_qvals <= level1)
    sel_eGenes_simes <- eGene_pvals[ind_sel_simes, ]
    rej_simes <- t(1 * apply(sel_eGenes_simes[, col_ind_pvals], 1, TreeQTL:::qsel_by_fam, q2_adj))
    colnames(rej_simes) <- tissue_names
    sel_eGenes_simes$n_sel_tissues <- rowSums(rej_simes)
    sel_eGenes_simes$n_tested_tissues <- rowSums(!is.na(sel_eGenes_simes[, col_ind_pvals]))
    eGene_xT_sel <- data.frame(gene = sel_eGenes_simes$gene)
    eGene_xT_sel <- cbind(eGene_xT_sel, rej_simes)
    names(eGene_xT_sel)[2:(n_tissue + 1)] <- tissue_names
    
    return(list(
      qvalue_table = eGene_xT_sel,
      raw_pval_table = eGene_pvals,
      qvals = eGene_xT_qvals
    ))
  }
  
  if (four_level) {
    shared_file <- list.files(path = data_dir, pattern = "shared_shared.cis_pairs.txt", full.names = TRUE)
    if (length(shared_file) == 0) stop("No shared_shared.cis_pairs.txt file found.")
    
    to_TreeBH_input(data_dir = data_dir, shared_file = shared_file, context_names = context_names, out_dir = out_dir)
    df <- read.table(file.path(out_dir, "treeBH_input.txt"), header = TRUE)
    treeBH_step(matrix = df, fdr_thres = fdr_thresh, out_dir = out_dir)
    
    message("TreeBH multiple testing finished. Output written to ", out_dir)
    return(invisible("TreeBH done"))
  } else {
    res = get_eGenes_multi_tissue_mod(
      m_eqtl_out_dir = data_dir, 
      treeQTL_dir = out_dir, 
      tissue_names = context_names,
      level1 = level1, level2 = level2, level3 = level3, 
      exp_suffix = "specific",
      four_level = four_level,
      shared_n_tests_per_gene = shared_n_tests_per_gene
    )
    
    specific_eGenes = res$qvalue_table
    
    if (nrow(specific_eGenes) == 0) {
      message("No specific eGenes discovered. Skipping downstream steps.")
    }
    
    write.table(x = specific_eGenes, file = paste0(out_dir, "specific_eGenes.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    
    pattern = ("shared.all_pairs.txt")
    shared_file <- list.files(path = data_dir, pattern = "shared_shared.cis_pairs.txt", full.names = TRUE)
    if (length(shared_file) != 1) stop("Expected one shared file but found: ", length(shared_file))
    
    shared_eGenes = get_eGenes(
      n_tests_per_gene = shared_n_tests_per_gene, 
      m_eqtl_out = shared_file, 
      method = "BH",
      level1 = level1, level2 = level2,
      slice_size = 1e+05,
      silent = FALSE
    )
    
    if (nrow(shared_eGenes) == 0) {
      message("No shared eGenes discovered. Skipping shared association analysis.")
      return(NULL)
    }
    
    write.table(x = shared_eGenes, file = paste0(out_dir, "shared_eGenes.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    
    eAssociations = get_eAssociations(
      eDiscoveries = shared_eGenes,
      n_tests = shared_n_tests_per_gene,
      m_eqtl_out = shared_file,
      out_file = paste0(out_dir, "eAssoc_by_gene.context_shared.txt"), 
      by_snp = FALSE, slice_size = 1e+05,
      silent = FALSE
    )
  }
}