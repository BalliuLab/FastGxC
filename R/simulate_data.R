#' Simulations
#'
#' Function to generate simulated SNP and expression data
#'
#' @param data_dir - full filepath of the directory where simulated data will be written to
#' @param N - number of individuals simulated (must be a whole number)
#' @param n_genes - number of genes simulated (must be a whole number)
#' @param n_snps_per_gene - number of cis-SNPs per gene (must be a whole number)
#' @param n_contexts - number of contexts to simulate (e.g. tissues or cell types) (must be a whole number)
#' @param maf - minor allele frequency for genotypes 
#' @param w_corr - error covariance between contexts
#' @param missing - decimal value signifying percentage of missingness in simulated expression data (e.g. parameter value of 0.3 would indicate 30% missing values in outputted expression matrix)
#' @param hsq - numeric vector with length of the number of contexts representing heritability explained by the eQTL in each context
#' @param mus - numeric vector with length of the number of contexts representing average expression in each context. 
#' @param cisDist - cis window used to create geneloc and snploc files. Note: this value must be larger than 1
#' @param seed - random seed for reproducibility. Default value is NULL
#' @return outputs an expression matrix file, a genotype matrix file, a SNP location file, and a gene location file all in the format needed for FastGxC's decomposition step and then subsequent eQTL mapping step with Matrix eQTL.
#'
#' @export
simulate_data = function(data_dir, N = 300, n_genes=100, n_snps_per_gene=1000, 
                         n_contexts=10, maf=0.2, w_corr=0.2, 
                         missing = 0.05, 
                         hsq = rep(0.2, n_contexts),
                         mus = rep(0, n_contexts),
                         cisDist = 1e6,
                         seed = NULL){

if(!dir.exists(data_dir)) dir.create(data_dir)
if(!is.null(seed)){
  set.seed(seed)
}

## parameter checks
stopifnot("N parameter must be greater than 0" = N > 0)
stopifnot("n_genes parameter must be greater than 0" = n_genes > 0)
stopifnot("n_snps_per_gene parameter must be greater than 0" = n_snps_per_gene > 0)
stopifnot("n_contexts parameter must be greater than 1" = n_contexts > 1)
stopifnot("length of hsq parameter must be equal to the number of contexts" = length(hsq)==n_contexts)
stopifnot("length of mus parameter must be equal to the number of contexts" = length(mus)==n_contexts)
stopifnot("values of hsq parameter must be [0,1)" = all(hsq >= 0 & hsq < 1))
stopifnot("w_corr parameter must be [-1,1]" = w_corr >= -1 & w_corr <= 1)


# Error variance-covariance matrix
v_e=1 # set error variance to 1 automatically
sigma = matrix(w_corr,nrow=n_contexts,ncol=n_contexts) # Error variance-covariance matrix
diag(sigma) = v_e
betas = sqrt((hsq * v_e)/((1 - hsq) * (2*maf*(1-maf)))) # set effect sizes

if(!is.positive.semi.definite(sigma)){
  stop("Variance covariance matrix is not positive semidefinite.")
}

# Simulate genotypes 
M <- n_genes * n_snps_per_gene
genos <- matrix(rbinom((N * M), 2,maf), nrow = N, ncol = M, byrow = FALSE)
storage.mode(genos) <- "integer"
dimnames(genos) <- list(paste0("ind", seq_len(N)), paste0("snp", seq_len(M)))
snp_out <- data.frame(snpid = paste0("snp", seq_len(M)), t(genos), check.names = FALSE)

if (requireNamespace("data.table", quietly = TRUE)) { 
  data.table::fwrite(snp_out, file = file.path(data_dir, "SNPs.txt"),
                     sep = "\t", row.names = F, col.names = T)
} else {
  fwrite(snp_out, file = file.path(data_dir, "SNPs.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


print("Finished simulating and saving genotypes")


# Save SNP location file
# Data frame with 3 initial columns (name, chrom, and position) that match standard SNP map file, followed by 1 column for each context with a 0/1 indicator of whether the given SNP passed QC in that tissue. 
snp_loc=data.frame(snpid=colnames(genos), chr="chr1",
                   pos=1, 
                   matrix(data = 1, nrow = ncol(genos), ncol = n_contexts, dimnames = list(NULL, paste0("context",1:n_contexts))))

# Save gene location file
# data frame with 4 initial columns (name, chrom, and start and end position) that match standard gene map file, followed by 1 column for each context with a 0/1 indicator of whether the given gene passed QC in that context 
gene_loc=data.frame(geneid=paste0("gene",1:n_genes), 
                    chr="chr1",
                    s1=1, 
                    s2=1,
                    matrix(data = 1, nrow = n_genes, ncol = n_contexts, dimnames = list(NULL, paste0("context",1:n_contexts))))

## populate snploc and geneloc dataframes
cur_chr = 1
max_position = 1e8
snp_pos = 1
for(gene in 1:n_genes){
  cur_position = (gene)*cisDist*10
  if(cur_position < max_position){
    gene_loc$chr[gene] = paste0("chr", cur_chr)
    gene_loc$s1[gene] = cur_position
    gene_loc$s2[gene] = cur_position + 1
    
    snp_loc$chr[snp_pos:(snp_pos + n_snps_per_gene -1)] = paste0("chr", cur_chr)
    snp_loc$pos[snp_pos:(snp_pos + n_snps_per_gene -1)] = cur_position
    snp_pos = snp_pos + n_snps_per_gene
  }
  else{
    cur_chr = cur_chr + 1
    cur_position = 1
    gene_loc$chr[gene] = paste0("chr", cur_chr)
    gene_loc$s1[gene] = cur_position
    gene_loc$s2[gene] = cur_position + 1
    
    snp_loc$chr[snp_pos:(snp_pos + n_snps_per_gene -1)] = paste0("chr", cur_chr)
    snp_loc$pos[snp_pos:(snp_pos + n_snps_per_gene -1)] = cur_position
    snp_pos = snp_pos + n_snps_per_gene
  }
}

fwrite(x = snp_loc, file = paste0(data_dir,"snpsloc.txt"), quote = F, sep = "\t", row.names = F,
            col.names = T)
fwrite(x = gene_loc, file = paste0(data_dir, "geneloc.txt"), quote = F, sep = "\t", row.names = F,
            col.names = T)


print("Finished saving snp and gene location files")

# Use only one snp per gene to generate expression
genos_with_effect = genos[,seq(from = 1, to = (n_snps_per_gene*n_genes), by = n_snps_per_gene)]

# Generate expression matrix
exp_mat=expand.grid(iid=paste0("ind",1:N),context=paste0("context",1:n_contexts))

for (i in 1:n_genes) {
  Y = matrix(0, nrow = N, ncol = n_contexts, dimnames = list(paste0("ind", 1:N), paste0("context", 1:n_contexts)))
  for (j in 1:n_contexts) Y[, j] = mus[j] + genos_with_effect[, i] * betas[j]
  for (j in 1:N) Y[j, ] = Y[j, ] + rmvnorm(1, rep(0, n_contexts), sigma)
  data_mat = melt(data = data.table(Y, keep.rownames = T),
                  id.vars = "rn")
  colnames(data_mat) = c("iid", "context", paste0("gene", i))
  exp_mat = merge(x = exp_mat, y = data_mat)
}


rownames(exp_mat) = paste(exp_mat$iid, exp_mat$context, sep = " - ")
exp_mat = cbind(data.frame(design = paste(exp_mat$iid, exp_mat$context,  sep = " - ")), exp_mat)
total_elements = prod(dim(exp_mat[, -c(1:3)]))
missing_elements = floor(missing * total_elements)
missing_indices = sample(1:total_elements, missing_elements)
identifiers = exp_mat[, 1:3]
exp_mat_missing = as.vector(as.matrix(exp_mat[, -c(1:3)]))
exp_mat_missing[missing_indices] = NA
exp_mat_missing = matrix(exp_mat_missing, nrow = nrow(exp_mat),  ncol = ncol(exp_mat) - 3)
exp_mat_missing = as.data.frame(exp_mat_missing)
exp_mat_missing = cbind(identifiers, exp_mat_missing)
colnames(exp_mat_missing) = colnames(exp_mat)
fwrite(x = exp_mat_missing[, -c(2:3)], file = paste0(data_dir, "expression.txt"), quote = F,
            sep = "\t", row.names = F, col.names = T)
print("Finished simulating and saving expression file")
}
