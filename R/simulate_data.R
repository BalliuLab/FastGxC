#' Simulate Gene Expression and Genotype Data
#'
#' Generates simulated SNP genotypes, gene expression data across multiple contexts (e.g., tissues or cell types)
#'
#' @param data_dir String. Full file path to the directory where simulated data will be saved.
#' @param N Integer. Number of individuals to simulate.
#' @param n_genes Integer. Number of genes to simulate.
#' @param n_snps_per_gene Integer. Number of cis-SNPs per gene.
#' @param n_contexts Integer. Number of contexts to simulate (e.g., tissues or cell types).
#' @param maf Numeric. Minor allele frequency for simulated SNPs (value between 0 and 1).
#' @param w_corr Numeric. Off-diagonal values for the error covariance matrix (correlation between contexts).
#' @param v_e Numeric. Error variance for expression in each context.
#' @param missing Numeric. Proportion of missing values in the simulated expression data (e.g., 0.3 = 30% missing).
#' @param hsq Numeric vector of length `n_contexts`. Heritability values in each context (defaults to 0 = null model).
#' @param mus Numeric vector of length `n_contexts`. Mean expression in each context.
#'
#' @return Writes the following files to `data_dir`:
#' \itemize{
#'   \item \code{SNPs.txt} – genotype matrix (individuals × SNPs)
#'   \item \code{snpsloc.txt} – SNP location file for cis-window definitions
#'   \item \code{geneloc.txt} – gene location file for cis-window definitions
#'   \item \code{simulated_expression.txt} – simulated gene expression matrix with optional missingness
#' }
#'
#' @export
simulate_data = function(data_dir, N = 300, n_genes = 100, n_snps_per_gene = 1000,
                         n_contexts = 10, maf = 0.2, w_corr = 0.2,
                         v_e = 1, missing = 0,
                         hsq = rep(0,n_contexts),
                         mus = rep(0, n_contexts)) {

  clean_expr_file <- function(f) {
    lines <- readLines(f)
    header <- strsplit(lines[1], "\t")[[1]]
    if (header[1] == "") header <- header[-1]
    header <- header[header != "V301"]
    correct_header <- c("geneID", paste0("ind", 1:N))
    
    rows <- strsplit(lines[-1], "\t")
    maxlen <- max(sapply(rows, length))
    rows <- lapply(rows, function(x) { length(x) <- maxlen; x })
    mat <- do.call(rbind, rows)
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
    colnames(df) <- correct_header
    for (j in 2:ncol(df)) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
    rownames(df) <- df$geneID
    df <- df[, -1]
    write.table(df, file = f, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
  }

  if (!dir.exists(data_dir)) dir.create(data_dir)
  sigma = matrix(w_corr, nrow = n_contexts, ncol = n_contexts)
  diag(sigma) = v_e
  genos = sapply(1:(n_snps_per_gene * n_genes), function(X) rbinom(N, 2, maf))
  colnames(genos) = paste0("snp", 1:(n_snps_per_gene * n_genes))
  rownames(genos) = paste0("ind", 1:N)
  write.table(t(data.table(genos, keep.rownames = TRUE) %>%
                  rename(snpid = rn)),
              file = paste0(data_dir, "SNPs.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
  print("Finished simulating and saving genotypes")
  
  snp_loc = data.frame(snpid = colnames(genos), chr = "chr1",
                       pos = rep(seq(1, n_genes * 1e9, by = 1e9), each = n_snps_per_gene),
                       matrix(1, nrow = ncol(genos), ncol = n_contexts,
                              dimnames = list(NULL, paste0("context", 1:n_contexts))))
  write.table(snp_loc, file = paste0(data_dir, "snpsloc.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  print("Finished saving snp location file")
  
  gene_loc = data.frame(geneid = paste0("gene", 1:n_genes),
                        chr = "chr1", s1 = seq(1, n_genes * 1e9, by = 1e9),
                        s2 = seq(1, n_genes * 1e9, by = 1e9) + 1000,
                        matrix(1, nrow = n_genes, ncol = n_contexts,
                               dimnames = list(NULL, paste0("context", 1:n_contexts))))
  write.table(gene_loc, file = paste0(data_dir, "geneloc.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  print("Finished saving gene location file")
  
  genos_with_effect = genos[, seq(from = 1, to = (n_snps_per_gene * n_genes), by = n_snps_per_gene)]
  exp_mat = expand.grid(iid = paste0("ind", 1:N), context = paste0("context", 1:n_contexts))
  which_context = rep_len(1:n_contexts, length.out = n_genes)
  
  for (i in 1:n_genes) {
    betas = sqrt((hsq * v_e) / ((1 - hsq) * var(genos_with_effect[, i])))
    Y = matrix(0, nrow = N, ncol = n_contexts,
               dimnames = list(paste0("ind", 1:N), paste0("context", 1:n_contexts)))
    for (j in 1:n_contexts) {
      Y[, j] = mus[j] + genos_with_effect[, i] * betas[j]
    }
    for (j in 1:N) {
      Y[j, ] = Y[j, ] + rmvnorm(1, rep(0, n_contexts), sigma)
    }
    data_mat = melt(data.table(Y, keep.rownames = TRUE), id.vars = "rn")
    colnames(data_mat) = c("iid", "context", paste0("gene", i))
    exp_mat = merge(exp_mat, data_mat)
  }
  
  rownames(exp_mat) = paste(exp_mat$iid, exp_mat$context, sep = " - ")
  exp_mat = cbind(data.frame(design = paste(exp_mat$iid, exp_mat$context, sep = " - ")), exp_mat)
  total_elements = prod(dim(exp_mat[, -c(1:3)]))
  missing_elements = floor(missing * total_elements)
  missing_indices = sample(1:total_elements, missing_elements)
  identifiers = exp_mat[, 1:3]
  exp_mat_missing = as.vector(as.matrix(exp_mat[, -c(1:3)]))
  exp_mat_missing[missing_indices] = NA
  exp_mat_missing = matrix(exp_mat_missing, nrow = nrow(exp_mat), ncol = ncol(exp_mat) - 3)
  exp_mat_missing = as.data.frame(exp_mat_missing)
  exp_mat_missing = cbind(identifiers, exp_mat_missing)
  colnames(exp_mat_missing) = colnames(exp_mat)

  expr_file_path <- file.path(data_dir, "expression.txt")
  write.table(exp_mat_missing[, -c(2:3)], file = expr_file_path,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  print("Finished simulating and saving expression file")

  clean_expr_file(expr_file_path)
  print("Cleaned expression file for downstream compatibility")
}