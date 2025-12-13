# Patched versions of TreeQTL functions to fix bugs with empty data handling
# These functions fix the "missing value where TRUE/FALSE needed" error
# that occurs when processing genes/components with no significant associations

#' Patched version of get_eGenes that uses patched get_selected_families
#' @noRd
get_eGenes <- function (n_tests_per_gene, m_eqtl_out, method = "BY", level1 = 0.05, 
    level2 = 0.05, slice_size = 1e+05, silent = FALSE, gene_pvals = NA) 
{
    get_selected_families(n_tests_per_gene, m_eqtl_out, method, 
        level1, level2, by_snp = FALSE, slice_size, silent, gene_pvals)
}

#' Helper function to get Matrix eQTL column structure
#' @noRd
get_m_eqtl_cols <- function (line1, by_snp, all = FALSE) 
{
    if ((all || by_snp) && line1[1] != "SNP") {
        stop("Matrix eQTL output does not have expected column names: \n\n         First column should be 'SNP'")
    }
    if ((all || !by_snp) && line1[2] != "gene" && line1[2] != 
        "trait") {
        stop("Matrix eQTL output does not have expected column names: \n\n         Second column should be 'gene' or 'trait'")
    }
    if (length(line1) == 5) {
        if (all) {
            if (line1[4] == "p-value" || line1[4] == "p.value") {
                cols_to_read <- list(SNP = character(), gene = character(), 
                  t.stat = 0, p.value = 0, FDR = 0)
            }
            else if (line1[5] == "p-value" || line1[5] == "p.value") {
                cols_to_read <- list(SNP = character(), gene = character(), 
                  beta = 0, t.stat = 0, p.value = 0)
            }
            else {
                stop("Matrix eQTL output does not have expected column names")
            }
        }
        else {
            if (by_snp) {
                if (line1[4] == "p-value" || line1[4] == "p.value") {
                  cols_to_read <- list(family = character(), 
                    NULL, NULL, p.value = 0, NULL)
                }
                else if (line1[5] == "p-value" || line1[5] == 
                  "p.value") {
                  cols_to_read <- list(family = character(), 
                    NULL, NULL, NULL, p.value = 0)
                }
                else {
                  stop("Matrix eQTL output does not have expected column names")
                }
            }
            else {
                if (line1[4] == "p-value" || line1[4] == "p.value") {
                  cols_to_read <- list(NULL, family = character(), 
                    NULL, p.value = 0, NULL)
                }
                else if (line1[5] == "p-value" || line1[5] == 
                  "p.value") {
                  cols_to_read <- list(NULL, family = character(), 
                    NULL, NULL, p.value = 0)
                }
                else {
                  stop("Matrix eQTL output does not have expected column names")
                }
            }
        }
    }
    else if (length(line1) == 6) {
        if (all) {
            cols_to_read <- list(SNP = character(), gene = character(), 
                beta = 0, t.stat = 0, p.value = 0, FDR = 0)
        }
        else {
            if (by_snp) {
                if (line1[5] == "p-value" || line1[5] == "p.value") {
                  cols_to_read <- list(family = character(), 
                    NULL, NULL, NULL, p.value = 0, NULL)
                }
                else {
                  stop("Matrix eQTL output does not have expected column names")
                }
            }
            else {
                if (line1[5] == "p-value" || line1[5] == "p.value") {
                  cols_to_read <- list(NULL, family = character(), 
                    NULL, NULL, p.value = 0, NULL)
                }
                else {
                  stop("Matrix eQTL output does not have expected column names")
                }
            }
        }
    }
    else {
        stop("Matrix eQTL output does not have expected number of columns")
    }
    cols_to_read
}

#' Patched version of get_selected_families with fixes for empty data handling
#' @noRd
get_selected_families <- function (fam_sizes, m_eqtl_out, method = "BY", level1 = 0.05, 
    level2 = 0.05, by_snp = TRUE, slice_size = 1e+05, 
    silent = FALSE, fam_pvals = NA) 
{
    names(fam_sizes) <- c("family", "n_tests")
    fam_sizes <- data.table(fam_sizes)
    setkey(fam_sizes, family)
    if (!is.na(fam_pvals)) {
        names(fam_pvals) <- c("family", "fam_p")
        family_table <- data.table(fam_pvals)
        setkey(fam_sizes, family)
    }
    else {
        family_table <- data.table(data.frame(family = fam_sizes$family), 
            j = 0, fam_p = 1)
        setkey(family_table, family)
        n_lines_done <- 0
        if (!silent) {
            print("Computing Simes p-values from Matrix eQTL output")
        }
        file_con <- file(m_eqtl_out, open = "r")
        line1 <- scan(file_con, what = "", nlines = 1, quiet = TRUE)
        cols_to_read <- get_m_eqtl_cols(line1, by_snp)
        more_file <- TRUE
        while (more_file) {
            if (!silent & n_lines_done > 0 & n_lines_done%%1e+06 == 
                0) {
                print(paste("Lines done ", n_lines_done))
            }
            cur_slice <- scan(file_con, what = cols_to_read, 
                nmax = slice_size, multi.line = FALSE, quiet = TRUE)
            n_lines_read <- length(cur_slice$family)
            if (n_lines_read > 0) {
              cur_slice <- data.frame(family = cur_slice$family, 
                  p.value = cur_slice$p.value, stringsAsFactors = FALSE)
              cur_slice <- data.table(cur_slice)
              setkey(cur_slice, family)
              n_lines_done <- n_lines_done + n_lines_read
              family_table[cur_slice, `:=`(c("j", "fam_p"), list(j + 
                  .N, min(fam_p, p.value/(j + 1)))), by = .EACHI]
            } else {
              n_lines_done <- n_lines_done + n_lines_read
            }
            if (n_lines_read < slice_size) {
                more_file <- FALSE
            }
        }
        close(file_con)
        family_table[fam_sizes, `:=`(c("fam_p"), list(fam_p * 
            max(n_tests, 1))), by = .EACHI]
    }
    n_fam <- nrow(family_table)
    if (method == "Bonf") {
        p_max <- level1/n_fam
        R <- sum(family_table$fam_p <= p_max)
    }
    else {
        if (method == "BY") {
            threshold_across_fams <- level1/sum(1/seq(1:n_fam))
        }
        else {
            threshold_across_fams <- level1
        }
        family_pvals <- family_table$fam_p
        family_pvals[which(family_pvals > 1)] <- 1
        family_pvals <- sort(family_pvals)
        ind_reject <- which((family_pvals/(1:n_fam)) <= threshold_across_fams/n_fam)
        if (length(ind_reject) > 0) {
            p_max <- max(family_pvals[ind_reject])
        }
        else {
            p_max <- 0
        }
        R <- ifelse(length(ind_reject) > 0, max(ind_reject), 
            0)
    }
    q_fam <- R * level2/n_fam
    ind_sel <- which(family_table$fam_p <= p_max)
    sel_families <- family_table[ind_sel, ]
    if (nrow(sel_families) == 0) {
        if (!silent) {
            print("No families were selected")
        }
        sel_families$n_sel <- numeric()
    }
    else {
        if (!silent) {
            print("Computing number of selections within each family")
        }
        sel_families$n_sel <- 0
        sel_families$j_fam <- 0
        sel_families <- data.table(sel_families)
        setkey(sel_families, family)
        sel_families[fam_sizes, `:=`(c("n_tests"), list(n_tests)), 
            by = .EACHI]
        n_lines_done <- 0
        file_con <- file(m_eqtl_out, open = "r")
        line1 <- scan(file_con, what = "", nlines = 1, quiet = TRUE)
        cols_to_read <- get_m_eqtl_cols(line1, by_snp)
        p_max_to_read <- q_fam
        cur_max_p <- 0
        more_file <- TRUE
        
        # FIX: Add defensive checks for while loop condition
        while (more_file) {
            # Check if we should continue based on p-value threshold
            if (!is.finite(cur_max_p) || !is.finite(p_max_to_read)) {
                break
            }
            if (cur_max_p >= p_max_to_read) {
                break
            }
            
            if (!silent & n_lines_done > 0 & n_lines_done%%1e+05 == 
                0) {
                print(paste("Lines done ", n_lines_done))
            }
            file_slice <- scan(file_con, what = cols_to_read, 
                nmax = slice_size, multi.line = FALSE, quiet = TRUE)
            n_lines_read <- length(file_slice$family)
            
            if (n_lines_read > 0) {
              file_slice <- data.frame(family = file_slice$family, 
                  p.value = file_slice$p.value, stringsAsFactors = FALSE)
              n_lines_done <- n_lines_done + n_lines_read
              sel_families[file_slice, `:=`(c("n_sel", "j_fam"), 
                  list(ifelse(((n_tests * p.value/(j_fam + 1)) <= 
                    q_fam) & (fam_p <= p_max), j_fam + 1, n_sel), 
                    j_fam + 1)), by = .EACHI]
              
              # FIX: Handle empty p.value vector properly
              valid_pvals <- file_slice$p.value[!is.na(file_slice$p.value) & is.finite(file_slice$p.value)]
              if (length(valid_pvals) > 0) {
                cur_max_p <- max(valid_pvals)
              } else {
                cur_max_p <- 0
              }
            } else {
              n_lines_done <- n_lines_done + n_lines_read
              cur_max_p <- 0
            }
            
            if (n_lines_read < slice_size) {
                more_file <- FALSE
            }
        }
        close(file_con)
        if (!silent && is.finite(cur_max_p) && is.finite(p_max_to_read) && 
            cur_max_p < p_max_to_read && more_file == FALSE) {
            warning("Matrix eQTL output threshold may be too small for given levels")
        }
    }
    as.data.frame(sel_families[, c("family", "fam_p", "n_sel"), 
        with = FALSE], stringsAsFactors = FALSE)
}

#' Patched version of get_nsel_SNPs_per_gene_tissue_pair with fixes for empty data handling
#' @noRd
get_nsel_SNPs_per_gene_tissue_pair <- function (sel_genes, tissue_name, m_eqtl_out, R_G, M, level3 = 0.05, 
    slice_size = 1e+05, silent = FALSE) 
{
    sel_genes <- data.table(sel_genes)
    setkey(sel_genes, family)
    sel_genes$q_adj <- R_G * sel_genes$n_sel_tissues * level3/M/sel_genes$n_tested_tissues
    if (!silent) {
        print("Computing number of selections for each gene x tissue selection")
    }
    sel_genes$n_sel_snp <- 0
    sel_genes$j_snp <- 0
    n_lines_done <- 0
    file_con <- file(m_eqtl_out, open = "r")
    line1 <- scan(file_con, what = "", nlines = 1, quiet = TRUE)
    cols_to_read <- get_m_eqtl_cols(line1, by_snp = FALSE)
    p_max_to_read <- max(sel_genes$q_adj)
    
    # FIX: Ensure p_max_to_read is finite
    if (!is.finite(p_max_to_read) || is.na(p_max_to_read)) {
      p_max_to_read <- 0
    }
    
    cur_max_p <- 0
    more_file <- TRUE
    
    # FIX: Add defensive checks for while loop condition
    while (more_file) {
        # Check if we should continue based on p-value threshold
        if (!is.finite(cur_max_p) || !is.finite(p_max_to_read)) {
            break
        }
        if (cur_max_p >= p_max_to_read) {
            break
        }
        
        if (!silent & n_lines_done > 0 & n_lines_done%%1e+05 == 
            0) {
            print(paste("Lines done ", n_lines_done))
        }
        file_slice <- scan(file_con, what = cols_to_read, nmax = slice_size, 
            multi.line = FALSE, quiet = TRUE)
        n_lines_read <- length(file_slice$family)
        
        if (n_lines_read > 0) {
          file_slice <- data.frame(family = file_slice$family, 
              p.value = file_slice$p.value, stringsAsFactors = FALSE)
          n_lines_done <- n_lines_done + n_lines_read
          sel_genes[file_slice, `:=`(c("n_sel_snp", "j_snp"), list(ifelse(((n_tests * 
              p.value/(j_snp + 1)) <= q_adj), j_snp + 1, n_sel_snp), 
              j_snp + 1)), by = .EACHI]
          
          # FIX: Handle empty p.value vector properly
          valid_pvals <- file_slice$p.value[!is.na(file_slice$p.value) & is.finite(file_slice$p.value)]
          if (length(valid_pvals) > 0) {
            cur_max_p <- max(valid_pvals)
          } else {
            cur_max_p <- 0
          }
        } else {
          n_lines_done <- n_lines_done + n_lines_read
          cur_max_p <- 0
        }
        
        if (n_lines_read < slice_size) {
            more_file <- FALSE
        }
    }
    close(file_con)
    
    # FIX: Only warn if values are finite
    if (is.finite(cur_max_p) && is.finite(p_max_to_read) && 
        cur_max_p < p_max_to_read && more_file == FALSE) {
        warning("Matrix eQTL output threshold may be too small for given levels")
    }
    
    as.data.frame(sel_genes[, c("family", "n_tests", "n_sel_snp"), 
        with = FALSE], stringsAsFactors = FALSE)
}

