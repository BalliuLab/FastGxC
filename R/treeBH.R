#' Converting eQTL Mapping Output to TreeBH Input
#' 
#' @param data_dir - full filepath of the directory where eQTL output files are stored. This function assumes that files are named as outputted by FastGxC's eQTL mapping function
#' @param shared_file - full filepath to the shared all_pairs file from the output of the mapping function
#' @param context_names - vector of all context names in the format c("tissue1", "tissue2", ..., etc.)
#' @param out_dir - full filepath of the output directory where the input of TreeBH can be stored
#' @return A data.frame with columns:
#'   - gene_name     : gene name
#'   - snp_name      : snp name
#'   - context_name  : context name (or "shared")
#'   - gene_snp      : concatenated "gene_SNP"
#'   - component_tag : "gene_SNP_specific" or "gene_SNP_shared"
#'   - context_tag   : "gene_SNP_component_context"
#'   - p_value       : numeric p-value
#' @export
to_TreeBH_input <- function(data_dir, shared_file, context_names, out_dir) {
  # Prepare list to hold each context-specific + shared df
  num_contexts <- length(context_names)
  df_list      <- vector("list", num_contexts + 1)

  # Inline file-reading and p-value/column validation per context
  for (i in seq_len(num_contexts)) {
    context <- context_names[i]
    pattern <- paste0(context, "_specific.all_pairs.txt$")
    files   <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) stop("No file matching pattern ", pattern)
    df      <- readr::read_tsv(files[[1]], col_types = readr::cols())

    # Ensure columns exist and types are correct
    if (!"p-value" %in% names(df)) stop("Missing 'p-value' in ", files[[1]])
    p <- df[["p-value"]]
    if (!is.numeric(p) || any(is.na(p))) stop("`p-value` must be numeric without NAs")
    if (any(p < 0 | p > 1)) stop("`p-value` entries must lie between 0 and 1")

    if (!"gene" %in% names(df)) stop("Missing 'gene' in ", files[[1]])
    if (!is.character(df[["gene"]])) stop("`gene` must be character in ", files[[1]])

    if (!"SNP" %in% names(df)) stop("Missing 'SNP' in ", files[[1]])
    if (!is.character(df[["SNP"]])) stop("`SNP` must be character in ", files[[1]])

    df_list[[i]] <- df
  }

  # Shared file: same inline checks
  df_shared <- readr::read_tsv(shared_file, col_types = readr::cols())
  for (col in c("p-value","gene","SNP")) {
    if (!col %in% names(df_shared)) stop("Missing '", col, "' in shared_file")
    col_data <- df_shared[[col]]
    if (col == "p-value") {
      if (!is.numeric(col_data) || any(is.na(col_data))) stop("`p-value` invalid in shared_file")
      if (any(col_data < 0 | col_data > 1)) stop("`p-value` entries must lie between 0 and 1 in shared_file")
    } else {
      if (!is.character(col_data)) stop("`", col, "` must be character in shared_file")
    }
  }
  df_list[[num_contexts + 1]] <- df_shared

  # Build grouping labels
  rows_per_df     <- vapply(df_list, nrow, integer(1))
  context_index   <- rep(seq_along(df_list), times = rows_per_df)
  context_label   <- c(context_names, "shared")[context_index]
  component_label <- ifelse(context_label == "shared", "shared", "specific")

  # Extract columns across all data.frames
  p_values   <- unlist(lapply(df_list, `[[`, "p-value"),   use.names = FALSE)
  gene_names <- unlist(lapply(df_list, `[[`, "gene"),      use.names = FALSE)
  snp_names  <- unlist(lapply(df_list, `[[`, "SNP"),       use.names = FALSE)

  # Build hierarchical identifiers
  gene_snp_id  <- paste0(gene_names, "_", snp_names)
  component_id <- paste0(gene_snp_id, "_", component_label)
  context_id   <- paste0(component_id, "_", context_label)

  # Assemble final data.frame
  out_df <- data.frame(
    gene_name          = gene_names,
    snp_name           = snp_names, 
    context_name       = context_label, 
    gene_snp           = gene_snp_id,
    component_tag      = component_id,
    context_tag        = context_id,
    p_value            = p_values,
    stringsAsFactors   = FALSE,
    check.names        = FALSE
  )

  # Write file
  write.table(
    out_df,
    file      = file.path(out_dir, "treeBH_input.txt"),
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
}

#' Multiple Testing Correction with TreeBH
#' Function to adjust for hierarchical multiple testing correction using TreeBH. 
#' @param matrix - output of to_TreeBH_Input, matrix with the groups and pvalues
#' @param fdr_thres -  value between 0 and 1 that signifies what FDR threshold for multiple testing correction. The same value will be used across all hierarchical levels.
#' @param out_dir - full filepath of the output directory where the output of TreeBH can be stored
#' @return A matrix with identification columns and output columns:
#' - gene                    : gene name
#' - snp                     : snp name
#' - context                 : context name or "specific"
#' - treeBH Output Columns   : eGene, eQTL, Component, Context Specific
#' @export
treeBH_step <- function(matrix, fdr_thres, out_dir) {
  
  id_mat <- matrix[, 1:3]
  groups <- matrix[, c(1, 4, 5, 6)]
  p_values <- matrix[, 7] |> as.vector() |> as.numeric()
  
  num_levels <- 4
  
  treeBH_results <- TreeBH::get_TreeBH_selections(pvals = p_values, 
                                                  groups = groups, 
                                                  q = rep(fdr_thres, num_levels)) 
  
  colnames(treeBH_results) <- c("eGene", "eQTL", "Component", "ContextSpecific")
  
  results <- cbind(id_mat, treeBH_results)
  
  write.table(
    results, 
    file       = file.path(out_dir, "treeBH_output.txt"),
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
}