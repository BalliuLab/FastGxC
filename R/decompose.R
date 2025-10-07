#' Decomposition Step
#'
#' Function to decompose expression into one shared component and specific components per context
#'
#' @param  exp_mat_filename - full input filepath where expression matrix is stored. This file should be in the same format as the expression data file outputted by FastGxC's simulate_data function
#' @param  data_dir - full filepath where decomposed output files will be written out to
#' @return outputs one file with the shared component of expression per individual and C files for each specific expression component for each of the C contexts
#'
#' @export
decomposition_step = function(exp_mat_filename, data_dir){
  if(!dir.exists(data_dir)) dir.create(data_dir)
  
  exp_mat = data.table::fread(file = exp_mat_filename, sep = '\t', data.table = FALSE)
  
  design = sapply(1:nrow(exp_mat), function(i) unlist(strsplit(exp_mat[,1][i], split = " - "))[1])
  context_names = sapply(1:nrow(exp_mat), function(i) unlist(strsplit(exp_mat[,1][i], split = " - "))[2])
  contexts = unique(context_names)
  
  print("Decomposing data")
  rownames(exp_mat) = exp_mat[,1]
  exp_mat = exp_mat[, -1]
  
  cat(sprintf("There are %s samples and %s genes. The max number of missing samples for a gene is %s. The max number of missing genes for a sample is %s.\n",
              nrow(exp_mat), ncol(exp_mat), max(colSums(is.na(exp_mat))), max(rowSums(is.na(exp_mat)))))
  
  dec_exp_all = decompose(X = exp_mat, design = design)
  bexp_all = dec_exp_all$Xb
  wexp_all = dec_exp_all$Xw
  bexp_all[is.nan(bexp_all)] = NA
  wexp_all[is.nan(wexp_all)] = NA
  
  cat(sprintf("Between-individual matrix: %s individuals, %s genes. Max missing per gene: %s. Max missing per sample: %s.\n",
              nrow(bexp_all), ncol(bexp_all), max(colSums(is.na(bexp_all))), max(rowSums(is.na(bexp_all)))))
  cat(sprintf("Within-individual matrix: %s samples, %s genes. Max missing per gene: %s. Max missing per sample: %s.\n",
              nrow(wexp_all), ncol(wexp_all), max(colSums(is.na(wexp_all))), max(rowSums(is.na(wexp_all)))))
  
  # Function to clean expression files
  clean_expr_file <- function(f, n_individuals = 300) {
    lines <- readLines(f)
    header <- strsplit(lines[1], "\t")[[1]]
    if (header[1] == "") header <- header[-1]
    header <- header[header != "V301"]
    correct_header <- c("geneID", paste0("ind", 1:n_individuals))
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
  
  # Save shared expression matrix
  print("Saving between-individuals variation matrix")
  shared_file <- paste0(data_dir, "context_shared_expression.txt")
  fwrite(
    x = data.table::data.table(t(bexp_all), keep.rownames = TRUE) %>%
      { setnames(., old = "rn", new = "geneID")[] },
    file = shared_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
  )
  clean_expr_file(shared_file)
  
  # Save context-specific expression matrices
  print("Saving within-individuals variation matrix for context:")
  for (i in seq_along(contexts)) {
    print(contexts[i])
    wexp_t = wexp_all[grep(pattern = paste0(contexts[i], "$"), rownames(wexp_all)), ]
    rownames(wexp_t) = gsub(pattern = paste0(" - ", contexts[i]), replacement = "", x = rownames(wexp_t))
    context_file <- paste0(data_dir, contexts[i], "_specific_expression.txt")
    fwrite(
      x = data.table::data.table(t(wexp_t), keep.rownames = TRUE) %>%
        { setnames(., old = "rn", new = "geneID")[] },
      file = context_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
    )
    clean_expr_file(context_file)
  }
}