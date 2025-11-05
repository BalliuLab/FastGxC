#' Decomposition Step
#'
#' Function to decompose expression into one shared component and specific components per context
#'
#' @param  exp_mat_filename - full input filepath where expression matrix is stored. This file should be in the same format as the expression data file outputted by FastGxC's simulate_data function
#' @param  data_dir - full filepath where decomposed output files will be written out to
#' @return outputs one file with the shared component of expression per individual and C files for each specific expression component for each of the C contexts
#'
#' @export
decomposition_step = function(exp_mat_filename, data_dir, out_dir = data_dir){
if(!dir.exists(data_dir)) dir.create(data_dir)
#%%%%%%%%%%%%%%% Read expression matrix, genes in columns, samples in rows.
#exp_mat=read.table(file = paste0(data_dir,exp_mat_filename), sep = '\t')
exp_mat=data.table::fread(file = exp_mat_filename, sep = '\t', data.table = F)

#%%%%%%%%%%%%%%% Sample and context names
design=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(exp_mat[,1][i], split = " - "))[1])
context_names=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(exp_mat[,1][i], split = " - "))[2])
contexts=unique(context_names)

#%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
print("Decomposing data")
rownames(exp_mat) = exp_mat[,1]
exp_mat = exp_mat[,-1]

#%%%%%%%%%%%%%%% Print number of genes and samples
string1 = sprintf("There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s. \n", nrow(exp_mat), ncol(exp_mat),max(colSums(is.na(exp_mat))),max(rowSums(is.na(exp_mat))))
cat(string1)

dec_exp_all=decompose(X = exp_mat, design = design)
bexp_all=dec_exp_all$Xb
wexp_all=dec_exp_all$Xw
bexp_all[is.nan(bexp_all)]=NA
wexp_all[is.nan(wexp_all)]=NA

string2 = sprintf("Between individual matrix: There are %s individuals and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s. \n", nrow(bexp_all), ncol(bexp_all),max(colSums(is.na(bexp_all))),max(rowSums(is.na(bexp_all))))
cat(string2)

string3 = sprintf("Within individual matrix: There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s. \n", nrow(wexp_all), ncol(wexp_all),max(colSums(is.na(wexp_all))),max(rowSums(is.na(wexp_all))))
cat(string3)

#%%%%%%%%%%%%%%% Save decomposed expression files 
print("Finished decomposition, saving files")

print("Saving between-individuals variation matrix")
fwrite(x = data.table::data.table(t(bexp_all),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]},
       file = paste0(out_dir,"context_shared_expression.txt"), quote = F, row.names = F,
       col.names = T, append = F, sep = '\t')

print("Saving within-individuals variation matrix for context: ")
for(i in 1:length(contexts)){
  print(contexts[i])
  wexp_t = wexp_all[grep(pattern = paste0(contexts[i],"$"), rownames(wexp_all)),]
  rownames(wexp_t)=gsub(pattern = paste0(" - ",contexts[i]), replacement = "", x = rownames(wexp_t))
  fwrite(x = data.table::data.table(t(wexp_t),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]},
         file = paste0(out_dir,contexts[i],"_specific_expression.txt"),quote = F, row.names = F,
         col.names = T, append = F, sep = '\t')
}
}
