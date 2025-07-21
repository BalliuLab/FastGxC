# FastGxC: A powerful and computationally efficient software for context-specific eQTL mapping in single-cell genomics data

FastGxC was originally developed for single-cell data, where each individual contributes gene expression measurements across multiple cell types.

However, it can also be applied to bulk RNA-seq data when the same individuals are profiled across multiple tissues or conditions.

In both settings, FastGxC models **repeated samples** from each individual, removing shared noise and enabling more accurate detection of context-specific genetic effects.

FastGxC is also **robust to missing data** —for example, when certain individuals or genes are missing in some cell types or tissues.

Please read the [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.17.448889v2) preprint for more details.

<!-- Extended data with FastGxC results on GTEx, OneK1K, and CLUES cohorts can be found [here](https://zenodo.org/record/5015123#.YNJ1WpNKjOR) -->

# Package Installation and Dependencies

FastGxC is an R package that can be loaded and used in any R environment. In order for FastGxC to run properly, the following packages are required and should be installed in R prior to using any FastGxC functions.

```         
library(qvalue)
library(devtools)
library(dplyr)
library(data.table)
library(MatrixEQTL)
library(mvtnorm)
library(reshape2)
library(magrittr)
library(TreeQTL)
```

**Note** : To install TreeQTL, qvalue must be installed first.

```         
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
```

To install TreeBH, download the following package: [TreeBH_1.0.tar.gz](https://github.com/user-attachments/files/21328935/TreeBH_1.0.tar.gz)

Then, run the following command:

```         
install.packages("~/Downloads/TreeBH_1.0.tar.gz", repos = NULL, type = "source")
library(TreeBH)
```

Once all dependencies are installed and loaded you can install FastGxC using:

```         
devtools::install_github("BalliuLab/FastGxC")
```

If devtools is not installed:

```         
install.packages("devtools")
```

Once FastGxC is installed, load all functions using

```         
library(FastGxC)
```

Optional: Set up Python dependencies for TensorQTL (used with method = "tensor")\
Requires: reticulate + compatible Python environment

To check your Python setup, run in your terminal:

```         
which python3
```

Install required Python packages (Terminal only)

```         
/your/python/path/here -m pip install numpy pandas pyarrow rpy2 torch tensorqtl
```

Next, start RStudio and set Python path

```         
library(reticulate)
use_python("/your/python/path/here", required = TRUE)
py <- import_builtins()
```

# Simulate toy data

To run a toy example, generate simulated data by running the following code in R:

```         
n_contexts = 10
  data_dir_sim= "~/simulated_example/"    
  simulate_data(data_dir = data_dir_sim, # Path to directory where output files will be saved   
              N = 300, # Number of individuals
              n_contexts = n_contexts, # Number of contexts
              n_genes = 100, # Number of Genes
              n_snps_per_gene = 1000, # Number of SNPs per gene
              maf = 0.2, # Minor allele frequency for SNPs
              w_corr = 0.2, # Intra-individual residual correlation between contexts
              v_e = 1, # Error variance per context
              missing = 0.05, # Fraction of missing values in the simulated expression matrix
              hsq = rep(0, n_contexts), # Vector with proportion of expression heritability                   explained by eQTL in each context
              mus = rep(0, n_contexts)) # Vector with average expression in each context 
```

**Note**: running the code above simulates data with default parameters (300 individuals, 1,000 SNPs, 100 genes, 10 contexts, etc.), but this function can be run with any combination of parameter values. See all possible parameters for `simulate_data()` by running `?simulate_data` in R.

**Output Files**\
The simulation will generate the following files in `data_dir`:

**SNPs.txt**

```         
snpid   ind1  ind2  ind3  ind4  ind5 ...
snp1     0     2     1     1     0   ...
snp2     1     2     0     1     1   ...
snp3     0     1     0     1     1   ...
snp4     1     2     0     1     1   ...
snp5     0     0     0     1     1   ...
```

**snpsloc.txt**

```         
snpid   chr   pos   context1  context2  ...  context10
snp1    chr1   1        1         1           1
snp2    chr1   1        1         1           1
snp3    chr1   1        1         1           1
snp4    chr1   1        1         1           1
snp5    chr1   1        1         1           1
```

**expression.txt**

```         
design             gene1     gene2     gene3    ...   gene100
ind1 - context1   0.4369    -1.7087    -1.8113  ...   1.1721
ind1 - context10  0.6437    -0.4357    -0.8296  ...  -0.3944
ind1 - context2   0.1092    -0.1294    -0.3222  ...  -0.4939
ind1 - context3   0.1234    -0.1947    -0.1111  ...  -0.9999
ind1 - context4   0.8347    -0.9876    -0.1490  ...  -0.3833
```

**geneloc.txt**

```         
geneid  chr   s1          s2          context1  context2  context3  context4  context5  context6
gene1   chr1  0           1001            1         1         1         1         1         1
gene2   chr1  100000001   100001001       1         1         1         1         1         1
gene3   chr1  200000001   200001001       1         1         1         1         1         1
gene4   chr1  300000001   300001001       1         1         1         1         1         1
gene5   chr1  400000001   400001001       1         1         1         1         1         1
gene6   chr1  500000001   500001001       1         1         1         1         1         1
```

# Running FastGxC

FastGxC works in two steps. In the first step, expression is decomposed into shared and context-specific components. In the second step, eQTLs are separately mapped on these components.

*Step 1 - Decomposition:* For each individual, decompose the phenotype of interest (e.g. gene expression) across C contexts (e.g. tissues or cell-types) into one context-shared and C context-specific components using the `decomposition_step()` function. This function takes as imput a file with gene expression data for all individuals, genes, and contexts (see output of `simulate_data()` for the right format) and outputs one file with context-shared expression (context_shared_expression.txt) and C files with expression specific to each context (CONTEXT_NAME_specific_expression.txt).

The following code example demonstrates how to use this function with the data we just simulated above.

```         
data_dir_decomp = "~/simulated_example/"
exp_mat_filename = paste0(data_dir_decomp, "expression.txt")
decomposition_step(exp_mat_filename, data_dir_decomp)
```

**Output Files**

context_shared_expression.txt / contextX_specific_expression.txt

```         
        ind1        ind2        ind3        ind4        ind5
gene1  0.1834      0.0456     -0.0912      0.2314     -0.0873
gene2 -0.1345      0.0543      0.1123     -0.0784      0.0912
gene3  0.0034     -0.0312      0.0402      0.0105     -0.0067
gene4 -0.1423     -0.1092     -0.0876     -0.1904     -0.1523
gene5  0.2145      0.1432      0.1784      0.1984      0.2013
```

*Step 2 - eQTL mapping:* FastGxC estimates genetic effects on the context-shared component and each of the C context-specific components separately using simple linear models. Note: Here we use the R package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) but any other software that can perform quick linear regression can be used (e.g. [FastQTL](http://fastqtl.sourceforge.net/) or [tensorqtl](https://github.com/broadinstitute/tensorqtl)). FastGxC implements eQTL mapping using its `eQTL_mapping_step()` function.

This function take as input data needed to run MatrixEQTL and outputs eQTL summary statistics in the MatrixEQTL format. In the end, you should have one file with summary statistics for shared eQTL and C files with summary statistics for each context C.

Below is a code example to map context-specific eQTLs and shared eQTLs using the decomposed simulated data from above.

```         
input_dir = "~/simulated_example/"
out_dir = input_dir
expr_files <- list.files(input_dir, pattern = "_specific_expression.txt$")
nC <- length(expr_files)
context_names <- paste0("context", seq(1, nC))

all_contexts <- c(context_names, "shared")
all_shared_specific <- c(rep("specific", nC), "shared")

SNP_file_name <- file.path(input_dir, "SNPs.txt")
snps_location_file_name <- file.path(input_dir, "snpsloc.txt")
gene_location_file_name <- file.path(input_dir, "geneloc.txt")

for (i in seq_along(all_contexts)) {
  context <- all_contexts[i]
  shared_specific <- all_shared_specific[i]
  
  if (shared_specific == "shared") {
    expression_file_name <- file.path(input_dir, "context_shared_expression.txt")
  } else {
    expression_file_name <- file.path(input_dir, paste0(context, "_specific_expression.txt"))
  }

  eQTL_mapping_step(
    SNP_file_name = SNP_file_name,
    snps_location_file_name = snps_location_file_name,
    expression_file_name = expression_file_name,
    gene_location_file_name = gene_location_file_name,
    context = context,
    shared_specific = shared_specific,
    out_dir = out_dir,
    output_file_name_cis = file.path(out_dir, paste0(context, "_", shared_specific, ".cis_pairs.txt")),
    output_file_name_tra = file.path(out_dir, paste0(context, "_", shared_specific, ".trans_pairs.txt")),
    method = "MatrixEQTL",
    use_model = modelLINEAR,
    cis_dist = 1e6,
    pv_threshold_cis = 1,
    pv_threshold_tra = 0,
    error_covariance = numeric()
  )
}
```

*Step 2 - eQTL mapping with TensorQTL (optional alternative):* FastGxC also supports cis-eQTL mapping using TensorQTL, a GPU-accelerated implementation for fast eQTL analysis. This is enabled through the R package reticulate, which allows Python functions to be called from within R.

To use this functionality, specify method = "tensorQTL" when calling eQTL_mapping(). All input formats remain the same as for MatrixEQTL. The function returns results in a format consistent with MatrixEQTL outputs for compatibility with downstream steps. TensorQTL is particularly useful when working with large datasets due to its computational efficiency.

**Note:** TensorQTL requires Python version 3.7 or higher to run properly.

To run TensorQTL, simply pass method = TensorQTL. The inputs remain unchanged.

**Output Files**

contextX_specific.cis_pairs.txt / shared_shared.cis_pairs.txt

```         
SNP         gene    beta                    p.value                 FDR
snp64928    gene65  0.418661087751389   8.4807042218939e-07     0.084807042218939
snp42566    gene43  -0.355810463428497  7.53491371859939e-06    0.306773379105091
snp22672    gene23  -0.394770801067352  9.20320137315272e-06    0.306773379105091
snp73350    gene74  0.344478219747543   1.59910672991937e-05    0.399776682479843
snp86064    gene87  -0.371489137411118  3.10237595856035e-05    0.550542627746407
```

**Note:** TensorQTL writes intermediate results as .parquet files by default, but FastGxC saves the final output in .txt format for consistency with MatrixEQTL outputs.

# Multiple testing adjustment

## With `TreeQTL`

To adjust for multiple testing across all contexts, genes, and genetic variants tested, FastGxC uses the hierarchical FDR procedures implemented in the R package [TreeQTL](http://bioinformatics.org/treeqtl/) via the `treeQTL_step()` function.

This function requires that you run MatrixEQTL to do eQTL mapping (see step 2 above). If you used another eQTL mapping softwares, please make sure the output is in the format required by TreeQTL. You can also replace TreeQTL with other methods, e.g. [mashR](https://github.com/stephenslab/mashr), which can also lead to a considerable increase in power.

The following code example demonstrates how to use this function with the data outputted from the eQTL mapping we just ran. This script take as input data needed to run TreeQTL and outputs shared and specific eGenes (two files) and eAssociation (C+1 files) summary statistics in the TreeQTL format.

Map specific-eGenes, i.e., genes with at least one context-specific eQTL

```         
## directory with all matrixeQTL files for all contexts 
## (note that this function expects files to be named in the same was as the 
output files from FastGxC's eQTL mapping function)
data_dir = "~/simulated_example/"
snps_location_file_name = "~/simulated_example/snpsloc.txt"
gene_location_file_name = "~/simulated_example/geneloc.txt"
out_dir = "~/simulated_example/"
fdr_thresh = 0.05

## Run multiple testing correction without the four level hierarchy:
## this will only output specfic eGenes and eAssociations
treeQTL_step(
       data_dir,
       snps_location_file_name,
       gene_location_file_name,
       context_names,
       out_dir,
       fdr_thresh = fdr_thresh
     )

## Run multiple testing correction with the four level hierarchy: 
## this will output specific and shared eGenes and eAssociations as well as 
global eGenes
treeQTL_step(
       data_dir,
       snps_location_file_name,
       gene_location_file_name,
       context_names,
       out_dir,
       fdr_thresh = fdr_thresh,
       four_level = T
     )
```

## With `TreeBH`

[TreeBH](https://github.com/cbpeterson/TreeBH) is the updated version of TreeQTL and supports the implementation of custom and more complex hierarchies. FastGxC provides the functions `to_TreeBH_input()` and `treeBH_step()` to leverage TreeBH's capibilities.

The following code example demonstrates how to use these functions with the data outputted from the eQTL mapping step. In the example, the output of TreeBH will be saved to the path `~/simulated_example/treeBH_input.txt`.

``` r
context_names <- paste0("context", 1:10)
shared_file <- "~/simulated_example/shared_shared.cis_pairs.txt"
data_dir <- "~/simulated_example/"
out_dir <- "~/simulated_example"

FastGxC::to_TreeBH_input(
  data_dir = data_dir,
  shared_file = shared_file, 
  context_names = context_names,
  out_dir = out_dir
)

out_dir <- "~/simulated_example/"
df <- read.table(out_dir)
fdr_thres <- 0.05

FastGxC::treeBH_step(
  matrix = df,
  fdr_thres = 0.05,
  out_dir
)
```
