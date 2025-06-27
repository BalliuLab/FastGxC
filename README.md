# FastGxC: A powerful and computationally efficient software for context-specific eQTL mapping in single-cell omics data

FastGxC was originally developed for single-cell data, where each individual contributes gene expression measurements across multiple cell types.

However, it can also be applied to bulk RNA-seq data when the same individuals are profiled across multiple tissues or conditions.

In both settings, FastGxC models **repeated samples** from each individual, removing shared noise and enabling more accurate detection of context-specific genetic effects.

FastGxC is also **robust to missing data** —for example, when certain individuals or genes are missing in some cell types or tissues.

Please read the [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.17.448889v2) preprint for more details.

<!-- Extended data with FastGxC results on GTEx, OneK1K, and CLUES cohorts can be found [here](https://zenodo.org/record/5015123#.YNJ1WpNKjOR) -->

# Package Installation and Dependencies

FastGxC is an R package that can be loaded and used in any R environment. In order for FastGxC to run properly, the following packages are required and should be installed in R prior to using any FastGxC functions.

```         
library(reticulate)
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

\*\* Note \*\* : to install TreeQTL, qvalue must be installed first.

Once all dependencies are installed and loaded you can install FastGxC using:

```         
  devtools::install_github("BalliuLab/FastGxC")
```

Once FastGxC is installed, load all functions using

```         
library(FastGxC)
```

Optional: Set up Python dependencies for TensorQTL (used with method = "tensor")\
Requires: reticulate + compatible Python environment

```
library(reticulate)
# Replace '/your/python/path/here' with the full path to your Python 3 executable (e.g., from `which python3` in terminal)
use_python("/your/python/path/here", required = TRUE)

reticulate::py_install(c(
  "tensorqtl",
  "pandas",
  "numpy",
  "pyarrow",
  "torch",
  "rpy2"
))

if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
```

# Simulate toy data

To run a toy example, generate simulated data by running the following code in R:

```         
  data_dir_sim = "~/simulations/"
  sim_scenario = "single_context_het"
  simulate_data(data_dir = data_dir_sim, sim_scenario = sim_scenario)
```

\*\* Note: running the code above simulates data with default parameters (300 individuals, 1,000 SNPs, 100 genes, and 10 contexts without missing data), but this function can be run with any combination of parameter values. See all possible parameters for `simulate_data()` by running `?simulate_data` in R.

Running the above code will generate and save the following files in the data_dir: (1) {sim_scenario}\_SNPs.txt: SNP genotype data for 1,000 SNPs and 300 individuals (individual IDs as columns and SNP IDs as rows)\

(2) {sim_scenario}\_snpsloc.txt: location information of the 1,000 simulated SNPs (MatrixEQTL input format)\

(3) {sim_scenario}\_geneloc.txt: location information of the 100 simulated genes (MatrixEQTL input format)\

(4) {sun_scenario}\_simulated_expression.txt: gene expression data for the 300 simulated individuals across 100 genes and 10 contexts\

# Running FastGxC

FastGxC works in two steps. In the first step, expression is decomposed into shared and context-specific components. In the second step, eQTLs are separately mapped on these components.

*Step 1 - Decomposition:* For each individual, decompose the phenotype of interest (e.g. gene expression) across C contexts (e.g. tissues or cell-types) into one context-shared and C context-specific components using the `decomposition_step()` function. This function takes as imput a file with gene expression data for all individuals, genes, and contexts (see output of `simulate_data()` for the right format) and outputs one file with context-shared expression (context_shared_expression.txt) and C files with expression specific to each context (CONTEXT_NAME_specific_expression.txt).

The following code example demonstrates how to use this function with the data we just simulated above.

```         
exp_mat_filename = "~/simulations/single_context_het_simulated_expression.txt"
data_dir_decomp = "~/example_output_single_context_het/"
decomposition_step(exp_mat_filename, data_dir_decomp)
```

*Step 2 - eQTL mapping:* FastGxC estimates genetic effects on the context-shared component and each of the C context-specific components separately using simple linear models. Note: Here we use the R package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) but any other software that can perform quick linear regression can be used (e.g. [FastQTL](http://fastqtl.sourceforge.net/) or [tensorqtl](https://github.com/broadinstitute/tensorqtl)). FastGxC implements eQTL mapping using its `eQTL_mapping_step()` function.

This function take as input data needed to run MatrixEQTL and outputs eQTL summary statistics in the MatrixEQTL format. In the end, you should have one file with summary statistics for shared eQTL and C files with summary statistics for each context C.

Below is a code example to map context-specific eQTLs and shared eQTLs using the decomposed simulated data from above.

```         
out_dir = "~/example_output_single_context_het/"
input_dir = "~/simulations/"
expr_files <- list.files(out_dir, pattern = "_specific_expression.txt$")
nC <- length(expr_files)
context_names <- paste0("context", seq(1, nC))

SNP_file_name <- file.path(input_dir, "single_context_het_SNPs.txt")
snps_location_file_name <- file.path(input_dir, "single_context_het_snpsloc.txt")
gene_location_file_name <- file.path(input_dir, "single_context_het_geneloc.txt")

for (context in context_names) {
  expression_file_name <- file.path(out_dir, paste0(context, "_specific_expression.txt"))
  shared_specific <- "specific"
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  eQTL_mapping_step(
    SNP_file_name = SNP_file_name,
    snps_location_file_name = snps_location_file_name,
    expression_file_name = expression_file_name,
    gene_location_file_name = gene_location_file_name,
    context = context,
    shared_specific = shared_specific,
    out_dir = out_dir,
    output_file_name_cis = paste0(out_dir, "/", context, "_", shared_specific, ".all_pairs.txt"),
    output_file_name_tra = paste0(out_dir, "/", context, "_", shared_specific, ".trans_pairs.txt"),
    method = "MatrixEQTL",
    python_dir,
    use_model = use_model,
    cis_dist = cis_dist,
    pv_threshold_cis = pv_threshold_cis,
    pv_threshold_tra = pv_threshold_tra,
    error_covariance = error_covariance
  )
}

expression_file_name <- file.path(out_dir, "context_shared_expression.txt")
shared_specific <- "shared"
context <- "shared"

eQTL_mapping_step(
  SNP_file_name = SNP_file_name,
  snps_location_file_name = snps_location_file_name,
  expression_file_name = expression_file_name,
  gene_location_file_name = gene_location_file_name,
  context = context,
  shared_specific = shared_specific,
  out_dir = out_dir,
  output_file_name_cis = paste0(out_dir, "/", context, "_", shared_specific, ".all_pairs.txt"),
  output_file_name_tra = paste0(out_dir, "/", context, "_", shared_specific, ".trans_pairs.txt"),
  method = method,
  python_dir = python_dir,
  use_model = use_model,
  cis_dist = cis_dist,
  pv_threshold_cis = pv_threshold_cis,
  pv_threshold_tra = pv_threshold_tra,
  error_covariance = error_covariance
)
```

*Step 2 - eQTL mapping with TensorQTL (optional alternative):* FastGxC also supports cis-eQTL mapping using TensorQTL, a GPU-accelerated implementation for fast eQTL analysis. This is enabled through the R package reticulate, which allows Python functions to be called from within R.

To use this functionality, specify method = "tensorQTL" when calling eQTL_mapping(). All input formats remain the same as for MatrixEQTL. The function returns results in a format consistent with MatrixEQTL outputs for compatibility with downstream steps. TensorQTL is particularly useful when working with large datasets due to its computational efficiency.

**Note:** TensorQTL requires Python version 3.7 or higher to run properly.

Below is an example of how to perform eQTL mapping using TensorQTL on context-specific and shared components:

```         
# Store the path in a variable for downstream use (e.g., inside function calls)
python_dir <- "/your/python/path/here"

out_dir = "~/example_output_single_context_het/"
input_dir = "~/simulations/"
expr_files <- list.files(out_dir, pattern = "_specific_expression.txt$")
nC <- length(expr_files)
context_names <- paste0("context", seq(1, nC))

SNP_file_name <- file.path(input_dir, "single_context_het_SNPs.txt")
snps_location_file_name <- file.path(input_dir, "single_context_het_snpsloc.txt")
gene_location_file_name <- file.path(input_dir, "single_context_het_geneloc.txt")

for (context in context_names) {
  expression_file_name <- file.path(out_dir, paste0(context, "_specific_expression.txt"))
  shared_specific <- "specific"
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  eQTL_mapping_step(
    SNP_file_name = SNP_file_name,
    snps_location_file_name = snps_location_file_name,
    expression_file_name = expression_file_name,
    gene_location_file_name = gene_location_file_name,
    context = context,
    shared_specific = shared_specific,
    out_dir = out_dir,
    output_file_name_cis = paste0(out_dir, "/", context, "_", shared_specific, ".all_pairs.txt"),
    output_file_name_tra = paste0(out_dir, "/", context, "_", shared_specific, ".trans_pairs.txt"),
    method = "tensorQTL",
    python_dir = python_dir,
    use_model = use_model,
    cis_dist = cis_dist,
    pv_threshold_cis = pv_threshold_cis,
    pv_threshold_tra = pv_threshold_tra,
    error_covariance = error_covariance
  )
}

expression_file_name <- file.path(out_dir, "context_shared_expression.txt")
shared_specific <- "shared"
context <- "shared"

eQTL_mapping_step(
  SNP_file_name = SNP_file_name,
  snps_location_file_name = snps_location_file_name,
  expression_file_name = expression_file_name,
  gene_location_file_name = gene_location_file_name,
  context = context,
  shared_specific = shared_specific,
  out_dir = out_dir,
  output_file_name_cis = paste0(out_dir, "/", context, "_", shared_specific, ".all_pairs.txt"),
  output_file_name_tra = paste0(out_dir, "/", context, "_", shared_specific, ".trans_pairs.txt"),
  method = method,
  python_dir = python_dir,
  use_model = use_model,
  cis_dist = cis_dist,
  pv_threshold_cis = pv_threshold_cis,
  pv_threshold_tra = pv_threshold_tra,
  error_covariance = error_covariance
)
```

Note: TensorQTL writes intermediate results as .parquet files by default, but FastGxC saves the final output in .txt format for consistency with MatrixEQTL outputs.

# Multiple testing adjustment

To adjust for multiple testing across all contexts, genes, and genetic variants tested, FastGxC uses the hierarchical FDR procedures implemented in the R package [TreeQTL](http://bioinformatics.org/treeqtl/) via the `treeQTL_step()` function.

This function requires that you run MatrixEQTL to do eQTL mapping (see step 2 above). If you used another eQTL mapping softwares, please make sure the output is in the format required by TreeQTL. You can also replace TreeQTL with other methods, e.g. [mashR](https://github.com/stephenslab/mashr), which can also lead to a considerable increase in power.

The following code example demonstrates how to use this function with the data outputted from the eQTL mapping we just ran. This script take as input data needed to run TreeQTL and outputs shared and specific eGenes (two files) and eAssociation (C+1 files) summary statistics in the TreeQTL format.

Map specific-eGenes, i.e., genes with at least one context-specific eQTL

```         
## directory with all matrixeQTL files for all contexts 
## (note that this function expects files to be named in the same was as the 
output files from FastGxC's eQTL mapping function)
data_dir = "~/example_output_single_context_het/" 
snps_location_file_name = "~/simulations/single_context_het_snpsloc.txt"
gene_location_file_name = "~/simulations/single_context_het_geneloc.txt"
out_dir = "~/example_output_single_context_het/"
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
