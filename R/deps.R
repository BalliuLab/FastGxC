# R/deps.R

#' Install FastGxC dependencies
#'
#' Installs CRAN, Bioconductor, and source packages required by FastGxC.
#' This includes optional packages like qvalue, TreeQTL, and TreeBH.
#'
#' @param ask Logical, whether to ask before installing (default: TRUE if interactive).
#' @export
install_deps <- function(ask = interactive()) {
  # --- Core CRAN packages (should already install via Imports:) ---
  cran_pkgs <- c("devtools","dplyr","data.table","MatrixEQTL",
                 "mvtnorm","reshape2","magrittr","reticulate")
  to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
  if (length(to_install)) {
    msg <- paste("Installing missing CRAN packages:", paste(to_install, collapse = ", "))
    if (!ask || utils::menu(c("Yes","No"), title = msg) == 1) {
      install.packages(to_install)
    }
  }
  
  # --- Bioconductor package: qvalue ---
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    message("Installing Bioconductor package: qvalue")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("qvalue")
  }
  
  # --- Source-only package: TreeQTL ---
  if (!requireNamespace("TreeQTL", quietly = TRUE)) {
    message("Installing TreeQTL from source tarball")
    utils::install.packages("http://bioinformatics.org/treeqtl/TreeQTL_2.0.tar.gz",
                            repos = NULL, type = "source")
  }
  
  # --- Source-only package: TreeBH ---
  if (!requireNamespace("TreeBH", quietly = TRUE)) {
    message("Installing TreeBH from source tarball")
    utils::install.packages("https://github.com/user-attachments/files/21328935/TreeBH_1.0.tar.gz",
                            repos = NULL, type = "source")
  }
  
  message("All FastGxC dependencies are installed or already present.")
}
