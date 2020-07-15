install.packages("Rcpp", repos = "http://cran.us.r-project.org")
install.packages("rstan", repos = "http://cran.us.r-project.org")
install.packages("Seurat", repos = "http://cran.us.r-project.org")
install.packages("devtools", repos = "http://cran.us.r-project.org")
install.packages("tydiverse", repos = "http://cran.us.r-project.org")
install.packages("here", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
install.packages("Matrix", repos = "http://cran.us.r-project.org")
install.packages("optparse", repos = "http://cran.us.r-project.org")
install.packages("bayesplot", repos = "http://cran.us.r-project.org")
install.packages("pROC", repos = "http://cran.us.r-project.org")

devtools::install_github("klmr/modules")


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("TCGAbiolinks")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("EDASeq")



