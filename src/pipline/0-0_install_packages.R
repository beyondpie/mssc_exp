install.packages("Rcpp")
install.packages("rstan")
install.packages("Seurat")
install.packages("devtools")
install.packages("tydiverse")
install.pacakges("here")
install.pacakges("data.table")
install.packages("ggplot2")
install.packages("Matrix")
install.packages("optparse")
install.packages("bayesplot")
install.packages("pROC")

devtools::install_github("klmr/modules")


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("TCGAbiolinks")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("EDASeq")



