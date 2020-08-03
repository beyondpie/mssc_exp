if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }
BiocManager::install("splatter")


## * SymSim
devtools::install_github(repo = "beyondpie/SymSim",
                         ref = "szu")
