library(TCGAWorkflowData)
library(data.table)
library(here)
library(tidyverse)

cancer_project <- "TCGA-UVM"

query <- GDCquery(project = cancer_project,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Queantification")
