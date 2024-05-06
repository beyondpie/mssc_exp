import os
import sys
import scanpy as sc
import anndata as ad

import pyprojroot

proj_root = pyprojroot.here()

covid19_cohort_mtx = os.path.join(proj_root, "data",
                                  "COVID19_large_cohort", "mtx")
cvdlgann = sc.read_10x_mtx(path = covid19_cohort_mtx,
                           prefix = "GSE158055_covid19_")

