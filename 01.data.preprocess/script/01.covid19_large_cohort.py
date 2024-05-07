import os
import sys
from scipy import sparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


import pyprojroot

proj_root = pyprojroot.here()

# construct anndata
covid19_cohort_mtx = os.path.join(proj_root, "data",
                                  "COVID19_large_cohort", "mtx")

cvdlgann = ad.read_mtx(
    filename = os.path.join(proj_root, "data",
                       "COVID19_large_cohort", "mtx",
                       "GSE158055_covid19_matrix.mtx"),
    dtype = np.uint)

barcodes = pd.read_csv(
    os.path.join(proj_root, "data", "COVID19_large_cohort", "mtx",
                 "GSE158055_covid19_barcodes.tsv"),
    sep = "\t", header = None
)

genes = pd.read_csv(
    os.path.join(proj_root, "data", "COVID19_large_cohort", "mtx",
                 "GSE158055_covid19_genes.tsv"),
    sep = "\t", header = None
)

cvd = cvdlgann.T
del cvdlgann

cvd.obs_names = barcodes.iloc[:, 0]
cvd.var_names = genes.iloc[:, 0]

# add metadata
cellmeta = pd.read_csv(
    os.path.join(proj_root, "data", "COVID19_large_cohort",
                 "GSE158055_cell_annotation.csv"),
    sep = ",", header = 0
)
cellmeta.set_index("cellName", drop = False, inplace = True)

cvd.obs_names.equals(cellmeta.index)
cvd.obs = cellmeta

# add sample-level metadata
samplemeta = pd.read_csv(
   os.path.join(proj_root, "data", "COVID19_large_cohort",
                "GSE158055_sample_metadata.csv"),
   sep = ",", header=0
)

cvd.uns['sample_info'] = samplemeta

cvd.write_h5ad(filename = os.path.join(proj_root, "data", "COVID19_large_cohort",
                                       "covid19.large.cohort.ann.h5ad"))
