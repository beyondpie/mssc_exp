import os
import sys
from scipy import sparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


import pyprojroot

proj_root = pyprojroot.here()

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

cvd.write_h5ad(filename = os.path.join(proj_root, "data", "COVID19_large_cohort",
                                       "covid19.large.cohort.ann.h5ad"))

# * to Seurat, organized by major types
R: str = os.path.join("/home/szu/miniforge3/envs/r/bin/Rscript")
scipyann2seurat_script = os.path.join
ann = ad.read_h5ad(
    filename = os.path.join(
        proj_root, "data", "COVID19_large_cohort",
        "covid19.large.cohort.ann.h5ad")
)

majorTypes = ann.obs.majorType.unique().to_list()

