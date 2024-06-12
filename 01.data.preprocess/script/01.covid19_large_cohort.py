import os
import sys
from subprocess import Popen
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
projd = pyprojroot.here()
R: str = "/home/szu/miniforge3/envs/r/bin/Rscript"
scanpyann2seurat_script = os.path.join(projd, "tools", "singlecell",
                                      "ScanpyAnnData2Seurat.R")
covid19_cohort_dir = os.path.join(projd,
                                  "data", "COVID19_large_cohort")
covid19_mt_dir = os.path.join(covid19_cohort_dir, "majortype")
os.makedirs(covid19_mt_dir, exist_ok = True)
conda = "/home/szu/miniforge3/bin/conda"


ann = ad.read_h5ad(
    filename = os.path.join(covid19_cohort_dir,
        "covid19.large.cohort.ann.h5ad")
)
majorTypes = ann.obs.majorType.unique().to_list()

def subset_by_mt(mt:str) -> ad.AnnData:
    r: ad.AnnData = ann[ann.obs.majorType == mt].copy()
    r.write_h5ad(filename = os.path.join(covid19_mt_dir,
                                         f"covid19.large.{mt}.h5ad"))
    return(r)
for mt in majorTypes:
    print(mt)
    subset_by_mt(mt)

def get_ann2seurat_exp(mt: str) -> str:
    annfnm = os.path.join(covid19_mt_dir, f"covid19.large.{mt}.h5ad")
    outfnm = os.path.join(covid19_mt_dir, f"covid19.large.{mt}.seu.rds")
    r = [R, scanpyann2seurat_script, f"-f {annfnm}",
         f"-o {outfnm}", f"--conda {conda}", f"--condaenv sa2"]
    return ' '.join(r)
    
for mt in majorTypes:
    print(mt)
    cmd = get_ann2seurat_exp(mt = mt)
    p = Popen(cmd, shell = True)
    p.wait()

