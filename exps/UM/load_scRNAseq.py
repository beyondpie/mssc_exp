import os
import re
from pathlib import Path, PurePath
from typing import List, Dict

import GEOparse
import harmonypy as hm
import matplotlib.pyplot as plt
import scipy
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pandas import DataFrame, Series
from pyprojroot import here

# * configs
datas: str = "data"
tumor: str = "UM"
mygeo: str = "GSE139829"

# * get scRNAseq metadata.
mygse = GEOparse.get_GEO(geo=mygeo, destdir=here(PurePath(datas, tumor)))
sample_id: List[str] = mygse.metadata['sample_id']
genders: Series = mygse.phenotype_data["characteristics_ch1.4.Sex"]

# * get the count matrix.
mydatat = here(PurePath(datas, tumor, f"{mygeo}_RAW"))
suffix: str = "*matrix.mtx.gz"
scRNAseqs: Dict[str, Path] = {
    re.sub("(_.*)$", "", x.name): x
    for x in mydatat.glob(suffix)
}

p = sample_id[1]
cnt = sc.read_10x_mtx(scRNAseqs[p])

dtype="float32"
with open(scRNAseqs[p], 'r') as f:
    x = scipy.io.mmread(f)
    a = csr_matrix(x)
