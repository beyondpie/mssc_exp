import scanpy as sc
from os import PathLike, fspath
import gzip

def read_mtx(filename: PathLike, dtype: str = "float32") -> AnnData:
    """\
    Read `.mtx` file.

    Parameters
    ----------
    filename
        The filename.
    dtype
        Numpy data type.
    """
    from scipy.io import mmread

    # could be rewritten accounting for dtype to be more performant
    X = mmread(fspath(filename)).astype(dtype)
    from scipy.sparse import csr_matrix

    X = csr_matrix(X)
    return AnnData(X, dtype=dtype)

def _read_v3_10x_mtx(
    path,
    var_names='gene_symbols',
    make_unique=True,
    cache=False,
    cache_compression=_empty,
):
    """
    Read mex from output from Cell Ranger v3 or later versions
    """
    path = Path(path)
    adata = read(
        path / 'matrix.mtx.gz', cache=cache, cache_compression=cache_compression
    ).T  # transpose the data
    genes = pd.read_csv(path / 'features.tsv.gz', header=None, sep='\t')
    if var_names == 'gene_symbols':
        var_names = genes[1].values
        if make_unique:
            var_names = anndata.utils.make_index_unique(pd.Index(var_names))
        adata.var_names = var_names
        adata.var['gene_ids'] = genes[0].values
    elif var_names == 'gene_ids':
        adata.var_names = genes[0].values
        adata.var['gene_symbols'] = genes[1].values
    else:
        raise ValueError("`var_names` needs to be 'gene_symbols' or 'gene_ids'")
    adata.var['feature_types'] = genes[2].values
    adata.obs_names = pd.read_csv(path / 'barcodes.tsv.gz', header=None)[0].values
    return adata
