#!/usr/bin/env python
import pickle
import os
import anndata as ad
import logging
import pandas as pd
from importlib.metadata import version

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def attach_squidpy_ligrec(adata, ligrec_path):
    with open(ligrec_path, "rb") as f:
        ligrec = pickle.load(f)

    log.info(f"got ligrec with keys: {list(ligrec.keys())}")

    # flatten indices to store in adata
    # https://github.com/scverse/squidpy/issues/533
    for k, v in list(ligrec.items()):
        v.index = v.index.to_flat_index().str.join("-")
        if k != "metadata":
            v.columns = v.columns.to_flat_index().str.join("-")
    adata.uns['ligrec_means'] = ligrec['means'].sparse.to_dense()
    adata.uns['ligrec_pvalues'] = ligrec['pvalues'].sparse.to_dense()

    return adata

def attach_spacemarkers_ligrec(adata, ligrec_path):
    log.info(f"reading spacemarkers ligrec from {ligrec_path}")
    ligrec_name = os.path.basename(ligrec_path).split(".")[0]

    if( ligrec_path != "null"):
        ligrec = pd.read_csv(ligrec_path, index_col=0)
    else:
        ligrec = pd.DataFrame() #empty dataframe to pass through adata

    #drop all na rows
    ligrec.dropna(how='all', inplace=True)

    if ligrec.empty:
        log.warning(f"ligrec dataframe is empty after dropping all-NA rows, skipping attachment")
    else:
        adata.uns[ligrec_name] = ligrec

    return adata

if __name__ == "__main__":
    adata_path = "${adata}"
    ligrec_path = "${ligrec}"
    sample = "${meta.id}"

    #output directory
    os.makedirs(sample, exist_ok=True)

    log.info(f"attaching ligand-receptor results to adata.uns")
    adata = ad.read_h5ad(adata_path)
    
    if "pickle" in ligrec_path:
        adata = attach_squidpy_ligrec(adata, ligrec_path)
    else:
        adata = attach_spacemarkers_ligrec(adata, ligrec_path)

    log.info(f"saving adata with attached ligand-receptor results to {sample}/staple.h5ad")
    adata.write_h5ad(f"{sample}/staple.h5ad", compression='gzip')

    #versions
    with open("versions.yml", "w") as f:
        f.write("${task.process}:\\n")
        f.write("    anndata: {}\\n".format(version("anndata")))
        f.write("    pandas: {}\\n".format(version("pandas")))
