#!/usr/bin/env python
import anndata as ad
import logging
import os
from importlib.metadata import version

def adata_preprocess(adata, drop_prefix):
    if drop_prefix:
        adata = adata[:, ~adata.var_names.str.startswith(drop_prefix)]
    adata.file.close()
    return adata

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    log = logging.getLogger()

    drop_prefix = "${params.drop_genes_prefix}"
    adata_path = "${adata}"
    sample = "${meta.id}"

    #parse drop prefix in case some provided
    if drop_prefix:
        drop_prefix = tuple(drop_prefix.split(","))

    #output directory
    os.makedirs(sample, exist_ok=True)
    
    adata = ad.read_h5ad(adata_path)
    log.info(f"Before preprocessing adata: {adata}")

    #preprocess adata
    adata = adata_preprocess(adata, drop_prefix)
    
    #save
    log.info(f"After preprocessing adata: {adata}")
    adata.write_h5ad(f"{sample}/adata.h5ad", compression='gzip')

    #record versions
    with open("versions.yml", "w") as f:
        f.write("${task.process}:\\n")
        f.write("    anndata: {}\\n".format(version("anndata")))