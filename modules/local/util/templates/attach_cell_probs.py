#!/usr/bin/env python
import anndata as ad
import pandas as pd
import os
import logging
from importlib.metadata import version

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

adata_path = "${adata}"
sample = "${sample}"
probs_path = "${cell_probs}"
out_name = "${out_name}"

#output directory
os.makedirs(sample, exist_ok=True)

#read inputs
adata = ad.read_h5ad(adata_path)
cell_probs = pd.read_csv(probs_path)
if('barcode' in cell_probs.columns):
    log.info("found barcode, setting as index")
    cell_probs.set_index('barcode', inplace=True)
else:
    log.warning("no barcode column found in cell_probs, assuming first column is index")
    cell_probs.set_index(cell_probs.columns[0], inplace=True)
    cell_probs.index.name = 'index'
    cell_probs.index = cell_probs.index.astype(str)

#put cell types to uns
adata.uns['cell_type_composition'] = cell_probs

#save most abundant cell type to obs
most_abundant = cell_probs.idxmax(axis=1)
adata.obs['cell_type'] = most_abundant
adata.obs['cell_type_prob'] = cell_probs.max(axis=1)
adata.obs['cell_type_zscore'] = adata.uns['cell_type_composition']\
    .apply(lambda x: (x-x.mean())/x.std(), axis=1).max(axis=1)

#save
adata.write_h5ad(f"{sample}/{out_name}.h5ad", compression='gzip')

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(version("anndata")))
    f.write("    pandas: {}\\n".format(version("pandas")))


