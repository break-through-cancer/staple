#!/usr/bin/env python
import anndata as ad
import os
import logging
from importlib.metadata import version

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

adata_path = "${adata}"
sample = "${meta.id}"
metadata = "${meta}"

#output directory
os.makedirs(sample, exist_ok=True)

#read inputs
adata = ad.read_h5ad(adata_path)

#parse groovy meta map to dict, example map: metadata="[response:, id:sample1]"
flist = metadata.replace("[","").replace("]","").strip()
fdict = {x.split(":",1)[0].strip():x.split(":",1)[1].strip() for x in flist.split(",")}

#attach metadata to obs if not blank
for key in fdict:
    if key in adata.obs.columns:
        log.warning(f"metadata {key} already exists in adata.obs, skipping")
        continue
    else:
        if(fdict[key] != ""):
            log.info(f"attaching metadata {key}: {fdict[key]} to adata.obs")
            adata.obs[key] = fdict[key]
        else:
            log.warning(f"metadata {key} is blank, skipping")

#save added fields to uns
field_list = list(fdict.keys())
adata.uns['staple_meta_fields'] = field_list
log.info(f"added metadata fields: {field_list} to '.obs' and recorded in '.uns[staple_meta_fields]'")

#save
log.info(f"saving adata with attached metadata to {sample}/adata.h5ad")
adata.write_h5ad(f"{sample}/adata.h5ad", compression='gzip')

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(version("anndata")))
