#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import logging
from importlib.metadata import version

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

sample = "${prefix}"
data = "${data}"
table = "${params.visium_hd}"

os.makedirs(sample, exist_ok=True)

#read visium_hd dataset
log.info(f"loading Visium HD table {table} for sample {sample}")
ds = sd.visium_hd(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
log.info(f"converting to anndata")
adata = to_legacy_anndata(ds, coordinate_system=sample, 
                          table_name = table, include_images=True)
adata.var_names_make_unique()

#grab cell areas from the geopandas series
areas = ds.shapes[f"{sample}_cell_segmentations"].area
if len(areas) == len(adata.obs_names):
    adata.obs['cell_area'] = areas.values
    log.info(f"added cell_area to adata.obs")
else:
    log.warning(f"Length mismatch: {len(areas)} areas vs {len(adata.obs_names)} observations. cell_area not added.")

log.info(f"adata is {adata}")

#save
outname = os.path.join(sample, "adata.h5ad")
adata.write_h5ad(filename=outname)

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(version("spatialdata_io")))
    f.write("    anndata: {}\\n".format(version("anndata")))