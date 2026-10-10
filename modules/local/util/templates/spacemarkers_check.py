#!/usr/bin/env python
import sys
import logging
import anndata as ad
import pandas as pd
from importlib.metadata import version

logging.basicConfig(level=logging.INFO, stream=sys.stderr)
log = logging.getLogger()

adata_path = "${adata}"
sample = "${meta.id}"
key = "${params.sm_patterns_uns}"

adata = ad.read_h5ad(adata_path, backed="r")

# nested keys are separated by chr(36), same as in SpaceMarkers
node = adata.uns
for part in key.split(chr(36)):
    node = node.get(part) if hasattr(node, "get") else None
    if node is None:
        break

eligible = (
    isinstance(node, pd.DataFrame)
    and node.shape[0] > 0
    and node.shape[1] >= 2
    and node.shape[1] == node.select_dtypes("number").shape[1]
)

if not eligible:
    log.warning(
        f"SpaceMarkers skipped for {sample}: no usable latent features in uns['{key}'] "
        "(need a numeric spots x patterns table with at least 2 patterns)"
    )

# stdout is consumed by the workflow
sys.stdout.write("true" if eligible else "false")

with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(version("anndata")))
    f.write("    pandas: {}\\n".format(version("pandas")))
