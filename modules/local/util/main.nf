process ATLAS_MATCH {
    //return adata_sc with gene index matching adata_st by gene name or gene id
    tag "$meta.id"
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    //adata_sc adata_sc
    //adata_st adata_st
    input:
    tuple val(meta), path(adata_sc), path(adata_st)
    output:
    tuple val(meta), path("${prefix}/adata_matched.h5ad"), emit: adata_matched
    path "versions.yml",                                   emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3
import os
import anndata as ad
import numpy as np

with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    numpy: {}\\n".format(np.__version__))

print("Reading adata_sc in the backed mode")
adata_sc = ad.read_h5ad("$adata_sc", backed='r')
print("adata_sc:")
print(adata_sc)

print("Reading adata_st")
adata_st = ad.read_h5ad("$adata_st")
print("adata_st:")
print(adata_st)

os.makedirs("${prefix}", exist_ok=True)

#look for matching indices
matching_index = adata_sc.var.index.intersection(adata_st.var.index)
print(f"Found {len(matching_index)} matching genes in var.index")

#look for adata_sc.index in var["gene_ids"] of adata_st
if 'gene_ids' in adata_st.var.columns:
    matching_gene_ids = adata_sc.var.index.intersection(adata_st.var["gene_ids"])
    print(f"Found {len(matching_gene_ids)} matching genes in var[gene_ids]")
else:
    matching_gene_ids = []

#look for adata_st.index in adata_sc.var["feature_names"]
if 'feature_name' in adata_sc.var.columns:
    matching_feature_names = adata_st.var.index.intersection(adata_sc.var["feature_name"])
    print(f"Found {len(matching_feature_names)} matching genes in var[feature_name]")
else:
    matching_feature_names = []

#find largest matching case
matching_lengths = [len(x) for x in [matching_index, matching_gene_ids, matching_feature_names]]
which_matching = np.argmax(matching_lengths)

if matching_lengths[which_matching] == 0:
    raise RuntimeError("no matching genes found")

if which_matching == 0:
    print("Matching by index")
    matching = matching_index
    adata_st[:, matching].write_h5ad("${prefix}/adata_matched.h5ad", compression='gzip')
    print(f"Saved adata_st with {len(matching)} matching genes")
elif which_matching == 1:
    print("Matching by gene_ids")
    matching = matching_gene_ids
    adata_st.var.reset_index(drop=False, inplace=True)
    adata_st.var.set_index("gene_ids", inplace=True)
    adata_st.var.index = adata_st.var.index.astype('object')
    adata_st[:, matching].write_h5ad("${prefix}/adata_matched.h5ad", compression='gzip')
    print(f"Saved adata_st with {len(matching)} matching genes")
elif which_matching == 2:
    print("Matching by feature_name")
    matching = matching_feature_names
    m = {value: key for key, value in zip(adata_sc.var.index, adata_sc.var["feature_name"])}
    adata_st.var["name_matched"] = adata_st.var.index.map(m)
    adata_st.var.dropna(subset=["name_matched"], inplace=True)
    adata_st.var.reset_index(drop=False, inplace=True)
    adata_st.var.set_index("name_matched", inplace=True)
    adata_st.var.index = adata_st.var.index.astype('object')
    adata_st[:, adata_st.var.index].write_h5ad("${prefix}/adata_matched.h5ad", compression='gzip')
    print(f"Saved adata_st with {len(matching)} matching genes")
else:
    raise RuntimeError("More cases than expected")

adata_sc.file.close()
adata_st.file.close()
"""
}

process ATLAS_GET {
    //download an atlas anndata file from a url
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        val(url)
    output:
        path("*.h5ad"),         emit: atlas
        path("versions.yml"),   emit: versions

    script:
    prefix = task.ext.prefix
"""
#!/usr/bin/env python3
import os
import requests
from urllib.parse import urlparse
import boto3

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    requests: {}\\n".format(requests.__version__))
    f.write("    boto3: {}\\n".format(boto3.__version__))

myurl = "${url}"

if not(myurl.endswith(".h5ad")):
    raise ValueError("URL must end with .h5ad")

parsed_url = urlparse(myurl)
file_key = parsed_url.path.lstrip('/')

if myurl.startswith("s3://"):
    print("Downloading from S3")
    bucket_name = parsed_url.netloc
    s3 = boto3.client('s3')
    s3.download_file(bucket_name, file_key, os.path.basename(file_key))
elif myurl.startswith("https://"):
    print("Downloading from https")
    r = requests.get(myurl)
    r.raise_for_status()
    with open(os.path.basename(file_key), "wb") as f:
        f.write(r.content)
elif myurl.startswith("http://"):
    raise ValueError("Insecure HTTP URLs are not allowed. Please use HTTPS for remote atlas files.")
else:
    print("Local file path specified")
    if not os.path.isfile(myurl):
        raise FileNotFoundError(f"File {myurl} not found")
    os.symlink(myurl, os.path.basename(myurl))

print(f"Got atlas from {myurl}")
"""
}

process QC {
    //generate a simple report of the atlas adata
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(adata), val(report_name)
    output:
        path("*report.csv"),                     emit: report,   optional: true
        path("versions.yml"),                    emit: versions, optional: true

    script:
"""
#!/usr/bin/env python3
import os
import anndata as ad
import pandas as pd
import scanpy as sc

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    pandas: {}\\n".format(pd.__version__))

adata_path = "$adata"
outname = "$report_name"
sample = "${meta.id}"

if outname in ["atlas_input", "adata_input", "adata_output"]:
    adata = ad.read_h5ad(adata_path)
    #basic scanpy qc metrics
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    qc = sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, inplace=False)
    report = pd.DataFrame({
        "Sample": [sample],
        "n_genes": adata.shape[1],
        "n_cells": adata.shape[0],
        "mean_genes_by_counts": qc[0]["n_genes_by_counts"].mean(),
        "mean_cells_by_counts": qc[1]["n_cells_by_counts"].mean(),
        "mean_total_nnz_counts": adata.X[adata.X.nonzero()].mean(),
        "mean_percent_mito": qc[0]["pct_counts_mito"].mean()
    })
    report.to_csv(f"{outname}_report.csv", index=False)
    adata.file.close()

if outname in ["atlas_counts", "adata_counts"]:
    #cell type statistics
    adata = ad.read_h5ad(adata_path, backed='r')
    if "cell_type" in adata.obs.columns:
        cell_type_col = "cell_type"
    else:
        cell_type_col = "${params.ref_scrna_type_col}"
    if cell_type_col in adata.obs.columns:
        ct_counts = adata.obs[cell_type_col].value_counts(dropna=False)
        ct_counts.columns = ["cell_type", "n_cells"]
        ct_counts = pd.DataFrame(ct_counts).transpose()
        ct_counts.insert(0, "Sample", sample)
        pd.DataFrame(ct_counts).to_csv(f"{outname}_report.csv", index=False)
    else:
        print(f"Cell type column {cell_type_col} not found in adata.obs")
    adata.file.close()

if outname in ["cell_probs"]:
    #cell type qc metrics
    adata = ad.read_h5ad(adata_path, backed='r')
    if "cell_type_prob" in adata.obs.columns:
        mean_probs = adata.obs.groupby('cell_type').agg({'cell_type_prob':'mean'}).transpose()
        mean_probs.insert(0, "Sample", sample)
        mean_probs.to_csv(f"{outname}_report.csv", index=False)
    else:
        print("cell_type_prob not found in adata.obs")
    adata.file.close()
"""
}

process ADATA_FROM_VISIUM_HD {
    //convert vhd file to h5ad
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("${prefix}/adata.h5ad"),   emit: adata
        path("versions.yml"),                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))

sample = "${prefix}"
data = "${data}"
table = "${params.visium_hd}"
os.makedirs(sample, exist_ok=True)

#read visium_hd dataset
ds = sd.visium_hd(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
adata = to_legacy_anndata(ds, coordinate_system=sample,
                          table_name=table, include_images=True)
adata.var_names_make_unique()

#make compatible with BayesTME (uses an older, scanpy notation)
adata.X = adata.X.astype(int)
adata.uns['layout'] = 'IRREGULAR'
sq.gr.spatial_neighbors(adata)
adata.obsp['connectivities'] = adata.obsp['spatial_connectivities'].astype(bool)

#save
outname = os.path.join(sample, "adata.h5ad")
adata.write_h5ad(filename=outname, compression='gzip')
"""
}

process ADATA_FROM_VISIUM {
    //convert visium dir to h5ad
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("${prefix}/adata.h5ad"),          emit: adata
        path("versions.yml"),                                   emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq

sample = "${prefix}"
data = "${data}"

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))

os.makedirs(sample, exist_ok=True)

#read visium dataset
ds = sd.visium(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
adata = to_legacy_anndata(ds, coordinate_system=sample,
                          include_images=True)
adata.var_names_make_unique()

#make compatible with BayesTME (uses an older, scanpy notation)
adata.X = adata.X.astype(int)
adata.uns['layout'] = 'IRREGULAR'
sq.gr.spatial_neighbors(adata)
adata.obsp['connectivities'] = adata.obsp['spatial_connectivities'].astype(bool)

#save
outname = os.path.join(sample, "adata.h5ad")
adata.write_h5ad(filename=outname, compression='gzip')
"""
}

process ADATA_FROM_SEGMENTED_VISIUM {
    //convert visium dir to h5ad with cells instead of spots 
    //and cell centers as coordinates
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("${prefix}/adata.h5ad"),   emit: adata
        path("versions.yml"),                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'adata_from_segmented_visium.py'
}

process ADATA_FROM_XENIUM {
    //convert xenium dir to h5ad
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("${prefix}/adata.h5ad"),   emit: adata
        path("versions.yml"),                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))

sample = "${prefix}"
data = "${data}"

os.makedirs(sample, exist_ok=True)

#read xenium dataset
ds = sd.xenium(data)

#convert to anndata
adata = to_legacy_anndata(ds, include_images=True)
adata.var_names_make_unique()

#save
outname = os.path.join(sample, "adata.h5ad")
adata.write_h5ad(filename=outname, compression='gzip')
    """
}

process ADATA_PREPROCESS {
    //filter genes from an adata file by dropping genes whose names match a specified prefix (via drop_genes_prefix)
    tag "$meta.id"
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(adata)
    output:
        tuple val(meta), path("${prefix}/adata.h5ad"),   emit: adata
        path("versions.yml"),                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'adata_preprocess.py'
}

process ATTACH_CELL_PROBS {
    //attach cell type probabilities to anndata obsm
    tag "$meta.id"
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(cell_probs), path(adata), val(out_name)
    output:
        tuple val(meta), path("${prefix}/${out_name}.h5ad"),       emit: adata
        path("versions.yml"),                                      emit: versions

    script:
    sample = "${meta.id}"
    prefix = task.ext.prefix ?: "${sample}"
    template 'attach_cell_probs.py'
}

process CELL_TYPES_FROM_COGAPS {
    //extract cell types from a cogaps object
    tag "$meta.id"
    label "process_low"
    container "ghcr.io/fertiglab/cogaps@sha256:15dc4d443d927a7876b0b0f18291055fe0b3be63f1f040c71db8b6002b73e5de"

    input:
        tuple val(meta), path(cogaps_obj)
    output:
        tuple val(meta), path("${prefix}/cogaps_cell_types.csv"), emit: cogaps_cell_types
        path("versions.yml"),                                     emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'cell_types_from_cogaps.r'
}

 process ADATA_ADD_METADATA {
    //add metadata columns to anndata obs from samplesheet
    //this overwrites previously created andata.h5ad files
    tag "$meta.id"
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(adata)

    output:
        tuple val(meta), path("${prefix}/adata.h5ad"), emit: adata
        path("versions.yml"),                          emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'attach_metadata.py'
}

process STAPLE_ATTACH_LIGREC {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520'

    input:
        tuple val(meta), path(adata), path(ligrec)
    output:
        tuple val(meta), path("${prefix}/staple.h5ad"), emit: adata
        path "versions.yml",                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'attach_ligrec.py'
}


process SPACEMARKERS_CHECK {
    tag "$meta.id"
    label 'process_single'
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:ed44380c6e6e73fc575b743eba864941c26880053e50c3d70b5c9bfc526c0520"

    input:
        tuple val(meta), path(adata)

    output:
        tuple val(meta), path(adata), stdout, emit: checked
        path "versions.yml",                  emit: versions

    script:
    template 'spacemarkers_check.py'
}


process SPACEMARKERS_HARMONIZE {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/deshpandelab/spacemarkers@sha256:e13854a27622a04293fd8c26e8829a0407ab08a91d06259ece02eb440eab9ae2'


    input:
        tuple val(meta), path(scores), val(source)

    output:
        tuple val(meta), path("${prefix}/${source}/${scores}.gz") , emit: spacemarkers 

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'harmonize_spacemarkers.r'
}