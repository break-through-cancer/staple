#!/usr/bin/env python
import logging
import os
from pathlib import Path

import pandas as pd
import anndata as ad
from importlib.metadata import version

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

if __name__ == "__main__":
    adata_path = "${adata}"
    mapping_path = "${harmonization_mapping}"
    sample_id = """${meta.id}"""
    source = """${source}"""
    prefix = "${prefix}"

    # output directory
    os.makedirs(prefix, exist_ok=True)

    log.info(f"attaching harmonized cell type names to {sample_id}")
    log.info(f"method: {source}")

    # load data
    log.info(f"loading adata from {adata_path}")
    adata = ad.read_h5ad(adata_path)

    log.info(f"loading mapping from {mapping_path}")
    mapping = pd.read_csv(mapping_path)
    log.info(f"mapping contains {len(mapping)} entries")

    # extract sample-specific mapping
    sample_mapping = mapping[mapping['original_cell_type'].str.startswith(f"{sample_id}.")]
    log.info(f"found {len(sample_mapping)} entries for {sample_id}")

    if len(sample_mapping) == 0:
        log.warning(f"no mapping found for {sample_id}")
        unique_samples = set(
            orig.split('.')[0] for orig in mapping['original_cell_type']
            if '.' in orig
        )
        if unique_samples:
            log.warning(f"available samples: {sorted(unique_samples)}")

    # build mapping dictionaries
    orig_to_harmonized = {}
    confidence_dict = {}

    for _, row in sample_mapping.iterrows():
        orig_full = row['original_cell_type']
        pattern_name = orig_full.split('.', 1)[1] if '.' in orig_full else orig_full

        orig_to_harmonized[pattern_name] = row['harmonized_cell_type']

        if 'confidence' in row and pd.notna(row['confidence']):
            confidence_dict[pattern_name] = row['confidence']
        elif 'correlation_to_consensus' in row and pd.notna(row['correlation_to_consensus']):
            confidence_dict[pattern_name] = row['correlation_to_consensus']

    # apply mapping
    # The column is always named 'cell_type' regardless of the deconvolution method
    cell_type_col = 'cell_type'
    log.info(f"looking for column: {cell_type_col}")

    if cell_type_col not in adata.obs.columns:
        log.warning(f"column {cell_type_col} not found in adata.obs")
        log.warning(f"available columns: {list(adata.obs.columns)}")
    else:
        # backup original names (convert to string to preserve values)
        adata.obs[f'{cell_type_col}_original'] = adata.obs[cell_type_col].astype(str)

        # convert cell_type column to string to allow new harmonized values
        original_values = adata.obs[cell_type_col].astype(str)

        # apply harmonized mapping
        harmonized = original_values.map(orig_to_harmonized)

        n_unmapped = harmonized.isna().sum()
        if n_unmapped > 0:
            log.warning(f"{n_unmapped} cells with unmapped cell types")
            # Keep original values for unmapped cells
            harmonized = harmonized.fillna(original_values)

        adata.obs[cell_type_col] = harmonized
        log.info(f"updated {cell_type_col} with {adata.obs[cell_type_col].nunique()} harmonized types")

        # add confidence scores
        adata.obs[f'{cell_type_col}_confidence'] = (
            adata.obs[f'{cell_type_col}_original'].map(confidence_dict)
        )
        n_conf = (~adata.obs[f'{cell_type_col}_confidence'].isna()).sum()
        log.info(f"added confidence scores to {n_conf} cells")

    # harmonize cell_type_composition columns if present
    if 'cell_type_composition' in adata.uns:
        comp_df = adata.uns['cell_type_composition']
        log.info(f"harmonizing cell_type_composition columns: {list(comp_df.columns)}")
        
        # rename columns based on the mapping
        new_columns = {}
        for col in comp_df.columns:
            if col in orig_to_harmonized:
                new_columns[col] = orig_to_harmonized[col]
            else:
                log.warning(f"column {col} not in mapping, keeping original name")
                new_columns[col] = col
        
        adata.uns['cell_type_composition'] = comp_df.rename(columns=new_columns)
        log.info(f"updated cell_type_composition columns to: {list(adata.uns['cell_type_composition'].columns)}")

    # save output
    log.info(f"saving to {prefix}/staple.h5ad")
    adata.write_h5ad(f"{prefix}/staple.h5ad", compression='gzip')

    # versions
    with open("versions.yml", "w") as f:
        f.write("${task.process}:\\n")
        f.write("    anndata: {}\\n".format(version("anndata")))
        f.write("    pandas: {}\\n".format(version("pandas")))
