#!/usr/bin/env python3
"""
Apply harmonized cell type names to an adata object.

Input: adata file, harmonization mapping CSV
Output: adata with harmonized cell type names in .obs
"""

import pandas as pd
import anndata as ad
from pathlib import Path
import numpy as np

# ============================================================================
# Parse inputs
# ============================================================================

adata_path = "${adata}"
mapping_path = "${harmonization_mapping}"
sample_id = """${meta.id}"""
source = """${source}"""

print(f"[AttachHarmonized] Processing sample: {sample_id}")
print(f"  Method: {source}")

# ============================================================================
# Load data
# ============================================================================

print(f"[AttachHarmonized] Loading adata from: {adata_path}")
adata = ad.read_h5ad(adata_path)
print(f"  Shape: {adata.shape}")

print(f"[AttachHarmonized] Loading mapping from: {mapping_path}")
mapping = pd.read_csv(mapping_path)
print(f"  Mapping rows: {len(mapping)}")

# ============================================================================
# Filter mapping for this sample
# ============================================================================

# The mapping CSV has 'original_cell_type' column with format "sample_id.pattern_name"
# Extract mappings for this sample based on the prefix
sample_mapping = mapping[mapping['original_cell_type'].str.startswith(f"{sample_id}.")]
print(f"[AttachHarmonized] Entries for {sample_id}: {len(sample_mapping)}")

if len(sample_mapping) == 0:
    print(f"[AttachHarmonized] WARNING: No mapping found for sample {sample_id}")
    # Extract sample_ids from the mapping
    unique_samples = set()
    for orig_type in mapping['original_cell_type']:
        if '.' in orig_type:
            unique_samples.add(orig_type.split('.')[0])
    print(f"  Available samples in mapping: {sorted(unique_samples)}")
    # Continue anyway, just won't update labels

# Create mapping dictionaries
# For the adata, we need to map from just the pattern name (without sample prefix)
# to the harmonized type
orig_to_harmonized = {}
confidence_dict = {}

for _, row in sample_mapping.iterrows():
    # Extract pattern name without sample prefix
    orig_full = row['original_cell_type']  # e.g., "sample1.P1"
    pattern_name = orig_full.split('.', 1)[1] if '.' in orig_full else orig_full
    
    orig_to_harmonized[pattern_name] = row['harmonized_cell_type']
    # Get confidence if it exists
    if 'confidence' in row and pd.notna(row['confidence']):
        confidence_dict[pattern_name] = row['confidence']
    elif 'correlation_to_consensus' in row and pd.notna(row['correlation_to_consensus']):
        confidence_dict[pattern_name] = row['correlation_to_consensus']

print(f"[AttachHarmonized] Mapping created with {len(orig_to_harmonized)} entries")

# ============================================================================
# Apply mapping to adata
# ============================================================================

cell_type_col = f'{source}_cell_type'
print(f"[AttachHarmonized] Looking for column: {cell_type_col}")

if cell_type_col in adata.obs.columns:
    print(f"[AttachHarmonized] Found {cell_type_col} in adata.obs")
    
    # Backup original names
    adata.obs[f'{cell_type_col}_original'] = adata.obs[cell_type_col].copy()
    print(f"[AttachHarmonized] Backed up original to: {cell_type_col}_original")
    
    # Apply harmonized mapping
    harmonized = adata.obs[cell_type_col].map(orig_to_harmonized)
    
    # Count unmapped values
    n_unmapped = harmonized.isna().sum()
    if n_unmapped > 0:
        print(f"[AttachHarmonized] WARNING: {n_unmapped} cells with unmapped cell types")
        print(f"  Original values: {adata.obs[cell_type_col][harmonized.isna()].unique()}")
        # Keep original names for unmapped cells
        harmonized = harmonized.fillna(adata.obs[cell_type_col])
    
    adata.obs[cell_type_col] = harmonized
    print(f"[AttachHarmonized] Updated {cell_type_col} with harmonized names")
    print(f"  Unique harmonized types: {adata.obs[cell_type_col].nunique()}")
    print(f"  Values: {sorted(adata.obs[cell_type_col].unique())}")
    
    # Add confidence scores
    adata.obs[f'{cell_type_col}_confidence'] = adata.obs[f'{cell_type_col}_original'].map(confidence_dict)
    n_conf = (~adata.obs[f'{cell_type_col}_confidence'].isna()).sum()
    print(f"[AttachHarmonized] Added confidence scores to {n_conf} cells")
    
else:
    print(f"[AttachHarmonized] WARNING: Column {cell_type_col} not found in adata.obs")
    print(f"  Available columns: {list(adata.obs.columns)}")

# ============================================================================
# Save output
# ============================================================================

output_dir = Path('${prefix}')
output_dir.mkdir(parents=True, exist_ok=True)

output_file = output_dir / 'adata_harmonized.h5ad'
print(f"[AttachHarmonized] Saving to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')
print(f"[AttachHarmonized] Saved successfully")

# ============================================================================
# Save versions
# ============================================================================

with open('versions.yml', 'w') as f:
    f.write("${task.process}:\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
    f.write(f"    anndata: {ad.__version__}\\n")

print("[AttachHarmonized] Saved: versions.yml")
print("=" * 80)
print("ATTACH HARMONIZED LABELS COMPLETE")
print("=" * 80)
