#!/usr/bin/env python
import base64
import csv
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from importlib.metadata import version
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

if __name__ == "__main__":
    cut_height = float('${cut_height}')
    min_samples = int('${min_samples}')
    confidence_threshold = float('${confidence_threshold}')

# load patterns
matrices = {}
all_patterns = []
sample_pattern_ranges = {}
pattern_names = []
pattern_sample_ids = []
pattern_columns = []

for entry in json.loads(base64.b64decode('${manifest}')):
    sample_id = entry['sample_id']
    filepath = Path(entry['path'])
    if not isinstance(sample_id, str) or not sample_id:
        raise ValueError(f"{filepath}: sample ID must be a non-empty string")
    if sample_id in matrices:
        raise ValueError(f"duplicate sample ID: {sample_id}")

    with filepath.open(newline='') as f:
        header = next(csv.reader(f))
    if len(header) < 2 or header[0] != '':
        raise ValueError(f"{filepath}: expected an unnamed feature index and at least one pattern column")
    if any(not name for name in header[1:]) or len(set(header[1:])) != len(header) - 1:
        raise ValueError(f"{filepath}: pattern names must be non-empty and unique")

    df = pd.read_csv(filepath, header=0, names=range(len(header)), dtype=str, keep_default_na=False)
    df = df.set_index(0)
    df.index.name = None
    df.columns = header[1:]
    if not df.index.is_unique or (df.index == '').any():
        raise ValueError(f"{filepath}: feature IDs must be non-empty and unique")
    df = df.apply(pd.to_numeric, errors='raise')
    if not np.isfinite(df.to_numpy()).all():
        raise ValueError(f"{filepath}: pattern values must be finite")
    matrices[sample_id] = df
    log.info(f"  {sample_id}: {df.shape[0]} features x {df.shape[1]} patterns")

if not matrices:
    raise ValueError("no pattern matrices provided")

sample_ids = sorted(matrices)
features = matrices[sample_ids[0]].index
for sample_id in sample_ids[1:]:
    features = features.intersection(matrices[sample_id].index, sort=False)
if len(features) < 2:
    raise ValueError(f"at least two shared features are required, found {len(features)}")
log.info(f"using {len(features)} shared features across {len(sample_ids)} samples")

for sample_id in sample_ids:
    df = matrices[sample_id].loc[features]

    start_idx = len(all_patterns)
    for col_name in df.columns:
        all_patterns.append(df[col_name].values)
        pattern_names.append(f"{sample_id}.{col_name}")
        pattern_sample_ids.append(sample_id)
        pattern_columns.append(col_name)
    end_idx = len(all_patterns)

    sample_pattern_ranges[sample_id] = (start_idx, end_idx)

log.info(f"total patterns: {len(pattern_names)}")
if len(pattern_names) < 2:
    raise ValueError("at least two patterns are required for clustering")

# compute distances
log.info("computing pattern correlations")
all_patterns_array = np.array(all_patterns).T
n_patterns = all_patterns_array.shape[1]

correlations = np.corrcoef(all_patterns_array.T)
correlations = np.nan_to_num(correlations, nan=0.0)
np.fill_diagonal(correlations, 1.0)

distances = 1 - correlations
distances = np.clip(distances, 0, 2)
distances = (distances + distances.T) / 2  # ensure symmetry

# hierarchical clustering
log.info("performing hierarchical clustering")
condensed_distances = squareform(distances)
Z = linkage(condensed_distances, method='complete')
clusters = fcluster(Z, cut_height, criterion='distance')

n_clusters = len(np.unique(clusters))
log.info(f"found {n_clusters} clusters at cut_height={cut_height}")

# cluster composition
cluster_composition = {}
for cluster_id in np.unique(clusters):
    cluster_composition[int(cluster_id)] = {sample_id: 0 for sample_id in sample_ids}

    for pattern_idx, assigned_cluster in enumerate(clusters):
        if assigned_cluster == cluster_id:
            for sample_id, (start_idx, end_idx) in sample_pattern_ranges.items():
                if start_idx <= pattern_idx < end_idx:
                    cluster_composition[cluster_id][sample_id] += 1
                    break

valid_clusters = []
for cluster_id, composition in cluster_composition.items():
    n_samples_with_pattern = sum(1 for count in composition.values() if count > 0)
    if n_samples_with_pattern >= min_samples:
        valid_clusters.append(cluster_id)

log.info(f"after filtering: {len(valid_clusters)} valid clusters")
consensus_names = {
    cluster_id: f"C{index:03d}"
    for index, cluster_id in enumerate(sorted(valid_clusters), start=1)
}

# generate mappings
log.info("generating harmonized cell type names")
mappings = []
mapped_pattern_indices = set()

# First, add mappings for valid clusters only
for cluster_id in sorted(valid_clusters):
    cluster_patterns = [i for i, c in enumerate(clusters) if c == cluster_id]
    cluster_pattern_names = [pattern_names[i] for i in cluster_patterns]

    cluster_pattern_vectors = [all_patterns[i] for i in cluster_patterns]
    consensus = np.mean(cluster_pattern_vectors, axis=0)

    correlations_to_consensus = [
        np.corrcoef(pattern_vec, consensus)[0, 1]
        for pattern_vec in cluster_pattern_vectors
    ]

    confidences = [
        min(1.0, max(0.0, corr))
        for corr in correlations_to_consensus
    ]

    composition = cluster_composition[cluster_id]
    composition_str = ','.join(
        [f"{sid}:{composition[sid]}" for sid in sample_ids]
    )

    harmonized_name = consensus_names[cluster_id]

    for pattern_idx, (pattern_name, confidence) in enumerate(zip(cluster_pattern_names, confidences)):
        if confidence >= confidence_threshold:
            mappings.append({
                'original_cell_type': pattern_name,
                'harmonized_cell_type': harmonized_name,
                'cluster': cluster_id,
                'confidence': confidence,
                'composition': composition_str,
                'sample_id': pattern_sample_ids[cluster_patterns[pattern_idx]],
                'original_pattern': pattern_columns[cluster_patterns[pattern_idx]]
            })
        
        # Mark this pattern as mapped (belongs to a valid cluster)
        mapped_pattern_indices.add(cluster_patterns[pattern_idx])

# Add patterns from invalid clusters as unique types (not meeting min_samples threshold)
next_unique_id = 1000
all_cluster_ids = np.unique(clusters)
invalid_clusters = [cid for cid in all_cluster_ids if int(cid) not in valid_clusters]

for cluster_id in sorted(invalid_clusters):
    cluster_patterns = [i for i, c in enumerate(clusters) if c == cluster_id]
    
    for pattern_idx in cluster_patterns:
        if pattern_idx not in mapped_pattern_indices:
            pattern_name = pattern_names[pattern_idx]
            sample_name = pattern_sample_ids[pattern_idx]
            
            mappings.append({
                'original_cell_type': pattern_name,
                'harmonized_cell_type': f"U{next_unique_id:03d}",
                'cluster': None,
                'confidence': 1.0,
                'composition': f"{sample_name}:1",
                'sample_id': sample_name,
                'original_pattern': pattern_columns[pattern_idx]
            })
            next_unique_id += 1
            mapped_pattern_indices.add(pattern_idx)

n_unique = len([m for m in mappings if m['cluster'] is None])
if n_unique > 0:
    log.info(f"added {n_unique} unique pattern mappings from invalid clusters")

mapping_df = pd.DataFrame(mappings, columns=[
    'original_cell_type', 'harmonized_cell_type', 'cluster', 'confidence',
    'composition', 'sample_id', 'original_pattern'
])
mapping_df.to_csv('cell_type_mapping.csv', index=False)
log.info(f"mapped {len(mapping_df)} patterns")

# consensus patterns
log.info("computing consensus patterns")
cluster_consensus = {}
for cluster_id in sorted(valid_clusters):
    cluster_patterns = [i for i, c in enumerate(clusters) if c == cluster_id]
    cluster_pattern_vectors = [all_patterns[i] for i in cluster_patterns]
    consensus = np.mean(cluster_pattern_vectors, axis=0)
    cluster_consensus[cluster_id] = consensus

consensus_cols = []
consensus_data = []

# Add consensus for valid clusters
for cluster_id in sorted(valid_clusters):
    consensus_cols.append(consensus_names[cluster_id])
    consensus_data.append(cluster_consensus[cluster_id])

# Add individual patterns from invalid clusters as unique consensus types
next_unique_id = 1000
all_cluster_ids = np.unique(clusters)
invalid_clusters = [cid for cid in all_cluster_ids if int(cid) not in valid_clusters]

for cluster_id in sorted(invalid_clusters):
    cluster_patterns = [i for i, c in enumerate(clusters) if c == cluster_id]
    for pattern_idx in cluster_patterns:
        consensus_cols.append(f"U{next_unique_id:03d}")
        consensus_data.append(all_patterns[pattern_idx])
        next_unique_id += 1

if consensus_data:
    consensus_df = pd.DataFrame(
        np.column_stack(consensus_data),
        columns=consensus_cols,
        index=features
    )
else:
    # Create empty dataframe when no valid clusters
    consensus_df = pd.DataFrame()
    
consensus_df.to_csv('consensus_patterns.csv')
if len(consensus_df) > 0:
    log.info(f"saved {len(consensus_df.columns)} consensus patterns")
else:
    log.info("no consensus patterns (no valid clusters)")

# diagnostics
log.info("creating diagnostics")
diagnostics = {
    'n_samples': len(sample_ids),
    'n_total_patterns': n_patterns,
    'n_clusters': n_clusters,
    'n_valid_clusters': len(valid_clusters),
    'cut_height': cut_height,
    'min_samples_per_cluster': min_samples,
    'confidence_threshold': confidence_threshold,
    'cluster_composition': cluster_composition,
    'algorithm': 'hierarchical_complete_linkage'
}

with open('harmonization_diagnostics.json', 'w') as f:
    json.dump(diagnostics, f, indent=2, default=str)

# dendrogram
log.info("plotting dendrogram")
fig, ax = plt.subplots(figsize=(14, 8))
dendrogram(Z, labels=pattern_names, ax=ax, leaf_font_size=8)
ax.axhline(y=cut_height, color='r', linestyle='--', label=f'cut_height={cut_height}')
ax.set_ylabel('Distance')
ax.set_title('Hierarchical Clustering')
ax.legend()
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('dendrogram.png', dpi=150, bbox_inches='tight')
plt.close()

# versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    pandas: {}\\n".format(version("pandas")))
    f.write("    numpy: {}\\n".format(version("numpy")))
    f.write("    scipy: {}\\n".format(version("scipy")))
    f.write("    matplotlib: {}\\n".format(version("matplotlib")))
