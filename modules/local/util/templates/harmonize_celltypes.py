#!/usr/bin/env python
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
    sample_ids_str = '''${sample_ids}'''
    if sample_ids_str.startswith('[') and sample_ids_str.endswith(']'):
        sample_ids_str = sample_ids_str[1:-1]
    sample_ids = [x.strip() for x in sample_ids_str.split(',') if x.strip()]

    file_paths_str = '''${file_paths_str}'''
    file_paths = [x.strip() for x in file_paths_str.split(chr(10)) if x.strip()]

    cut_height = float('${cut_height}')
    min_samples = int('${min_samples}')
    confidence_threshold = float('${confidence_threshold}')

    log.info(f"loading {len(sample_ids)} samples from {len(file_paths)} files")

# load patterns
matrices = {}
all_patterns = []
sample_pattern_ranges = {}
pattern_names = []

for sample_id, filepath in zip(sample_ids, file_paths):
    p = Path(filepath)
    if not p.exists():
        log.error(f"file not found: {filepath}")
        exit(1)

    df = pd.read_csv(filepath, index_col=0)
    matrices[sample_id] = df
    log.info(f"  {sample_id}: {df.shape[0]} cells x {df.shape[1]} patterns")

    start_idx = len(all_patterns)
    for col_name in df.columns:
        all_patterns.append(df[col_name].values)
        pattern_names.append(f"{sample_id}.{col_name}")
    end_idx = len(all_patterns)

    sample_pattern_ranges[sample_id] = (start_idx, end_idx)

log.info(f"total patterns: {len(pattern_names)}")

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

# generate mappings
log.info("generating harmonized cell type names")
mappings = []

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

    harmonized_name = f"C{cluster_id:03d}"

    for pattern_name, confidence in zip(cluster_pattern_names, confidences):
        if confidence >= confidence_threshold:
            mappings.append({
                'original_cell_type': pattern_name,
                'harmonized_cell_type': harmonized_name,
                'cluster': cluster_id,
                'confidence': confidence,
                'composition': composition_str
            })

mapping_df = pd.DataFrame(mappings)
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

if cluster_consensus:
    consensus_df = pd.DataFrame(
        np.column_stack([cluster_consensus[cid] for cid in sorted(valid_clusters)]),
        columns=[f"C{cid:03d}" for cid in sorted(valid_clusters)]
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

