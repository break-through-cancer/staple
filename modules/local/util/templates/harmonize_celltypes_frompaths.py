#!/usr/bin/env python3
import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata as ad
import sys
import os
from importlib.metadata import version

sample_ids_str = '''${sample_ids}'''
if sample_ids_str.startswith('[') and sample_ids_str.endswith(']'):
    sample_ids_str = sample_ids_str[1:-1]
sample_ids = [x.strip() for x in sample_ids_str.split(',') if x.strip()]

file_paths_str = '''${file_paths_str}'''
file_paths = [x.strip() for x in file_paths_str.split(chr(10)) if x.strip()]

cut_height = float('${cut_height}')
min_samples = int('${min_samples}')
confidence_threshold = float('${confidence_threshold}')

print(f"[Harmonize] Loading {len(sample_ids)} samples from {len(file_paths)} files")

matrices = {}
all_pattern_labels = []

for sample_id, filepath in zip(sample_ids, file_paths):
    p = Path(filepath)
    if not p.exists():
        print(f"ERROR: File not found: {filepath}")
        sys.exit(1)
    
    df = pd.read_csv(filepath, index_col=0)
    matrices[sample_id] = df
    print(f"  {sample_id}: {df.shape[0]} cells x {df.shape[1]} patterns")
    
    for col in df.columns:
        all_pattern_labels.append(f"{sample_id}.{col}")

print(f"[Harmonize] Total patterns: {len(all_pattern_labels)}")
print("[Harmonize] Computing pattern correlations...")

all_patterns = []
sample_pattern_ranges = {}
pattern_names = []

for sample_id in sample_ids:
    start_idx = len(all_patterns)
    df = matrices[sample_id]
    
    for col_name in df.columns:
        all_patterns.append(df[col_name].values)
        pattern_names.append(f"{sample_id}.{col_name}")
    
    end_idx = len(all_patterns)
    sample_pattern_ranges[sample_id] = (start_idx, end_idx)

all_patterns_array = np.array(all_patterns).T
n_patterns = all_patterns_array.shape[1]
print(f"  Computing {n_patterns} patterns")

correlations = np.corrcoef(all_patterns_array.T)
correlations = np.nan_to_num(correlations, nan=0.0)
np.fill_diagonal(correlations, 1.0)

distances = 1 - correlations
distances = np.clip(distances, 0, 2)

# Ensure symmetry (fix numerical precision issues)
distances = (distances + distances.T) / 2

condensed_distances = squareform(distances)

print(f"[Harmonize] Performing hierarchical clustering...")
Z = linkage(condensed_distances, method='complete')
clusters = fcluster(Z, cut_height, criterion='distance')

n_clusters = len(np.unique(clusters))
print(f"  Found {n_clusters} clusters at cut_height={cut_height}")

pattern_to_cluster = {i: int(clusters[i]) for i in range(len(pattern_names))}

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

print(f"  After filtering: {len(valid_clusters)} valid clusters")

print(f"[Harmonize] Generating harmonized cell type names...")

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
    avg_correlation = np.nanmean(correlations_to_consensus)
    
    confidences = []
    for pattern_vec, corr_to_consensus in zip(cluster_pattern_vectors, correlations_to_consensus):
        confidence = min(1.0, max(0.0, corr_to_consensus))
        confidences.append(confidence)
    
    composition = cluster_composition[cluster_id]
    composition_str = ','.join([f"{sample_id}:{composition[sample_id]}" for sample_id in sample_ids])
    
    harmonized_name = f"C{cluster_id:03d}"
    
    for pattern_idx, pattern_name, confidence in zip(cluster_patterns, cluster_pattern_names, confidences):
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
print(f"  Mapped {len(mapping_df)} patterns")

print(f"[Harmonize] Computing consensus patterns...")

consensus_patterns = []
for cluster_id in sorted(valid_clusters):
    cluster_patterns = [i for i, c in enumerate(clusters) if c == cluster_id]
    cluster_pattern_vectors = [all_patterns[i] for i in cluster_patterns]
    consensus = np.mean(cluster_pattern_vectors, axis=0)
    
    consensus_patterns.append(consensus)

if consensus_patterns:
    consensus_df = pd.DataFrame(
        np.array(consensus_patterns).T,
        columns=[f"C{cluster_id:03d}" for cluster_id in sorted(valid_clusters)]
    )
    consensus_df.to_csv('consensus_patterns.csv')
    print(f"  Saved {len(consensus_df.columns)} consensus patterns")

print(f"[Harmonize] Creating diagnostics...")

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

print(f"  Saved diagnostics")

print(f"[Harmonize] Plotting dendrogram...")

fig, ax = plt.subplots(figsize=(14, 8))
dendrogram(Z, labels=pattern_names, ax=ax, leaf_font_size=8)
ax.axhline(y=cut_height, color='r', linestyle='--', label=f'cut_height={cut_height}')
ax.set_ylabel('Distance')
ax.set_title(f'Hierarchical Clustering')
ax.legend()
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('dendrogram.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"  Saved dendrogram plot")

versions_txt = f'pandas: {version("pandas")}{chr(10)}numpy: {version("numpy")}{chr(10)}scipy: {version("scipy")}{chr(10)}anndata: {version("anndata")}{chr(10)}matplotlib: {version("matplotlib")}{chr(10)}scikit-learn: {version("scikit-learn")}'

with open('versions.yml', 'w') as f:
    f.write(versions_txt)

print(f"[Harmonize] COMPLETE!")
