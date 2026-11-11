#!/usr/bin/env python
"""Generate test pattern matrices for harmonization tests."""

import numpy as np
import pandas as pd
from pathlib import Path


def create_cross_sample_patterns(n_genes=50, seed=42):
    """Create patterns where sample1.Pattern_1 == sample2.Pattern_2."""
    np.random.seed(seed)
    
    shared_pattern = np.random.rand(n_genes)
    
    patterns = {
        'sample1': {
            'Pattern_1': shared_pattern,
            'Pattern_2': np.random.rand(n_genes) * 0.5,
            'Pattern_3': np.random.rand(n_genes) * 0.5,
        },
        'sample2': {
            'Pattern_1': np.random.rand(n_genes) * 0.5,
            'Pattern_2': shared_pattern,
            'Pattern_3': np.random.rand(n_genes) * 0.5,
        }
    }
    return patterns


def create_random_patterns(n_samples=2, n_patterns=3, n_genes=50, seed=42):
    """Create completely random uncorrelated patterns."""
    np.random.seed(seed)
    
    patterns = {}
    for sample_idx in range(n_samples):
        patterns[f'sample{sample_idx+1}'] = {
            f'Pattern_{p+1}': np.random.rand(n_genes)
            for p in range(n_patterns)
        }
    return patterns


def create_unequal_patterns(n_genes=50, seed=42):
    """Create 3 samples with different numbers of patterns."""
    np.random.seed(seed)
    
    shared_pattern1 = 0.5 + np.random.rand(n_genes) * 0.5
    shared_pattern2 = 0.3 + np.random.rand(n_genes) * 0.4
    
    patterns = {
        'sample1': {
            'Pattern_1': shared_pattern1,
            'Pattern_2': np.random.rand(n_genes) * 0.5,
        },
        'sample2': {
            'Pattern_1': np.random.rand(n_genes) * 0.5,
            'Pattern_2': shared_pattern1,
            'Pattern_3': np.random.rand(n_genes) * 0.5,
            'Pattern_4': shared_pattern2,
        },
        'sample3': {
            'Pattern_1': np.random.rand(n_genes) * 0.5,
            'Pattern_2': shared_pattern2,
        }
    }
    return patterns


def create_threshold_zero_patterns(n_genes=50, seed=42):
    """Create patterns with varying similarity levels."""
    np.random.seed(seed)
    
    similar1 = 0.7 + np.random.rand(n_genes) * 0.3
    similar2 = similar1 + (np.random.rand(n_genes) - 0.5) * 0.1
    
    patterns = {
        'sample1': {
            'Pattern_1': similar1,
            'Pattern_2': np.random.rand(n_genes),
        },
        'sample2': {
            'Pattern_1': np.random.rand(n_genes),
            'Pattern_2': similar2,
        }
    }
    return patterns


def patterns_to_csv(patterns_dict, n_genes=50):
    """Convert pattern dict to CSV format string."""
    header = "Gene," + ",".join(patterns_dict.keys())
    rows = []
    for gene_idx in range(n_genes):
        row = f"GENE_{gene_idx:04d}"
        for pattern_values in patterns_dict.values():
            row += f",{pattern_values[gene_idx]:.4f}"
        rows.append(row)
    return header + "\\n" + "\\n".join(rows)


def save_patterns(patterns_dict, output_dir):
    """Save patterns to CSV files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    files = {}
    for sample_name, patterns in sorted(patterns_dict.items()):
        csv_content = patterns_to_csv(patterns)
        csv_path = output_dir / f"{sample_name}_patterns.csv"
        csv_path.write_text(csv_content)
        files[sample_name] = str(csv_path)
    
    return files


if __name__ == "__main__":
    # Nextflow template variables
    test_type = "${test_type}"
    
    if test_type == "cross_sample":
        patterns = create_cross_sample_patterns()
    elif test_type == "random":
        patterns = create_random_patterns()
    elif test_type == "unequal":
        patterns = create_unequal_patterns()
    elif test_type == "threshold_zero":
        patterns = create_threshold_zero_patterns()
    else:
        raise ValueError(f"Unknown test type: {test_type}")
    
    save_patterns(patterns, "patterns")
