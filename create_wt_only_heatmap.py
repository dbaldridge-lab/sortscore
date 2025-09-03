#!/usr/bin/env python3
"""
Creates a WT-only heatmap showing wild-type sequence activity.

Usage:
    python create_wt_only_heatmap.py /path/to/experiment/directory

The script automatically finds the config.json and latest scores file in the directory.
It uses the synonymous WT score and creates a heatmap showing WT activity across positions.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from pathlib import Path
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.visualization.heatmap_matrix import make_dms_matrix


def create_wt_only_heatmap(experiment_dir: str):
    """Create a WT-only heatmap from an experiment directory."""
    
    exp_path = Path(experiment_dir)
    
    # Find config file
    config_file = exp_path / 'config.json'
    if not config_file.exists():
        raise FileNotFoundError(f"No config.json found in {exp_path}")
    
    # Load experiment configuration
    config = ExperimentConfig.from_json(str(config_file))
    
    if config.variant_type != 'aa':
        raise ValueError("This script only supports amino acid variant data")
    
    # Find scores directory and files
    scores_dir = exp_path / 'scores'
    if not scores_dir.exists():
        raise FileNotFoundError(f"No scores directory found in {exp_path}")
    
    # Find and load scores file for the DMS matrix
    scores_files = list(scores_dir.glob('*aa_scores*.csv'))
    if not scores_files:
        raise FileNotFoundError(f"No AA scores CSV files found in {scores_dir}")
    
    scores_file = max(scores_files, key=lambda p: p.stat().st_mtime)
    scores_df = pd.read_csv(scores_file)
    print(f"Using scores file: {scores_file.name}")
    
    # Get WT score from stats file
    stats_files = list(scores_dir.glob('*stats*.json'))
    if not stats_files:
        raise FileNotFoundError(f"No stats JSON files found in {scores_dir}")
    
    # Use the most recent stats file
    stats_file = max(stats_files, key=lambda p: p.stat().st_mtime)
    
    import json
    with open(stats_file, 'r') as f:
        stats = json.load(f)
    
    # Get synonymous (WT) average score
    if 'syn_avg' not in stats:
        raise ValueError("Could not find 'syn_avg' key in stats file")
    
    wt_score = stats['syn_avg']
    print(f"Found WT score from stats file (syn_avg): {wt_score:.1f}")
    
    # Create DMS matrix to get WT positions
    dms_matrix = make_dms_matrix(
        scores_df, 
        'avgscore', 
        config.num_aa, 
        config.wt_seq,
        variant_type=config.variant_type,
        mutagenesis_variants=config.mutagenesis_variants,
        position_type=config.position_type
    )
    wt_mask = dms_matrix == 'WT'
    
    # Create a matrix showing only WT positions with the WT score
    wt_only_matrix = pd.DataFrame(
        np.nan, 
        index=dms_matrix.index, 
        columns=dms_matrix.columns
    )
    
    # Fill WT positions with the WT score
    wt_only_matrix[wt_mask] = wt_score
    
    # Use the package's plot_heatmap function instead of manual plotting
    from sortscore.visualization.heatmaps import plot_heatmap
    
    # Create minimal scores dataframe for the heatmap function with correct columns
    # We only need WT positions, so create dummy data that matches the WT mask
    dummy_scores = []
    for i in range(len(dms_matrix.index)):
        for j in range(len(dms_matrix.columns)):
            if wt_mask.iloc[i, j]:  # Only add WT positions
                aa_var = dms_matrix.index[i]
                pos = int(dms_matrix.columns[j])
                dummy_scores.append({
                    'aa_seq_diff': f'{aa_var}.{pos}.{aa_var}',  # WT format: same AA before and after
                    'annotate_aa': 'synonymous',
                    'avgscore': wt_score,
                    'aa_pos': pos,
                    'aa_var': aa_var
                })
    
    dummy_df = pd.DataFrame(dummy_scores)
    
    # Save to figures directory
    figures_dir = exp_path / 'figures'
    figures_dir.mkdir(exist_ok=True)
    
    output_file = figures_dir / f'{config.experiment_name}_wt_only_heatmap.png'
    
    # Use the package's heatmap function
    plot_heatmap(
        dummy_df,
        'avgscore',
        config,
        wt_score=wt_score,
        export_heatmap=True,
        output=str(figures_dir),
        title=f'{config.experiment_name} - Wild Type Sequence Only',
        transparent=True
    )
    
    print(f"WT-only heatmap saved to figures directory")
    print(f"WT score used: {wt_score:.1f}")
    print(f"WT positions found: {wt_mask.sum().sum()}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python create_wt_only_heatmap.py /path/to/experiment/directory")
        sys.exit(1)
    
    create_wt_only_heatmap(sys.argv[1])