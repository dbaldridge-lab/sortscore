"""
Main entry point for Sort-seq variant analysis.

This script loads the experiment configuration, ensures output directories exist, and orchestrates the analysis workflow.

Usage:
    python -m sortscore.run_analysis --config path/to/your_config.json
"""
import argparse
import json
import logging
import sys
import os
import pandas as pd
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.utils import ensure_output_subdirs
from sortscore.analysis.score import calculate_full_activity_scores
from sortscore.visualization.plots import plot_heatmap


def main():
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    parser = argparse.ArgumentParser(description="Run Sort-seq variant analysis.")
    parser.add_argument('-c', '--config', type=str, required=True, help='Path to experiment config JSON file')
    args = parser.parse_args()

    # Load experiment config using dataclass
    try:
        experiment = ExperimentConfig.from_json(args.config)
    except Exception as e:
        logging.error(f"Failed to load config: {e}")
        sys.exit(1)

    output_dir = experiment.output_dir or '.'
    # Ensure output subdirectories exist
    try:
        ensure_output_subdirs(output_dir)
    except PermissionError as e:
        logging.error(f"Permission denied: cannot create output directories in '{output_dir}'. {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Failed to create output directories in '{output_dir}': {e}")
        sys.exit(1)
    logging.info(f"Output directories ensured in {output_dir}")

    # Load all counts using the dataclass method
    try:
        experiment.load_counts()
    except Exception as e:
        logging.error(f"Failed to load counts: {e}")
        sys.exit(1)
    logging.info(f"Loaded counts for {len(experiment.counts)} replicates.")
    # experiment.counts[rep][bin] gives the DataFrame for each replicate/bin

    print("Sort-seq analysis setup complete. Counts loaded for all replicates and bins.")
    
    # Calculate activity scores using dataclass parameters
    try:
        logging.info("Calculating activity scores...")
        scores_df = calculate_full_activity_scores(
            counts=experiment.counts,
            median_gfp=experiment.median_gfp,
            min_bins=experiment.bins_required,
            min_reps=experiment.reps_required,
            minread_threshold=experiment.minread_threshold,
            avg_method=experiment.avg_method
        )
        logging.info(f"Calculated scores for {len(scores_df)} variants.")
        logging.info(f"Score columns: {list(scores_df.columns)}")
    except Exception as e:
        logging.error(f"Failed to calculate activity scores: {e}")
        sys.exit(1)
    
    # Annotate sequences if needed
    try:
        experiment.annotate_counts(experiment.wt_seq)
        logging.info("Added sequence annotations.")
        
        # Add codon_diff to scores DataFrame for DNA variant types
        if experiment.variant_type == 'dna' and experiment.counts:
            from sortscore.sequence_parsing import compare_codon_lists
            scores_df['codon_diff'] = scores_df['variant_seq'].apply(
                lambda x: compare_codon_lists(experiment.wt_seq, x)
            )
            scores_df['codon_diff'] = scores_df['codon_diff'].fillna('')
            logging.info("Added codon_diff column to scores DataFrame.")
    except Exception as e:
        logging.warning(f"Failed to add annotations: {e}")
    
    # Save results
    timestamp = pd.Timestamp.now().strftime('%Y%m%d')
    
    # Save DNA scores (full data)
    dna_scores_file = os.path.join(output_dir, f"dna-scores_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.csv")
    scores_df.to_csv(dna_scores_file, index=False)
    logging.info(f"Saved DNA scores to {dna_scores_file}")
    
    # Save AA scores (simplified) if available
    if 'aa_seq_diff' in scores_df.columns:
        aa_cols = ['aa_seq_diff', 'annotate_aa'] + [col for col in scores_df.columns if 'score' in col.lower()]
        aa_scores = scores_df[aa_cols].copy()
        aa_scores_file = os.path.join(output_dir, f"aa-scores_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.csv")
        aa_scores.to_csv(aa_scores_file, index=False)
        logging.info(f"Saved AA scores to {aa_scores_file}")
    
    # Generate visualizations
    try:
        logging.info("Generating visualizations...")
        
        # Convert avg_method to column name format (replace hyphens with underscores)
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        
        # AA heatmap
        if 'aa_seq_diff' in scores_df.columns:
            aa_heatmap_file = os.path.join(output_dir, f"aa_heatmap_{experiment.avg_method}_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.png")
            plot_heatmap(scores_df, score_col, experiment, 
                        export=True, output=aa_heatmap_file)
            logging.info(f"Saved AA heatmap to {aa_heatmap_file}")
        
        # Codon heatmap  
        if experiment.variant_type == 'dna':
            codon_heatmap_file = os.path.join(output_dir, f"codon_heatmap_{experiment.avg_method}_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.png")
            plot_heatmap(scores_df, score_col, experiment,
                        export=True, output=codon_heatmap_file)
            logging.info(f"Saved codon heatmap to {codon_heatmap_file}")
            
    except Exception as e:
        logging.error(f"Failed to generate visualizations: {e}")
        sys.exit(1)
    
    # Save summary statistics
    try:
        stats = {}
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        if score_col in scores_df.columns:
            stats['all_avg'] = float(scores_df[score_col].mean())
            stats['all_min'] = float(scores_df[score_col].min())
            stats['all_max'] = float(scores_df[score_col].max())
            
            # Add annotation-based stats if available
            if 'annotate_aa' in scores_df.columns:
                for annot in ['synonymous', 'nonsense', 'missense_aa']:
                    subset = scores_df[scores_df['annotate_aa'] == annot]
                    if len(subset) > 0:
                        stats[f'{annot.replace("_aa", "")}_avg'] = float(subset[score_col].mean())
                        stats[f'{annot.replace("_aa", "")}_min'] = float(subset[score_col].min())
                        stats[f'{annot.replace("_aa", "")}_max'] = float(subset[score_col].max())
        
        stats_file = os.path.join(output_dir, f"stats_{experiment.avg_method}_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.json")
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        logging.info(f"Saved statistics to {stats_file}")
        
    except Exception as e:
        logging.warning(f"Failed to save statistics: {e}")
    
    print(f"Analysis complete! Results saved to {output_dir}")

if __name__ == "__main__":
    main()
