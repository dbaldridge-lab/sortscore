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
import numpy as np
from scipy import stats as scipy_stats
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.utils import ensure_output_subdirs
from sortscore.analysis.score import calculate_full_activity_scores
from sortscore.analysis.data_processing import aggregate_aa_data
from sortscore.analysis.annotation import annotate_scores_dataframe, aggregate_synonymous_variants
from sortscore.sequence_parsing import translate_dna
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
        
        # Get merged counts DataFrame
        merged_df = experiment.get_merged_counts()
        
        scores_df = calculate_full_activity_scores(
            counts=experiment.counts,
            median_gfp=experiment.median_gfp,
            min_bins=experiment.bins_required,
            min_reps=experiment.reps_required,
            minread_threshold=experiment.minread_threshold,
            avg_method=experiment.avg_method,
            total_reads=experiment.total_reads,
            merged_df=merged_df
        )
        logging.info(f"Calculated scores for {len(scores_df)} variants.")
        logging.info(f"Score columns: {list(scores_df.columns)}")
    except Exception as e:
        logging.error(f"Failed to calculate activity scores: {e}")
        sys.exit(1)
    
    # Annotate sequences
    try:
        # Annotate the scores DataFrame with sequence information
        scores_df = annotate_scores_dataframe(scores_df, experiment.wt_seq, experiment.variant_type)
        logging.info("Added sequence annotations to scores DataFrame.")
    except Exception as e:
        logging.warning(f"Failed to add annotations: {e}")
    
    # Save results
    timestamp = pd.Timestamp.now().strftime('%Y%m%d')
    scores_dir = os.path.join(output_dir, 'scores')
    
    # Save DNA scores
    scores_df_rounded = scores_df.copy()
    
    # Calculate standard deviation, coefficient of variation, and 95% CI of replicate scores
    rep_score_columns = [col for col in scores_df_rounded.columns if col.startswith('Rep') and col.endswith('.score')]
    if len(rep_score_columns) >= 2:
        rep_mean = scores_df_rounded[rep_score_columns].mean(axis=1)
        rep_std = scores_df_rounded[rep_score_columns].std(axis=1, ddof=1)
        
        # Calculate n_measurements dynamically based on non-empty replicate values
        # Each DNA row represents 1 codon, so n_measurements = number of non-empty replicates
        n_measurements = scores_df_rounded[rep_score_columns].notna().sum(axis=1)
        
        # Calculate 95% CI using t-distribution with actual degrees of freedom
        df_actual = n_measurements - 1
        t_critical = scipy_stats.t.ppf(0.975, df_actual)
        margin_of_error = t_critical * (rep_std / np.sqrt(n_measurements))
        
        scores_df_rounded['SD_rep'] = rep_std.round().astype('Int64')
        scores_df_rounded['CV%'] = (rep_std / rep_mean * 100).round().astype('Int64')
        scores_df_rounded['n_measurements'] = n_measurements.astype('Int64')
        scores_df_rounded['CI_lower'] = (rep_mean - margin_of_error).round().astype('Int64')
        scores_df_rounded['CI_upper'] = (rep_mean + margin_of_error).round().astype('Int64')
    
    # Round other score columns to integers
    score_columns = [col for col in scores_df_rounded.columns if 'score' in col.lower()]
    for col in score_columns:
        if scores_df_rounded[col].dtype in ['float64', 'float32']:
            scores_df_rounded[col] = scores_df_rounded[col].round().astype('Int64')
    
    # Save DNA scores for DNA variant type, skip for AA variant type
    if experiment.variant_type == 'dna':
        dna_scores_file = os.path.join(scores_dir, f"dna-scores_{experiment.experiment_name}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.csv")
        scores_df_rounded.to_csv(dna_scores_file, index=False)
        logging.info(f"Saved DNA scores to {dna_scores_file}")
    
    # Save AA scores 
    if 'aa_seq_diff' in scores_df.columns:
        # Filter out rows with NaN values first
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        scores_df_drop_nan = scores_df.dropna(subset=[score_col])
        
        # Find replicate score columns
        rep_score_columns = [col for col in scores_df_drop_nan.columns if col.startswith('Rep') and col.endswith('.score')]
        
        # Determine if we need to aggregate DNA->AA variants or if we have AA-only variants
        # Check if there are multiple codons per AA variant (DNA->AA case)
        aa_variant_counts = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa']).size()
        needs_aggregation = (aa_variant_counts > 1).any()
        
        if needs_aggregation:
            # DNA->AA aggregation case: aggregate synonymous variants
            columns_to_average = ['avgscore', 'avgscore_rep_weighted', 'avgscore_codon_weighted'] + rep_score_columns
            
            # Calculate standard deviation and count of codon-level scores before AA aggregation
            aa_scores_std = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[score_col].agg(['std', 'count']).reset_index()
            aa_scores_std.columns = ['aa_seq_diff', 'annotate_aa', 'SD_codon', 'n_codons']
            
            # Calculate mean scores for aggregation
            aa_scores = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[columns_to_average].mean().reset_index()
            
            # Merge the standard deviation and count of codon scores
            aa_scores = aa_scores.merge(aa_scores_std, on=['aa_seq_diff', 'annotate_aa'], how='left')
            
            # Calculate statistics with codon and replicate variance decomposition
            if len(rep_score_columns) >= 2:
                aa_rep_mean = aa_scores[rep_score_columns].mean(axis=1)
                aa_rep_std = aa_scores[rep_score_columns].std(axis=1, ddof=1)
                
                # Calculate variances
                var_codon = aa_scores['SD_codon'] ** 2
                var_rep = aa_rep_std ** 2
                
                # Calculate total variance (sum of codon and replicate variances)
                var_total = var_codon.fillna(0) + var_rep
                
                # Calculate n_measurements dynamically
                n_measurements = aa_scores['n_codons'] * aa_scores[rep_score_columns].notna().sum(axis=1)
                
                # Calculate SEM considering total sample size
                sem = np.sqrt(var_total / n_measurements)
                
                # Calculate 95% CI using t-distribution
                df_total = n_measurements - 1
                t_critical = scipy_stats.t.ppf(0.975, df_total)
                aa_margin_of_error = t_critical * sem
                
                aa_scores['SD_rep'] = aa_rep_std.round().astype('Int64')
                aa_scores['CV%_rep'] = (aa_rep_std / aa_rep_mean * 100).round().astype('Int64')
                aa_scores['CV%_codon'] = (aa_scores['SD_codon'] / aa_rep_mean * 100).round().astype('Int64')
                aa_scores['n_measurements'] = n_measurements.astype('Int64')
                aa_scores['SEM'] = sem.round().astype('Int64')
                aa_scores['CI_lower'] = (aa_rep_mean - aa_margin_of_error).round().astype('Int64')
                aa_scores['CI_upper'] = (aa_rep_mean + aa_margin_of_error).round().astype('Int64')
            
            # Round SD_codon column to integers
            if 'SD_codon' in aa_scores.columns:
                aa_scores['SD_codon'] = aa_scores['SD_codon'].round().astype('Int64')
                
        else:
            # AA-only case: no aggregation needed, just copy the data
            columns_to_include = ['aa_seq_diff', 'annotate_aa', 'avgscore', 'avgscore_rep_weighted'] + rep_score_columns
            # Only include avgscore_codon_weighted if it exists (it won't for AA-only)
            if 'avgscore_codon_weighted' in scores_df_drop_nan.columns:
                columns_to_include.append('avgscore_codon_weighted')
            
            aa_scores = scores_df_drop_nan[columns_to_include].copy()
            
            # Calculate simple replicate statistics (no codon variance)
            if len(rep_score_columns) >= 2:
                aa_rep_mean = aa_scores[rep_score_columns].mean(axis=1)
                aa_rep_std = aa_scores[rep_score_columns].std(axis=1, ddof=1)
                
                # Calculate n_measurements (just number of non-empty replicates)
                n_measurements = aa_scores[rep_score_columns].notna().sum(axis=1)
                
                # Calculate SEM using only replicate variance
                sem = aa_rep_std / np.sqrt(n_measurements)
                
                # Calculate 95% CI using t-distribution
                df_actual = n_measurements - 1
                t_critical = scipy_stats.t.ppf(0.975, df_actual)
                aa_margin_of_error = t_critical * sem
                
                aa_scores['SD_rep'] = aa_rep_std.round().astype('Int64')
                aa_scores['CV%_rep'] = (aa_rep_std / aa_rep_mean * 100).round().astype('Int64')
                aa_scores['n_measurements'] = n_measurements.astype('Int64')
                aa_scores['SEM'] = sem.round().astype('Int64')
                aa_scores['CI_lower'] = (aa_rep_mean - aa_margin_of_error).round().astype('Int64')
                aa_scores['CI_upper'] = (aa_rep_mean + aa_margin_of_error).round().astype('Int64')
        
        # Round score columns to integers
        score_columns = [col for col in aa_scores.columns if 'score' in col.lower()]
        for col in score_columns:
            if aa_scores[col].dtype in ['float64', 'float32']:
                aa_scores[col] = aa_scores[col].round().astype('Int64')
        
        aa_scores_file = os.path.join(scores_dir, f"aa-scores_{experiment.experiment_name}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.csv")
        aa_scores.to_csv(aa_scores_file, index=False)
        logging.info(f"Saved AA scores to {aa_scores_file} ({len(aa_scores)} unique AA variants)")
    
    # Generate visualizations
    try:
        logging.info("Generating visualizations...")
        
        # Convert avg_method to column name format (replace hyphens with underscores)
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        
        figures_dir = os.path.join(output_dir, 'figures')
        
        # AA heatmap
        if 'aa_seq_diff' in scores_df.columns:
            # Use AA data directly for AA-only experiments, or aggregate for DNA experiments
            if experiment.variant_type == 'aa':
                # For AA-only data, use the data directly
                aa_data = scores_df[['aa_seq_diff', 'annotate_aa', score_col]].copy()
            else:
                # For DNA data, aggregate to AA level
                aa_data = aggregate_aa_data(scores_df, score_col)
            
            # Get WT score from the aggregated AA data for AA heatmap
            wt_score = pd.NA  # default
            if 'annotate_aa' in aa_data.columns:
                if experiment.variant_type == 'dna':
                    # For DNA mode: average all wt_dna entries
                    wt_subset = aa_data[aa_data['annotate_aa'] == 'wt_dna']
                    if len(wt_subset) > 0 and score_col in wt_subset.columns:
                        score_val = wt_subset[score_col].mean()
                        wt_score = float(score_val) if pd.notna(score_val) else pd.NA
                        logging.info(f"Found WT score from DNA->AA data (averaged from {len(wt_subset)} wt_dna): {wt_score}")
                else:
                    # For AA mode: average all synonymous variants
                    wt_subset = aa_data[aa_data['annotate_aa'] == 'synonymous']
                    if len(wt_subset) > 0 and score_col in wt_subset.columns:
                        score_val = wt_subset[score_col].mean()
                        wt_score = float(score_val) if pd.notna(score_val) else pd.NA
                        logging.info(f"Found WT score from AA data (averaged from {len(wt_subset)} synonymous variants): {wt_score}")
            
            # Create a simple config object for AA heatmap
            from types import SimpleNamespace
            aa_config = SimpleNamespace(
                num_aa=experiment.num_aa,
                min_pos=experiment.min_pos,
                wt_seq=translate_dna(experiment.wt_seq),
                variant_type='aa'
            )
            
            # Set up colorbar with WT score if available
            tick_values = None
            tick_labels = None
            if not pd.isna(wt_score):
                # Get min/max values from data for colorbar range
                data_min = aa_data[score_col].min()
                data_max = aa_data[score_col].max()
                # Add WT score to tick marks if min/max are not NaN
                if pd.notna(data_min) and pd.notna(data_max):
                    tick_values = [data_min, wt_score, data_max]
                    tick_labels = [f'{data_min:.0f}', f'WT={wt_score:.0f}', f'{data_max:.0f}']
            
            aa_heatmap_file = os.path.join(figures_dir, f"aa_heatmap_{experiment.avg_method}_{experiment.experiment_name}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.png")
            plot_heatmap(aa_data, score_col, aa_config, wt_score=wt_score,
                        tick_values=tick_values, tick_labels=tick_labels,
                        export=True, output=aa_heatmap_file, format='png', export_matrix=True)
            logging.info(f"Saved AA heatmap to {aa_heatmap_file}")
        
        # Codon heatmap  
        if experiment.variant_type == 'dna':
            # Get WT score for codon heatmap tick marks
            wt_score_codon = pd.NA  # default
            if 'annotate_dna' in scores_df.columns:
                wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
                if len(wt_subset) > 0 and score_col in wt_subset.columns:
                    # Average all wt_dna entries for codon heatmap
                    score_val = wt_subset[score_col].mean()
                    wt_score_codon = float(score_val) if pd.notna(score_val) else pd.NA
                    logging.info(f"Found WT score for codon heatmap (averaged from {len(wt_subset)} wt_dna): {wt_score_codon}")
            
            # Set up colorbar with WT score if available
            tick_values_codon = None
            tick_labels_codon = None
            if not pd.isna(wt_score_codon):
                # Get min/max values from data for colorbar range
                data_min = scores_df[score_col].min()
                data_max = scores_df[score_col].max()
                # Add WT score to tick marks if min/max are not NaN
                if pd.notna(data_min) and pd.notna(data_max):
                    tick_values_codon = [data_min, wt_score_codon, data_max]
                    tick_labels_codon = [f'{data_min:.0f}', f'WT={wt_score_codon:.0f}', f'{data_max:.0f}']
            
            codon_heatmap_file = os.path.join(figures_dir, f"codon_heatmap_{experiment.avg_method}_{experiment.experiment_name}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.png")
            plot_heatmap(scores_df, score_col, experiment, wt_score=wt_score_codon,
                        tick_values=tick_values_codon, tick_labels=tick_labels_codon,
                        export=True, output=codon_heatmap_file, format='png', export_matrix=True)
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
            mean_val = scores_df[score_col].mean()
            min_val = scores_df[score_col].min()
            max_val = scores_df[score_col].max()
            stats['all_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
            stats['all_min'] = round(float(min_val)) if pd.notna(min_val) else None
            stats['all_max'] = round(float(max_val)) if pd.notna(max_val) else None
            
            # Add annotation-based stats
            if 'annotate_dna' in scores_df.columns:
                # WT stats from DNA level
                wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
                if len(wt_subset) > 0:
                    if experiment.barcoded:
                        # For barcoded experiments, include avg, min, max
                        mean_val = wt_subset[score_col].mean()
                        min_val = wt_subset[score_col].min()
                        max_val = wt_subset[score_col].max()
                        stats['wt_dna_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                        stats['wt_dna_min'] = round(float(min_val)) if pd.notna(min_val) else None
                        stats['wt_dna_max'] = round(float(max_val)) if pd.notna(max_val) else None
                    else:
                        # For non-barcoded experiments, include only avg
                        mean_val = wt_subset[score_col].mean()
                        stats['wt_dna'] = round(float(mean_val)) if pd.notna(mean_val) else None
                
                # Synonymous (WT) stats from DNA level
                syn_subset = scores_df[scores_df['annotate_dna'] == 'synonymous']
                if len(syn_subset) > 0:
                    mean_val = syn_subset[score_col].mean()
                    min_val = syn_subset[score_col].min()
                    max_val = syn_subset[score_col].max()
                    stats['synonymous_wt_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                    stats['synonymous_wt_min'] = round(float(min_val)) if pd.notna(min_val) else None
                    stats['synonymous_wt_max'] = round(float(max_val)) if pd.notna(max_val) else None
            
            # Missense and nonsense: use AA-aggregated scores with annotate_aa
            if 'aa_seq_diff' in scores_df.columns:
                aa_scores = aggregate_synonymous_variants(scores_df)
                if score_col in aa_scores.columns:
                    # Nonsense stats from AA level
                    nonsense_subset = aa_scores[aa_scores['annotate_aa'] == 'nonsense']
                    if len(nonsense_subset) > 0:
                        mean_val = nonsense_subset[score_col].mean()
                        min_val = nonsense_subset[score_col].min()
                        max_val = nonsense_subset[score_col].max()
                        stats['nonsense_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                        stats['nonsense_min'] = round(float(min_val)) if pd.notna(min_val) else None
                        stats['nonsense_max'] = round(float(max_val)) if pd.notna(max_val) else None
                    
                    # Missense stats from AA level
                    missense_subset = aa_scores[aa_scores['annotate_aa'] == 'missense_aa']
                    if len(missense_subset) > 0:
                        mean_val = missense_subset[score_col].mean()
                        min_val = missense_subset[score_col].min()
                        max_val = missense_subset[score_col].max()
                        stats['missense_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                        stats['missense_min'] = round(float(min_val)) if pd.notna(min_val) else None
                        stats['missense_max'] = round(float(max_val)) if pd.notna(max_val) else None
        
        stats_file = os.path.join(scores_dir, f"stats_{experiment.avg_method}_{experiment.experiment_name}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.json")
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        logging.info(f"Saved statistics to {stats_file}")
        
    except Exception as e:
        logging.warning(f"Failed to save statistics: {e}")
    
    print(f"Analysis complete! Results saved to {output_dir}")

if __name__ == "__main__":
    main()
