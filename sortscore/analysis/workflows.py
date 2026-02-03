"""
Variant analysis workflows for different experiment types.

This module provides separate workflow functions for DNA, AA, and SNV analysis,
allowing for clear separation of concerns and type-specific processing logic.
"""
import os
import logging
import pandas as pd
from typing import Optional, Tuple
from sortscore.analysis.score import calculate_full_activity_scores
from sortscore.analysis.annotation import annotate_scores_dataframe
from sortscore.analysis.data_processing import aggregate_synonymous_variants
from sortscore.analysis.statistics import calculate_replicate_statistics, round_score_columns, get_replicate_score_columns
from sortscore.analysis.summary_stats import calculate_summary_stats, save_summary_stats
from sortscore.analysis.aa_scores import process_and_save_aa_scores


def calculate_variant_scores(experiment, merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Common scoring function used by all workflow types.
    
    Parameters
    ----------
    experiment : ExperimentConfig
        Experiment configuration
    merged_df : pd.DataFrame
        Merged counts DataFrame
        
    Returns
    -------
    pd.DataFrame
        Scored and annotated DataFrame
    """
    # Calculate activity scores
    scores_df = calculate_full_activity_scores(
        counts=experiment.counts,
        mfi=experiment.mfi,
        min_bins=experiment.bins_required,
        min_reps=experiment.reps_required,
        minread_threshold=experiment.minread_threshold,
        avg_method=experiment.avg_method,
        total_reads=experiment.total_reads,
        cell_prop=experiment.cell_prop,
        merged_df=merged_df,
        max_cv=experiment.max_cv
    )
    
    # Annotate sequences
    scores_df = annotate_scores_dataframe(
        scores_df, 
        experiment.wt_seq, 
        experiment.variant_type, 
        experiment.transcript_id
    )
    
    return scores_df


def process_dna_workflow(experiment, output_dir: str, output_suffix: str, analysis_logger) -> Optional[str]:
    """
    Process DNA-level analysis workflow (codon or snv analysis types).
    
    Parameters
    ----------
    experiment : ExperimentConfig
        Experiment configuration
    output_dir : str
        Output directory
    output_suffix : str
        Output file suffix
    analysis_logger : AnalysisLogger
        Analysis logger
        
    Returns
    -------
    Optional[str]
        Path to DNA scores file if created, None if skipped
    """
    if experiment.variant_type != 'dna':
        logging.info("Skipping DNA workflow (auto-detected variant type is not 'dna')")
        return None
    
    if experiment.analysis_type not in ['codon', 'snv']:
        logging.info(
            "Skipping DNA workflow (analysis_type '%s' is not 'codon' or 'snv')",
            experiment.analysis_type,
        )
        return None
    
    logging.info(f"Processing DNA workflow for analysis_type '{experiment.analysis_type}'...")
    
    # Get merged counts DataFrame
    merged_df = experiment.get_merged_counts()
    
    # Calculate scores using common function
    scores_df = calculate_variant_scores(experiment, merged_df)
    logging.info(f"Calculated DNA scores for {len(scores_df)} variants.")
    analysis_logger.set_processing_stats(len(scores_df))
    
    # Prepare for saving
    scores_dir = os.path.join(output_dir, 'scores')
    scores_df_rounded = scores_df.copy()
    rep_score_columns = get_replicate_score_columns(scores_df_rounded)
    scores_df_rounded = calculate_replicate_statistics(scores_df_rounded, rep_score_columns)
    scores_df_rounded = round_score_columns(scores_df_rounded)
    
    # Save DNA scores
    analysis_type_suffix = f"_{experiment.analysis_type}" if experiment.analysis_type != 'codon' else ""
    dna_scores_file = os.path.join(scores_dir, f"{experiment.experiment_name}_dna_scores{analysis_type_suffix}_{output_suffix}.csv")
    scores_df_rounded.to_csv(dna_scores_file, index=False)
    logging.info(f"Saved DNA scores to {dna_scores_file}")
    
    # Log output
    analysis_logger.log_output_file(
        'dna_scores', 
        f"{experiment.experiment_name}_dna_scores{analysis_type_suffix}_{output_suffix}.csv",
        dna_scores_file,
        variant_count=len(scores_df_rounded)
    )
    
    # Calculate and save statistics on the final processed DNA data
    stats = calculate_summary_stats(scores_df_rounded, experiment)
    save_summary_stats(
        stats,
        experiment,
        scores_dir,
        output_suffix,
        analysis_logger,
        stats_basename="dna_stats",
        output_field="dna_statistics",
    )
    logging.info("Calculated statistics on final DNA scores")
    
    return dna_scores_file


def process_aa_workflow(experiment, output_dir: str, output_suffix: str, analysis_logger, 
                       dna_scores_file: Optional[str] = None) -> str:
    """
    Process AA-level analysis workflow.
    
    This handles both:
    - AA analysis from DNA data (aggregate synonymous variants)
    - AA analysis from AA data (direct processing)
    
    Parameters
    ----------
    experiment : ExperimentConfig
        Experiment configuration
    output_dir : str
        Output directory
    output_suffix : str
        Output file suffix
    analysis_logger : AnalysisLogger
        Analysis logger
    dna_scores_file : Optional[str]
        Path to DNA scores file if available
        
    Returns
    -------
    str
        Path to AA scores file
    """
    should_run = (experiment.analysis_type == 'aa') or (experiment.variant_type == 'dna')
    if not should_run:
        logging.info(
            "Skipping AA workflow (analysis_type '%s' with variant_type '%s')",
            experiment.analysis_type,
            experiment.variant_type,
        )
        return ""
    
    logging.info("Processing AA workflow...")
    scores_dir = os.path.join(output_dir, 'scores')
    
    if dna_scores_file and os.path.exists(dna_scores_file):
        # Load DNA scores and aggregate to AA level
        logging.info(f"Loading DNA scores from {dna_scores_file} for AA aggregation")
        dna_scores_df = pd.read_csv(dna_scores_file)
        
        if 'aa_seq_diff' not in dna_scores_df.columns:
            raise ValueError("DNA scores file missing aa_seq_diff column needed for AA aggregation")
        
        # Aggregate to AA level
        scores_df = aggregate_synonymous_variants(dna_scores_df)
        logging.info(f"Aggregated to {len(scores_df)} unique AA variants from DNA data")
        
    else:
        # Build AA table directly from counts (AA-only experiment or no DNA scores)
        if experiment.variant_type == 'aa':
            logging.info("Building AA scores directly from AA input data")
        else:
            logging.info("Building AA scores directly from DNA input data (no pre-existing DNA scores)")
        
        merged_df = experiment.get_merged_counts()
        scores_df = calculate_variant_scores(experiment, merged_df)
        logging.info(f"Calculated AA scores for {len(scores_df)} variants.")
        analysis_logger.set_processing_stats(len(scores_df))
    
    # Process and save AA scores using existing function
    process_and_save_aa_scores(scores_df, experiment, scores_dir, output_suffix, analysis_logger)
    
    # Calculate and save statistics on the final processed AA data
    aa_scores_file = os.path.join(scores_dir, f"{experiment.experiment_name}_aa_scores_{output_suffix}.csv")
    
    if os.path.exists(aa_scores_file):
        final_scores_df = pd.read_csv(aa_scores_file)
        stats = calculate_summary_stats(final_scores_df, experiment)
        save_summary_stats(
            stats,
            experiment,
            scores_dir,
            output_suffix,
            analysis_logger,
            stats_basename="aa_stats",
            output_field="aa_statistics",
        )
        logging.info("Calculated statistics on final AA scores")
    
    return aa_scores_file


def run_variant_analysis_workflow(experiment, output_dir: str, output_suffix: str, analysis_logger) -> Tuple[Optional[str], str]:
    """
    Run the complete variant analysis workflow based on analysis_type.
    
    This orchestrates the DNA and AA workflows in the correct order:
    1. DNA workflow (if needed for codon/snv analysis)
    2. AA workflow (if needed for aa analysis)
    
    Parameters
    ----------
    experiment : ExperimentConfig
        Experiment configuration
    output_dir : str
        Output directory
    output_suffix : str
        Output file suffix
    analysis_logger : AnalysisLogger
        Analysis logger
        
    Returns
    -------
    Tuple[Optional[str], str]
        (dna_scores_file_path, aa_scores_file_path)
        Either path can be empty string if that workflow was skipped
    """
    logging.info(f"Running variant analysis workflow for analysis_type: '{experiment.analysis_type}'")
    logging.info(f"Auto-detected variant_type: '{experiment.variant_type}'")
    
    dna_scores_file = process_dna_workflow(experiment, output_dir, output_suffix, analysis_logger)
    aa_scores_file = process_aa_workflow(experiment, output_dir, output_suffix, analysis_logger, dna_scores_file)
    
    return dna_scores_file, aa_scores_file
