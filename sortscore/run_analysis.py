"""
Main entry point for Sort-seq variant scoring.

This script loads the experiment configuration and orchestrates the analysis workflow.

Usage:
sortscore score -n EXPERIMENT -e path/to/experiment_setup.csv -c path/to/config.json
"""
import logging
import sys
import os
from typing import Optional, Dict, Any
import pandas as pd
from sortscore.utils.load_experiment import ExperimentConfig
from sortscore.utils.file_utils import ensure_output_subdirs
from sortscore.utils.analysis_logger import AnalysisLogger, generate_date_suffix
from sortscore.analysis.workflows import run_variant_analysis_workflow
from sortscore.visualization.heatmap_workflow import generate_heatmap_visualizations
from sortscore.utils.console_utils import parse_analysis_args, build_merged_analysis_config
from sortscore.utils.experiment_setup import load_experiment_setup
from sortscore.utils.tile_configs import (
    find_tile_output_dir_column,
    get_tile_output_dirs,
)


def _select_experiment_entry(merged: Dict[str, Any], tile: Optional[int] = None) -> Optional[Dict[str, Any]]:
    """Return per-experiment config entry from merged['experiments']."""
    entries = merged.get("experiments")
    if not isinstance(entries, list) or not entries:
        return None
    if tile is None:
        return entries[0] if len(entries) == 1 else None
    for entry in entries:
        try:
            if int(entry.get("tile")) == int(tile):
                return entry
        except Exception:
            continue
    return None


def _apply_experiment_overrides(base: Dict[str, Any], entry: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Apply supported per-experiment keys over a base config dict."""
    if not entry:
        return base
    out = dict(base)
    for key in ("wt_seq", "min_pos", "max_pos", "output_dir", "mutagenesis_type", "mutagenesis_variants"):
        if key in entry and entry[key] not in (None, ""):
            out[key] = entry[key]
    return out


def _run_single_analysis(
    experiment: ExperimentConfig,
    args,
    output_suffix: str,
    *,
    fail_on_viz_error: bool = True,
) -> None:
    """Run scoring + visualization for a single experiment config."""
    output_dir = experiment.output_dir or '.'
    ensure_output_subdirs(output_dir)

    logging.info(f"Loaded counts for {len(experiment.counts)} replicates.")
    print("Counts loaded.")

    fig_format = getattr(args, 'fig_format', None) or getattr(experiment, 'fig_format', None) or 'png'
    analysis_logger = AnalysisLogger(experiment, args, output_suffix, output_dir)

    try:
        dna_scores_file, aa_scores_file = run_variant_analysis_workflow(
            experiment, output_dir, output_suffix, analysis_logger
        )
        dna_scores_df = pd.read_csv(dna_scores_file) if dna_scores_file and os.path.exists(dna_scores_file) else None
        aa_scores_df = pd.read_csv(aa_scores_file) if aa_scores_file and os.path.exists(aa_scores_file) else None
        if dna_scores_df is not None:
            logging.info(f"Loaded DNA scores for visualization: {len(dna_scores_df)} variants")
        if aa_scores_df is not None:
            logging.info(f"Loaded AA scores for visualization: {len(aa_scores_df)} variants")
        if dna_scores_df is None and aa_scores_df is None:
            raise ValueError("No scores files were generated")
    except Exception as e:
        analysis_logger.add_error(f"Failed to run variant analysis workflow: {e}")
        logging.error(f"Failed to run variant analysis workflow: {e}")
        raise

    try:
        if experiment.mutagenesis_type == 'aa':
            if aa_scores_df is None:
                raise ValueError("mutagenesis_type='aa' requires AA scores for plotting, but none were generated.")
            generate_heatmap_visualizations(
                scores_df=aa_scores_df,
                experiment=experiment,
                output_dir=output_dir,
                output_suffix=output_suffix,
                fig_format=fig_format,
                export_positional_averages=args.pos_color,
            )
        elif experiment.mutagenesis_type == 'codon':
            if dna_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=dna_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
            if aa_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=aa_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
        else:
            if dna_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=dna_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
            if aa_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=aa_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
    except Exception as e:
        if fail_on_viz_error:
            logging.error(f"Failed to generate visualizations: {e}")
            raise
        logging.warning(f"Visualization step failed but scoring outputs were generated: {e}")

    log_file = analysis_logger.finish()
    print(f"Analysis complete! Results saved to {output_dir}")
    print(f"Analysis log saved to {log_file}")


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    args = parse_analysis_args()
    
    if not args.experiment_name:
        logging.error("Missing -n/--experiment-name (required).")
        sys.exit(1)
    if not args.experiment_setup_file:
        logging.error("Missing -e/--experiment-setup-file (required).")
        sys.exit(1)

    try:
        merged, config_dir = build_merged_analysis_config(args)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(1)

    # Generate string that will be used to name output files
    if args.suffix:
        output_suffix = args.suffix
        logging.info(f"Using custom output suffix: {output_suffix}")
    else:
        output_suffix = generate_date_suffix()
        logging.info(f"Using date-based output suffix: {output_suffix}")

    try:
        setup_df, setup_cols = load_experiment_setup(merged["experiment_setup_file"], require_tile=False)
    except Exception as e:
        logging.error(f"Failed to load experiment setup CSV: {e}")
        sys.exit(1)

    has_multitile = (
        setup_cols.tile is not None and
        setup_df[setup_cols.tile].dropna().astype(int).nunique() > 1
    )

    if not has_multitile:
        try:
            if isinstance(merged.get("experiments"), list) and len(merged["experiments"]) > 1:
                raise ValueError("Single-experiment scoring requires exactly one entry in config.experiments.")
            entry = _select_experiment_entry(merged, tile=None)
            merged_single = _apply_experiment_overrides(merged, entry)
            experiment = ExperimentConfig.from_dict(merged_single, config_file_dir=config_dir)
            _run_single_analysis(experiment, args, output_suffix, fail_on_viz_error=True)
        except Exception as e:
            logging.error(f"Failed to run analysis: {e}")
            sys.exit(1)
        return

    base_output_dir = merged.get("output_dir", ".")
    existing_col = find_tile_output_dir_column(setup_df.columns)
    tile_output_dirs = get_tile_output_dirs(
        setup_df,
        setup_cols.tile,
        base_output_dir,
        existing_output_dir_col=existing_col,
    )

    logging.info(
        f"Detected {len(tile_output_dirs)} tiles. Running scoring independently per tile using one shared setup file."
    )
    failures = 0
    for tile in sorted(tile_output_dirs.keys()):
        tile_entry = _select_experiment_entry(merged, tile=tile)
        if isinstance(merged.get("experiments"), list) and merged["experiments"] and tile_entry is None:
            failures += 1
            logging.error(f"Tile {tile} failed: no matching config.experiments entry for tile={tile}")
            continue
        tile_merged = _apply_experiment_overrides(dict(merged), tile_entry)
        tile_merged["output_dir"] = tile_merged.get("output_dir") or tile_output_dirs[tile]
        tile_merged["experiment_name"] = f"{merged['experiment_name']}_tile{tile}"
        tile_merged["tile_id"] = int(tile)

        logging.info(f"Running tile {tile}: output_dir={tile_merged['output_dir']}")
        try:
            experiment = ExperimentConfig.from_dict(tile_merged, config_file_dir=config_dir)
            _run_single_analysis(experiment, args, output_suffix, fail_on_viz_error=True)
        except Exception as e:
            failures += 1
            logging.error(f"Tile {tile} failed: {e}")

    if failures > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
