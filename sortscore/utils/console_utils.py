"""
Console utility functions for command-line argument parsing.

This module provides reusable argument parsing functions for sortscore CLI tools.
"""
import argparse
import json
from pathlib import Path
from typing import Dict, Any, Optional, Tuple


def create_analysis_parser() -> argparse.ArgumentParser:
    """
    Create argument parser for Sort-seq variant analysis.
    
    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for the main analysis CLI.
    """
    parser = argparse.ArgumentParser(
        description="Run Sort-seq variant analysis (single-experiment scoring).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-n", "--experiment-name",
                        type=str,
                        required=True,
                        help="Experiment name used for output file naming."
    )
    parser.add_argument("-e", "--experiment-setup-file",
                        type=str,
                        required=True,
                        help="Path to experiment setup CSV."
    )               
    parser.add_argument("-c", "--config",
                        type=str,
                        required=False,
                        help="Optional JSON config file. CLI options take precedence over config values."
    )
    parser.add_argument("-w","--wt-seq",
                        type=str,
                        required=False,
                        help="Wild-type sequence. Required unless provided in --config."
    )
    parser.add_argument("-o", "--output-dir",
                        type=str,
                        required=False,
                        help="Output directory. If not provided, uses config value or current directory."
    )
    parser.add_argument("-B", "--bins-required", 
                        type=int, 
                        required=False, 
                        help="Minimum number of bins required. If not provided, uses config value or default (1)."
                            )
    parser.add_argument("-R", "--reps-required", 
                        type=int, 
                        required=False, 
                        help="Minimum number of replicates required. If not provided, uses config value or default (1)."
                            )
    parser.add_argument("-a", "--avg-method",
                        type=str,
                        choices=["rep-weighted", "simple-avg"],
                        required=False,
                        help="Averaging method across replicates. If not provided, uses config value or default ('rep-weighted')."
                            )
    parser.add_argument("-m", "--minread-threshold", 
                        type=int, 
                        required=False, 
                        help="Minimum read threshold. If not provided, uses config value or default (0)."
                            )
    parser.add_argument("--max-cv", 
                        type=float, 
                        required=False, 
                        help="Maximum coefficient of variation allowed. If not provided, uses config value or default (None)."
                            )
    parser.add_argument("-v", "--mutagenesis-variants",
                        type=str,
                        required=False,
                        help="Comma-separated list of mutagenesis variants. If not provided, uses config value or default (all 20 AAs + *)."
                            )
    parser.add_argument("--mutagenesis-type",
                        type=str,
                        choices=["aa", "codon", "snv"],
                        required=False,
                        help="Mutagenesis type. Default is 'aa'. Set to 'codon' or 'snv' for DNA-based analyses."
                            )
    parser.add_argument("-O", "--position-offset",
                        type=int, 
                        required=False, 
                        help="Offset for position numbering. If not provided, uses config value or default (0)."
                            )
    parser.add_argument("--biophysical-prop",
                        action="store_true",
                        default=None,
                        required=False,
                        help="Show biophysical properties panel in heatmaps. If not provided, uses config value or default (False)."
                            )
    parser.add_argument("--min-pos", 
                        type=int, 
                        required=False, 
                        help="Minimum position (1-based). If not provided, uses config value or default (1)."
                            )
    parser.add_argument("--max-pos", 
                        type=int, 
                        required=False, 
                        help="Maximum position (1-based). If not provided, uses config value or default (None)."
                            )
    parser.add_argument('-s', '--suffix', 
                        type=str, 
                        required=False,
                        help='Custom suffix for output files. If not provided, uses date-based suffix (YYYYMMDD).'
                            )
    parser.add_argument('-p', '--pos-color', 
                        action='store_true', 
                        help='Export positional averages with colors for protein structure visualization. If not provided, uses config value or default (False).'
                            )
    parser.add_argument('--fig-format', 
                        choices=['png', 'svg', 'pdf'], 
                        required=False,
                        help='Output format for figures. If not provided, uses config value or default (png).')
    return parser


def parse_analysis_args(argv=None) -> argparse.Namespace:
    """
    Parse CLI arguments for the main analysis entrypoint.
    """
    parser = create_analysis_parser()
    return parser.parse_args(argv)


def build_merged_analysis_config(
    args: argparse.Namespace,
) -> Tuple[Dict[str, Any], Optional[Path]]:
    """
    Merge optional JSON config with CLI overrides.

    Returns
    -------
    Tuple[Dict[str, Any], Optional[pathlib.Path]]
        (merged_config_mapping, config_dir)
    """
    merged: Dict[str, Any] = {}
    config_dir: Optional[Path] = None

    if args.config:
        try:
            config_path = Path(args.config).expanduser().resolve()
            config_dir = config_path.parent
            with open(config_path, "r") as f:
                merged.update(json.load(f))
        except Exception as e:
            raise ValueError(f"Failed to load config JSON '{args.config}': {e}") from e

    merged["experiment_name"] = args.experiment_name
    merged["experiment_setup_file"] = str(Path(args.experiment_setup_file).expanduser().resolve())
    if args.wt_seq is not None:
        merged["wt_seq"] = args.wt_seq
    if args.output_dir is not None:
        merged["output_dir"] = str(Path(args.output_dir).expanduser().resolve())
    if args.bins_required is not None:
        merged["bins_required"] = args.bins_required
    if args.reps_required is not None:
        merged["reps_required"] = args.reps_required
    if args.avg_method is not None:
        merged["avg_method"] = args.avg_method
    if args.minread_threshold is not None:
        merged["minread_threshold"] = args.minread_threshold
    if args.max_cv is not None:
        merged["max_cv"] = args.max_cv
    if args.mutagenesis_variants is not None:
        if isinstance(args.mutagenesis_variants, str):
            merged["mutagenesis_variants"] = [v.strip() for v in args.mutagenesis_variants.split(",") if v.strip()]
        elif isinstance(args.mutagenesis_variants, list):
            merged["mutagenesis_variants"] = args.mutagenesis_variants
        else:
            merged["mutagenesis_variants"] = list(args.mutagenesis_variants)
    if args.mutagenesis_type is not None:
        merged["mutagenesis_type"] = args.mutagenesis_type
    if args.biophysical_prop is not None:
        merged["biophysical_prop"] = bool(args.biophysical_prop)
    if args.min_pos is not None:
        merged["min_pos"] = args.min_pos
    if args.max_pos is not None:
        merged["max_pos"] = args.max_pos

    return merged, config_dir
