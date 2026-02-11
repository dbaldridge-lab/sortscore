"""
Console utility functions for command-line argument parsing.

This module provides reusable argument parsing functions for sortscore CLI tools.
"""
import argparse


def create_analysis_parser() -> argparse.ArgumentParser:
    """
    Create argument parser for Sort-seq variant analysis.
    
    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for the main analysis CLI.
    """
    parser = argparse.ArgumentParser(
        description="Run Sort-seq variant analysis. Supports scoring and a batch normalization mode with -b (for combining tiled experiments).",
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
    parser.add_argument('-b', '--batch', 
                        action='store_true', 
                        required=False,
                        help='Run batch processing mode for combining multiple experiments. If not provided, uses config value or default (False).'
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
    parser.add_argument("-O", "--position-offset",
                        type=int, 
                        required=False, 
                        help="Offset for position numbering. If not provided, uses config value or default (0)."
                            )
    parser.add_argument("--biophysical-prop",
                        action="store_true",
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
