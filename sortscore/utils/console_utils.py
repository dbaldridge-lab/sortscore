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
    parser = argparse.ArgumentParser(description="Run Sort-seq variant analysis.")
    parser.add_argument(
        "-n",
        "--experiment-name",
        type=str,
        required=False,
        help="Experiment name used for output file naming (required unless --batch).",
    )
    parser.add_argument(
        "-e",
        "--experiment-setup-file",
        type=str,
        required=False,
        help="Path to experiment setup CSV (required unless --batch).",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=False,
        help="Optional JSON config file. CLI options take precedence over config values.",
    )
    parser.add_argument(
        "-w",
        "--wt-seq",
        type=str,
        required=False,
        help="Wild-type sequence. Required unless provided in --config.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=False,
        help="Output directory. Defaults to config value or current directory.",
    )
    parser.add_argument("--bins-required", type=int, required=False, help="Minimum number of bins required.")
    parser.add_argument("--reps-required", type=int, required=False, help="Minimum number of replicates required.")
    parser.add_argument(
        "--avg-method",
        type=str,
        choices=["rep-weighted", "simple-avg"],
        required=False,
        help="Averaging method across replicates.",
    )
    parser.add_argument("--minread-threshold", type=int, required=False, help="Minimum read threshold.")
    parser.add_argument("--max-cv", type=float, required=False, help="Maximum coefficient of variation allowed.")
    parser.add_argument(
        "--mutagenesis-variants",
        type=str,
        required=False,
        help="Comma-separated list of mutagenesis variants (e.g. 'G,C,T,A').",
    )
    parser.add_argument("--position-offset", type=int, required=False, help="Offset for position numbering.")
    parser.add_argument(
        "--biophysical-prop",
        action="store_true",
        default=False,
        help="Show biophysical properties panel in heatmaps.",
    )
    parser.add_argument(
        "--position-type",
        choices=["aa", "dna"],
        required=False,
        help="Position axis for plots ('aa' or 'dna'). Defaults based on detected variant type.",
    )
    parser.add_argument("--min-pos", type=int, required=False, help="Minimum position (1-based).")
    parser.add_argument("--max-pos", type=int, required=False, help="Maximum position (1-based).")
    parser.add_argument('-s', '--suffix', type=str, 
                       help='Custom suffix for output files (default: auto-generated from config)')
    parser.add_argument('-b', '--batch', action='store_true', 
                       help='Run batch processing mode for combining multiple experiments')
    parser.add_argument('-p', '--pos-color', action='store_true', 
                       help='Export positional averages with colors for protein structure visualization')
    parser.add_argument('--fig-format', choices=['png', 'svg', 'pdf'], default='png',
                       help='Output format for figures (default: png)')
    return parser
