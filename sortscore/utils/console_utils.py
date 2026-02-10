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
        description="Run Sort-seq variant analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-n",
        "--experiment-name",
        type=str,
        required=True,
        help="Experiment name used for output file naming.",
    )
    parser.add_argument(
        "-e",
        "--experiment-setup-file",
        type=str,
        required=True,
        help="Path to experiment setup CSV.",
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
        default="."
    )
    parser.add_argument("--bins-required", 
                        type=int, 
                        required=False, 
                        help="Minimum number of bins required.",
                        default=1)
    parser.add_argument("--reps-required", 
                        type=int, 
                        required=False, 
                        help="Minimum number of replicates required.",
                        default=1)
    parser.add_argument("--avg-method",
                        type=str,
                        choices=["rep-weighted", "simple-avg"],
                        required=False,
                        help="Averaging method across replicates.",
                        default="rep-weighted"
    )
    parser.add_argument("--minread-threshold", 
                        type=int, 
                        required=False, 
                        help="Minimum read threshold.",
                        default=0)
    parser.add_argument("--max-cv", 
                        type=float, 
                        required=False, 
                        help="Maximum coefficient of variation allowed.",
                        default=None)
    parser.add_argument("--mutagenesis-variants",
                        type=str,
                        required=False,
                        help="Comma-separated list of mutagenesis variants.",
                        default=['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']
    )
    parser.add_argument("--position-offset", 
                        type=int, 
                        required=False, 
                        help="Offset for position numbering.",
                        default=0)
    parser.add_argument("--biophysical-prop",
                        action="store_true",
                        default=False,
                        help="Show biophysical properties panel in heatmaps.",
    )
    parser.add_argument("--position-type",
                        choices=["aa", "dna"],
                        required=False,
                        help="Position axis for plots ('aa' or 'dna'). Defaults based on detected variant type.",
                        default='aa'
    )
    parser.add_argument("--min-pos", 
                        type=int, 
                        required=False, 
                        help="Minimum position (1-based).",
                        default=1)
    parser.add_argument("--max-pos", 
                        type=int, 
                        required=False, 
                        help="Maximum position (1-based).",
                        default=None)
    parser.add_argument('-s', '--suffix', 
                        type=str, 
                        help='Custom suffix for output files.',
                        default="YYYYMMDD (current date)")
    parser.add_argument('-b', '--batch', 
                        action='store_true', 
                        help='Run batch processing mode for combining multiple experiments',
                        default=False)
    parser.add_argument('-p', '--pos-color', action='store_true', 
                        help='Export positional averages with colors for protein structure visualization.')
    parser.add_argument('--fig-format', choices=['png', 'svg', 'pdf'], default='png',
                        help='Output format for figures')
    return parser
