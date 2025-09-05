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
    parser.add_argument('-c', '--config', type=str, required=True, 
                       help='Path to experiment config JSON file (or batch config with --batch)')
    parser.add_argument('-s', '--suffix', type=str, 
                       help='Custom suffix for output files (default: auto-generated from config)')
    parser.add_argument('-b', '--batch', action='store_true', 
                       help='Run batch processing mode for combining multiple experiments')
    parser.add_argument('-p', '--pos-color', action='store_true', 
                       help='Export positional averages with colors for protein structure visualization')
    parser.add_argument('--fig-format', choices=['png', 'svg', 'pdf'], default='png',
                       help='Output format for figures (default: png)')
    return parser