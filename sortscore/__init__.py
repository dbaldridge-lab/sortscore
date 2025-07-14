"""
sortscore: A modular Python package for Sort-seq variant analysis.

This package provides tools for analyzing Sort-seq experimental data,
calculating activity scores, and generating visualizations.
"""

__version__ = "0.1.0"
__author__ = "Caitlyn Chitwood"
__email__ = "c.chitwood@wustl.edu"

# Import main classes for convenience
from .analysis.load_experiment import ExperimentConfig
from .analysis.score import calculate_full_activity_scores

__all__ = [
    "ExperimentConfig", 
    "calculate_full_activity_scores"
]