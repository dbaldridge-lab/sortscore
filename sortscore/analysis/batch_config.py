"""
Batch configuration module for Sort-seq variant analysis.

This module provides framework for combining and normalizing multiple 
Sort-seq experiments in batch processing workflows.
"""

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Any


@dataclass
class BatchConfig:
    """
    Dataclass for batch processing configuration.
    
    This class manages configuration for combining and normalizing multiple Sort-seq experiments.
    It supports two normalization methods:
    
    1. **2-pole normalization**: Uses synonymous and pathogenic variants as reference points
    2. **Z-score scaled three-step normalization** (default): Creates standardized scale where synonymous 
       variants center around 0 with unit variance, making cross-experiment comparisons meaningful
    
    Attributes
    ----------
    experiment_configs : List[str]
        List of paths to individual experiment configuration files
    batch_normalization_method : str, optional
        Normalization method to use ('zscore_2pole', '2pole', or 'zscore_center'), default 'zscore_2pole'
    pathogenic_control_type : str, optional
        Type of pathogenic control ('nonsense' or 'custom'), default 'nonsense'  
    pathogenic_variants : Optional[List[str]], optional
        Custom pathogenic variants if using 'custom' pathogenic_control_type
    combined_output_dir : str, optional
        Directory for final combined results, default current directory
    global_min_pos : Optional[int], optional
        Overall minimum position across all experiments (for tiled heatmaps)
    global_max_pos : Optional[int], optional
        Overall maximum position across all experiments (for tiled heatmaps)
    allow_position_breaks : bool, optional
        Whether to allow gaps/breaks in tiled position display, default True
    cleanup_individual_files : bool, optional
        Whether to delete individual experiment output files after combination, default True
    """
    experiment_configs: List[str]
    batch_normalization_method: str = 'zscore_2pole'
    pathogenic_control_type: str = 'nonsense'
    pathogenic_variants: Optional[List[str]] = None
    combined_output_dir: str = '.'
    global_min_pos: Optional[int] = None
    global_max_pos: Optional[int] = None
    allow_position_breaks: bool = True
    cleanup_individual_files: bool = True
    
    @staticmethod
    def from_json(json_path: str) -> 'BatchConfig':
        """
        Load batch configuration from a JSON file.
        
        Parameters
        ----------
        json_path : str
            Path to the JSON batch config file
            
        Returns
        -------
        config : BatchConfig
            Loaded batch configuration
            
        Examples
        --------
        Example JSON structure:
        
        .. code-block:: json
        
            {
                "experiment_configs": [
                    "/path/to/experiment1/config.json",
                    "/path/to/experiment2/config.json",
                    "/path/to/experiment3/config.json"
                ],
                "batch_normalization_method": "zscore_2pole",
                "pathogenic_control_type": "nonsense",
                "combined_output_dir": "/path/to/combined/results",
                "global_min_pos": 1,
                "global_max_pos": 500,
                "allow_position_breaks": true,
                "cleanup_individual_files": true
            }
        """
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Required field
        if 'experiment_configs' not in data:
            raise ValueError("experiment_configs is required in batch configuration")
        
        args = {'experiment_configs': data['experiment_configs']}
        
        # Optional fields (only add if present to preserve defaults)
        optional_fields = [
            'batch_normalization_method', 'pathogenic_control_type', 
            'pathogenic_variants', 'combined_output_dir', 'global_min_pos', 
            'global_max_pos', 'allow_position_breaks', 'cleanup_individual_files'
        ]
        
        for field in optional_fields:
            if field in data:
                args[field] = data[field]
        
        return BatchConfig(**args)
    
    def validate_config(self) -> None:
        """
        Validate batch configuration parameters.
        
        Raises
        ------
        ValueError
            If configuration parameters are invalid
        FileNotFoundError
            If experiment config files don't exist
        """
        # Check that experiment config files exist
        for config_path in self.experiment_configs:
            if not os.path.exists(config_path):
                raise FileNotFoundError(f"Experiment config file not found: {config_path}")
        
        # Validate normalization method
        valid_methods = ['zscore_2pole', '2pole', 'zscore_center']
        if self.batch_normalization_method not in valid_methods:
            raise ValueError(f"Invalid normalization method: {self.batch_normalization_method}. "
                           f"Must be one of: {valid_methods}")
        
        # Validate pathogenic control type
        valid_control_types = ['nonsense', 'custom']
        if self.pathogenic_control_type not in valid_control_types:
            raise ValueError(f"Invalid pathogenic control type: {self.pathogenic_control_type}. "
                           f"Must be one of: {valid_control_types}")
        
        # If using custom pathogenic controls, ensure variants are specified
        if self.pathogenic_control_type == 'custom' and not self.pathogenic_variants:
            raise ValueError("pathogenic_variants must be specified when using 'custom' pathogenic_control_type")
        
        # Validate position parameters
        if self.global_min_pos is not None and self.global_max_pos is not None:
            if self.global_min_pos >= self.global_max_pos:
                raise ValueError("global_min_pos must be less than global_max_pos")
    
    def get_batch_config_dict(self) -> Dict[str, Any]:
        """
        Convert batch configuration to dictionary format for processing.
        
        Returns
        -------
        Dict[str, Any]
            Configuration dictionary for batch processing functions
        """
        return {
            'experiment_configs': self.experiment_configs,
            'batch_normalization_method': self.batch_normalization_method,
            'pathogenic_control_type': self.pathogenic_control_type,
            'pathogenic_variants': self.pathogenic_variants,
            'combined_output_dir': self.combined_output_dir,
            'global_min_pos': self.global_min_pos,
            'global_max_pos': self.global_max_pos,
            'allow_position_breaks': self.allow_position_breaks,
            'cleanup_individual_files': self.cleanup_individual_files
        }