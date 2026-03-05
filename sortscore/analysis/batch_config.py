"""
Batch configuration module for Sort-seq variant analysis.

This module provides framework for combining and normalizing multiple 
Sort-seq experiments in batch processing workflows.
"""

import json
from pathlib import Path
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
    batch_normalization_method : str, optional
        Normalization method to use ('zscore_2pole', '2pole', or 'zscore_center'), default 'zscore_2pole'
    pathogenic_control_type : str, optional
        Type of pathogenic control ('nonsense' or 'custom'), default 'nonsense'  
    pathogenic_variants : Optional[List[str]], optional
        Custom pathogenic variants if using 'custom' pathogenic_control_type
    combined_output_dir : str, optional
        Directory for final combined results, default current directory
    """
    experiments: List[Dict[str, Any]]
    batch_normalization_method: str = 'zscore_2pole'
    pathogenic_control_type: str = 'nonsense'
    pathogenic_variants: Optional[List[str]] = None
    combined_output_dir: str = './normalized'
    
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
        """
        json_path_obj = Path(json_path).expanduser().resolve()
        with open(json_path_obj, 'r') as f:
            data = json.load(f)
        if 'experiments' not in data:
            raise ValueError("Batch configuration requires 'experiments'")
        args: Dict[str, Any] = {'experiments': data['experiments']}
        
        # Optional fields (only add if present to preserve defaults)
        optional_fields = [
            'batch_normalization_method', 'pathogenic_control_type', 
            'pathogenic_variants', 'combined_output_dir'
        ]
        
        for field in optional_fields:
            if field in data:
                args[field] = data[field]

        # Resolve output paths relative to the config file directory.
        config_dir = json_path_obj.parent
        experiments = []
        for entry in args['experiments']:
            entry_cfg = dict(entry)
            if 'output_dir' in entry_cfg and entry_cfg['output_dir'] is not None:
                entry_cfg['output_dir'] = str(
                    (config_dir / Path(str(entry_cfg['output_dir'])).expanduser()).resolve()
                )
            experiments.append(entry_cfg)
        args['experiments'] = experiments

        if 'combined_output_dir' in args and args['combined_output_dir'] is not None:
            args['combined_output_dir'] = str(
                (config_dir / Path(str(args['combined_output_dir'])).expanduser()).resolve()
            )
        
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
        if not self.experiments:
            raise ValueError("experiments must contain at least one entry")

        seen_tiles = set()
        for idx, cfg in enumerate(self.experiments):
            if 'tile' not in cfg:
                raise ValueError(f"experiments[{idx}] missing 'tile'")
            if 'output_dir' not in cfg:
                raise ValueError(f"experiments[{idx}] missing 'output_dir'")
            if 'wt_seq' not in cfg:
                raise ValueError(f"experiments[{idx}] missing 'wt_seq'")
            if 'min_pos' not in cfg:
                raise ValueError(f"experiments[{idx}] missing 'min_pos'")
            if 'max_pos' not in cfg:
                raise ValueError(f"experiments[{idx}] missing 'max_pos'")
            try:
                tile = int(cfg['tile'])
            except Exception as e:
                raise ValueError(f"experiments[{idx}].tile must be int-like") from e
            try:
                min_pos = int(cfg['min_pos'])
                max_pos = int(cfg['max_pos'])
            except Exception as e:
                raise ValueError(f"experiments[{idx}].min_pos/max_pos must be int-like") from e
            if min_pos >= max_pos:
                raise ValueError(f"experiments[{idx}] min_pos must be less than max_pos")
            if tile in seen_tiles:
                raise ValueError(f"Duplicate tile value in experiments: {tile}")
            seen_tiles.add(tile)

            out_dir = Path(str(cfg['output_dir'])).expanduser().resolve()
            if not out_dir.exists():
                raise FileNotFoundError(f"experiments[{idx}] output_dir does not exist: {out_dir}")
        
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
        
    def get_batch_config_dict(self) -> Dict[str, Any]:
        """
        Convert batch configuration to dictionary format for processing.
        
        Returns
        -------
        Dict[str, Any]
            Configuration dictionary for batch processing functions
        """
        return {
            'experiments': self.experiments,
            'batch_normalization_method': self.batch_normalization_method,
            'pathogenic_control_type': self.pathogenic_control_type,
            'pathogenic_variants': self.pathogenic_variants,
            'combined_output_dir': self.combined_output_dir
        }
