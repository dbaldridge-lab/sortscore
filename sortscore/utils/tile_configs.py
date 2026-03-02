"""Helpers for tile-aware scoring and normalization workflows."""

from pathlib import Path
from typing import Dict, Optional, Any

import pandas as pd
import json


def find_tile_output_dir_column(columns) -> Optional[str]:
    """Find a column that appears to store tile output directories."""
    for col in columns:
        norm = str(col).strip().lower()
        if "tile" in norm and "output" in norm and "dir" in norm:
            return col
    return None


def get_tile_output_dirs(
    setup_df: pd.DataFrame,
    tile_col: str,
    base_output_dir: str,
    existing_output_dir_col: Optional[str] = None,
) -> Dict[int, str]:
    """
    Resolve output directory per tile.

    If `existing_output_dir_col` is provided and values are present, those are used.
    Otherwise defaults to `<base_output_dir>/tile_<tile>`.
    """
    base_dir = Path(base_output_dir).expanduser().resolve()
    tile_dirs: Dict[int, str] = {}

    grouped = setup_df.groupby(setup_df[tile_col].astype(int), sort=True)
    for tile, tile_df in grouped:
        resolved = None
        if existing_output_dir_col is not None and existing_output_dir_col in tile_df.columns:
            values = (
                tile_df[existing_output_dir_col]
                .dropna()
                .astype(str)
                .map(str.strip)
            )
            values = values[values != ""]
            if not values.empty:
                resolved = str(Path(values.iloc[0]).expanduser().resolve())

        if resolved is None:
            resolved = str((base_dir / f"tile_{int(tile)}").resolve())

        tile_dirs[int(tile)] = resolved

    return tile_dirs


def write_tile_output_dirs_to_setup(
    experiment_setup_file: str,
    setup_df: pd.DataFrame,
    tile_col: str,
    tile_output_dirs: Dict[int, str],
    *,
    output_col: str = "Tile Output Dir",
) -> None:
    """Write/update tile output directory column in experiment setup CSV."""
    updated = setup_df.copy()
    updated[output_col] = (
        updated[tile_col]
        .astype(int)
        .map(lambda tile: tile_output_dirs.get(int(tile)))
    )
    updated.to_csv(experiment_setup_file, index=False)


def build_tile_experiments(
    base_config: dict,
    shared_setup_file: str,
    tile_output_dirs: Dict[int, str],
) -> list[dict]:
    """
    Build one batch config entry per tile.
    """
    tile_configs: list[dict] = []

    for tile in sorted(tile_output_dirs.keys()):
        tile_cfg = {
            "tile": int(tile),
            "output_dir": str(Path(tile_output_dirs[tile]).expanduser().resolve()),
            "wt_seq": base_config.get("wt_seq"),
            "min_pos": base_config.get("min_pos"),
            "max_pos": base_config.get("max_pos"),
        }
        tile_configs.append(tile_cfg)

    return tile_configs


def write_batch_config_json(
    base_config: Dict[str, Any],
    tile_experiments: list[dict],
    *,
    output_path: str,
) -> str:
    """Write a batch config JSON and return its path."""
    payload = {
        "experiments": tile_experiments,
        "batch_normalization_method": base_config.get("batch_normalization_method", "zscore_2pole"),
        "pathogenic_control_type": base_config.get("pathogenic_control_type", "nonsense"),
        "pathogenic_variants": base_config.get("pathogenic_variants"),
        "combined_output_dir": base_config.get("combined_output_dir", "./normalized"),
    }
    out = Path(output_path).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    return str(out)
