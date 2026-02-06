"""
Experiment setup CSV loading and validation.

This module centralizes validation and column handling for the experiment setup
CSV (often named experiment_setup.csv). Downstream code relies on a consistent
set of required fields (replicate, bin, count file path, and median fluorescence 
intensity), but users may provide slightly different column names depending on 
their pipeline producing the count file inputs.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple

import pandas as pd

@dataclass(frozen=True)
class ExperimentSetupColumns:
    replicate: str
    bin: str
    count_file: str
    mfi: str
    read_count: Optional[str] = None
    proportion_of_cells: Optional[str] = None
    tile: Optional[str] = None


def load_experiment_setup(
    experiment_setup_file: str,
    *,
    validate_count_files: bool = False,
) -> Tuple[pd.DataFrame, ExperimentSetupColumns]:
    """
    Load and validate an experiment setup CSV.

    Parameters
    ----------
    experiment_setup_file : str
        Path to experiment setup CSV file.
    validate_count_files : bool
        If True, checks that every referenced count file exists on disk (after
        resolving paths relative to the setup file directory).

    Returns
    -------
    (pd.DataFrame, ExperimentSetupColumns)
        The raw DataFrame and the resolved column names to use.

    Raises
    ------
    FileNotFoundError
        If the setup file does not exist.
    ValueError
        If the setup file is unreadable, empty, missing required columns, or
        contains invalid values.
    """
    setup_path = Path(experiment_setup_file)
    if not setup_path.exists():
        raise FileNotFoundError(f"Experiment setup file not found: {experiment_setup_file}")

    try:
        setup_df = pd.read_csv(setup_path)
    except Exception as e:
        raise ValueError(f"Failed to read experiment setup CSV '{experiment_setup_file}': {e}") from e

    if setup_df.empty:
        raise ValueError(f"Experiment setup CSV is empty: {experiment_setup_file}")

    cols = _resolve_setup_columns(setup_df.columns)

    _validate_required_values(setup_df, cols, experiment_setup_file)

    if validate_count_files:
        _validate_count_files_exist(setup_df, cols, experiment_setup_file)

    return setup_df, cols


def _resolve_setup_columns(df_columns: Iterable[str]) -> ExperimentSetupColumns:
    columns = list(df_columns)

    def find_one(
        candidates: Iterable[str],
        *,
        required: bool = True,
        exclude_if_contains: Iterable[str] = (),
    ) -> Optional[str]:
        """
        Find a column name in `columns` that matches any of the provided candidates.

        This function normalizes both the actual column names and the candidate names
        (by stripping whitespace and converting to lowercase) to allow for flexible matching.
        It matches if any candidate appears as a substring in the column name.

        Parameters
        ----------
        candidates : Iterable[str]
            Possible column names to match against the DataFrame columns.
        required : bool, default True
            If True, raises an error if no match is found. If False, returns None.

        Returns
        -------
        Optional[str]
            The actual column name found in the DataFrame, or None if not found and not required.

        Raises
        ------
        ValueError
            If more than one column matches the candidates (ambiguous).
        """
        normalized_to_actual = {}
        for actual in columns:
            normalized = str(actual).strip().lower()
            normalized_to_actual.setdefault(normalized, []).append(actual)

        matches = []
        for candidate in candidates:
            norm_candidate = str(candidate).strip().lower()
            for norm_col, actuals in normalized_to_actual.items():
                if norm_candidate in norm_col:
                    if exclude_if_contains:
                        excluded = False
                        for token in exclude_if_contains:
                            if str(token).strip().lower() in norm_col:
                                excluded = True
                                break
                        if excluded:
                            continue
                    matches.extend(actuals)

        unique_matches = list(dict.fromkeys(matches))
        if not unique_matches:
            return None if not required else None
        if len(unique_matches) > 1:
            raise ValueError(f"Ambiguous columns matched {list(candidates)}: {unique_matches}")
        return unique_matches[0]

    replicate = find_one(["rep"])
    bin = find_one(["bin", "gate", "tube"])
    mfi = find_one(["mfi", 
                    "median fluorescence", 
                    "median_fluorescence",
                    "mean fluorescence"
                    "mean_fluorescence",
                    "afu",
                    "gfp"])

    count_file = find_one(
        [
            "path",
            "file",
            "counts",
            "tsv",
            "parquet"
        ],
        required=False,
    )

    if count_file is None:
        raise ValueError(
            "Experiment setup CSV must include a count file path column. "
            "Accepted column names must contain at least one of: 'path', 'file'."
        )

    read_count = find_one(
        [
            "read count",
            "read_count",
            "total reads",
            "total_reads",
        ],
        required=False,
        # Avoid treating file-path columns as numeric read-count columns.
        exclude_if_contains=("path", "file", "csv", "tsv", "parquet"),
    )

    proportion_of_cells = find_one(
        ["proportion of cells", 
        "proportion_of_cells",
        "cell prop",
        "cell_prop",
        "cell proportion", 
        "cell_proportion",
        "cell fraction",
        "cell_fraction"], required=False
    )

    tile = find_one(["tile", "oligo"], required=False)
    
    # replicate/bin/mfi are required
    missing = []
    if replicate is None:
        missing.append("Replicate")
    if bin is None:
        missing.append("Bin")
    if mfi is None:
        missing.append("MFI")
    if missing: 
        raise ValueError(f"Experiment setup CSV missing required column(s): {', '.join(missing)}")

    return ExperimentSetupColumns(
        replicate=replicate,
        bin=bin,
        count_file=count_file,
        mfi=mfi,
        read_count=read_count,
        proportion_of_cells=proportion_of_cells,
        tile=tile,
    )


def _validate_required_values(
    setup_df: pd.DataFrame,
    cols: ExperimentSetupColumns,
    experiment_setup_file: str,
) -> None:
    # Replicate must be int-like
    replicate_series = setup_df[cols.replicate]
    if replicate_series.isna().any():
        raise ValueError(f"Missing replicate values in experiment setup CSV: {experiment_setup_file}")
    try:
        replicate_series.astype(int)
    except Exception as e:
        raise ValueError(
            f"Replicate column '{cols.replicate}' must contain integers in {experiment_setup_file}"
        ) from e

    # Bin must be int-like
    bin_series = setup_df[cols.bin]
    if bin_series.isna().any():
        raise ValueError(f"Missing bin values in experiment setup CSV: {experiment_setup_file}")
    try:
        bin_series.astype(int)
    except Exception as e:
        raise ValueError(
            f"Bin column '{cols.bin}' must contain integers in {experiment_setup_file}"
        ) from e
        
    # Count file paths must be present
    count_series = setup_df[cols.count_file]
    if count_series.isna().any():
        raise ValueError(
            f"Missing count file path values in column '{cols.count_file}' for {experiment_setup_file}"
        )
    if (count_series.astype(str).str.strip() == "").any():
        raise ValueError(
            f"Empty count file path values found in column '{cols.count_file}' for {experiment_setup_file}"
        )

    # MFI must be float-like
    mfi_series = setup_df[cols.mfi]
    if mfi_series.isna().any():
        raise ValueError(f"Missing MFI values in experiment setup CSV: {experiment_setup_file}")
    try:
        mfi_series.astype(float)
    except Exception as e:
        raise ValueError(
            f"MFI column '{cols.mfi}' must contain numeric values in {experiment_setup_file}"
        ) from e


def _validate_count_files_exist(
    setup_df: pd.DataFrame,
    cols: ExperimentSetupColumns,
    experiment_setup_file: str,
) -> None:
    setup_dir = Path(experiment_setup_file).expanduser().resolve().parent

    missing = []
    for path in setup_df[cols.count_file].dropna().astype(str):
        user_expanded = Path(path.strip()).expanduser()
        if user_expanded.is_absolute():
            resolved = user_expanded
        else:
            resolved = setup_dir / user_expanded
        if not Path(resolved).exists():
            missing.append(path)

    if missing:
        preview = ", ".join(missing[:5])
        suffix = "" if len(missing) <= 5 else f" (and {len(missing) - 5} more)"
        raise FileNotFoundError(
            f"{len(missing)} count file(s) referenced in {experiment_setup_file} do not exist: {preview}{suffix}"
        )
