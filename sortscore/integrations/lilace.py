"""
Lilace export helpers.
"""
import re
from typing import Literal
from pathlib import Path
import pandas as pd

COUNT_COLUMN_PATTERN = re.compile(r"^count\.r(?P<replicate>\d+)b(?P<bin>\d+)$")
POSITION_PATTERN = re.compile(r"(?P<position>\d+)")


def _get_mutagenesis_columns(mutagenesis_type: Literal["aa", "codon"]) -> tuple[str, str]:
    """Return the variant and annotation column names for one Lilace export type."""
    if mutagenesis_type == "aa":
        return "aa_seq_diff", "annotate_aa"
    if mutagenesis_type == "codon":
        return "dna_seq_diff", "annotate_dna"
    raise ValueError(f"Unsupported Lilace export mutagenesis_type: {mutagenesis_type}")


def _extract_count_layout(df: pd.DataFrame) -> dict[int, dict[int, str]]:
    """
    Inspect count columns and map each replicate/bin pair to its source column.

    Parameters
    ----------
    df : pd.DataFrame
        Input table containing sortscore count columns named like
        ``count.r1b1``, ``count.r2b4``, and similar.

    Returns
    -------
    dict[int, dict[int, str]]
        Nested mapping of replicate number to bin number to original
        column name, for example ``{1: {1: "count.r1b1", 2: "count.r1b2"}}``.

    Raises
    ------
    ValueError
        If no Lilace-compatible count columns are present or if replicates do
        not all share the same bin layout.
    """
    replicate_bins: dict[int, dict[int, str]] = {}

    for column in df.columns:
        match = COUNT_COLUMN_PATTERN.match(column)
        if match is None:
            continue
        replicate = int(match.group("replicate"))
        bin_number = int(match.group("bin"))
        replicate_bins.setdefault(replicate, {})[bin_number] = column

    if not replicate_bins:
        raise ValueError("No count columns matching 'count.r<replicate>b<bin>' were found.")

    expected_bins = None
    for replicate, bins in replicate_bins.items():
        current_bins = sorted(bins)
        if expected_bins is None:
            expected_bins = current_bins
            continue
        if current_bins != expected_bins:
            raise ValueError(
                "Replicates do not share the same bin layout: "
                f"replicate {replicate} has bins {current_bins}, expected {expected_bins}."
            )

    return replicate_bins


def _get_position_values(
    score_table: pd.DataFrame,
    *,
    variant_column: str,
) -> pd.Series:
    """
    Return per-row amino-acid positions for Lilace export.

    Parameters
    ----------
    score_table : pd.DataFrame
        Input score table that may already contain a ``position`` column.
    variant_column : str
        Column containing variant identifiers from which positions can be parsed.

    Returns
    -------
    pd.Series
        Position values aligned to ``score_table`` rows. Rows without a parsed
        position are returned as missing values.
    """
    if "position" in score_table.columns:
        return pd.to_numeric(score_table["position"], errors="coerce")

    return pd.to_numeric(
        score_table[variant_column].astype(str).str.extract(POSITION_PATTERN)["position"],
        errors="coerce",
    )


def _build_lilace_input_table(
    score_table: pd.DataFrame,
    *,
    mutagenesis_type: Literal["aa", "codon"] = "aa",
    batch: str | None = None,
    batch_column: str = "batch",
    metadata_columns: dict[str, str] | None = None,
    metadata_constants: dict[str, object] | None = None,
) -> pd.DataFrame:
    """
    Convert a sortscore score table into Lilace format.

    Parameters
    ----------
    score_table : pd.DataFrame
        sortscore score table containing a variant identifier, an annotation
        column, count columns named like ``count.r<replicate>b<bin>``, and
        optionally a batch label plus extra metadata columns to preserve.
    mutagenesis_type : {"aa", "codon"}, default "aa"
        Whether to export amino-acid level or codon-level variant identifiers and
        annotations into the Lilace table.
    batch : str | None, default None
        Optional batch or tile label to export, e.g. ``tile4``. When not
        provided, the full input table is used.
    batch_column : str, default "batch"
        Column containing batch or tile labels used to filter the input table
        when ``batch`` is provided.
    metadata_columns : dict[str, str] | None, default None
        Optional mapping from Lilace output metadata column names to source
        columns in ``score_table``.
    metadata_constants : dict[str, object] | None, default None
        Optional constant metadata values to add to every exported row.

    Returns
    -------
    pd.DataFrame
        Lilace-compatible table with one row per variant per replicate and
        columns for the variant id, annotation, replicate number, and
        ``c_0``, ``c_1``, ... count values.

    Raises
    ------
    ValueError
        If required columns are missing, the requested batch is absent, count
        columns are malformed, or the filtered table contains duplicate
        variant identifiers.
    """
    variant_column, annotation_column = _get_mutagenesis_columns(mutagenesis_type)
    required_columns = {variant_column, annotation_column}
    if batch is not None:
        required_columns.add(batch_column)
    if metadata_columns:
        required_columns.update(metadata_columns.values())
    missing_columns = sorted(required_columns - set(score_table.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    filtered = score_table.copy()
    if batch is not None:
        filtered = filtered.loc[filtered[batch_column] == batch].copy()
        if filtered.empty:
            raise ValueError(f"No rows found where {batch_column!r} == {batch!r}.")
    filtered["position"] = _get_position_values(filtered, variant_column=variant_column)
    filtered = filtered.loc[filtered["position"].notna()].copy()
    if filtered.empty:
        raise ValueError(
            f"No rows with usable position values remained after filtering {variant_column}."
        )
    filtered["position"] = filtered["position"].astype(int)

    replicate_bins = _extract_count_layout(filtered)
    duplicate_mask = filtered[variant_column].duplicated(keep=False)
    if duplicate_mask.any():
        duplicated_variants = filtered.loc[duplicate_mask, variant_column].astype(str).drop_duplicates().tolist()
        variants = ", ".join(duplicated_variants[:5])
        raise ValueError(
            f"Duplicate {variant_column!r} values found in Lilace input source rows: {variants}"
        )

    output_rows: list[dict[str, object]] = []
    sorted_replicates = sorted(replicate_bins)
    sorted_bins = sorted(next(iter(replicate_bins.values())))
    output_bin_columns = [f"c_{index}" for index, _ in enumerate(sorted_bins)]
    metadata_columns = metadata_columns or {}
    metadata_output_columns = list(metadata_columns)
    metadata_constants = metadata_constants or {}

    for _, row in filtered.sort_values(variant_column, kind="stable").iterrows():
        for replicate in sorted_replicates:
            output_row: dict[str, object] = {
                "variant_id": row[variant_column],
                "mutation_type": row[annotation_column],
                "position": row["position"],
                "replicate": replicate,
            }
            for index, bin_number in enumerate(sorted_bins):
                output_row[f"c_{index}"] = row[replicate_bins[replicate][bin_number]]
            for metadata_output_column in metadata_output_columns:
                output_row[metadata_output_column] = row[metadata_columns[metadata_output_column]]
            output_row.update(metadata_constants)
            output_rows.append(output_row)

    output = pd.DataFrame(
        output_rows,
        columns=[
            "variant_id",
            "mutation_type",
            "position",
            "replicate",
            *output_bin_columns,
            *metadata_output_columns,
            *metadata_constants.keys(),
        ],
    )
    output.sort_values(["variant_id", "replicate"], inplace=True, kind="stable")
    output.reset_index(drop=True, inplace=True)
    return output


def score_table_to_lilace_csv(
    input_path: str | Path,
    output_path: str | Path,
    *,
    mutagenesis_type: Literal["aa", "codon"] = "aa",
    batch: str | None = None,
    batch_column: str = "batch",
    metadata_columns: dict[str, str] | None = None,
    metadata_constants: dict[str, object] | None = None,
) -> Path:
    """
    Read a sortscore score table CSV and write one Lilace CSV for a selected
    batch.

    Parameters
    ----------
    input_path : str | Path
        Path to the source score table CSV file.
    output_path : str | Path
        Destination path for the Lilace-compatible CSV file.
    mutagenesis_type : {"aa", "codon"}, default "aa"
        Whether to export amino-acid level or codon-level variant identifiers and
        annotations into the Lilace table.
    batch : str | None, default None
        Optional batch or tile label to export from the input file.
    batch_column : str, default "batch"
        Column containing batch or tile labels in the input CSV.
    metadata_columns : dict[str, str] | None, default None
        Optional mapping from Lilace output metadata columns to source columns.
    metadata_constants : dict[str, object] | None, default None
        Optional constant metadata values to add to every output row.

    Returns
    -------
    Path
        The resolved output file path that was written.

    Raises
    ------
    ValueError
        Propagated from :func:`build_lilace_input` when the input data cannot
        be converted into Lilace format.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    batch_scores = pd.read_csv(input_path)
    lilace_input = _build_lilace_input_table(
        batch_scores,
        mutagenesis_type=mutagenesis_type,
        batch=batch,
        batch_column=batch_column,
        metadata_columns=metadata_columns,
        metadata_constants=metadata_constants,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    lilace_input.to_csv(output_path, index=False)
    return output_path
