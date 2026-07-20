# Command Line Arguments Reference for sortscore

This document lists command line arguments for the scoring command:

```bash
sortscore score [options]
```

| Argument | Short | Type | Description | Default |
|----------|-------|------|-------------|---------|
| `--experiment-name` | `-n` | str | Experiment name used for output file naming **(required)** | - |
| `--experiment-setup-file` | `-e` | str | Path to experiment setup CSV **(required)** | - |
| `--config` | `-c` | str | Optional experiment config JSON file (used as fallback defaults; CLI takes precedence) | - |
| `--wt-seq` | `-w` | str | Wild-type sequence (required unless provided in `--config`). Use DNA WT for `--mutagenesis-type codon`; use protein WT for `--mutagenesis-type aa`. | - |
| `--output-dir` | `-o` | str | Output directory. Relative values are resolved from the current working directory. | `.` |
| `--bins-required` | - | int | Minimum number of bins required | 1 |
| `--reps-required` | - | int | Minimum number of replicates required | 1 |
| `--avg-method` | - | str | Averaging method: `rep-weighted` or `simple-avg` | rep-weighted |
| `--minread-threshold` | - | int | Minimum read threshold | 0 |
| `--max-cv` | - | float | Maximum coefficient of variation allowed | None |
| `--mutagenesis-type` | - | str | Mutagenesis type: `aa` or `codon` | aa |
| `--mutagenesis-variants` | - | str | Comma-separated list (e.g. `G,C,T,A`) | W,F,Y,P,M,I,L,V,A,G,C,S,T,Q,N,D,E,H,R,K,* |
| `--position-offset` | - | int | Offset for position numbering | 0 |
| `--biophysical-prop` | - | bool | Show biophysical properties panel in heatmaps | False |
| `--min-pos` | - | int | Minimum position (1-based) | 1 |
| `--max-pos` | - | int | Maximum position (1-based) | None |
| `--relative-path-base` | - | str | Base directory for resolving relative paths when this setting is applied: `setup` or `cwd` | setup |
| `--pos-color` | `-p` | flag | See [visualization.md](visualization.md) for exporting positional averages with colors for protein structure visualization | False |
| `--fig-format` | - | str | Output format for figures: png, svg, pdf | png |
Run cross-tile normalization with:

```bash
sortscore norm -c batch_config.json
sortscore norm -c batch_config.json --method zscore_2pole
sortscore norm -c batch_config.json --output-dir /path/to/combined
```

Use the helper command for related scoring tools:

```bash
sortscore integrate lilace --input path/to/aa_scores.csv --output Lilace_input_counts.csv
sortscore integrate lilace --input path/to/batch_aa_scores.csv --output Lilace_input_counts.csv --batch tile4
sortscore integrate lilace --input path/to/batch_dna_scores.csv --output Lilace_input_counts.csv --mutagenesis-type codon --batch tile4
```

`lilace` writes `variant_id`, `mutation_type`, `position`, `replicate`, and bin count columns. `--batch` is optional, and `--mutagenesis-type` defaults to `aa`.

Path resolution note:
- `--output-dir` (CLI) uses the current working directory as the base for relative paths.
- `output_dir` in JSON config files is resolved relative to the config file directory.
