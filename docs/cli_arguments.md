# Command Line Arguments Reference for sortscore

This document lists all command line arguments available for the `sortscore` tool.

| Argument | Short | Type | Description | Default |
|----------|-------|------|-------------|---------|
| `--experiment-name` | `-n` | str | Experiment name used for output file naming **(required)** | - |
| `--experiment-setup-file` | `-e` | str | Path to experiment setup CSV **(required)** | - |
| `--config` | `-c` | str | Optional experiment config JSON file (used as fallback defaults; CLI takes precedence) | - |
| `--wt-seq` | `-w` | str | Wild-type sequence (required unless provided in `--config`) | - |
| `--output-dir` | `-o` | str | Output directory | `.` |
| `--bins-required` | - | int | Minimum number of bins required | 1 |
| `--reps-required` | - | int | Minimum number of replicates required | 1 |
| `--avg-method` | - | str | Averaging method: `rep-weighted` or `simple-avg` | rep-weighted |
| `--minread-threshold` | - | int | Minimum read threshold | 0 |
| `--max-cv` | - | float | Maximum coefficient of variation allowed | None |
| `--mutagenesis-variants` | - | str | Comma-separated list (e.g. `G,C,T,A`) | W,F,Y,P,M,I,L,V,A,G,C,S,T,Q,N,D,E,H,R,K,* |
| `--position-offset` | - | int | Offset for position numbering | 0 |
| `--biophysical-prop` | - | bool | Show biophysical properties panel in heatmaps | False |
| `--position-type` | - | str | Position axis for plots: `aa` or `dna` | aa |
| `--min-pos` | - | int | Minimum position (1-based) | 1 |
| `--max-pos` | - | int | Maximum position (1-based) | None |
| `--suffix` | `-s` | str | Custom suffix for all output files | (auto: current date) |
| `--batch` | `-b` | flag | Enable batch processing mode | False |
| `--pos-color` | `-p` | flag | See [visualization.md](visualization.md) for exporting positional averages with colors for protein structure visualization | False |
| `--fig-format` | - | str | Output format for figures: png, svg, pdf | png |
