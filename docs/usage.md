# Usage Guide for sortscore

This guide provides detailed instructions for running Sort-seq variant analysis using the `sortscore` package.

For visualization options, see [visualization.md](visualization.md).


## 1. Confirm Count File Formatting (TSV or Parquet)

The first column of input variant count files must contain the variant sequences. 

The second column must contain the unique counts for each variant.

Column names and the presence of a header row are flexible. Only this ordering of columns is strictly required. Any additional annotations may be listed in the remaining columns, which will be ignored by `sortscore`.

Example:

| seq      | count |
|----------|-------|
| ATGCGT...|  123  |
| GCTTAA...|   45  |

Note: `sortscore` expects full variant sequences (DNA or protein) in count files. Work to include additional variant nomenclatures is ongoing.

## 2. Prepare Experiment Configuration File
- Create your (see `config/experimental_setup.csv`) experiment setup CSV.
- Edit this file to match your experiment's parameters and data file locations. Relative paths resolve relative to `experimental_setup.csv`.
- 
### Experiment Setup Expected Columns

The experiment setup CSV must contain the following required columns:
- `Replicate`: Technical replicate number. Set to 1 if your experiment doesn't have technical replicates.
- `Bin`: Bin number. Label from 1 up to the number of tubes that cells were sorted into. The order does not matter, so long as MFI is mapped to the correct bin in this file.
- `Read Counts (CSV)`: Path to the variant count file for this replicate/bin.
- `MFI`: Median fluorescence value for this replicate/bin.

Optional setup columns:
- `Tile`: Used for normalization across experiments in tiled/batch workflows. Not required for standard single-experiment scoring.

## 3. Running the Scoring Workflow

The entry point for running the Sort-seq scoring workflow is the `sortscore` command. The required arguments are listed below and must be provided in the bash shell (terminal) command line.

**Required CLI Arguments**

| Key                   | Type  | Description |
|------------------------|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| `experiment_name`      | str   | Name/ID of the experiment or submission. |
| `experiment_setup_file`| str   | Path to the experiment setup CSV file (see below). |
| `wt_seq`               | str   | Wild-type reference sequence for the region analyzed. Use **DNA** sequence when `mutagenesis_type` is `codon` or `snv`; use **protein** sequence when `mutagenesis_type` is `aa`. |

Additional optional fields can be used to customize the analysis. These can be selected by providing a JSON configuration with the `-c` option, or through CLI flags. If a parameter is provided both in the CLI and the config file, the CLI value takes precedence.

**Optional Fields (CLI or Config File)**

| Key                   | Type  | Description |
|------------------------|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| `bins_required`        | int   | Minimum number of bins per replicate a variant must appear in to be scored. **Default:** `1`. |
| `reps_required`        | int   | Minimum number of replicates a variant must appear in to be scored. **Default:** `1`. |
| `avg_method`           | str   | Method for averaging scores (e.g., `rep-weighted`, `simple-avg`). **Default:** `rep-weighted`. |
| `minread_threshold`    | int   | Minimum reads per bin for a variant to be scored. **Default:** `0`. |
| `max_cv`               | float | Maximum coefficient of variation (CV) allowed across replicates. Variants exceeding this are filtered out. |
| `mutagenesis_type`     | str   | Mutagenesis type: `aa`, `codon`, or `snv`. **Default:** `aa`. Set this in config or CLI when running DNA-based analysis (`codon`/`snv`). |
| `read_count`           | list  | List of demultiplexed read counts for each sample/bin. |
| `output_dir`           | str   | Directory where all results and figures will be saved. **Default:** `.`. |
| `mutagenesis_variants` | list  | Custom list of variants for heatmap y-axis. **Default:** `all 20 AAs plus stop codon`. |
| `min_pos`              | int   | Minimum position (1-based). **Default:** `1`. Interpreted as amino acid or DNA position depending on `wt_seq` and the variant sequences provided in the count files. |
| `max_pos`              | int   | Maximum position (1-based). **Default:** `1 + length of wild-type sequence`. Interpreted as amino acid or DNA position depending on `wt_seq` and the variant sequences provided in the count files. |

Note:
Any relative file paths specified in the experiment setup file are resolved relative to the location of the setup file itself, not the current working directory.

`wt_seq` format requirement by mutagenesis type:
- `mutagenesis_type: codon` or `snv` -> provide DNA `wt_seq`
- `mutagenesis_type: aa` -> provide protein `wt_seq`

### Basic Usage

```bash
# After installation (recommended)
sortscore score -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json

# Without installation (from project root)
python -m sortscore score -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
```


### CLI Arguments Reference

See [docs/cli_arguments.md](cli_arguments.md) for a complete list of command line arguments.

### Parameter Organization

**CLI (highest precedence):**
- Required: `-n/--experiment-name`, `-e/--experiment-setup-file`
- Optional: pass any analysis parameters to override the config file.

**Optional JSON config (fallback defaults):**
- Provide `-c/--config` to set defaults like `wt_seq`, thresholds, and plotting options.

**Export Options:**


### File Naming Conventions

**Current Implementation:**
```
# Score files
scores/{experiment_name}_dna_scores_{suffix}.csv
scores/{experiment_name}_aa_scores_{suffix}.csv
# TODO: #31 add SNV functionality and test
scores/{experiment_name}_dna_scores_snv_{suffix}.csv

# Summary statistics
scores/{experiment_name}_dna_stats_{suffix}.json         # when DNA scores are produced
scores/{experiment_name}_aa_stats_{suffix}.json

# For a list of all visualization output files, see [visualization.md](visualization.md).
```

**Auto-generated suffix format:** `YYYYMMDD` (current date)

### Examples

```bash
# Basic analysis
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json
```
```bash
# DNA codon-level analysis
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json --mutagenesis-type codon
```
```bash
# With custom output suffix
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json -s "final_analysis"
```
```bash
# Generate SVG figures
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json --fig-format svg
```
```bash
# Batch processing multiple experiments
sortscore norm -c batch_config.json
```

## Tiled Mutagenesis Batch Processing

For tiled experimental designs where different sequencing datasets cover different regions of the same protein, sortscore supports automatic batch processing with cross-tile normalization:

```bash
# Run tiled batch analysis (requires a `tile` column in the setup CSV)
sortscore norm -c batch_config.json


```

### Tiled Experimental Setup

Add a `tile` column to your experiment setup CSV to indicate which sequence region each count file represents:

```csv
Replicate,Bin,Read Counts (CSV),MFI,tile
1,1,tile1_rep1_bin1.tsv,1000,N_terminal_1-50
1,2,tile1_rep1_bin2.tsv,2000,N_terminal_1-50
1,3,tile1_rep1_bin3.tsv,3000,N_terminal_1-50
2,1,tile2_rep1_bin1.tsv,1500,middle_domain_51-100
2,2,tile2_rep1_bin2.tsv,2500,middle_domain_51-100
2,3,tile2_rep1_bin3.tsv,3500,middle_domain_51-100
3,1,tile3_rep1_bin1.tsv,1800,C_terminal_101-150
3,2,tile3_rep1_bin2.tsv,2800,C_terminal_101-150
3,3,tile3_rep1_bin3.tsv,3800,C_terminal_101-150
```

### Automatic Tiled Workflow

When sortscore detects a `tile` column in your experiment setup CSV, it automatically processes each sequence tile separately and then combines them using the selected normalization method and controls. Each tile is analyzed independently to generate raw activity scores, then all tiles are combined to allow comparisons of scores across the entire protein sequence.

The system generates unified tiled heatmaps that properly map each tile's positions to the global protein coordinate system, with visual indicators showing tile boundaries and handling any gaps between non-contiguous sequence regions. Individual tile outputs are automatically cleaned up after successful combination, leaving only the final combined results. Intermediate files for individual experiments can be optionally retained.

### Batch Configuration Parameters

| Key                        | Type      | Description                                                              |
|---------------------------|-----------|--------------------------------------------------------------------------|
| batch_normalization_method| str       | "zscore_2pole" (default), "2pole", or "zscore_center"                |
| pathogenic_control_type   | str       | "nonsense" (default) or "custom"                                       |
| pathogenic_variants       | list      | Custom pathogenic variants (required when using "custom")              |
| combined_output_dir       | str       | Directory for final combined results                                     |
| global_min_pos            | int       | Overall minimum position across all experiments (updates tiled heatmaps axis)    |
| global_max_pos            | int       | Overall maximum position across all experiments (updates tiled heatmaps axis)    |

### Batch Normalization Methods

**1. Z-score scaled 2-pole normalization** (recommended, default)
- Creates standardized scale where synonymous variants center around 0 with unit variance
- Enables meaningful cross-experiment comparisons
- Process:
  1. WT normalization: `norm1 = raw_score * (global_wt / experiment_wt)`
  2. Z-score transformation: `norm2 = (norm1 - syn_mean) / syn_std_dev`  
  3. Pathogenic control normalization: `final = norm2 * (global_pathogenic / experiment_pathogenic)`

**2. 2-pole normalization**
- Uses synonymous and pathogenic variants as reference points
- Formula: `(b/(a-c))*(A-C)` where:
  - `b` = individual variant score
  - `a` = experiment synonymous median, `c` = experiment pathogenic median
  - `A` = global synonymous median, `C` = global pathogenic median

**3. Z-score centering normalization** 
- WT-only normalization
- Process:
  1. WT normalization: `norm1 = raw_score * (global_reference / experiment_reference)`
  2. Z-score transformation: `final = (norm1 - syn_mean) / syn_std_dev`
- For DNA variants: uses `wt_dna` scores as reference (fallback to synonymous)
- For AA variants: uses synonymous variants as reference
- Use when pathogenic controls are unavailable

### Batch Processing Workflow

1. Individual Analysis: Each experiment analyzed separately to generate raw scores
2. Data Combination: Raw scores and statistics combined across experiments  
3. Global Calculations: Reference values computed from combined data
4. Normalization: Selected method applied to standardize scores
5. Statistics Recalculation: Final statistics computed from normalized data
6. Visualization: Combined tiled heatmaps generated with position mapping
7. Output: Combined results saved, individual files cleaned up (if requested)

### Mutagenesis Type Selection

`sortscore` uses `mutagenesis_type` to determine analysis mode:
- `aa` (default)
- `codon`
- `snv` (not yet implemented; see open issues)

Set it either in your config JSON (`"mutagenesis_type": "codon"`) or on the CLI with `--mutagenesis-type codon`.


### Position Numbering Convention
All positions are relative to the provided `wt_seq` unless otherwise specified:

- **Default**: Position 1 corresponds to the first character of `wt_seq`
- **Pre-annotated data**: If count files contain position annotations (e.g., "K.2.E"),
  those positions are used as-is, with optional offset adjustment applied
- **Sequence flexibility**: `wt_seq` can be DNA or amino acid sequence, but analysis mode is controlled by `mutagenesis_type`

### Common Mutagenesis Libraries

#### Amino Acid Analysis
- **Input**: DNA sequences (aggregates synonymous variants) OR AA sequences (direct processing)
- **Output**: AA substitution heatmap, AA-level statistics, AA scores file
- **Use for**: Deep mutational scanning, protein function studies

#### Codon-Level Analysis  
- **Input**: DNA sequences (required)
- **Output**: Codon heatmap, DNA scores file, codon variance quantification, synonymous vs non-synonymous analysis
- **Use for**: Codon optimization studies, synonymous variant effects, quantifying codon-level variance
