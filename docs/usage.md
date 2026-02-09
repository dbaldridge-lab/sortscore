# Usage Guide for sortscore

This guide provides detailed instructions for running Sort-seq variant analysis using the `sortscore` package.

## 1. Prepare Your Experiment Configuration
- Create your experiment configuration JSON (see `config/example_experiment.json`) and experiment setup CSV.
- Edit these files to match your experiment's parameters and data file locations. You can place them anywhere; just provide the correct path when running the analysis.
  - `variant_type` is not a config option; it is auto-detected from your count files. Current supported mutagenesis schemes are AA, codon, and SNV.

### Variant Type Auto-Detection (Accepted Nomenclature)

`sortscore` auto-detects whether your inputs are **DNA variants** or **amino-acid (AA) variants** by sampling the **first column** of the count files referenced by your experiment setup CSV.

**Count file expectations**
- The **first column** is treated as the variant identifier/sequence.
- The detector samples variants for classification.

**Accepted DNA variant formats**
# TODO: test this
- Full DNA sequences using `A`, `T`, `C`, `G` (ambiguous IUPAC bases like `N`, `R`, `Y`, etc. are also treated as DNA).
- Simple substitutions: `A123T` (also accepts `A.123.T`, `A-123-T`, `A_123_T`).
- HGVS DNA notation like `c.123A>T`, `g.123A>T`, `n.123A>T` (supported when `mavehgvs` is installed).

**Accepted AA variant formats**
- Full AA sequences using the 20 one-letter amino acids plus `*` (stop).
- One-letter substitutions: `M1V`, `R98*` (also accepts `M.1.V`, `M-1-V`, `M_1_V`).
- Three-letter substitutions: `Ala123Val`, and `Ter` for stop.
- HGVS protein notation like `p.Arg97Cys` (supported when `mavehgvs` is installed).

## 2. Command Line Interface (CLI)

### Basic Usage

```bash
# After installation (recommended)
sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json

# Without installation (from project root)
python -m sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
```

### CLI Arguments Reference

| Argument | Short | Type | Description | Default |
|----------|-------|------|-------------|---------|
| `--experiment-name` | `-n` | str | Experiment name used for output file naming **(required unless `--batch`)** | - |
| `--experiment-setup-file` | `-e` | str | Path to experiment setup CSV **(required unless `--batch`)** | - |
| `--config` | `-c` | str | Optional experiment config JSON file (used as fallback defaults; CLI takes precedence) | - |
| `--wt-seq` | `-w` | str | Wild-type sequence (required unless provided in `--config`) | - |
| `--output-dir` | `-o` | str | Output directory | `.` |
| `--bins-required` | - | int | Minimum number of bins required | from config |
| `--reps-required` | - | int | Minimum number of replicates required | from config |
| `--avg-method` | - | str | Averaging method: `rep-weighted` or `simple-avg` | from config |
| `--minread-threshold` | - | int | Minimum read threshold | from config |
| `--max-cv` | - | float | Maximum coefficient of variation allowed | from config |
| `--mutagenesis-variants` | - | str | Comma-separated list (e.g. `G,C,T,A`) | from config |
| `--position-offset` | - | int | Offset for position numbering | from config |
| `--biophysical-prop` / `--no-biophysical-prop` | - | bool | Show biophysical properties panel in heatmaps | from config |
| `--position-type` | - | str | Position axis for plots: `aa` or `dna` | auto |
| `--min-pos` | - | int | Minimum position (1-based) | from config / auto |
| `--max-pos` | - | int | Maximum position (1-based) | from config / auto |
| `--suffix` | `-s` | str | Custom suffix for all output files | Current date (YYYYMMDD) |
| `--batch` | `-b` | flag | Enable batch processing mode | False |
| `--pos-color` | `-p` | flag | Export positional averages with colors for protein structure visualization | False |
| `--fig-format` | - | str | Output format for figures: png, svg, pdf | png |

### Parameter Organization

**CLI (highest precedence):**
- Required: `-n/--experiment-name`, `-e/--experiment-setup-file`
- Optional: pass any analysis parameters to override the config file.

**Optional JSON config (fallback defaults):**
- Provide `-c/--config` to set defaults like `wt_seq`, thresholds, and plotting options.

**Export Options:**
- `--pos-color` / `-p`: Generates `{experiment_name}_positional_averages.csv` with position, average score, and hex color columns for protein structure visualization

### File Naming Conventions

**Current Implementation:**
```
# Score files
scores/{experiment_name}_dna_scores_{suffix}.csv
scores/{experiment_name}_aa_scores_{suffix}.csv
# TODO: (later) add SNV functionality and test
scores/{experiment_name}_dna_scores_snv_{suffix}.csv

# Summary statistics
scores/{experiment_name}_dna_stats_{suffix}.json         # when DNA scores are produced
scores/{experiment_name}_aa_stats_{suffix}.json

# Visualization files
figures/{experiment_name}_aa_heatmap_{suffix}.{png|svg|pdf}
figures/{experiment_name}_aa_heatmap_matrix_{suffix}.csv
figures/{experiment_name}_codon_heatmap_{suffix}.{png|svg|pdf}         # when plotting DNA-level scores
figures/{experiment_name}_codon_heatmap_matrix_{suffix}.csv            # when plotting DNA-level scores
figures/{experiment_name}_positional_averages_{suffix}.csv   # when using --pos-color
```

**Auto-generated suffix format:** `YYYYMMDD` (current date)

### Examples

```bash
# Basic analysis
sortscore -n my_experiment -e experiment_setup.csv -c my_experiment.json

# With custom output suffix
sortscore -n my_experiment -e experiment_setup.csv -c my_experiment.json -s "final_analysis"

# Export positional averages for protein structure visualization
sortscore -n my_experiment -e experiment_setup.csv -c my_experiment.json -p

# Generate SVG figures
sortscore -n my_experiment -e experiment_setup.csv -c my_experiment.json --fig-format svg

# Batch processing multiple experiments
sortscore -b -c batch_config.json
```

## Tiled Mutagenesis Batch Processing

For tiled experimental designs where different sequencing datasets cover different regions of the same protein, sortscore supports automatic batch processing with cross-tile normalization:

```bash
# Run tiled analysis with automatic tile detection
sortscore --config experiment_config.json

# With custom suffix for combined results
sortscore --config experiment_config.json --suffix "tiled_analysis"
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
| experiment_configs        | list      | Paths to individual experiment JSON configuration files                   |
| batch_normalization_method| str       | "zscore_2pole" (default), "2pole", or "zscore_center"                |
| pathogenic_control_type   | str       | "nonsense" (default) or "custom"                                       |
| pathogenic_variants       | list      | Custom pathogenic variants (required when using "custom")              |
| combined_output_dir       | str       | Directory for final combined results                                     |
| global_min_pos            | int       | Overall minimum position across all experiments (for tiled heatmaps)    |
| global_max_pos            | int       | Overall maximum position across all experiments (for tiled heatmaps)    |
| allow_position_breaks     | bool      | Whether to allow gaps/breaks in tiled position display (default: true) |
| cleanup_individual_files  | bool      | Remove individual experiment outputs after combination (default: true)  |

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

1. **Individual Analysis**: Each experiment analyzed separately to generate raw scores
2. **Data Combination**: Raw scores and statistics combined across experiments  
3. **Global Calculations**: Reference values computed from combined data
4. **Normalization**: Selected method applied to standardize scores
5. **Statistics Recalculation**: Final statistics computed from normalized data
6. **Visualization**: Combined tiled heatmaps generated with position mapping
7. **Output**: Combined results saved, individual files cleaned up (if requested)

### Tiled Heatmap Features

Batch processing automatically generates tiled heatmaps with:
- **Position mapping**: Experiments mapped to global coordinate system
- **Gap handling**: Visualization of non-contiguous coverage
- **Boundary markers**: Visual indicators of experiment transitions  
- **Unified scaling**: All data on same normalized scale
- **Automatic sizing**: Figure dimensions adjust to position range

## 3. Python API Usage

### Basic Analysis
You can import and use the package directly in your own scripts or notebooks:

```python
from sortscore.analysis.load_experiment import ExperimentConfig
config = ExperimentConfig.from_json('config.json')
# ...proceed with analysis using the loaded config...
```

### Heatmap Visualization Customization

The package provides direct control over heatmap background transparency through the Python API:

```python
from sortscore.visualization.heatmaps import plot_heatmap

# Load your data and experiment config
config = ExperimentConfig.from_json('config.json')
# ... load and process data ...

# Generate heatmap with transparent background (default)
plot_heatmap(
    data=scores_df,
    score_col='avgscore', 
    experiment=config,
    transparent=True,          # Transparent background
    fig_format='png',
    export_heatmap=True,
    output='output/figures'
)

# Generate heatmap with white background
plot_heatmap(
    data=scores_df,
    score_col='avgscore',
    experiment=config, 
    transparent=False,         # White background
    fig_format='pdf',          # PDF works better with white backgrounds
    export_heatmap=True,
    output='output/figures'
)
```

**Transparency Options:**
- `transparent=True` (default): Creates transparent background, ideal for presentations and overlays
- `transparent=False`: Creates white background, better for PDF output and traditional publications

**Format Recommendations:**
- PNG and SVG: Work excellently with transparent backgrounds
- PDF: Use `transparent=False` for better compatibility with PDF viewers

## 4. Output
- Results and plots will be saved to the `output_dir` specified in your config.

## 5. Input Requirements
- All count files listed in your experiment setup CSV must have barcodes already mapped to their correct sequences. The Sort-seq pipeline does not perform barcode-to-sequence mapping; it assumes all input files are pre-processed in this way.
- See `config/example_experiment.json` and `config/experiment_setup.csv` for configuration and file manifest examples.

## Experiment Configuration JSON Reference
Note:
Any relative file paths specified in the experiment setup file are resolved relative to the location of the setup file itself, not the current working directory.

The main configuration file (JSON) defines all parameters for your Sort-seq analysis. Below are the standard keys and their meanings:
# TODO: Update this with latest config updates
| Key                   | Type    | Required | Description                                                                                 |
|-----------------------|---------|----------|---------------------------------------------------------------------------------------------|
| **Required Parameters** | | | |
| experiment_name       | str     | Yes      | Name/ID of the experiment or submission.                                                    |
| experiment_setup_file | str     | Yes      | Path to the experiment setup CSV file (see below).                                          |
| wt_seq                | str     | Yes      | Wild-type reference sequence (DNA or amino acid) for the region analyzed.                                     |
| **Optional Parameters** | | | |
| bins_required         | int     | No       | Minimum number of bins per replicate a variant must appear in to be scored. Default: 1.                               |
| reps_required         | int     | No       | Minimum number of replicates a variant must appear in to be scored. Default: 1.                         |
| avg_method            | str     | No       | Method for averaging scores (e.g., 'rep-weighted', 'simple-avg'). Default: 'rep-weighted'.         |
| minread_threshold     | int     | No       | Minimum reads per bin for a variant to be scored. Default: 0.                                |
| max_cv                | float   | No       | Maximum coefficient of variation (CV) allowed across replicates. Variants exceeding this are filtered out. |
| read_count            | list    | No       | List of demultiplexed read counts for each sample/bin.                                      |
| output_dir            | str     | No       | Directory where all results and figures will be saved. Default: current directory.                                       |
| mutagenesis_variants  | list    | No       | Custom list of variants for heatmap y-axis. Default: all 20 AAs + stop codon.               |
| position_offset       | int     | No       | What position 1 of wt_seq becomes (e.g., position_offset=50 means position 1 → 50). Default: 0. |

See also the [experiment setup CSV reference](#experiment-setup-csv-reference) for details on the CSV file format.

### Automatic Variant Format Detection

The system automatically detects the input variant format from your count files.

**DNA Formats Detected:**
- Full DNA sequences: `ATGCGTAAC...`
- Nucleotide changes: `A123T`, `G.252.C`, `C_45_T`
- HGVS DNA notation: `c.123A>T`

**Amino Acid Formats Detected:**
- Full AA sequences: `MKILVAGD...`
- Single-letter changes: `M1V`, `R.98.C`, `P171*`
- Three-letter changes: `Met1Val`, `Arg98Cys`, `Pro171Ter`
- HGVS protein notation: `p.Met1Val`, `p.Arg98Cys`


### Position Numbering Convention

All positions are relative to the provided `wt_seq` unless otherwise specified:

- **Default**: Position 1 corresponds to the first character of `wt_seq`
- **With position_offset**: Position 1 of `wt_seq` becomes the offset value
  (e.g., `position_offset=50` means position 1 → position 50, position 2 → position 51)
- **Pre-annotated data**: If count files contain position annotations (e.g., "K.2.E"),
  those positions are used as-is, with optional offset adjustment applied
- **Sequence flexibility**: `wt_seq` can be DNA or amino acid sequence - the system auto-detects the type and handles accordingly

**Examples**:
```json
// Default: positions relative to wt_seq
{"wt_seq": "MKVLIVAG", "position_offset": 0}
// Position 1 = M, Position 2 = K, etc.

// With offset: position 1 of wt_seq becomes position 50  
{"wt_seq": "MKVLIVAG", "position_offset": 50}
// Position 1 → 50, Position 2 → 51, Position 8 → 57
```

### Common Mutagenesis Libraries

#### Amino Acid Analysis
- **Input**: DNA sequences (aggregates synonymous variants) OR AA sequences (direct processing)
- **Output**: AA substitution heatmap, AA-level statistics, AA scores file
- **Use for**: Deep mutational scanning, protein function studies

#### Codon-Level Analysis  
- **Input**: DNA sequences (required)
- **Output**: Codon heatmap, DNA scores file, codon variance quantification, synonymous vs non-synonymous analysis
- **Use for**: Codon optimization studies, synonymous variant effects, quantifying codon-level variance

#### Single Nucleotide Variant Analysis (in development)
- **Input**: DNA sequences (required)  
- **Output**: Position-by-nucleotide heatmap, SNV-specific statistics
- **Use for**: Saturation genome editing (SGE), base editing, nucleotide-level screens


## Experiment Setup CSV Reference

The experiment setup CSV must contain the following columns:
- `Replicate`: Replicate number (integer)
- `Bin`: Bin number (integer)
- `Read Counts (CSV)`: Path to the variant count file for this replicate/bin
- `MFI`: Median fluorescence value for this replicate/bin

### Input Variant Count File Format
Each input variant count file must:
- Have the first column named `seq` containing the variant sequences.
- Have the second column containing the unique counts for each variant (column name can be anything, but must be the second column).

Example (TSV):
```
seq	count
ATGCGT...	123
GCTTAA...	45
...
```

The pipeline assumes all input files are pre-processed and formatted as above.

## Heatmap Customization

The package supports flexible heatmap generation with customizable axes for different experimental designs.

### Controlling Heatmap Axes

# TODO: check if this is still valid after updating config docs
**X-Axis (Positions)** - Controlled by `position_type`:
- `"aa"` (default): Amino acid positions using `min_pos` to `max_pos` range
# TODO: make position offset apply to DNA plots
- `"dna"`: DNA nucleotide positions (1 to length of `wt_seq`)

**Y-Axis (Variants)** - Controlled by `mutagenesis_variants`:
# TODO: test that changing this doesn't break the codon heatmap
- Default: All 20 amino acids + stop codon `["W", "F", "Y", "P", "M", "I", "L", "V", "A", "G", "C", "S", "T", "Q", "N", "D", "E", "H", "R", "K", "*"]`
- Custom: Any subset or reordering of variants

### Example Configurations

#### Standard Deep Mutational Scanning (Default)
```json
{
  "experiment_name": "MyProtein_DMS",
  "wt_seq": "MKVLIVAG...",
  "mutagenesis_variants": ["W", "F", "Y", "P", "M", "I", "L", "V", "A", "G", "C", "S", "T", "Q", "N", "D", "E", "H", "R", "K", "*"]
}
```
- **X-axis**: Amino acid positions (auto-detected from data)
- **Y-axis**: All 20 amino acids + stop codon
- **Matrix size**: 21 × (detected position range)

#### GCTA (Single Nucleotide Scanning)
```json
{
  "experiment_name": "MyGene_GCTA",
  "mutagenesis_variants": ["G", "C", "T", "A"],
  "wt_seq": "ATGCGTAAC..."
}
```
- **X-axis**: DNA positions (1, 2, 3... up to DNA sequence length)
- **Y-axis**: Four DNA bases (G, C, T, A)
- **Matrix size**: 4 × DNA sequence length

#### Hydrophobic Amino Acid Screen
```json
{
  "experiment_name": "Hydrophobic_Screen",
  "wt_seq": "MKVLIVAG...",
  "position_offset": 50,
  "mutagenesis_variants": ["M", "I", "L", "V", "F", "W", "Y", "A"]
}
```
- **X-axis**: Amino acid positions (auto-detected, with offset so position 1 becomes 50)
- **Y-axis**: Only hydrophobic amino acids
- **Matrix size**: 8 × (detected position range)

#### Custom DNA Base Subset
```json
{
  "experiment_name": "GC_Content_Study",
  "mutagenesis_variants": ["G", "C"],
  "wt_seq": "ATGCGTAAC..."
}
```
- **X-axis**: DNA positions
- **Y-axis**: Only G and C bases
- **Matrix size**: 2 × DNA sequence length

### Impact on Analysis

- **Dropout Calculation**: Automatically adjusts based on actual variants studied

---
For more details, see the docstrings in each module and the example configuration files in the `config/` directory.
