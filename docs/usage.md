# Usage Guide for sortscore

This guide provides detailed instructions for running Sort-seq variant analysis using the `sortscore` package.

## 1. Prepare Your Experiment Configuration
- Create your experiment configuration JSON (see `config/example_experiment.json`) and experiment setup CSV.
- Edit these files to match your experiment's parameters and data file locations. You can place them anywhere; just provide the correct path when running the analysis.

## 2. Command Line Interface (CLI)

### Basic Usage

```bash
# Standard analysis
python -m sortscore.run_analysis --config path/to/config.json

# Or using the console (after pip install)
sortscore --config path/to/config.json
```

### CLI Arguments Reference

| Argument | Short | Type | Description | Default |
|----------|-------|------|-------------|---------|
| `--config` | `-c` | str | Path to experiment config JSON file **(required)** | - |
| `--suffix` | `-s` | str | Custom suffix for all output files | Current date (YYYYMMDD) |
| `--batch` | `-b` | flag | Enable batch processing mode | False |
| `--pos-color` | `-p` | flag | Export positional averages with colors for protein structure visualization | False |
| `--fig-format` | - | str | Output format for figures: png, svg, pdf | png |

### Parameter Organization

**Experiment Configuration (JSON file):**
- Core parameters that define the experiments
- Examples: `experiment_name`, `wt_seq`, `avg_method`, `minread_threshold`, `bins_required`

**Runtime Control (CLI arguments):**
- Parameters that control how analysis runs and output formatting
- Examples: `--suffix`, `--pos-color`, `--batch`

**Export Options:**
- `--pos-color` / `-p`: Generates `{experiment_name}_positional_averages.csv` with position, average score, and hex color columns for protein structure visualization

### File Naming Conventions

**Current Implementation:** 
```
# Score files (consistent with suffix)
{experiment_name}_scores_{suffix}.csv
stats_{avg_method}_{suffix}.json

# Visualization files (mixed patterns - being standardized)
{experiment_name}_heatmap.png              # New API (no suffix yet)
{experiment_name}_heatmap_matrix.csv       # New API (no suffix yet) 
{experiment_name}_positional_averages.csv # New API (no suffix yet)
aa_heatmap_{avg_method}_{suffix}.png       # Legacy naming
```

**Auto-generated suffix format:** `{experiment_name}_{bins}bins_{minreads}minreads_{max_cv}cv_{date}`

### Examples

```bash
# Basic analysis
sortscore -c my_experiment.json

# With custom output suffix
sortscore -c my_experiment.json -s "final_analysis" 

# Export positional averages for protein structure visualization
sortscore -c my_experiment.json -p

# Generate SVG figures
sortscore -c my_experiment.json --fig-format svg

# Batch processing multiple experiments
sortscore -b -c batch_config.json
```

## Batch Processing (Multiple Experiments)

For tiled experimental designs or combining multiple experiments, use batch processing mode:

```bash
# Run batch analysis
sortscore --batch --config path/to/batch_config.json

# With custom suffix
sortscore --batch --config path/to/batch_config.json --suffix "combined_analysis"
```

### Batch Configuration File

Create a separate JSON configuration file for batch processing:

```json
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
```

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

The main configuration file (JSON) defines all parameters for your Sort-seq analysis. Below are the standard keys and their meanings:

| Key                   | Type    | Description                                                                                 |
|-----------------------|---------|---------------------------------------------------------------------------------------------|
| experiment_name       | str     | Name/ID of the experiment or submission.                                                    |
| bins_required         | int     | Minimum number of bins per replicate a variant must appear in to be scored.                               |
| reps_required         | int     | Minimum number of replicates a variant must appear in to be scored.                         |
| avg_method            | str     | Method for averaging scores (e.g., 'rep-weighted', 'simple-avg').         |
| minread_threshold     | int     | Minimum reads per bin for a variant to be scored.                                |
| max_cv                | float   | Maximum coefficient of variation (CV) allowed across replicates. Variants exceeding this are filtered out. |
| read_count            | list    | List of demultiplexed read counts for each sample/bin.                                      |
| output_dir            | str     | Directory where all results and figures will be saved. Default value is the current directory.                                       |
| experiment_setup_file | str     | Path to the experiment setup CSV file (see below).                                          |
| wt_seq                | str     | Wild-type sequence (DNA) for the region analyzed.                                     |
| variant_type          | str     | Type of variant ('aa' for amino acid, 'dna' for nucleotide/codon).                           |
| min_pos               | int     | Starting amino acid sequence position for the region analyzed (used for plot labels).            |
| max_pos               | int     | Ending amino acid sequence position for the region analyzed (used for plot labels).              |
| mutagenesis_variants  | list    | Custom list of variants for heatmap y-axis. Default: all 20 AAs + stop codon.               |
| position_type         | str     | Position type for heatmap x-axis ('aa' for amino acid positions, 'dna' for DNA positions). Default: 'aa'. |

See also the [experiment setup CSV reference](#experiment-setup-csv-reference) for details on the CSV file format.

## Experiment Setup CSV Reference

The experiment setup CSV must contain the following columns:
- `Replicate`: Replicate number (integer)
- `Bin`: Bin number (integer)
- `Read Counts (CSV)`: Path to the variant count file for this replicate/bin
- `Median GFP`: Median GFP value for this replicate/bin

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

**X-Axis (Positions)** - Controlled by `position_type`:
- `"aa"` (default): Amino acid positions using `min_pos` to `max_pos` range
- `"dna"`: DNA nucleotide positions (1 to length of `wt_seq`)

**Y-Axis (Variants)** - Controlled by `mutagenesis_variants`:
- Default: All 20 amino acids + stop codon `["W", "F", "Y", "P", "M", "I", "L", "V", "A", "G", "C", "S", "T", "Q", "N", "D", "E", "H", "R", "K", "*"]`
- Custom: Any subset or reordering of variants

### Example Configurations

#### Standard Deep Mutational Scanning (Default)
```json
{
  "experiment_name": "MyProtein_DMS",
  "position_type": "aa",
  "variant_type": "aa",
  "min_pos": 1,
  "max_pos": 150,
  "mutagenesis_variants": ["W", "F", "Y", "P", "M", "I", "L", "V", "A", "G", "C", "S", "T", "Q", "N", "D", "E", "H", "R", "K", "*"]
}
```
- **X-axis**: Amino acid positions (1, 2, 3... up to 150)
- **Y-axis**: All 20 amino acids + stop codon
- **Matrix size**: 21 × 150

#### GCTA (Single Nucleotide Scanning)
```json
{
  "experiment_name": "MyGene_GCTA",
  "position_type": "dna",
  "variant_type": "dna", 
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
  "position_type": "aa",
  "variant_type": "aa",
  "min_pos": 50,
  "max_pos": 100,
  "mutagenesis_variants": ["M", "I", "L", "V", "F", "W", "Y", "A"]
}
```
- **X-axis**: Amino acid positions (50 to 100)
- **Y-axis**: Only hydrophobic amino acids
- **Matrix size**: 8 × 51

#### Custom DNA Base Subset
```json
{
  "experiment_name": "GC_Content_Study",
  "position_type": "dna",
  "variant_type": "dna",
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
