# Usage Guide for sortscore

This guide provides detailed instructions for running Sort-seq variant analysis using the `sortscore` package.

## 1. Prepare Your Experiment Configuration
- Create your experiment configuration JSON (see `config/example_experiment.json`) and experiment setup CSV.
- Edit these files to match your experiment's parameters and data file locations. You can place them anywhere; just provide the correct path when running the analysis.

## 2. Run Analysis from the Command Line
- Use the provided Python API or command-line interface (CLI) to run your analysis.
- **Example CLI usage:**
  ```bash
  python -m sortscore.run_analysis --config path/to/your_config.json
  ```
  - Replace `path/to/your_config.json` with the path to your config file.
  - The config JSON should reference your experiment setup CSV via the `experiment_setup_file` key (with its path).

- **Custom output file naming:**
  ```bash
  python -m sortscore.run_analysis --config path/to/your_config.json --suffix "custom_name"
  ```
  - Use `--suffix` or `-s` to specify a custom suffix for all output files
  - Default: auto-generated from experiment parameters (experiment_name, bins_required, etc.)

## 3. Python API Usage
- You can also import and use the package directly in your own scripts or notebooks:
  ```python
  from sortscore.analysis.load_experiment import ExperimentConfig
  config = ExperimentConfig.from_json('config.json')
  # ...proceed with analysis using the loaded config...
  ```

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
| avg_method            | str     | Method for averaging scores (e.g., 'rep-weighted', 'simple-avg', 'codon-weighted').         |
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
