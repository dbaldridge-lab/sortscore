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
| read_count            | list    | List of demultiplexed read counts for each sample/bin.                                      |
| output_dir            | str     | Directory where all results and figures will be saved. Default value is the current directory.                                       |
| experiment_setup_file | str     | Path to the experiment setup CSV file (see below).                                          |
| wt_seq                | str     | Wild-type sequence (DNA) for the region analyzed.                                     |
| variant_type          | str     | Type of variant ('aa' for amino acid, 'dna' for nucleotide/codon).                           |
| min_pos               | int     | Starting amino acid sequence position for the region analyzed (used for plot labels).            |
| max_pos               | int     | Ending amino acid sequence position for the region analyzed (used for plot labels).              |

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

---
For more details, see the docstrings in each module and the example configuration files in the `config/` directory.
