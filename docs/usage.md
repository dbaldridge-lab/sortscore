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
  from sortscore.analysis.config import get_submission_config
  config = get_submission_config('your_submission_name')
  # ...proceed with analysis using the loaded config...
  ```

## 4. Output
- Results and plots will be saved to the `output_dir` specified in your config.

## 5. Input Requirements
- All count files listed in your experiment setup CSV must have barcodes already mapped to their correct sequences. The Sort-seq pipeline does not perform barcode-to-sequence mapping; it assumes all input files are pre-processed in this way.
- See `config/example_experiment.json` and `config/experiment_setup.csv` for configuration and file manifest examples.

---

For more details, see the docstrings in each module and the example configuration files in the `config/` directory.
