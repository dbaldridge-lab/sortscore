# sortscore

## Quickstart

1. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Example usage**

   ```python
   import pandas as pd
   from sortscore.analysis.score import calculate_activity_scores
   from sortscore.visualization.plots import plot_activity_score_distribution, plot_beeswarm, plot_histogram, plot_heatmap

   # Load your data
   df = pd.read_csv('your_data.csv')

   # Calculate activity scores
   scores = calculate_activity_scores([df], method='simple-avg')

   # Plot distribution
   plot_activity_score_distribution(scores, score_col='avgscore')

   # Beeswarm plot
   plot_beeswarm(scores, x='annotate_aa', y='avgscore')

   # Histogram by annotation
   plot_histogram(scores, score_col='avgscore', group_col='annotate_dna')

   # Heatmap (requires additional columns and config)
   plot_heatmap(
       data=scores,
       score_col='avgscore',
       num_aa=10,
       wt_seq='ABCDEFGHIJ',
       min_pos=1,
       mutant_type='aa',
       wt_score=1.0,
       fig_size='small',
       export=False,
       output=None,
       tick_values=[0.5, 1.0, 2.0],
       tick_labels=['Low', 'WT', 'High'],
       motif_indices=None,
       row_avg=False,
       title='Example Heatmap'
   )
   ```

3. **Configuration**
   - Place experiment parameters in a JSON or YAML file in the `config/` directory.
   - See `config/example_experiment.json` for structure.

4. **Testing**
   ```bash
   pytest sortscore/sortscore/visualization/tests/
   ```

## Documentation
- All public functions use NumPy-style docstrings.
- See module docstrings and examples for details.

## System Requirements
- Python 3.11+
- See `requirements.txt` for dependencies.

## License
MIT

## Input Requirements

- All count files listed in your experiment setup CSV must have barcodes already mapped to their correct sequences. The Sort-seq pipeline does not perform barcode-to-sequence mapping; it assumes all input files are pre-processed in this way.
- See `config/example_experiment.json` and `config/experiment_setup.csv` for configuration and file manifest examples.

## How to Use sortscore

### 1. Prepare Your Experiment Configuration
- Create an experiment configuration file (see `sortscore/config/example_experiment.json`) and experiment setup file (see `sortscore/config/experiment_setup.csv`).
- Edit these files to match your experiment's parameters and data file locations. 

### 2. Run Analysis from the Command Line
- Use the provided Python API or command-line interface (CLI) to run your analysis.
- **Example CLI usage:**
  ```bash
  python -m sortscore.run_analysis --config path/to/your_config.json
  ```
  - Replace `path/to/your_config.json` with the path to your config file.
  - The config JSON should reference your experiment setup CSV via the `experiment_setup_file` key (provide full path).

### 3. Python API Usage
- You can also import and use the package directly in your own scripts or notebooks:
  ```python
  from sortscore.analysis.config import get_submission_config
  config = get_submission_config('your_submission_name')
  # ...proceed with analysis using the loaded config...
  ```

### 4. Output
- Results and plots will be saved to the `output_dir` specified in your config.

---

For more details, see the docstrings in each module and the example configuration files in the `config/` directory.