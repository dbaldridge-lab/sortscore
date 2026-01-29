# sortscore

## Quickstart

## Scoring Analysis Pipeline
### Configuration
- See `config.json` and `experiment_setup.csv` for configuration file templates.
### Usage
```bash
sortscore -c config.json
```
### Example usage

## Tile Normalization Pipeline
### Configuration
- See `batch_config.json` for configuration file template.
### Usage
```bash
sortscore -cb batch_config.json
```
### Example usage

## Python API
- All public functions use NumPy-style docstrings. See module docstrings and examples for API details.
### TODO: test these example workflows
### Loading Counts
   ```python
      import pandas as pd
      # Load your data
      df = pd.read_csv('your_data.csv')
   ```

### Calculating Activity scores
   ```python
   from sortscore.analysis.score import calculate_activity_scores

   # Calculate activity scores. Average over replicates using a simple average, instead of default averaging weighted by variant read count.
   scores = calculate_activity_scores([df], method='simple-avg')
```

### Plotting
#### Histogram
   ```python
   from sortscore.visualization.plots import plot_activity_score_distribution
   plot_activity_score_distribution(scores, score_col='avgscore')
   ```

## Further Documentation
- For detailed installation, configuration, and usage instructions, see the [docs/](docs/) folder:
  - [Installation Guide](docs/installation.md)
  - [Usage Guide](docs/usage.md)
  - [Documentation Index](docs/index.md)

## System Requirements
- Python 3.10+
- Bash shell

## License
MIT