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
   from sortscore.visualization.plots import plot_activity_score_distribution

   # Load your data
   df = pd.read_csv('your_data.csv')

   # Calculate activity scores
   scores = calculate_activity_scores([df], method='simple-avg')

   # Plot distribution
   plot_activity_score_distribution(scores, score_col='avgscore')
   ```

## Further Documentation
- For detailed installation, configuration, and usage instructions, see the [docs/](docs/) folder:
  - [Installation Guide](docs/installation.md)
  - [Usage Guide](docs/usage.md)
  - [Documentation Index](docs/index.md)
- See `config/example_experiment.json` and `config/experiment_setup.csv` for configuration file templates.
- All public functions use NumPy-style docstrings; see module docstrings and examples for API details.

## System Requirements
- Python 3.11+
- See `requirements.txt` for dependencies.

## License
MIT