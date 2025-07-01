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
- Python 3.8+
- See `requirements.txt` for dependencies.

## License
MIT