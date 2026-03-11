# sortscore

`sortscore` is a Python package for Sort-seq variant analysis, including scoring, normalization, and visualization.

## Quick Start

```bash
python -m venv .venv
source .venv/bin/activate
pip install sortscore
```
Run a standard variant scoring analysis:

```bash
sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
```
If the `sortscore` command is not on your `PATH`, run:

```bash
python -m sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
```

## Example Heatmaps

Single Experiment Scoring Heatmap

![Single Experiment Scoring Heatmap](https://raw.githubusercontent.com/dbaldridge-lab/sortscore/main/docs/assets/single_experiment_heatmap.png)

Batch Normalization Heatmap

![Batch Normalization Heatmap](https://raw.githubusercontent.com/dbaldridge-lab/sortscore/main/docs/assets/batch_normalization_heatmap.png)

## Documentation
- [Installation](https://github.com/dbaldridge-lab/sortscore/blob/main/docs/installation.md)
- [Usage](https://github.com/dbaldridge-lab/sortscore/blob/main/docs/usage.md)
- [CLI Arguments](https://github.com/dbaldridge-lab/sortscore/blob/main/docs/cli_arguments.md)
- [Visualization](https://github.com/dbaldridge-lab/sortscore/blob/main/docs/visualization.md)
- [Batch Processing](https://github.com/dbaldridge-lab/sortscore/blob/main/docs/batch_processing.md)
- [Troubleshooting](https://github.com/dbaldridge-lab/sortscore/blob/main/TROUBLESHOOTING.md)
- [Contributing](https://github.com/dbaldridge-lab/sortscore/blob/main/CONTRIBUTING.md)

## Demo

- [Single Experiment Scoring Notebook Demo](https://github.com/dbaldridge-lab/sortscore/blob/main/demo_data/single_experiment_demo.ipynb)
- [Batch Normalization Notebook Demo](https://github.com/dbaldridge-lab/sortscore/blob/main/demo_data/tiled_demo.ipynb)
- [Example Config](https://github.com/dbaldridge-lab/sortscore/blob/main/demo_data/GLI2_oPool5b/config.json)
- [Example Experiment Setup](https://github.com/dbaldridge-lab/sortscore/blob/main/demo_data/GLI2_oPool5b/experiment_setup.csv)
