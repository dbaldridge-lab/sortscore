# sortscore

`sortscore` is a Python package for Sort-seq variant analysis, including scoring, normalization, and visualization.

## Quick Start

`sortscore` is a Python package for Sort-seq variant analysis, including scoring, normalization, and visualization.

## Quick Start

```bash
git clone https://github.com/dbaldridge-lab/sortscore
git clone https://github.com/dbaldridge-lab/sortscore
cd sortscore
python -m venv .venv
source .venv/bin/activate
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Run a standard analysis:
Run a standard analysis:

```bash
sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
```

If you did not install the CLI entry point, run:
If you did not install the CLI entry point, run:

```bash
python -m sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
python -m sortscore -n EXPERIMENT_NAME -e path/to/experiment_setup.csv -c path/to/config.json
```

## Documentation
- [Installation](docs/installation.md)
- [Usage](docs/usage.md)
- [CLI Arguments](docs/cli_arguments.md)
- [Visualization](docs/visualization.md)
- [Batch Processing](docs/batch_processing.md)
- [Troubleshooting](TROUBLESHOOTING.md)
- [Contributing](CONTRIBUTING.md)

## Demo

- [Notebook demo](demo_data/GLI2_oPool5b/demo.ipynb)
- [Example config](demo_data/GLI2_oPool5b/config.json)
- [Example experiment setup](demo_data/GLI2_oPool5b/experiment_setup.csv)
