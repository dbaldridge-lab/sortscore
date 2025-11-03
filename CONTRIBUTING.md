# Contributing to sortscore

This guide provides instructions for testing the package with example fixtures.

## Setup for Testing

### 1. Clone the Repository

```bash
git clone https://github.com/dbaldridge-lab/sortscore.git
cd sortscore
```

**All paths below are relative to this `sortscore/` directory (the repository root).**

### 2. Create a Virtual Environment (Recommended)

```bash
python -m venv sortscore-env
source sortscore-env/bin/activate  # On Windows: sortscore-env\Scripts\activate
```

This ensures the `sortscore` command will be in your PATH.

### 3. Install the Package

```bash
# Install the package from the local directory
pip install .
```

This installs `sortscore` like any PyPI package, making the `sortscore` command available in your virtual environment.

### 4. Verify Installation

```bash
# Check that the command is available
sortscore --help

# Verify you can run as a module
python -m sortscore --help
```

**Troubleshooting:** If you get `command not found: sortscore`:

1. Check the package is installed:
   ```bash
   pip show sortscore
   ```

2. Find where scripts are installed:
   ```bash
   python -c "import sysconfig; print(sysconfig.get_path('scripts'))"
   ```

3. Make sure that directory is in your PATH, or use the full path to the script

4. **Recommended workaround:** Use `python -m sortscore` instead:
   ```bash
   python -m sortscore --config tests/fixtures/GLI2_oPool5b/config.json
   ```

## Running Tests with Fixtures

Test fixtures are included in the repository under `tests/fixtures/`. You can use these to verify the package works correctly.

### Command Line Testing

```bash
# From the repository root directory (sortscore/)
sortscore --config tests/fixtures/GLI2_oPool5b/config.json
```

This will run the analysis on the GLI2 example dataset and generate outputs in the fixture directory.

### Jupyter Notebook Testing

Create a new notebook in the repository root and use these commands:

```python
# Cell 1: Verify you're in the correct directory
import os
print("Current directory:", os.getcwd())
# Should show: /path/to/sortscore

# Cell 2: Check installation
```
!sortscore --help
```

# Cell 3: Run analysis on test fixtures
```
!sortscore --config tests/fixtures/GLI2_oPool5b/config.json
```

# Cell 4: Use the Python API
```
from pathlib import Path
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.score import calculate_full_activity_scores
```

# Load experiment configuration
```
config_path = Path("tests/fixtures/GLI2_oPool5b/config.json")
experiment = ExperimentConfig.from_json(config_path)
```

# Load count data
```
experiment.load_counts()
```

# Calculate activity scores
```
results = calculate_full_activity_scores(experiment)

print(f"Analyzed {len(results)} variants")
print("\nFirst few results:")
print(results.head())
```

**Console Commands:**
```bash
# Navigate to repository
cd /path/to/sortscore
```
```
# Verify current directory
pwd
```
```
# Check installation
sortscore --help
```

```
# Run analysis
sortscore --config tests/fixtures/GLI2_oPool5b/config.json
```

## Running Unit Tests

If you want to run the automated test suite:

```bash
# Install testing dependencies
pip install pytest
# Run all tests
pytest
```
```bash
# Run tests for a specific module
pytest sortscore/analysis/tests/
pytest sortscore/visualization/tests/
```

## Available Test Fixtures

The repository includes example datasets in `tests/fixtures/`:

- `GLI2_oPool5b/` - Example Sort-seq experiment for GLI2 TF DNA binding domain
  - Config file: `tests/fixtures/GLI2_oPool5b/config.json`
  - Includes compressed count files (.tsv.gz) and experiment setup CSV

### Making Changes to Code
If you need to modify the package code and test changes without reinstalling:

```bash
# Install in editable/development mode
pip install -e .
```

This creates a link to the source code, so any changes are immediately available without reinstalling.

## Getting Help

If you encounter issues:
1. Check that you're in the correct directory (`pwd` should show `.../sortscore`)
2. Verify the package is installed (`pip show sortscore`)
3. Make sure test fixtures exist (`ls tests/fixtures/GLI2_oPool5b/`)
4. Report issues at https://github.com/dbaldridge-lab/sortscore/issues