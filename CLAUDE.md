# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`sortscore` is a modular Python package for Sort-seq variant analysis that calculates activity scores from sequencing count data. The package has been refactored from project-specific "oPool" and "MAVE activity analysis" to generalized Sort-seq analysis. It processes variant count files across multiple replicates and bins, normalizes the data, and computes weighted activity scores using different averaging methods.

**Key Design Principles:**
- Modular architecture with logical separation (analysis, visualization, utilities)
- Self-contained plotting utilities independent of external notebooks
- User-configurable file paths (no hardcoded locations)
- Generalized for any Sort-seq experiment (not project-specific)

## Common Commands

### Running Analysis
```bash
# For testing pipeline, activate venv and use console script
source venv/bin/activate
sortscore --config path/to/config.json

# Alternative: Main analysis command
python -m sortscore.run_analysis --config path/to/config.json

# Run tests
pytest sortscore/analysis/tests/
pytest sortscore/visualization/tests/

# Install in development mode
pip install -e .

# Install dependencies
pip install -r requirements.txt
```

### Package Structure Commands
```bash
# Entry point via console script (after pip install)
sortscore --config config.json
```

## Code Architecture

### Core Components

**ExperimentConfig** (`sortscore/analysis/load_experiment.py`): Central dataclass that loads experiment configuration from JSON and manages count data loading. Key methods:
- `from_json()`: Loads config from JSON file
- `load_counts()`: Populates `counts[rep][bin]` and `median_gfp[rep][bin]` dictionaries
- `annotate_counts()`: Adds sequence difference annotations

**Activity Score Calculation** (`sortscore/analysis/score.py`): 
- `calculate_full_activity_scores()`: Main scoring function implementing Sort-seq logic
- Supports three averaging methods: 'simple-avg', 'rep-weighted', 'codon-weighted'
- Processes count normalization, bin proportions, and replicate aggregation

**Visualization** (`sortscore/visualization/plots.py`):
- Plotting functions for activity score distributions, heatmaps, and beeswarms
- Integrates with `heatmap_matrix.py` for MAVE matrix visualization

### Data Flow

1. **Configuration**: JSON config file references experiment setup CSV
2. **Data Loading**: ExperimentConfig loads count files specified in CSV  
3. **Scoring**: Count data → normalized counts → bin proportions → activity scores
4. **Visualization**: Score data → plots and heatmaps

### Configuration Files

**Experiment JSON**: Contains parameters like `experiment_name`, `bins_required`, `reps_required`, `avg_method`, `minread_threshold`, `wt_seq`, `variant_type`, `min_pos`, `max_pos`, `output_dir`, `experiment_setup_file`

**Experiment Setup CSV**: Maps replicates/bins to count files with columns: `Replicate`, `Bin`, `Read Counts (CSV)`, `Median GFP`

**Count Files**: TSV/CSV with columns `seq` (variant sequence) and count column (any name, must be second column)

### Key Conventions

- Count data stored as nested dictionaries: `counts[replicate][bin] = DataFrame`
- Median GFP values: `median_gfp[replicate][bin] = float`
- Activity scores calculated per replicate, then averaged across replicates
- Filtering requires minimum bins per replicate and minimum replicates per variant
- Sequence parsing supports both DNA ('dna') and amino acid ('aa') variant types

### Testing

- Test files in `sortscore/analysis/tests/` and `sortscore/visualization/tests/`
- Example/test data in `tests/data/` for automated testing
- No pytest configuration file - uses pytest defaults
- Test command: `pytest sortscore/analysis/tests/` for analysis tests

## Development Guidelines

**Package Requirements:**
- Python 3.11+ (enforced in setup.py, requirements.txt, and docs)
- All dependencies listed in requirements.txt

**Code Organization:**
- All package source code in `sortscore/` directory
- Non-package files (documentation, notebooks, papers) outside package source
- Configuration files user-defined (no hardcoded paths)
- Self-contained modules with minimal cross-dependencies

**Coding Standards:**
- Follow PEP 8 for code style
- Use type hints for all function arguments and return values
- Use NumPy-style docstrings for all public functions, classes, and modules
- Include Parameters, Returns, and Examples sections in docstrings where applicable
- Prefer explicit imports over wildcard imports
- Use logging instead of print statements for status messages

**Version Control:**
- Make atomic, logically separated Git commits for each significant change
- Write clear, descriptive commit messages summarizing purpose and scope
- Avoid combining unrelated changes in single commits
- Ignore unnecessary files in .gitignore (.ipynb_checkpoints, __pycache__, .env)
- Commit configuration templates, not user-specific configs

**Future Considerations:**
- Package is ready for PyPI distribution
- Modular design supports easy feature additions
- Generalized codebase eliminates project-specific references

## Checkpoint Management

**Automatic Checkpoint Updates:**
Claude should periodically update `sortscore_checkpoint_summary.md` (condensed version) when making significant changes to the codebase. This provides a context-window-friendly summary of project status.

**Update Triggers:**
- After completing major refactoring tasks
- When adding new modules or significant features  
- Before/after structural changes to package architecture
- When test coverage or configuration changes significantly

**Checkpoint Content:**
- Current module count and key capabilities
- Architecture overview (brief)
- Recent significant changes
- Readiness status for production/distribution

The full detailed checkpoint remains in `sortscore_refactor_checkpoint.md` for reference.

## Git Workflow

**Periodic Commits:**
When working on extended tasks or making multiple related changes, Claude should suggest periodic git commits to preserve progress. This workflow requires user approval:

1. **Commit Triggers:**
   - After completing logical units of work (single feature, bug fix, refactor)
   - Before starting major structural changes
   - When switching between different types of tasks
   - At natural breakpoints in multi-step processes

2. **Commit Process:**
   - Show `git status` and `git diff` to display all changes
   - Propose a descriptive commit message following project conventions
   - Wait for user approval of both the diff and commit message
   - Only proceed with commit after explicit user confirmation

3. **Commit Message Style:**
   - Follow existing repository patterns (check `git log` for examples)
   - Use imperative mood ("Add feature" not "Added feature")
   - Be concise but descriptive
   - Include scope when relevant (e.g., "analysis:", "viz:", "tests:")

**Never commit without explicit user approval of the changes and message.**

## Data Handling Guidelines

**File Operations:**
- Never unzip or decompress files without explicit user permission
- Preserve compressed data formats (.gz, .bz2, etc.) in their original state
- Use appropriate tools that can read compressed files directly when needed
- Always ask before modifying file formats or compression states

**Test Data:**
- Test data in `tests/` directory should remain unchanged unless explicitly requested
- Compressed count files (e.g., `.tsv.gz`) should stay compressed
- Use pandas or other tools that can read compressed files directly: `pd.read_csv('file.tsv.gz', sep='\t')`

## Memory Notes

- Don't unzip files, use zcat to read them as needed