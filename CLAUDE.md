0;276;0c# CLAUDE.md

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

# With positional averages export for protein structure visualization
sortscore --config path/to/config.json --pos-color

# Run tests
pytest sortscore/analysis/tests/
pytest sortscore/visualization/tests/

# Install in development mode
pip install -e .

# Install dependencies
pip install -r requirements.txt
```

### CLI Arguments vs Config Parameters

**Parameter Organization:**
- **JSON Config File**: Core experimental parameters that define the scientific analysis (reproducible)
  - `experiment_name`, `wt_seq`, `avg_method`, `bins_required`, `minread_threshold`, etc.
- **CLI Arguments**: Runtime control and output formatting (can vary between runs)
  - `--suffix`, `--pos-color`, `--batch`

**Current CLI Arguments:**
- `-c, --config`: Path to config JSON **(required)**
- `-s, --suffix`: Custom suffix for output files
- `-b, --batch`: Enable batch processing mode  
- `-p, --pos-color`: Export positional averages with colors
- `--fig-format`: Figure output format (png, svg, pdf)

### Package Structure Commands
```bash
# Entry point via console script (after pip install)  
sortscore --config config.json

# Batch processing multiple experiments
sortscore --batch --config batch_config.json

# With custom suffix
sortscore --batch --config batch_config.json --suffix custom_name
```

## Batch Normalization

The package supports batch processing to combine and normalize multiple Sort-seq experiments, enabling cross-experiment comparisons through systematic normalization. This is particularly useful for tiled experimental designs where different experiments cover different regions of a protein.

### Normalization Methods

**1. Z-score scaled 2-pole normalization** (default)
- **Step 1**: WT normalization to global reference: `norm1 = raw_score * (global_wt / experiment_wt)`
- **Step 2**: Z-score transformation using synonymous distribution: `norm2 = (norm1 - syn_mean) / syn_std_dev`
- **Step 3**: Pathogenic control normalization: `final = norm2 * (global_pathogenic / experiment_pathogenic)`

Creates standardized scale where synonymous variants center around 0 with unit variance, making cross-experiment comparisons meaningful.

**2. 2-pole normalization**
- Formula: `(b/(a-c))*(A-C)` where:
  - `b` = individual variant score
  - `a` = experiment synonymous median, `c` = experiment pathogenic median
  - `A` = global synonymous median, `C` = global pathogenic median

### Batch Configuration

Batch processing uses a dedicated configuration file format:

```json
{
    "experiment_configs": [
        "/path/to/experiment1/config.json",
        "/path/to/experiment2/config.json",
        "/path/to/experiment3/config.json"
    ],
    "batch_normalization_method": "zscore_threestep",
    "pathogenic_control_type": "nonsense",
    "combined_output_dir": "/path/to/combined/results",
    "global_min_pos": 1,
    "global_max_pos": 500,
    "allow_position_breaks": true,
    "cleanup_individual_files": true
}
```

**Configuration Parameters:**
- `experiment_configs`: List of paths to individual experiment JSON files
- `batch_normalization_method`: "zscore_2pole" (default) or "2pole"
- `pathogenic_control_type`: "nonsense" (default) or "custom" 
- `pathogenic_variants`: Custom pathogenic variants (when using "custom")
- `global_min_pos`/`global_max_pos`: Overall position range for tiled heatmaps
- `allow_position_breaks`: Handle gaps in tiled designs
- `cleanup_individual_files`: Remove individual experiment outputs after combination

### Batch Workflow

1. **Load individual experiments** with their specific configurations
2. **Run individual analyses** to generate raw scores and statistics  
3. **Combine raw data** from all experiments
4. **Calculate global references** (medians, averages) across combined data
5. **Apply normalization method** with all transformation steps
6. **Recalculate statistics** from fully normalized data
7. **Generate tiled heatmaps** with proper position mapping
8. **Save combined results** and cleanup individual files

### Tiled Heatmap Visualization

Batch processing generates combined heatmaps that properly handle:
- **Position mapping**: Each experiment's positions mapped to global coordinates
- **Gap handling**: Visualization of non-contiguous experimental coverage
- **Experiment boundaries**: Visual indicators of experiment transitions
- **Unified scaling**: All data on the same normalized scale

The tiled heatmap automatically adjusts figure size based on the global position range and includes experiment boundary markers for clarity.

## Code Architecture

### Core Components

**ExperimentConfig** (`sortscore/analysis/load_experiment.py`): Central dataclass that loads experiment configuration from JSON and manages count data loading. Key methods:
- `from_json()`: Loads config from JSON file
- `load_counts()`: Populates `counts[rep][bin]` and `median_gfp[rep][bin]` dictionaries
- `annotate_counts()`: Adds sequence difference annotations

**Activity Score Calculation** (`sortscore/analysis/score.py`): 
- `calculate_full_activity_scores()`: Main scoring function implementing Sort-seq logic
- Supports two averaging methods: 'simple-avg', 'rep-weighted'
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

**Experiment JSON**: Contains parameters like `experiment_name`, `bins_required`, `reps_required`, `avg_method`, `minread_threshold`, `wt_seq`, `variant_type`, `min_pos`, `max_pos`, `output_dir`, `experiment_setup_file`, `mutagenesis_variants`, `position_type`

**Experiment Setup CSV**: Maps replicates/bins to count files with columns: `Replicate`, `Bin`, `Path` (or `Read Counts (CSV)`), `Median GFP`, and optionally `Read Count` and `Proportion of Cells`

**Count Files**: TSV/CSV with columns `seq` (variant sequence) and count column (any name, must be second column)

#### Normalization Options

The package supports comprehensive normalization to account for technical variability:

**Base normalization**: All counts are normalized to reads per million to account for sequencing depth differences between samples.

**Optional cell proportion normalization**: If the `Proportion of Cells` column is provided in the experiment setup CSV, the package applies additional normalization to account for different numbers of cells sorted into each quantile gate:

This was added as an opportunity to specify unequal bins. It has not been tested.

- **Formula**: `normalized_reads = (variant_reads / total_reads / cell_proportion) * 1e6`
- **Values**: Cell proportions should be decimal values representing the fraction of viable, singlet cells that fall within each quantile gate
- **Purpose**: Corrects for cases where quantile gates capture different proportions of the sorted cell population (after dead/doublet removal)

**Cell sorter data**: Use the `%Gate` values from your cell sorter output. Convert percentages to decimal values for the CSV.

**Important**: The cell proportions do NOT need to add up to 100%. This is normal for Sort-seq experiments because quantile gates are exclusive and there are gaps between gates.

**Example from real cell sorter data**:
```
Gate    %Gate     Decimal Value for CSV
R4      25.06%    0.2506
R5      24.66%    0.2466  
R6      25.00%    0.2500
R7      23.90%    0.2390
Total:  98.62%    (Normal - doesn't need to sum to 100%)
```

This dual normalization ensures that activity scores reflect true biological differences rather than technical artifacts from sequencing depth or cell sorting variations.

#### Count File Formats

The package supports multiple formats for variant sequences in count files:

1. **Full Sequences**: Complete DNA or amino acid sequences (original format)
   - DNA: Full nucleotide sequences for `variant_type: "dna"`
   - AA: Full amino acid sequences for `variant_type: "aa"`

2. **Pre-annotated Amino Acid Changes** (Auto-detected): Individual amino acid changes in various formats
   - **Single-letter codes**: `M1M`, `R98C`, `P171*`
   - **Three-letter codes**: `Met1Met`, `Arg98Cys`, `Pro171Ter`
   - **HGVS p. notation**: `p.M1M`, `p.Arg98Cys`, `p.Pro171*`
   - **With separators**: `M.1.M`, `R-98-C`, `P_171_*`

**Stop Codon Notation**: Nonsense/stop-gained variants are represented with `*` (asterisk) in the output, regardless of input format (`X`, `Ter`, or `*`). For example, `P171X` becomes `P.171.*` in the `aa_seq_diff` column.

**Auto-detection Features**:
- **Headers**: System automatically detects if files have headers by checking when the second column becomes numeric
- **Variant format**: Automatically detects pre-annotated amino acid changes vs full sequences
- **No configuration needed**: Just use `variant_type: "aa"` in your config file for amino acid data

#### Example Count Files

**Full sequence format** (traditional):
```
seq,count
MGKLIVTAGHLYSLMNDQTDKEVNAKLRGFMCDVIVEVDQFQGVSGFDGMVDTLQDVTLVGAGDGVHQFVLKDGDLVLHFSGHVLSGSTYHLPLSRNVLPNVSVMGRKVVVLMGRNSDKGTLGDLPVHPFPFHGKGVMVTGFNGRDGAAILANLLSRMKGK,1523
MGKLIVTAGHLYSLMNDQTDKEVNAKLRGFMCDVIVEVDQFQGVSGFDGMVDTLQDVTLVGAGDGVHQFVLKDGDLVLHFSGHVLSGSTYHLPLSRNVLPNVSVMGRKVVVLMGRNSDKGTLGDLPVHPFPFHGKGVMVTGFNGRDGAAILANLLSRMKGK,892
```

**Pre-annotated format with headers** (auto-detected):
```
seq,count
M1M,1365835
R98C,137270
P171X,125092
W93Q,124464
```

**Pre-annotated format without headers** (auto-detected):
```
M1M,1365835
R98C,137270
P171X,125092
W93Q,124464
```

**Three-letter format** (auto-detected):
```
seq,count
Met1Met,1365835
Arg98Cys,137270
Pro171Ter,125092
Trp93Gln,124464
```

**HGVS format** (auto-detected):
```
seq,count
p.M1M,1365835
p.Arg98Cys,137270
p.Pro171Ter,125092
p.Trp93Gln,124464
```

## Heatmap Customization

The package supports flexible heatmap generation with customizable axes for different experimental designs.

### X-Axis (Positions)

Controlled by the `position_type` parameter:

- **`"aa"`** (default): Amino acid positions using `min_pos` to `max_pos` range
- **`"dna"`**: DNA nucleotide positions (1 to length of `wt_seq`)

### Y-Axis (Variants)

Controlled by the `mutagenesis_variants` parameter:

- **Default**: All 20 amino acids + stop codon `['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']`
- **Custom list**: Any subset or reordering of variants

### Experimental Designs

#### Standard Deep Mutational Scanning
```json
{
  "position_type": "aa",
  "variant_type": "aa",
  "mutagenesis_variants": ["W", "F", "Y", "P", "M", "I", "L", "V", "A", "G", "C", "S", "T", "Q", "N", "D", "E", "H", "R", "K", "*"]
}
```
- **X-axis**: Amino acid positions (e.g., 1, 2, 3... up to protein length)
- **Y-axis**: All 20 amino acids + stop codon
- **Matrix size**: ~21 × (max_pos - min_pos + 1)

#### GCTA (Single Nucleotide Scanning)
```json
{
  "position_type": "dna",
  "variant_type": "dna", 
  "mutagenesis_variants": ["G", "C", "T", "A"]
}
```
- **X-axis**: DNA positions (1, 2, 3... up to DNA sequence length)
- **Y-axis**: Four DNA bases
- **Matrix size**: 4 × DNA sequence length

#### Hydrophobic Amino Acid Screen
```json
{
  "position_type": "aa",
  "variant_type": "aa",
  "mutagenesis_variants": ["M", "I", "L", "V", "F", "W", "Y", "A"]
}
```
- **X-axis**: Amino acid positions
- **Y-axis**: Only hydrophobic amino acids
- **Matrix size**: 8 × (max_pos - min_pos + 1)

#### Custom DNA Base Subset
```json
{
  "position_type": "dna",
  "variant_type": "dna",
  "mutagenesis_variants": ["G", "A"]
}
```
- **X-axis**: DNA positions  
- **Y-axis**: Only G and A bases
- **Matrix size**: 2 × DNA sequence length

### Impact on Analysis

**Dropout Calculation**: Automatically adjusts based on the actual variants studied:
- Standard DMS: (positions × 21 variants) total possible
- GCTA: (DNA length × 4 bases) total possible  
- Custom: (positions × length of mutagenesis_variants) total possible

**Heatmap Size**: Matrix dimensions scale with your experiment design, improving visualization clarity for focused studies.

### Background Transparency

Heatmap functions support customizable background transparency through the `transparent` parameter:

**Python API**:
```python
from sortscore.visualization.heatmaps import plot_heatmap

# Transparent background (default)
plot_heatmap(data, 'avgscore', experiment, transparent=True)

# White background 
plot_heatmap(data, 'avgscore', experiment, transparent=False)
```

**Default Behavior**: All heatmaps use transparent backgrounds (`transparent=True`) by default, making them suitable for presentations and publications where you want to overlay on different background colors.

**Format Support**: Transparency works with PNG and SVG formats. For PDF output, use `transparent=False` for better compatibility.

### Key Conventions

- Count data stored as nested dictionaries: `counts[replicate][bin] = DataFrame`
- Median GFP values: `median_gfp[replicate][bin] = float`
- Activity scores calculated per replicate, then averaged across replicates
- Filtering requires minimum bins per replicate and minimum replicates per variant
- Sequence parsing supports both DNA ('dna') and amino acid ('aa') variant types
- Heatmap axes automatically adjust to experiment design via `position_type` and `mutagenesis_variants`

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
- Don't edit config files unless asked
- Create a new module if it is over 500 lines
- the python package name is sortscore. Do not use pascal case, only underscores.
- Add to package documentation means add to usage documentation in the python package.