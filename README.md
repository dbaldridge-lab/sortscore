# sortscore
[Documentation Index](docs/index.md)

## Installation
#TODO: Update this when released on PYPI #1
**Option 1: Using a virtual environment (recommended)**

The following are bash commands:
```bash
# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install sortscore
git clone git+https://github.com/dbaldridge-lab/sortscore.git
cd sortscore
pip install -e .
```

**Option 2: Using conda/anaconda**

```bash
# Activate your conda environment
conda activate your-env-name

# Verify you're using conda's pip
which pip  # Should show path to conda environment

# Install sortscore
git clone git+https://github.com/dbaldridge-lab/sortscore.git
cd sortscore
pip install -e .
```

**Option 3: Install directly (may require adding scripts directory to PATH)**

```bash
# Install sortscore
git clone git+https://github.com/dbaldridge-lab/sortscore.git
cd sortscore
pip install -e .
```

## Usage
### Scoring
With virtual environment or conda environment activated:
```bash
# only required arguements
sortscore -n EXPERIMENT_NAME -e experiment_config.csv -w WT_SEQ
```

Otherwise, add `python -m` before the sortscore command (or `python3 -m`, depending on your python installation).

```bash
python -m sortscore ...
```

### Tile Normalization

```bash
# Including optional -c option and batch config file
sortscore -b -n EXPERIMENT_NAME -e experiment_config.csv -c batch_config.json
```

## Config Templates

[config.json](https://github.com/dbaldridge-lab/sortscore/blob/main/tests/fixtures/GLI2_oPool5b/config.json)
- define additional parameters impacting the analysis execution and outputs

[experiment_setup.csv](https://github.com/dbaldridge-lab/sortscore/blob/main/tests/fixtures/GLI2_oPool5b/experiment_setup.csv)
- define replicates and bins for each sample
- provide parameters for each sample (MFI, cell counts, etc.)

`batch_config.json`

## System Requirements

- Python 3.10+
- See `requirements.txt` for dependencies

## Troubleshooting
If you encounter any issues during installation:

Ensure you're using Python 3.10 or higher: `python --version`
Try updating pip: `pip install --upgrade pip`
For dependency conflicts, consider using a virtual environment `python -m venv .venv && source .venv/bin/activate`.

```mermaid
flowchart TB
    %% Entry Points
    CLI(("CLI: run_analysis.py")):::entry
    API(("API: sortscore")):::entry

    %% Data Stores
    Config(("Config Files<br/>(JSON/CSV)")):::data
    RawData(("Counts Files<br/>(CSV/Parquet)")):::data
    
    %% Pipeline Stages
    subgraph "Pipeline Stages"
        direction TB
        EL["Experiment Loader"]:::pipelineStage
        DP["Data Processing"]:::pipelineStage
        FL["Filtering"]:::pipelineStage
        NM["Normalization"]:::pipelineStage
        AN["Annotation"]:::pipelineStage
        SC["Scoring"]:::pipelineStage
        EX["Export"]:::pipelineStage
        Utils["Utils"]:::pipelineStage
    end

    %% Visualization
    subgraph "Visualization" 
        direction TB
        VH["Heatmap generate_heatmap_matrix()"]:::viz
    end

    %% Outputs
    OF(("Output Figures<br/>(PNG/SVG)")):::entry
    OT(("Output Tables<br/>(CSV/DataFrame)")):::entry

    %% Main Flow
    CLI -->|"run_analysis.main()"| EL
    API -->|"sortscore API"| EL
    EL -->|"load_experiment_data"| DP
    DP -->|"filter_variants"| FL
    FL -->|"normalize_read_depth"| NM
    NM -->|"annotate_sequences"| AN
    AN --> SC
    SC -->|"calculate_full_activity_scores"| EX
    EX -->|"export_results"| VH
  
    VH -->|"generate_heatmap_matrix"| OF
    EX -->|"export_results"| OT

    %% Data Inputs
    Config --> EL
    RawData --> EL

    %% Click Events
    click CLI "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/run_analysis.py"
    click EL "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/load_experiment.py"
    click AN "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/annotation.py"
    click DP "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/data_processing.py"
    click FL "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/filtering.py"
    click NM "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/normalize_read_depth.py"
    click SC "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/score.py"
    click EX "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/export.py"
    click VD "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/visualization/plots.py"
    click VH "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/visualization/heatmap_matrix.py"
    click Utils "https://github.com/dbaldridge-lab/sortscore/blob/main/sortscore/analysis/utils.py"

    %% Styles
    classDef entry fill:#cccccc,stroke:#333,stroke-width:1px
    classDef pipeline fill:#bbeeff,stroke:#333,stroke-width:1px
    classDef pipelineStage fill:#bbeeff,stroke:#333,stroke-width:1px,font-size:22px
    classDef viz fill:#bbeecc,stroke:#333,stroke-width:1px
    classDef data fill:#ffeb99,stroke:#333,stroke-width:1px
    classDef ext fill:none,stroke-dasharray: 5 5,stroke:#666,stroke-width:1px
    classDef test fill:#f9d5e5,stroke:#333,stroke-width:1px
```

## License
MIT
