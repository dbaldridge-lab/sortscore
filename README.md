# sortscore
[Documentation Index](docs/index.md)

## Quickstart

## Installation
Installation can be done via pip to pull code directly from GitHub.
```
pip install git+https://git@github.com/idptools/sparrow.git
```

**Example usage**

[Add info on how to run standard pipeline.]

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

## Config Templates
- See `config/example_experiment.json` and `config/experiment_setup.csv` for configuration file templates.

## API
- All public functions use NumPy-style docstrings; see module docstrings and examples for API details.

## System Requirements
- Python 3.11+
- See `requirements.txt` for dependencies.

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
