# Batch Processing Guide

This guide provides detailed instructions for combining and normalizing multiple Sort-seq experiments using the batch processing functionality in sortscore.

## Overview

Batch processing enables you to:
- **Combine multiple experiments** covering different regions of the same protein (tiled design)
- **Apply cross-experiment normalization** to enable meaningful comparisons
- **Generate unified visualizations** with proper position mapping
- **Standardize scoring scales** across experimental batches

## When to Use Batch Processing

- **Tiled experiments**: Multiple experiments covering non-overlapping regions
- **Technical replicates**: Same experimental conditions run at different times
- **Cross-platform comparisons**: Normalizing experiments run on different instruments
- **Large protein analysis**: Experiments divided by protein domains or regions

## Quick Start

1. **Run individual experiments** using standard sortscore analysis
2. **Create batch configuration** JSON file listing all experiments
3. **Run batch analysis**:
   ```bash
   sortscore --batch --config batch_config.json
   ```

## Batch Configuration

### Basic Configuration

```json
{
    "experiment_configs": [
        "/path/to/experiment1/config.json",
        "/path/to/experiment2/config.json",
        "/path/to/experiment3/config.json"
    ],
    "combined_output_dir": "/path/to/combined_results"
}
```

### Full Configuration with All Options

```json
{
    "experiment_configs": [
        "/path/to/experiment1/config.json",
        "/path/to/experiment2/config.json", 
        "/path/to/experiment3/config.json"
    ],
    "batch_normalization_method": "zscore_2pole",
    "pathogenic_control_type": "nonsense",
    "pathogenic_variants": null,
    "combined_output_dir": "/path/to/combined_results",
    "global_min_pos": 1,
    "global_max_pos": 500,
    "allow_position_breaks": true,
    "cleanup_individual_files": true
}
```

### Configuration Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `experiment_configs` | list | **Required** | Paths to individual experiment JSON files |
| `batch_normalization_method` | string | `"zscore_2pole"` | Normalization method: `"zscore_2pole"`, `"2pole"`, or `"zscore_center"` |
| `pathogenic_control_type` | string | `"nonsense"` | Type of pathogenic control: `"nonsense"` or `"custom"` |
| `pathogenic_variants` | list | `null` | Custom pathogenic variants (required when `pathogenic_control_type` is `"custom"`) |
| `combined_output_dir` | string | `"."` | Output directory for combined results |
| `global_min_pos` | int | `null` | Global minimum position for tiled heatmaps (auto-calculated if not provided) |
| `global_max_pos` | int | `null` | Global maximum position for tiled heatmaps (auto-calculated if not provided) |
| `allow_position_breaks` | bool | `true` | Allow gaps in tiled heatmap visualization |
| `cleanup_individual_files` | bool | `true` | Remove individual experiment outputs after combining |

## Normalization Methods

### Z-score Scaled 2-Pole Normalization (Recommended)

This is the default method that creates a standardized scale making cross-experiment comparisons meaningful.

**Mathematical Steps:**
1. **WT Normalization**: `norm1 = raw_score × (global_wt_score / experiment_wt_score)`
2. **Z-score Transform**: `norm2 = (norm1 - synonymous_mean) / synonymous_std_dev`  
3. **Pathogenic Normalization**: `final_score = norm2 × (global_pathogenic_avg / experiment_pathogenic_avg)`

**Result**: Synonymous variants center around 0 with unit variance across all experiments.

**When to use**: 
- Default choice for most applications
- When you want standardized effect sizes
- Cross-experiment statistical comparisons needed

### 2-Pole Normalization from BRCA1 2018 DMS analysis

Uses synonymous and pathogenic variants as reference anchor points.

**Formula**: `(b/(a-c)) × (A-C)`

Where:
- `b` = individual variant score
- `a` = experiment synonymous median
- `c` = experiment pathogenic median  
- `A` = global synonymous median
- `C` = global pathogenic median

**When to use**:
- When you need to preserve relative score ranges
- Compatibility with published normalization approaches
- Direct interpretation of fold-changes relative to controls

### Z-score Centering Normalization

WT-only normalization without requiring pathogenic controls.

**Mathematical Steps:**
1. **WT Normalization**: `norm1 = raw_score × (global_reference_score / experiment_reference_score)`
2. **Z-score Transform**: `final_score = (norm1 - synonymous_mean) / synonymous_std_dev`

**Reference Selection**:
- For DNA variants: Uses `wt_dna` scores as reference (falls back to synonymous if unavailable)
- For AA variants: Uses synonymous variants as reference

**Result**: Creates standardized scale centered on synonymous variants with unit variance.

**When to use**:
- When pathogenic controls are unavailable or not needed
- For experiments focused on beneficial variants
- When you want z-score standardization without pathogenic anchoring

## Pathogenic Controls

### Nonsense Variants (Default)

Uses stop codon variants (`*`) as pathogenic controls.

```json
{
    "pathogenic_control_type": "nonsense"
}
```

**Requirements**: Your experiments must include stop codon variants in the data.

### Custom Pathogenic Variants

Specify your own list of pathogenic control variants.

```json
{
    "pathogenic_control_type": "custom",
    "pathogenic_variants": ["R98C", "P171X", "W93Q"]
}
```

**Use cases**:
- Known deleterious mutations
- Disease-associated variants
- Specific experimental controls

## Tiled Heatmap Visualization

Batch processing automatically generates combined heatmaps with advanced features:

### Position Mapping

Experiment position index are mapped to a global position:
```
Experiment 1: positions 1-100  → Positions 1-100 in combined heatmap
Experiment 2: positions 1-100 with 100 start position offset specified → Positions 101-200 in combined heatmap
```

### Gap Handling

When `allow_position_breaks: true` (default):
- Only positions with data are shown
- Gaps between experiments are visually clear
- Figure width adjusts to actual coverage

When `allow_position_breaks: false`:
- Full global range shown (including gaps)
- Missing positions displayed as empty/white
- Consistent position spacing maintained

### Visual Features

- **Experiment boundaries**: Vertical lines mark transitions between experiments
- **Unified scaling**: Single colorbar applies to all data
- **Automatic sizing**: Figure dimensions scale with position range
- **Biophysical properties**: Optional amino acid property panel

## Output Files

Batch processing generates:

### Score Files
- `batch_scores_[suffix].csv`: Combined normalized scores with batch tracking
- `batch_stats_[suffix].json`: Comprehensive statistics and normalization factors

### Visualizations  
- `tiled_heatmap_[suffix].png`: Combined heatmap with position mapping
- `tiled_heatmap_[suffix]_matrix.csv`: Underlying matrix data for heatmap

### Batch Tracking

The combined scores file includes a `batch` column identifying source experiments:
```csv
seq,avgscore,batch,Rep1.score,Rep2.score,...
M1A,150,experiment1,145,155,...
R98C,-200,experiment2,-195,-205,...
```

## Common Use Cases

### Tiled Protein Analysis

For large proteins analyzed in overlapping or non-overlapping segments:

```json
{
    "experiment_configs": [
        "domain1_analysis/config.json",    // positions 1-150
        "domain2_analysis/config.json",    // positions 151-300
        "domain3_analysis/config.json"     // positions 301-450
    ],
    "global_min_pos": 1,
    "global_max_pos": 450,
    "allow_position_breaks": false
}
```

### Multi-Batch Comparison

For experiments run at different times requiring normalization:

```json
{
    "experiment_configs": [
        "batch1/config.json",
        "batch2/config.json", 
        "batch3/config.json"
    ],
    "batch_normalization_method": "zscore_2pole",
    "cleanup_individual_files": false  // Keep individual results
}
```

### Custom Controls Analysis

Using specific variants as pathogenic controls:

```json
{
    "pathogenic_control_type": "custom",
    "pathogenic_variants": ["R175H", "R248Q", "R273H"],  // Known cancer mutations
    "batch_normalization_method": "2pole"
}
```

## Troubleshooting

### Common Issues

**Missing pathogenic controls**:
```
Error: No pathogenic variants found for normalization
```
- Solution: Ensure your experiments include the specified pathogenic controls
- Check `pathogenic_control_type` and `pathogenic_variants` settings

**Position range conflicts**:  
```
Warning: Experiment position ranges overlap
```
- Solution: Check `min_pos` and `max_pos` in individual experiment configs
- Set appropriate `global_min_pos` and `global_max_pos` in batch config

**Insufficient normalization references**:
```
Error: Cannot calculate global synonymous median - no synonymous variants found
```
- Solution: Ensure experiments include synonymous variants for normalization
- Check that `annotate_dna` or `annotate_aa` columns are properly generated

### Validation

Before running batch analysis, validate your configuration:

```python
from sortscore.analysis.batch_config import BatchConfig

# Load and validate config
config = BatchConfig.from_json('batch_config.json')
config.validate_config()
```

### Performance Tips

- **Large datasets**: Consider using `cleanup_individual_files: true` to save disk space
- **Memory usage**: Process experiments with similar position ranges together
- **Visualization**: Use `allow_position_breaks: true` for better figure clarity with sparse coverage

## API Usage

For programmatic access to batch processing:

```python
from sortscore.analysis.batch_config import BatchConfig
from sortscore.analysis.batch_normalization import run_batch_analysis

# Load configuration
config = BatchConfig.from_json('batch_config.json')
config_dict = config.get_batch_config_dict()

# Run batch analysis
results = run_batch_analysis(config_dict)

# Access results
normalized_scores = results['normalized_scores']
combined_stats = results['combined_stats']
```

For more examples and API details, see the docstrings in the batch processing modules.