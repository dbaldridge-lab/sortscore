# Batch Processing Guide

This guide provides detailed instructions for combining and normalizing multiple Sort-seq experiments using the batch processing functionality in `sortscore`.

## Overview

Batch processing enables you to:
- **Combine multiple experiments** covering different regions of the same protein (tiled design)
- **Apply cross-experiment normalization** to enable meaningful comparisons
- **Generate unified visualizations** with proper position mapping
- **Standardize scoring scales** across experimental batches

## When to Use Batch Processing

- Tiled experiments: Multiple experiments covering non-overlapping regions
- Technical replicates: Same experimental conditions run at different times
- Cross-platform comparisons: Normalizing experiments run on different instruments
- Large protein analysis: Experiments divided by protein domains or regions

## Quick Start

1. Run individual experiments using standard sortscore analysis
2. List 2 or more tiles in experimental_setup file
3. Run batch analysis:
   ```bash
   sortscore -b -e experimental_setup.csv
   ```

## Batch Configuration

### Example Config
Example batch config options are shown below. This template uses default options (besides global_min and global_max, which have no default). This may be added to the same config used for scoring or to a new, seperate config. The options will only be read during batch procesing initiated by the -b option.
```json
{
    "batch_normalization_method": "zscore_2pole",
    "pathogenic_control_type": "nonsense",
    "pathogenic_variants": ,
    "combined_output_dir": "/path/to/combined_results",
    "global_min_pos": 1,
    "global_max_pos": 500
}
```

### Configuration Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `batch_normalization_method` | string | `"zscore_2pole"` | Normalization method: `"zscore_2pole"`, `"2pole"`, or `"zscore_center"` |
| `pathogenic_control_type` | string | `"nonsense"` | Type of pathogenic control: `"nonsense"` or `"custom"` |
| `pathogenic_variants` | list | `null` | Custom pathogenic variants (required when `pathogenic_control_type` is `"custom"`) |
| `combined_output_dir` | string | `"./normalized"` | Output directory for combined results |
| `global_min_pos` | int | `null` | Global minimum position for tiled heatmaps (auto-calculated if not provided) |
| `global_max_pos` | int | `null` | Global maximum position for tiled heatmaps (auto-calculated if not provided) |


## Normalization Methods

### Z-score 2-Pole Normalization (Default)

This is the default method that creates a standardized scale making cross-experiment comparisons meaningful.

**Mathematical Steps:**
1. WT Normalization
2. Z-score Transform
3. Pathogenic Normalization

**Result**: Synonymous variants center around 0 with unit variance across all experiments.

### 2-Pole Normalization
**When to use**:
- When you need to preserve relative score ranges
- Compatibility with published normalization approaches
- Direct interpretation of fold-changes relative to controls

Uses synonymous and pathogenic variants as reference anchor points.

**Formula**: `(b/(a-c)) × (A-C)`

Where:
- `b` = individual variant score
- `a` = experiment synonymous median
- `c` = experiment pathogenic median  
- `A` = global synonymous median
- `C` = global pathogenic median

### Z-score Centering Normalization
**When to use**:
- When pathogenic controls are unavailable
- When you want z-score standardization without pathogenic anchoring

**Mathematical Steps:**
1. WT Normalization
2. Z-score Transform

**Reference Selection**:
- For DNA variants: Uses `wt_dna` scores as reference (falls back to synonymous if unavailable)
- For AA variants: Uses synonymous variants as reference

**Result**: Creates standardized scale centered on synonymous variants with unit variance.


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
- Experimental controls

## Output Files
Batch processing generates:

### Score Files
- `batch_scores_[suffix].csv`: Combined normalized scores with batch tracking
- `batch_stats_[suffix].json`: Comprehensive statistics and normalization factors

### Visualizations  
- Unified scaling: Single colorbar applies to all data  
- Automatic sizing: Figure dimensions scale with position range
`tiled_heatmap_[suffix].png`: Combined heatmap with position mapping
`tiled_heatmap_[suffix]_matrix.csv`: Underlying matrix data for heatmap