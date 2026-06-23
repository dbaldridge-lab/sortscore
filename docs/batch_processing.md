# Batch Processing Guide

This guide provides detailed instructions for cross-tile normalization across
multiple Sort-seq experiments.

## Overview

This workflow enables you to:
1. Standardize scoring scales across multiple tiled experiments
2. Generate unified visualizations with proper position mapping

## When to Use Cross-Tile Normalization

- Tiled experiments: Profiling a large region by splitting it across multiple short-read sequencing tiles
- Large protein analysis: Experiments divided by protein domains or regions, with the caveat that current heatmaps may be sparse

## Quick Start

1. Run tiled scoring to generate tile outputs
2. List 2 or more tiles in experimental_setup file
3. Create a `batch_config.json` with an `experiments` list of tile output directories
4. Run cross-tile normalization:
   ```bash
   sortscore score -n my_experiment -e experiment_setup.csv -c experiment_config.json
   sortscore norm -c batch_config.json
   ```

## Batch Configuration

### Example Config
Example batch config options are shown below. The `experiments` list is required and should include each tile plus the tile output directory where scoring outputs were written.
```json
{
    "experiments": [
        {
            "tile": 1,
            "output_dir": "/path/to/output/tile_1",
            "wt_seq": "ATG...",
            "min_pos": 1,
            "max_pos": 150
        },
        {
            "tile": 2,
            "output_dir": "/path/to/output/tile_2",
            "wt_seq": "ATG...",
            "min_pos": 151,
            "max_pos": 300
        }
    ],
    "batch_normalization_method": "zscore_2pole",
    "pathogenic_control_type": "nonsense",
    "pathogenic_variants": ["E501K"],
    "combined_output_dir": "/path/to/combined_results"
}
```

### Configuration Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `experiments` | list | required | Tile entries for normalization. Each entry must include `tile` (int), `output_dir`, `wt_seq`, `min_pos`, and `max_pos`. |
| `batch_normalization_method` | string | `"zscore_2pole"` | Normalization method: `"zscore_2pole"`, `"linear_range"`, or `"zscore_onepole"` |
| `pathogenic_control_type` | string | `"nonsense"` | Type of pathogenic control: `"nonsense"` or `"custom"` |
| `pathogenic_variants` | list | `null` | Custom pathogenic variants (required when `pathogenic_control_type` is `"custom"`) |
| `combined_output_dir` | string | `"./normalized"` | Output directory for combined results |

CLI override:
- `sortscore norm -c batch_config.json --method zscore_2pole`
- `sortscore norm -c batch_config.json --output-dir /path/to/combined_results`


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
- `<combined_output_dir>/normalized/<method>/scores/batch_dna_scores.csv`: Combined normalized DNA-level scores with batch tracking when DNA scores are available
- `<combined_output_dir>/normalized/<method>/scores/batch_aa_scores.csv`: Combined amino acid level normalized scores when AA aggregation is available
- `<combined_output_dir>/normalized/<method>/scores/batch_stats.json`: Statistics and normalization factors

### Visualizations  
- Unified scaling: Single colorbar applies to all data  
- Automatic sizing: Figure dimensions scale with position range
`tiled_heatmap.png`: Combined heatmap with position mapping
`tiled_heatmap_matrix.csv`: Underlying matrix data for heatmap
