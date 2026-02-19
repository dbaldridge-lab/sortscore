# Visualization files
figures/{experiment_name}_aa_heatmap_{suffix}.{png|svg|pdf}
figures/{experiment_name}_aa_heatmap_matrix_{suffix}.csv
figures/{experiment_name}_codon_heatmap_{suffix}.{png|svg|pdf}         # when plotting DNA-level scores
figures/{experiment_name}_codon_heatmap_matrix_{suffix}.csv            # when plotting DNA-level scores
### Exporting positional averages for structure visualization

To export positional averages with hex color codes for protein structure visualization, use the `--pos-color` (or `-p`) flag:

```bash
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json -p
```

This generates a file:

```
figures/{experiment_name}_positional_averages_{suffix}.csv
```

The file contains columns for position, average score, and a hex color suitable for mapping onto protein structures.
# Visualization Guide for sortscore

This guide covers heatmap generation, output formats, and visualization-specific settings.

## Generated visualization outputs

When figures are produced, they are written under `figures/` in your output directory.

Common files:

```text
figures/{experiment_name}_aa_heatmap_{suffix}.{png|svg|pdf}
figures/{experiment_name}_aa_heatmap_matrix_{suffix}.csv
figures/{experiment_name}_codon_heatmap_{suffix}.{png|svg|pdf}
figures/{experiment_name}_codon_heatmap_matrix_{suffix}.csv
figures/{experiment_name}_positional_averages_{suffix}.csv
```

Notes:

- Codon heatmap files are produced when DNA-level scores are plotted.
- `positional_averages` is produced when `--pos-color` is enabled.

## CLI options for visualization

| Option | Type | Description |
|---|---|---|
| `--fig-format` | `png` / `svg` / `pdf` | Output file format for figures. |
| `--pos-color` (`-p`) | flag | Export positional averages with hex colors for structure visualization workflows. |
| `--biophysical-prop` | flag | Include biophysical properties panel in heatmaps. |
| `--min-pos` | int | Minimum plotted position (1-based). |
| `--max-pos` | int | Maximum plotted position (1-based). |
| `--position-offset` (`-O`) | int | Offset applied to position labels in plots. |
| `--mutagenesis-variants` (`-v`) | comma-separated string | Order and subset of variant letters shown on heatmap y-axis. |

## Position-axis controls

Use these options to control the x-axis:

- `--min-pos` and `--max-pos` constrain plotted position range.
- `--position-offset` shifts axis labels to match external numbering (for example, mature protein numbering).

Position interpretation depends on detected variant type and parsed variant identifiers.

## Variant-axis controls

`--mutagenesis-variants` controls the y-axis variant order.

Example:

```bash
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json -v "W,F,Y,P,M,I,L,V,A,G,C,S,T,Q,N,D,E,H,R,K,*"
```

You can provide a subset or a different ordering as needed.

## Tiled heatmaps in batch mode

Batch workflows generate unified tiled heatmaps that:

- map each experiment to a global coordinate system
- show boundaries between tiles
- support non-contiguous regions (gaps)
- use one normalized scale across all tiles

Use batch config keys to control global axis range:

- `global_min_pos`
- `global_max_pos`

## Examples

Generate SVG figures:

```bash
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json --fig-format svg
```

Export positional averages with colors:

```bash
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json -p
```

Constrain plotted region and relabel positions:

```bash
sortscore score -n my_experiment -e experiment_setup.csv -c my_experiment.json --min-pos 10 --max-pos 120 -O 99
```
