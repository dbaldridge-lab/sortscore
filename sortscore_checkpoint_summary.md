# Sortscore Refactor Checkpoint Summary
*Generated: 2025-08-21*

## Package Status
- **Core modules**: 13+ (load_experiment, score, annotation, analysis_logger, batch_config, plots, heatmaps, heatmap_matrix, run_analysis, run_batch_analysis, ...)
- **Test coverage**: 2 test directories (analysis/tests, visualization/tests)
- **Config templates**: Production ready with dataclass-based configuration
- **Batch processing**: Full support for multi-experiment normalization and tiled visualizations

## Key Capabilities
- Sort-seq activity score calculation with two averaging methods (simple-avg, rep-weighted)
- Batch processing and cross-experiment normalization (z-score and 2-pole methods)
- Positional averages with hex colors for protein structure visualization
- RESTful-style analysis logging with complete audit trails
- Configurable experiment setup via JSON + CSV with flexible file naming
- Integrated visualization (heatmaps, tiled heatmaps, positional averages)
- Modular architecture (analysis/visualization/utilities)
- Comprehensive file naming consistency with optional suffixes

## Architecture
- Modular design: analysis, visualization, utilities
- User-configurable paths with consistent naming patterns
- Self-contained plotting independent of external notebooks
- Generalized for any Sort-seq experiment (not project-specific)
- RESTful logging approach for analysis audit trails

## Recent Updates (August 2025)
- **Positional averages**: Added hex color export for protein structure visualization
- **Analysis logging**: Implemented RESTful AnalysisLogger with complete audit trails
- **File naming consistency**: All files follow {experiment_name}_{file_type}_{suffix}.{ext} pattern
- **Batch normalization**: Full z-score and 2-pole normalization support for multi-experiment analysis
- **Codon weighting removal**: Simplified to two core averaging methods (simple-avg, rep-weighted)
- **Suffix integration**: Consistent suffix support across all output files with API flexibility
- **CLI improvements**: Added --pos-color, --fig-format flags with organized parameter structure

## Installation & Usage
```bash
# Setup (one time)
python3 -m venv venv
source venv/bin/activate
pip install -e .

# Single experiment analysis
sortscore --config path/to/config.json --suffix experiment_name --pos-color --fig-format svg

# Batch processing multiple experiments
sortscore --batch --config batch_config.json --suffix batch_analysis

# Alternative module invocation
python -m sortscore.run_analysis --config config.json
python -m sortscore.run_batch_analysis --config batch_config.json
```

## Ready For
- Production use with comprehensive logging and audit trails
- Single and batch experiment processing with cross-experiment normalization
- Protein structure visualization via positional averages with hex colors
- Complete file management with consistent naming and suffix support
- PyPI packaging and distribution with stable API
- Extension with new analysis methods, normalization approaches, and visualizations

*This is a condensed checkpoint. See full details in sortscore_refactor_checkpoint.md*