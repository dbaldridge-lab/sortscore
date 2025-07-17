# OTX2 Reanalysis with sortscore

## Overview
Reanalyzing OTX2 deep mutational scanning data using the refactored `sortscore` package. This analysis processes DNA-level variant counts across 4 bins and 2 biological replicates.

## Key Analysis Parameters
- **Experiment**: OTX2 transcription factor
- **Variant type**: DNA-level analysis
- **Positions**: 1-297 (full coding sequence)
- **Bins required**: 3 minimum
- **Replicates required**: 1 minimum  
- **Averaging method**: rep-weighted
- **Minimum read threshold**: 3 reads

## Data Processing Notes

### Read Count Source
**Important**: Using read counts **before filtering** for this reanalysis, whereas previous analysis used read counts **after filtering**. The difference between filtered and unfiltered counts is minimal:

| Sample | Before Filtering | After Filtering | Difference |
|--------|------------------|-----------------|------------|
| OTRA-1A | 57,295,729 | 57,079,872 | 0.38% |
| OTRA-1B | 39,086,796 | 38,997,625 | 0.23% |
| OTRA-1C | 44,059,489 | 43,970,598 | 0.20% |
| OTRA-1D | 51,319,515 | 51,214,537 | 0.20% |
| OTRA-2A | 44,754,126 | 44,644,058 | 0.25% |
| OTRA-2B | 50,712,006 | 50,590,118 | 0.24% |
| OTRA-2C | 47,667,963 | 47,513,811 | 0.32% |
| OTRA-2D | 35,165,584 | 35,232,194 |  |

### Data Structure
- **8 samples** across **4 bins** and **2 biological replicates**
- Count data stored as parquet files in `/tests/OTX2/Barcodes/`
- Sample naming: OTRA-{rep}{bin} (e.g., OTRA-1A, OTRA-1B, etc.)

## Analysis Workflow
1. Convert parquet count files to CSV format for sortscore compatibility
2. Configure experiment setup CSV mapping samples to bins/replicates
3. Run sortscore analysis pipeline
4. Generate quality metrics (SD, CV, SEM, CI)
5. Create DNA-level and amino acid heatmaps
6. Compare results with previous analysis

## Expected Outputs
- DNA-level activity scores with statistical measures
- Amino acid-level aggregated scores
- MAVE heatmaps with WT position markers
- Quality metrics for variant filtering and assessment

## Notes
- Full OTX2 coding sequence analysis (positions 1-297)
- Quality filtering may identify "misbehaving" variants with high CV
- Statistical measures include replicate and codon-level variance decomposition