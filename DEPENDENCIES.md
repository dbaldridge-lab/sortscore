# Sortscore Dependencies and References

## Core Dependencies

### Required Python Version
- **Python 3.11+**: The package requires Python 3.11 or higher

### Data Processing and Analysis

#### pandas (>=2.0.0)
- **Purpose**: Data manipulation, CSV/TSV file handling, DataFrame operations
- **Usage**: Primary data structure for scores, count data, and experiment results
- **Reference**: McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 51-56.
- **URL**: https://pandas.pydata.org/
- **Citation**: `pandas development team. (2020). pandas-dev/pandas: Pandas (v1.0.3). Zenodo. https://doi.org/10.5281/zenodo.3509134`

#### numpy (>=1.24.0)  
- **Purpose**: Numerical computing, array operations, statistical calculations
- **Usage**: Matrix operations, normalization, mathematical functions
- **Reference**: Harris, C.R., Millman, K.J., van der Walt, S.J. et al. (2020). Array programming with NumPy. *Nature* 585, 357–362.
- **URL**: https://numpy.org/
- **Citation**: `Harris, C.R., et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2`

### Visualization

#### matplotlib (>=3.6.0)
- **Purpose**: Core plotting library, figure generation, color management
- **Usage**: Heatmap creation, scatter plots, figure export with transparent backgrounds
- **Reference**: Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90-95.
- **URL**: https://matplotlib.org/
- **Citation**: `J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.`

#### seaborn (>=0.12.0)
- **Purpose**: Statistical visualization, enhanced heatmap functionality
- **Usage**: Heatmap generation with advanced styling, color normalization
- **Reference**: Waskom, M. L. (2021). seaborn: statistical data visualization. *Journal of Open Source Software*, 6(60), 3021.
- **URL**: https://seaborn.pydata.org/
- **Citation**: `Waskom, M. L., (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021, https://doi.org/10.21105/joss.03021.`

### Bioinformatics

#### biopython (>=1.81)
- **Purpose**: Sequence manipulation, DNA/protein sequence handling
- **Usage**: DNA to amino acid translation, sequence parsing
- **Reference**: Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.
- **URL**: https://biopython.org/
- **Citation**: `Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11):1422-3 pmid:19304878`

#### hgvs (>=1.5.0)
- **Purpose**: HGVS (Human Genome Variation Society) nomenclature parsing
- **Usage**: Processing variant annotations in standard HGVS format
- **Reference**: HGVS nomenclature for describing sequence variants
- **URL**: https://github.com/biocommons/hgvs
- **Citation**: `den Dunnen JT, et al. HGVS Recommendations for the Description of Sequence Variants: 2016 Update. Hum Mutat. 2016 Jun;37(6):564-9. pmid:26931183`

### Statistical Analysis

#### scipy (imported as needed)
- **Purpose**: Statistical functions, correlation analysis
- **Usage**: Statistical testing in correlation analysis modules
- **Reference**: Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., ... & SciPy 1.0 Contributors. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature methods*, 17(3), 261-272.
- **URL**: https://scipy.org/
- **Citation**: `Pauli Virtanen, et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.`

## Development Dependencies

### Standard Library Modules Used
- **argparse**: Command-line argument parsing
- **json**: JSON configuration file handling  
- **logging**: Comprehensive logging and audit trails
- **sys**: System-specific parameters and functions
- **os**: Operating system interface
- **pathlib**: Object-oriented filesystem paths
- **typing**: Type hints and annotations
- **dataclasses**: Configuration data structures

## Package Structure

### Console Scripts
- **sortscore**: Main entry point for single experiment analysis
- **sortscore.run_analysis**: Alternative module invocation
- **sortscore.run_batch_analysis**: Batch processing of multiple experiments

### Key Modules
- **sortscore.analysis**: Core analysis functions (scoring, normalization, filtering)
- **sortscore.visualization**: Plotting and visualization functions  
- **sortscore.utilities**: Helper functions and utilities

## Compatibility

### Python Version Support
- **Python 3.11+**: Required for modern type hints and performance improvements
- **Python 3.12**: Fully supported and tested

### Platform Support  
- **Cross-platform**: Works on Linux, macOS, and Windows
- **Architecture**: Supports both x86_64 and ARM64 (Apple Silicon)

## Installation

### From Source
```bash
pip install -e .
```

### Dependencies Installation
```bash  
pip install -r requirements.txt
```

## License and Attribution

This package builds upon the scientific Python ecosystem. When using sortscore in publications, please cite the relevant dependencies as listed above, particularly:

- **NumPy** and **pandas** for data processing
- **matplotlib** and **seaborn** for visualization  
- **BioPython** for sequence analysis
- **SciPy** for statistical functions

## Version History

- **v0.1.0**: Initial release with core Sort-seq analysis functionality
- **Dependencies locked**: Minimum versions specified for reproducibility
- **Active development**: Regular updates with backwards compatibility maintained