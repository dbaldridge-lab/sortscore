# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
with [PEP 440](https://peps.python.org/pep-0440/) pre-release identifiers for package versions.

## [Unreleased]

### Added
- 

### Changed
- 

### Fixed
- 

## [0.1.0] - 2026-05-12

### Added
- First stable release after the alpha/beta publishing cycle and successful testing by external users.
- Batch amino acid score outputs for combined normalization runs.
- Average positional score bar in batch heatmap outputs.
- GitHub Actions workflow for pull request pytest runs.
- Integration coverage for `run_batch_analysis(...)` batch-config loading.

### Changed
- Promoted package metadata from beta to stable release status.
- Simplified batch output naming by removing date/output suffixes from score file paths and related tests.
- Updated batch visualization layout, including grid refinements and removal of the batch figure title.
- Represent no sequence difference explicitly as `=` during annotation and parsing.
- Preserve synonymous/no-difference rows in amino acid score outputs.
- Write batch normalization outputs under `normalized/<method>/`.
- Simplify summary stats output to focus on the scored file and expose normalization inputs directly in batch stats.

### Fixed
- Preserved decimal precision in batch score exports, including `score.r#b#`, replicate score columns, and `avgscore`.
- Preserved decimal precision in AA score statistics exports, including `SD_codon`, `SD_rep`, `SEM`, and confidence interval columns.
- Stopped rounding scores before plotting and corrected global pathogenic heatmap tick handling.

## [0.1.0b3] - 2026-03-11

### Added
- README visual examples with generated heatmap figures for single-experiment scoring and batch normalization workflows.

### Changed
- Updated README Quick Start to use `pip install sortscore` as the default install path.
- Updated installation and usage docs to consistently prioritize PyPI installation.

### Fixed
- None.

## [0.1.0b2] - 2026-03-11

### Added
- Second beta release candidate for package index testing.

### Changed
- Bumped package version metadata from `0.1.0b1` to `0.1.0b2`.

### Fixed
- Converted README documentation/demo links to absolute GitHub URLs so they render correctly on TestPyPI/PyPI.

## [0.1.0b1] - 2026-03-11

### Added
- First beta release workflow for trusted publishing via GitHub Actions.
- TestPyPI-to-PyPI release pipeline with manual gate before production publish.

### Changed
- 

### Fixed
- 

## [0.1.0a1] - 2026-03-05

### Added
- Initial alpha pre-release.

### Changed
- 

### Fixed
- 

[Unreleased]: https://github.com/dbaldridge-lab/sortscore/compare/0.1.0...HEAD
[0.1.0]: https://github.com/dbaldridge-lab/sortscore/compare/0.1.0b3...0.1.0
[0.1.0b1]: https://github.com/dbaldridge-lab/sortscore/compare/0.1.0a1...0.1.0b1
[0.1.0a1]: https://github.com/dbaldridge-lab/sortscore/releases/tag/0.1.0a1
[0.1.0b2]: https://github.com/dbaldridge-lab/sortscore/compare/0.1.0b1...0.1.0b2
[0.1.0b3]: https://github.com/dbaldridge-lab/sortscore/compare/0.1.0b2...0.1.0b3
