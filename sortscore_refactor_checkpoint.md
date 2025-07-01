# Sortscore Refactor Progress Checkpoint

## Current Status
- All code and documentation references to "oPool" and "DMS activity analysis" have been generalized to "Sort-seq variant analysis".
- The package is modular, with logical separation into analysis, visualization, and utility modules.
- All plotting and matrix utilities are self-contained and do not depend on the notebook folder.
- Project-level files (README, LICENSE, requirements, setup, docs) are at the repo top level, as recommended.
- Example configuration and usage are provided in the README and config directory.
- All non-package files (e.g., paper outline) are outside the package source.
- Unit tests exist for all major modules and plotting functions.

## Outstanding/Next Steps
- Continue to generalize or refactor as needed for new features or user requests.
- Add more usage examples or tutorials if desired.
- Review and update dependencies in requirements.txt as needed.
- Prepare for packaging and distribution (PyPI, etc.) if desired.

## Reference
This checkpoint summarizes the state of the `sortscore` package refactor as of July 1, 2025. Reference this file in future chats to resume or review progress.
