# Sortscore Refactor Progress Checkpoint

## Current Status
- All code and documentation references to "oPool" and "DMS activity analysis" have been generalized to "Sort-seq variant analysis".
- The package is modular, with logical separation into analysis, visualization, and utility modules.
- All plotting and matrix utilities are self-contained and do not depend on the notebook folder.
- Project-level files (README, LICENSE, requirements, setup, docs) are at the repo top level, as recommended.
- Example configuration and usage are provided in the README and config directory.
- All non-package files (e.g., paper outline) are outside the package source.
- Unit tests exist for all major modules and plotting functions.
- Documentation is now split between a concise README and detailed guides in the docs/ folder (installation, usage, index).
- Example and test data folders are set up in tests/data/.
- Config and setup files can be placed anywhere; paths are user-defined.
- Python 3.11+ is now required and enforced in setup.py, requirements.txt, and documentation.

## Outstanding/Next Steps
- Continue to generalize or refactor as needed for new features or user requests.
- Add more usage examples or tutorials in docs/usage.md if desired.
- Review and update dependencies in requirements.txt as needed.
- Prepare for packaging and distribution (PyPI, etc.) if desired.
- Add or update minimal test configs in tests/data/ if needed for automated testing.

## Reference
This checkpoint summarizes the state of the `sortscore` package refactor as of July 7, 2025. Reference this file in future chats to resume or review progress.
