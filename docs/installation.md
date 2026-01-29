# Installation Guide

## Requirements

### Core Requirements
- Python 3.10 or higher
- Python package dependencies (see [requirements.txt](https://github.com/dbaldridge-lab/sortscore/blob/main/requirements.txt))

**Variant scoring pipeline:**
Bash is used to run the default scoring pipeline. This has not been tested with other shells (e.g., Zsh, Fish).
- **Linux:** Bash (default on most distributions)
- **macOS:** Bash available (run `bash` to switch from default Zsh)
- **Windows:** [Git Bash](https://git-scm.com/downloads) (included with Git for Windows)

**Python API for incorporating scoring functionality into your scripts:**
- You can use public functions directly in Python scripts or notebooks after importing sortscore
---

## Install Options 
### TODO: Do a test release on PyPI. Use option 2 for now.
### Option 1: Install from PyPI (Recommended)
The easiest way to install sortscore is directly from PyPI:
```bash
pip install sortscore
```

### Option 2: Install from Source
Clone the repository and install dependencies:

```bash
git clone https://github.com/dbaldridge-lab/sortscore
cd sortscore
pip install -e .
```

## Troubleshooting
If you encounter any issues during installation:

Ensure you're using Python 3.10 or higher: `python --version`
Try updating pip: `pip install --upgrade pip`
For dependency conflicts, consider using a virtual environment `python -m venv .venv && source .venv/bin/activate`

---
For usage instructions, see [here](https://github.com/dbaldridge-lab/sortscore/edit/main/docs/usage.md).

For contributing to sortscore development, see [here](https://github.com/dbaldridge-lab/sortscore/edit/main/docs/contributing.md).






