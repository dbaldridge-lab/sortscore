# Installation Guide

## Requirements

### Core Requirements
- Python 3.10 or higher
- Python package dependencies (see `requirements.txt`)

### Shell Environment Requirements
If you're only importing sortscore functions in Python code, any environment works. Bash is only needed for the command-line scoring pipeline.

**Command-line scoring pipeline:**
Bash shell environment is required to run the default scoring pipeline. This has not been tested with other shells (e.g., Zsh, Fish).
- **Linux:** Bash (default on most distributions)
- **macOS:** Bash available (run `bash` to switch from default Zsh)
- **Windows:** [Git Bash](https://git-scm.com/downloads) (included with Git for Windows)

**Python API for incorporating scoring functionality into your scripts:**
- No shell environment requirement
- Use sortscore public functions directly in any Python environment (scripts, notebooks, IDEs)
---

## Install Options
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
pip install -r requirements.txt
```

## Troubleshooting
If you encounter any issues during installation:

Ensure you're using Python 3.10 or higher: `python --version`
Try updating pip: `pip install --upgrade pip`
For dependency conflicts, consider using a virtual environment

---

For usage instructions, see [here](https://github.com/dbaldridge-lab/sortscore/edit/main/docs/usage.md).

For contributing to sortscore development, see [here](https://github.com/dbaldridge-lab/sortscore/edit/main/docs/contributing.md).





