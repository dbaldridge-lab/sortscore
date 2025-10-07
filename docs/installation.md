# Installation Guide

## Requirements
- Python 3.11 or higher
- See `requirements.txt` for Python package dependencies

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

### Optional 3: Install in Editable/Development Mode
If you want to make changes to the code and have them reflected immediately:

```bash
pip install -e .
```

## Testing Your Installation
Run the test suite to verify your setup:

```bash
pytest sortscore/sortscore/visualization/tests/
```

---

For usage instructions, see [docs/usage.md](usage.md).





