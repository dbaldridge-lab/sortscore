# Installation Guide

## Requirements

- Python 3.10 or higher
- Python package dependencies in [`requirements.txt`](../requirements.txt)
- Bash-compatible shell for the default CLI workflow

Platform notes:
- Linux: supported
- macOS: supported
- Windows: not officially supported; Git Bash is recommended if needed

## Install Options

### Option 1: Install from Source (current recommended path)

```bash
git clone https://github.com/dbaldridge-lab/sortscore
cd sortscore
pip install -e .
```

### Option 2: Install from PyPI

```bash
pip install sortscore
```

Use this option only if your required package version is available on PyPI.

## Verify Installation

```bash
python -m sortscore --help
sortscore --help
```

If `sortscore` is not found, use the module form (`python -m sortscore`) or confirm your environment is activated.

## Troubleshooting

See [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) for diagnostic commands and common environment issues.

## Next Step

Continue with the [Usage Guide](usage.md).
