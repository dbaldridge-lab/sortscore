# Troubleshooting Guide

This guide helps resolve common installation and usage issues with sortscore.


Ensure you're using Python 3.10 or higher.
```
 python --version 
```
 Try updating pip: 
 ```
 pip install --upgrade pip 
 ```
 For dependency conflicts, consider using a virtual environment:
 ```
 # Create virtual environment (do this once)
 python -m venv .venv

 # Activate prior to running sortscore in a terminal
 source .venv/bin/activate.
```

## Diagnostic Commands

If you're having issues, run these diagnostic commands:

```bash
echo "=== Python location ==="
which python

echo "=== Pip location ==="
which pip

echo "=== Is sortscore installed? ==="
python -m pip show sortscore

echo "=== Can Python import it? ==="
python -c "import sortscore; print('Import successful')"

echo "=== Can we run as module? ==="
python -m sortscore --help 2>&1 | head -3

echo "=== Where is sortscore command? ==="
which sortscore

echo "=== Try running sortscore ==="
sortscore --help 2>&1 | head -3
```

This will help identify where the installation is failing.

## Getting More Help
Report issues at https://github.com/dbaldridge-lab/sortscore/issues

Include in your issue report:
- Python version: `python --version` or `python3 --version`
- Full error message
- Operating system
- Environment
- Diagnostic commands output (see above)