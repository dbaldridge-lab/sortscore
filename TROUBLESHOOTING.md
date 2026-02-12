# Troubleshooting Guide

This guide helps resolve common installation and usage issues with sortscore.

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
- Python version: `python --version`
- Installation method used
- Full error message
- Operating system
- Diagnostic commands output

## Installation Issues

### Using Conda/Anaconda Environments

**If you're using Conda and nothing works:**

1. **Verify you're in the correct conda environment:**
   ```bash
   conda env list
   # The active environment has a * next to it
   ```

2. **Activate your environment:**
   ```bash
   conda activate your-env-name
   ```

3. **Verify pip is using conda's pip:**
   ```bash
   which pip
   # Should show: /path/to/anaconda3/envs/your-env-name/bin/pip

   which python
   # Should show: /path/to/anaconda3/envs/your-env-name/bin/python
   ```

4. **Install with conda's pip:**
   ```bash
   pip install git+https://github.com/dbaldridge-lab/sortscore.git
   ```

5. **Verify installation in the conda environment:**
   ```bash
   python -m pip show sortscore
   # Should show installation details
   ```

6. **Test it works:**
   ```bash
   python -m sortscore --help
   ```

**Common conda issue:** Installing with system pip instead of conda's pip. The package gets installed to the wrong Python!

**Solution:** Always activate your conda environment FIRST, then verify `which pip` points to your conda environment before installing.

### "command not found: sortscore"

**Problem:** After installing, the `sortscore` command isn't recognized.

**Why this happens:** The Python scripts directory isn't in your system PATH.

**Solutions:**

1. **Use a virtual environment (Recommended):**
   ```bash
   python -m venv sortscore-env
   source sortscore-env/bin/activate  # On Windows: sortscore-env\Scripts\activate
   pip install git+https://github.com/dbaldridge-lab/sortscore.git
   sortscore --help  # Should work now
   ```

2. **Use `python -m sortscore` instead:**
   ```bash
   python -m sortscore --config path/to/config.json
   ```
   This works regardless of PATH settings.

3. **Check where scripts are installed:**
   ```bash
   python -c "import sysconfig; print(sysconfig.get_path('scripts'))"
   ```
   Add that directory to your PATH, or use the full path to the script.

### "'sortscore' is a package and cannot be directly executed"

**Problem:** Running `python -m sortscore` gives this error.

**Causes:**
- Package not installed
- Wrong Python environment (installed with different Python than you're using)
- Working from inside the source directory before installing

**Solutions:**

1. **Verify package is installed:**
   ```bash
   pip show sortscore
   ```
   If not found, install it:
   ```bash
   pip install git+https://github.com/dbaldridge-lab/sortscore.git
   ```

2. **Check you're using the correct Python:**
   ```bash
   which python
   python -m pip show sortscore
   ```

3. **If in the source directory, either:**
   - Install the package: `pip install .`
   - Or navigate out: `cd ..` then try again

### Module import errors (scipy, pandas, etc.)

**Problem:** Errors like `ModuleNotFoundError: No module named 'scipy'`

**Solution:** Reinstall to ensure all dependencies are installed:
```bash
pip uninstall sortscore
pip install git+https://github.com/dbaldridge-lab/sortscore.git
```

## Virtual Environment Issues

### How to check if virtual environment is active

**Method 1: Check your prompt**
When active, you'll see the venv name in parentheses:
```bash
(sortscore-env) user@machine:~/sortscore$
```

**Method 2: Check which Python**
```bash
which python
# In venv: /path/to/sortscore-env/bin/python
# Not in venv: /usr/bin/python or similar
```

**Method 3: Check VIRTUAL_ENV variable**
```bash
echo $VIRTUAL_ENV
# In venv: /path/to/sortscore-env
# Not in venv: (empty/no output)
```

### How to deactivate virtual environment

```bash
deactivate
```

### How to delete virtual environment

```bash
# Just delete the directory
rm -rf sortscore-env
```

## Jupyter Notebook Issues

### Commands not working in Jupyter

**Problem:** Running `pip install` or `sortscore` in a Python notebook cell gives syntax errors.

**Solution:** Use bash/shell commands with `!` prefix:

```python
# In a Jupyter Python cell, prefix shell commands with !
!pip install git+https://github.com/dbaldridge-lab/sortscore.git
!sortscore --config path/to/config.json

# Or use Python API directly
from sortscore.analysis.load_experiment import ExperimentConfig
experiment = ExperimentConfig.from_json("path/to/config.json")
```

### Kernel doesn't see installed package

**Problem:** Installed sortscore but Jupyter kernel can't import it.

**Solution:** Make sure Jupyter is using the same Python where sortscore is installed:

```bash
# Activate your virtual environment first
source sortscore-env/bin/activate

# Install Jupyter in the same environment
pip install jupyter

# Start Jupyter from within the environment
jupyter notebook
```

## Runtime Errors

### "Could not find config file"

**Problem:** Error loading configuration file.

**Solution:**
- Use absolute paths: `/full/path/to/config.json`
- Or relative paths from your current working directory
- Check current directory: `pwd` or in Python: `import os; print(os.getcwd())`

### "No count file column found"

**Problem:** Error loading experiment setup CSV.

**Solution:** Ensure your experiment setup CSV has one of these column names:
- `Path`
- `Read Counts (CSV)`
- `Count File`
- `File Path`