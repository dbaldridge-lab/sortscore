---
applyTo: "**"
---
# Copilot General Coding Instructions

## Documentation
- Use NumPy-style docstrings for all public functions, classes, and modules. Reference: [NumPy docstring guide](https://numpydoc.readthedocs.io/en/latest/format.html)
- Include Parameters, Returns, and Examples sections where applicable.
- Add docstrings to all public-facing code.

## Coding Standards
- Follow [PEP 8](https://peps.python.org/pep-0008/) for code style.
- Use type hints for all function arguments and return values.
- Prefer explicit imports over wildcard imports.
- Use logging instead of print statements for status messages.

## Version Control
- Ignore unnecessary files in .gitignore (e.g., .ipynb_checkpoints, __pycache__, .env).
- Commit configuration templates, not user-specific configs.

## Atomic GitHub Commits
- Make atomic, logically separated Git commits for each significant change or feature.
- Write clear, descriptive commit messages summarizing the purpose and scope of each change.
- Avoid combining unrelated changes in a single commit to facilitate easier code review and rollback.
---
