# Instructions for Refactoring Jupyter Notebooks into a Python Package

## General Goals
- Refactor code from Jupyter notebooks into a clean, modular, and reusable Python package.
- Organize code into logical modules (e.g., data loading, preprocessing, analysis, visualization).
- Ensure all code is compatible with Python 3.8+.

## Notebook Refactoring
- Remove all notebook-specific code (e.g., `%matplotlib inline`, cell magics).
- Replace inline plotting with functions that can be called from scripts.
- Move parameter definitions to function arguments or configuration files.
- Avoid hardcoding file paths; use parameters or config files instead.

## Structure Guidelines
- Place package source code in a `src/` directory.
- Use submodules for different functional areas (e.g., `src/data/`, `src/analysis/`, `src/visualization/`).
- Include an `__init__.py` in each package directory.
- Move notebook code into functions and classes as appropriate.
- Place utility/helper functions in a `utils.py` or `utils/` submodule.
- Write tests for each module in a `tests/` directory.

## Dependency Management
- List all dependencies in a requirements.txt or pyproject.toml file.
- Avoid unnecessary dependencies; use standard libraries where possible.
- Document any system-level dependencies or installation steps in the README.md.

## Packaging & Distribution
- Include a setup.py or pyproject.toml for package installation.
- Add metadata such as author, version, and description.
- Ensure the package can be installed and imported with standard Python tooling.

## Large Notebook Refactoring
- When refactoring code from multiple cells, group related code into cohesive functions or classes, even if the code was originally spread across several cells.
- If a function becomes too large or complex, break it down into smaller helper functions and place them in a `utils.py` or appropriate submodule.
- Preserve the logical order and dependencies of code as they appeared in the notebook, but reorganize for clarity and maintainability.
- If code relies on variables or state defined in previous cells, refactor to use function arguments and return values instead of global state.
- When moving code that produces plots or figures, encapsulate plotting logic in functions that accept data and configuration as arguments.
- For code that loads or saves data, replace hardcoded file paths with function arguments or configuration parameters.
- If the notebook contains repeated code patterns, refactor them into reusable functions or methods.
- Add comments or docstrings to explain any non-obvious transformations or refactoring decisions, especially when merging code from multiple cells.
- If there are magic commands or notebook-specific constructs, remove or replace them with standard Python equivalents.
- When in doubt, prefer modularity and testability over preserving the exact cell structure of the original notebook.

## Example Prompt Guidance
- When asked to refactor a notebook cell, wrap the code in a function or class and place it in the appropriate module.
- When generating new modules, include an `__init__.py` and add docstrings.
- When moving code, ensure all dependencies are imported at the top of the module.


## Configuration Management
- Move all hardcoded variables and parameters from notebooks into configuration files that live outside of the package source code (e.g., in a `config/` directory at the project root).
- Use standard formats such as YAML, JSON, or `.env` for configuration files.
- Refactor code to load parameters from these config files using appropriate libraries.
- Pass configuration values to functions and classes via arguments or initialization, rather than relying on hardcoded values.

## User Customization
- Prioritize the ability for users to specify samples, replicates, and experiment-specific variables using a human-readable configuration file (such as YAML or JSON).
- Design code to read experiment parameters, sample lists, and replicate information from these config files, rather than hardcoding them.
- Document the expected structure and required fields of the configuration file in module docstrings and usage examples.
- Ensure that all functions and classes that require experiment-specific variables accept them as arguments or load them from the config file.
- Place example configuration files in a `config/` directory at the project root, outside the package source code.

## Example Usage & Tutorials
- Provide example scripts or notebooks demonstrating how to use the refactored package.
- Include a section in the README.md with quickstart instructions.

## Testing
- Write unit tests for all new modules and functions.
- Use `pytest` as the test framework.
- Place tests in the `tests/` directory, mirroring the package structure.

## Logging & Error Handling
- Use the logging module for status and error messages.
- Avoid using bare except: statements; catch specific exceptions.
- Provide helpful error messages for common user mistakes (e.g., missing config fields).