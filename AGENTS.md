# AGENTS.md - MacroscopicQED

Guidelines for agentic coding tasks in this repository.

## Project Overview

Python scientific computing package for calculating Dyadic Green's Functions in layered media. Uses hydra for configuration management, pytest for testing, and qutip for quantum dynamics.

## Build/Install

```bash
# Create environment from yaml
conda env create -f environment.yaml
conda activate mqed

# Install package in development mode
pip install -e .
```

## Running Tests

```bash
# Run all tests
pytest

# Run a single test file
pytest test/test_hydra.py

# Run a single test function
pytest test/test_hydra.py::test_config_composes

# Run tests matching a pattern
pytest -k "test_config"

# Run with verbose output
pytest -v
```

## Code Style Guidelines

### General
- Python >= 3.8
- Follow PEP 8 style conventions
- Use `snake_case` for function/variable names, `PascalCase` for classes
- Maximum line length: 100 characters (soft guideline)
- Add docstrings to all public functions and classes

### Imports
- Standard library imports first
- Third-party imports second
- Local imports third
- Sort imports alphabetically within each group
- Use absolute imports (e.g., `from mqed.Dyadic_GF.main import run_simulation`)

Example:
```python
import os
from pathlib import Path

import numpy as np
from hydra import compose
from loguru import logger

from mqed.Dyadic_GF.data_provider import DataProvider
from mqed.utils.SI_unit import eV_to_J, hbar
```

### Type Hints
- Use type hints for function signatures when beneficial
- Common types: `int`, `float`, `str`, `bool`, `List`, `Dict`, `Tuple`, `Optional`
- Use `np.ndarray` for numpy arrays, not generic `list`

### Naming
- Variables: descriptive names (`energy_J`, `target_lambdas_m`)
- Constants: uppercase with underscores (`MAX_ITERATIONS = 1000`)
- Private functions/variables: prefix with underscore (`_internal_func()`)

### Error Handling
- Use specific exception types
- Include informative error messages
- Let exceptions propagate when appropriate

### Docstrings
Use Google-style docstrings:
```python
def compute_gf_grid(energy_J, target_lambdas_m, rx_values_m, sim_params, data_provider):
    """Computes the Green's function grid over specified energies and Rx values.

    Args:
        energy_J: Array of energies in Joules.
        target_lambdas_m: Array of wavelengths in meters.
        rx_values_m: Array of Rx position values.
        sim_params: Simulation parameters from config.
        data_provider: Data provider for material properties.

    Returns:
        results_total: 4D numpy array [M,N,3,3] of total Green's functions.
        results_vacuum: 4D numpy array [M,N,3,3] of vacuum Green's functions.
    """
```

### Configuration (Hydra)
- All configuration via YAML files in `configs/` directory
- Use hydra's `@hydra.main` decorator for CLI applications
- Configs should be well-documented with comments

### Data Storage
- Use HDF5 (.h5/.hdf5) for large numerical data
- Follow existing patterns in `mqed/utils/save_hdf5.py`

### Logging
- Use `loguru` logger: `from loguru import logger`
- Use appropriate log levels: `logger.info()`, `logger.success()`, `logger.warning()`, `logger.error()`

### Testing
- Place tests in `test/` directory
- Use pytest fixtures for common setup
- Test hydra configs using `compose()` from hydra
- Use `tmp_path` fixture for temporary file operations

### Documentation
- Use Sphinx for API documentation (see `docs/` directory)
- Build docs: `make docs`
- Update README.md for user-facing changes
