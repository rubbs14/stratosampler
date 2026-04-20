# UV Setup Guide

This project uses **UV** for fast, reliable Python dependency management and package installation.

## What is UV?

[UV](https://docs.astral.sh/uv/) is a modern, fast Python package installer written in Rust. It provides:

- **Speed**: Much faster than pip (10-100x for lock resolution)
- **Reproducibility**: `uv.lock` ensures identical installations across environments
- **Reliability**: Deterministic dependency resolution
- **Single tool**: Replaces pip, venv, and poetry for most use cases

## Installation

### Option 1: Using pip (Recommended)

```bash
pip install uv
```

### Option 2: Using Homebrew (macOS/Linux)

```bash
brew install uv
```

### Option 3: Using Windows Package Manager

```bash
winget install astral-sh.uv
```

### Option 4: Using Cargo (if you have Rust)

```bash
cargo install uv
```

### Option 5: Direct Download

Visit [GitHub Releases](https://github.com/astral-sh/uv/releases) for pre-built binaries.

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/rubbs14/stratosampler.git
cd stratosampler
```

### 2. Install Dependencies

```bash
# Install core dependencies
uv sync

# Or install with documentation dependencies
uv sync --extra docs

# Or install all optional dependencies
uv sync --all-extras
```

This will:

- Create a `.venv` virtual environment (if one doesn't exist)
- Download and install all dependencies
- Create a `.venv/pyvenv.cfg` file with configuration

### 3. Activate the Virtual Environment

After `uv sync`, the environment is ready to use with `uv run`.

To manually activate:

**macOS/Linux:**

```bash
source .venv/bin/activate
```

**Windows (PowerShell):**

```powershell
.\.venv\Scripts\Activate.ps1
```

**Windows (Command Prompt):**

```cmd
.venv\Scripts\activate.bat
```

## Common Commands

### Running Python Code

```bash
# Run scripts with uv (no manual activation needed)
uv run python script.py

# Build documentation
uv run properdocs build

# Serve documentation locally
uv run properdocs serve

# Run tests
uv run pytest

# Run specific test file
uv run pytest tests/test_splitter.py
```

### Package Management

```bash
# Install additional package
uv pip install requests

# Add package to project (adds to pyproject.toml)
uv add requests

# Add development dependency
uv add --dev ruff

# Remove package
uv remove requests

# List installed packages
uv pip list
```

### Lock File Management

```bash
# Update lock file (respects version constraints in pyproject.toml)
uv lock

# Update all packages to latest versions
uv lock --upgrade

# Update specific package
uv lock --upgrade-package numpy

# Sync with lock file (default: gets latest versions matching constraints)
uv sync

# Force sync to exact lock file versions
uv sync --frozen
```

## Project Structure

The project uses the following configuration:

### `pyproject.toml`

Defines project metadata and dependencies:

```toml
[project]
name = "stratosampler"
version = "0.1.0"
requires-python = ">=3.9"
dependencies = [
    "numpy>=1.23",
    "pandas>=1.5",
    "scikit-learn>=1.1",
    "scipy>=1.9",
    "matplotlib>=3.6",
]

[project.optional-dependencies]
rdkit = ["rdkit>=2022.03"]
dev = ["pytest", "pytest-cov", "black", "ruff", "jupyter"]
docs = ["properdocs>=1.6.0", "properdocs-theme-mkdocs>=1.6.0"]
```

### `uv.lock`

Auto-generated lock file with pinned versions. **Always commit this file** to ensure reproducible builds.

## Dependency Groups

### Core Dependencies

Required for basic functionality:

- `numpy` - Numerical computing
- `pandas` - Data manipulation
- `scikit-learn` - Machine learning
- `scipy` - Scientific computing
- `matplotlib` - Plotting

Install: `uv sync`

### RDKit Extra

Optional molecular chemistry library:

- `rdkit` - Molecular property calculations from SMILES

Install: `uv sync --extra rdkit`

### Development Extra

Tools for development and testing:

- `pytest` - Testing framework
- `pytest-cov` - Coverage reporting
- `black` - Code formatting
- `ruff` - Linting
- `jupyter` - Interactive notebooks

Install: `uv sync --extra dev`

### Documentation Extra

Tools for building documentation:

- `properdocs` - Documentation generator
- `properdocs-theme-mkdocs` - Documentation theme

Install: `uv sync --extra docs`

### All Extras

```bash
uv sync --all-extras
```

## Development Workflow

### 1. Clone and Setup

```bash
git clone https://github.com/rubbs14/stratosampler.git
cd stratosampler
uv sync --all-extras
```

### 2. Run Tests

```bash
uv run pytest
uv run pytest -v  # Verbose
uv run pytest --cov=stratosampler  # With coverage
```

### 3. Format Code

```bash
uv run black stratosampler tests
uv run ruff check --fix stratosampler tests
```

### 4. Build Documentation

```bash
uv run properdocs build
uv run properdocs serve  # Preview locally
```

### 5. Commit Changes

```bash
git add .
git commit -m "Description of changes"
git push
```

## Troubleshooting

### Issue: "uv: command not found"

**Solution:**

```bash
# Reinstall UV
pip install --upgrade uv

# Or verify installation
python -m pip show uv
```

### Issue: Virtual environment not created

**Solution:**

```bash
# Explicitly create virtual environment
uv venv

# Then sync
uv sync
```

### Issue: Dependency conflicts

**Solution:**

```bash
# Check what's locked
uv lock --show

# Update lock file
uv lock --upgrade

# Verify constraints in pyproject.toml are compatible
```

### Issue: "ModuleNotFoundError" after sync

**Solution:**

```bash
# Ensure running with uv run
uv run python -c "import package_name"

# Or activate venv manually
source .venv/bin/activate  # macOS/Linux
.\.venv\Scripts\Activate.ps1  # Windows PowerShell
```

### Issue: Different Python versions needed

**Solution:**

```bash
# UV can manage Python versions (requires UV 0.1.0+)
uv python install 3.11
uv sync --python 3.11
```

## GitHub Integration

The repository includes GitHub Actions workflow (`.github/workflows/deploy-docs.yml`) that:

1. Installs UV
2. Runs `uv sync --extra docs`
3. Builds documentation
4. Deploys to GitHub Pages

See [DEPLOYMENT.md](DEPLOYMENT.md) for details.

## Best Practices

### 1. Always Commit `uv.lock`

```bash
git add uv.lock
git commit -m "Update dependencies"
```

### 2. Use Specific Version Constraints

In `pyproject.toml`:

```toml
dependencies = [
    "numpy>=1.23,<2.0",  # Good: specific range
    "pandas>=1.5",        # Good: minimum version
]
```

### 3. Separate Dev Dependencies

Use optional dependencies groups for development tools:

```toml
[project.optional-dependencies]
dev = ["pytest", "black", "ruff"]
```

### 4. Document Installation Steps

Include UV installation and sync in project README:

```bash
pip install uv
uv sync
```

### 5. Use `uv run` for CI/CD

Don't activate venv manually in scripts:

```bash
# Good: direct with uv run
uv run pytest

# Avoid: manual activation
source .venv/bin/activate && pytest
```

## References

- **UV Documentation**: <https://docs.astral.sh/uv/>
- **Python Packaging**: <https://packaging.python.org/>
- **PEP 508 - Dependency Specification**: <https://peps.python.org/pep-0508/>
- **PEP 517 - Build System**: <https://peps.python.org/pep-0517/>

## Support

For issues with UV:

- [UV GitHub Issues](https://github.com/astral-sh/uv/issues)
- [UV GitHub Discussions](https://github.com/astral-sh/uv/discussions)

For issues with this project:

- [Project Issues](https://github.com/rubbs14/stratosampler/issues)
