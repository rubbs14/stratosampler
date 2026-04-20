# UV Usage Examples

Complete examples of using UV for development, testing, and deployment.

## Development Setup

### Fresh Start

```bash
# Clone and setup
git clone https://github.com/rubbs14/stratosampler.git
cd stratosampler

# Install UV
pip install uv

# Install all dependencies (core + all extras)
uv sync --all-extras

# Verify installation
uv run python -c "import stratosampler; print(stratosampler.__version__)"
```

### Add a New Dependency

```bash
# Add to core dependencies
uv add numpy

# Add to dev dependencies
uv add --dev pytest-xdist

# View changes
uv lock
git add uv.lock pyproject.toml
git commit -m "Add new dependency"
```

## Testing

### Run Tests

```bash
# Install dev dependencies
uv sync --extra dev

# Run all tests
uv run pytest

# Run with verbose output
uv run pytest -v

# Run specific test file
uv run pytest tests/test_splitter.py

# Run with coverage
uv run pytest --cov=stratosampler --cov-report=html

# Run specific test
uv run pytest tests/test_splitter.py::test_function_name
```

### Code Quality

```bash
# Format code
uv run black stratosampler tests

# Lint code
uv run ruff check stratosampler tests

# Lint and fix
uv run ruff check --fix stratosampler tests

# Check type hints (if using mypy)
uv run mypy stratosampler  # (requires: uv add --dev mypy)
```

## Documentation

### Build and Preview

```bash
# Install docs dependencies
uv sync --extra docs

# Build documentation
uv run properdocs build

# Preview locally
uv run properdocs serve

# Access at http://localhost:8000
```

### Update Documentation

```bash
# Edit files in docs/
# Then rebuild
uv run properdocs build

# Verify it looks good locally
uv run properdocs serve

# Commit when ready
git add docs/ uv.lock
git commit -m "Update documentation"
git push  # Triggers GitHub Actions
```

## CI/CD Examples

### GitHub Actions Workflow

```yaml
name: Test and Build

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.9", "3.10", "3.11"]
    
    steps:
      - uses: actions/checkout@v4
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}
      
      - name: Install UV
        uses: astral-sh/setup-uv@v2
      
      - name: Install dependencies
        run: uv sync --all-extras
      
      - name: Run tests
        run: uv run pytest --cov
      
      - name: Lint
        run: |
          uv run black --check stratosampler tests
          uv run ruff check stratosampler tests
      
      - name: Build docs
        run: uv run properdocs build
```

### Local CI Simulation

```bash
# Simulate GitHub Actions locally
uv sync --all-extras

# Run all checks
uv run pytest --cov
uv run black --check stratosampler tests
uv run ruff check stratosampler tests
uv run properdocs build

echo "✓ All checks passed!"
```

## Docker/Container Usage

### Dockerfile Example

```dockerfile
FROM python:3.11-slim

# Install UV
RUN pip install --no-cache-dir uv

WORKDIR /app
COPY . .

# Install dependencies using UV
RUN uv sync

# Run application
CMD ["uv", "run", "python", "-m", "stratosampler"]
```

### .devcontainer/devcontainer.json

```json
{
  "image": "mcr.microsoft.com/devcontainers/python:3.11",
  "features": {
    "ghcr.io/astral-sh/uv:latest": {}
  },
  "postCreateCommand": "uv sync --all-extras",
  "customizations": {
    "vscode": {
      "extensions": ["ms-python.python", "charliermarsh.ruff"]
    }
  }
}
```

## Dependency Management

### Check What's Locked

```bash
# View all locked packages
uv lock --show

# View only direct dependencies
cat pyproject.toml | grep dependencies -A 20
```

### Update Dependencies

```bash
# Update everything (respecting version constraints)
uv lock --upgrade

# Update specific packages
uv lock --upgrade-package numpy pandas

# Update to latest (may break constraints)
uv lock --upgrade-package numpy --upgrade-all

# Check for updates
uv pip list --outdated
```

### Pin Specific Versions

Edit `pyproject.toml`:

```toml
dependencies = [
    "numpy>=1.23,<2.0",     # Specific range
    "pandas>=1.5,<3.0",     # Compatible release
    "scikit-learn>=1.1",    # Minimum version
]
```

Then run:

```bash
uv lock
```

## Virtual Environment

### Manual Virtual Environment

```bash
# Create venv
uv venv

# Activate (macOS/Linux)
source .venv/bin/activate

# Activate (Windows PowerShell)
.\.venv\Scripts\Activate.ps1

# Or use uv run directly (no activation needed)
uv run python script.py
```

### Multiple Python Versions

```bash
# Create with specific Python version
uv venv --python 3.10

# Install dependencies
uv sync

# Run with specific Python
uv run --python 3.10 python script.py
```

## Performance Tips

### Caching Dependencies

```bash
# UV caches by default - subsequent syncs are fast
uv sync
uv sync  # This is nearly instant

# Cache directory location
# Unix: ~/.cache/uv
# Windows: %LOCALAPPDATA%\uv\cache
```

### Parallel Installation

```bash
# UV installs in parallel by default
uv sync

# Check progress
uv sync --verbose
```

### Frozen Dependencies

```bash
# Install exact versions from lock file (faster, safer for prod)
uv sync --frozen

# Useful in Docker to avoid lock file changes
RUN uv sync --frozen
```

## Troubleshooting

### Issue: Module Not Found After Installation

```bash
# Make sure to use uv run
uv run python -c "import numpy"  # ✓ Works

# OR activate the venv
source .venv/bin/activate
python -c "import numpy"         # ✓ Also works
```

### Issue: Dependency Conflict

```bash
# Resolve conflicts
uv lock

# If that doesn't work, check constraints
cat pyproject.toml | grep dependencies

# May need to relax version constraints
# Then retry
uv lock
```

### Issue: Old Cached Version

```bash
# Clear cache
uv cache clean

# Reinstall
uv sync
```

## Best Practices

1. **Always commit `uv.lock`**

   ```bash
   git add uv.lock
   git commit -m "Update dependencies"
   ```

2. **Use `uv run` in CI/CD**

   ```yaml
   - run: uv run pytest
   ```

3. **Separate development dependencies**

   ```toml
   [project.optional-dependencies]
   dev = ["pytest", "black", "ruff"]
   ```

4. **Use `--frozen` in production**

   ```bash
   uv sync --frozen
   ```

5. **Keep Python version flexible but tested**

   ```toml
   requires-python = ">=3.9"
   ```

## Resources

- [UV Documentation](https://docs.astral.sh/uv/)
- [UV GitHub](https://github.com/astral-sh/uv)
- [PEP 508 - Dependency Specification](https://peps.python.org/pep-0508/)
- [PyPI](https://pypi.org/)

---

**Need help?** Check [UV_QUICK_REFERENCE.md](UV_QUICK_REFERENCE.md) or [UV_SETUP.md](UV_SETUP.md)
