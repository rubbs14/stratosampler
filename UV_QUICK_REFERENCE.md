# UV Quick Reference

**Fast Python package manager & installer** | [Docs](https://docs.astral.sh/uv/)

## Installation

```bash
pip install uv
```

## Essential Commands

### Setup

```bash
uv sync                    # Install all core dependencies
uv sync --all-extras       # Install all extras (rdkit, dev, docs)
uv sync --extra docs       # Install specific extra
uv venv                    # Create virtual environment
```

### Running Code

```bash
uv run python script.py    # Run Python script
uv run pytest              # Run tests
uv run properdocs build    # Build docs
uv run properdocs serve    # Preview docs (localhost:8000)
```

### Package Management

```bash
uv pip list                # List installed packages
uv add requests            # Add package (updates pyproject.toml)
uv add --dev pytest        # Add dev dependency
uv remove requests         # Remove package
```

### Lock File

```bash
uv lock                    # Update lock file (respects constraints)
uv lock --upgrade          # Update to latest versions
uv sync --frozen           # Sync to exact lock file versions
```

### Virtual Environment

```bash
source .venv/bin/activate  # Activate (macOS/Linux)
.\.venv\Scripts\Activate.ps1  # Activate (Windows PowerShell)
deactivate                 # Deactivate
```

## Dependency Groups (Optional)

- **rdkit**: `uv sync --extra rdkit`
- **dev**: `uv sync --extra dev` (testing, linting)
- **docs**: `uv sync --extra docs` (documentation)
- **all**: `uv sync --all-extras`

## Project Structure

- `pyproject.toml` - Project metadata & dependencies
- `uv.lock` - Locked dependency versions (commit this!)
- `.venv/` - Virtual environment (auto-created)

## Key Benefits

✓ **10-100x faster** than pip  
✓ **Reproducible** installs with lock file  
✓ **Single tool** replaces pip + venv  
✓ **Written in Rust** for performance  

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `uv: command not found` | `pip install --upgrade uv` |
| `ModuleNotFoundError` | Use `uv run` or activate `.venv` |
| Dependency conflicts | `uv lock --upgrade` |
| Need specific Python | `uv venv --python 3.11` |

## Why UV?

- **Speed**: Parallel dependency resolution in Rust
- **Reliability**: Deterministic `uv.lock` file
- **Developer experience**: Simpler command interface
- **Compatibility**: Works with pip-compatible packages
- **Zero config**: Works out of the box

## Common Workflows

### Development

```bash
uv sync --all-extras       # Install everything
uv run pytest              # Run tests
uv run black stratosampler # Format code
uv run ruff check stratosampler  # Lint
```

### Documentation

```bash
uv sync --extra docs       # Install docs tools
uv run properdocs serve    # Live preview
uv run properdocs build    # Build static site
```

### Production

```bash
uv sync                    # Install core only
uv run python -m stratosampler  # Run your code
```

---

📚 Full guide: [UV_SETUP.md](UV_SETUP.md)  
📖 Documentation: [DEPLOYMENT.md](DEPLOYMENT.md)
