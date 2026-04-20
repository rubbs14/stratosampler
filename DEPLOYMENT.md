# Documentation Deployment Guide

This project uses **ProperDocs** to generate documentation and deploy it to GitHub Pages.

## Setup

### 1. Initial Setup (One-Time)

The documentation setup is already configured in this repository:

- **Configuration**: `properdocs.yml` - Main documentation configuration
- **Source**: `docs/` - Markdown documentation files
- **Build Output**: `site/` - Generated HTML (committed to `gh-pages` branch)
- **GitHub Actions**: `.github/workflows/deploy-docs.yml` - Automatic deployment workflow

### 2. Enable GitHub Pages

1. Go to your repository on GitHub
2. Navigate to **Settings** → **Pages**
3. Under "Build and deployment":
   - **Source**: Select "GitHub Actions"
   - This will automatically use the workflow defined in `.github/workflows/deploy-docs.yml`

### 3. Deploy Documentation

#### Automatic Deployment

Documentation is automatically built and deployed to GitHub Pages when:

- You push changes to the `main` branch
- You modify files in the `docs/` directory
- You update `properdocs.yml`
- You manually trigger the workflow

#### Manual Deployment

To manually trigger the deployment:

1. Go to **Actions** tab in your repository
2. Select "Deploy Documentation" workflow
3. Click "Run workflow"

### 4. Access Your Documentation

Once deployed, your documentation will be available at:

```
https://rubbs14.github.io/stratosampler
```

(Replace `rubbs14` with your GitHub username)

## Local Development

### Install UV

First, install UV (a fast Python package installer):

```bash
# Using pip
pip install uv

# Or on macOS/Linux using Homebrew
brew install uv
```

### Install All Dependencies

```bash
# Install core and dev dependencies
uv sync

# Install with documentation dependencies
uv sync --extra docs

# Install with all optional dependencies (rdkit, dev, docs)
uv sync --all-extras
```

### Build Locally

```bash
# Using UV (recommended)
uv run properdocs build

# Or with uv sync
uv sync --extra docs
properdocs build

# View the built site
# Open: site/index.html in your browser
```

### Live Server

```bash
# Using UV
uv run properdocs serve

# Or after running uv sync
properdocs serve
```

The site will be available at `http://localhost:8000`

## Documentation Structure

```
docs/
├── index.md              # Home page
├── getting_started.md    # Installation and quick start
├── api.md               # API reference
└── examples.md          # Usage examples
```

## Updating Documentation

### Adding New Pages

1. Create a new `.md` file in the `docs/` directory
2. Add it to the navigation in `properdocs.yml` under the `nav:` section:

```yaml
nav:
  - Home: index.md
  - New Page: new_page.md
```

1. Commit and push to trigger automatic deployment

### Modifying Content

1. Edit the `.md` files in the `docs/` directory
2. Test locally with `properdocs serve`
3. Commit and push changes
4. GitHub Actions will automatically rebuild and deploy

## Configuration

### Theme

The documentation uses the `mkdocs` theme. To customize:

Edit `properdocs.yml`:

```yaml
theme: mkdocs
```

### Site Metadata

Update `properdocs.yml` with your project information:

```yaml
site_name: stratosampler
site_description: Stratified molecular dataset splitting for QSAR model development
repo_url: https://github.com/rubbs14/stratosampler
repo_name: stratosampler
```

### Navigation

Customize the navigation in `properdocs.yml`:

```yaml
nav:
  - Home: index.md
  - Section Name:
    - Page 1: section/page1.md
    - Page 2: section/page2.md
```

## Troubleshooting

### Workflow Failed to Deploy

1. Check the GitHub Actions log:
   - Go to **Actions** → "Deploy Documentation"
   - Click the failed workflow run
   - Review the error in the logs

2. Common issues:
   - Missing dependencies in workflow (update `deploy-docs.yml`)
   - Syntax errors in `properdocs.yml`
   - Missing files referenced in navigation

### Local Build Fails

1. Ensure all dependencies are installed:

   ```bash
   pip install -r requirements-docs.txt
   ```

2. Check for syntax errors in YAML:

   ```bash
   properdocs build --verbose
   ```

3. Verify Markdown syntax:
   - Use a Markdown linter or preview

### Site Not Updating

1. Force rebuild locally:

   ```bash
   properdocs build --clean
   ```

2. Manually trigger the workflow in GitHub Actions

3. Check that your default branch is set to `main` in repository settings

## UV Dependency Management

This project uses **UV** for fast, reliable Python dependency management.

### What is UV?

UV is a modern Python package installer written in Rust. It's much faster than pip and provides deterministic dependency resolution through `uv.lock`.

### Installation

```bash
# Using pip
pip install uv

# Using Homebrew (macOS/Linux)
brew install uv

# Using Windows Package Manager
winget install astral-sh.uv
```

### Quick Start

```bash
# Sync all dependencies (creates virtual environment if needed)
uv sync

# Install with specific extras
uv sync --extra docs          # Documentation dependencies
uv sync --extra rdkit         # RDKit for molecular property calculations
uv sync --all-extras          # All optional dependencies

# Run commands in the virtual environment
uv run properdocs build       # Build docs
uv run pytest                 # Run tests
uv run python script.py       # Run Python scripts
```

### Dependency Groups

The project has the following dependency groups defined in `pyproject.toml`:

- **Core**: numpy, pandas, scikit-learn, scipy, matplotlib
- **dev**: pytest, pytest-cov, black, ruff, jupyter
- **rdkit**: rdkit (for molecular property calculations)
- **docs**: properdocs, properdocs-theme-mkdocs

### Lock File (`uv.lock`)

The `uv.lock` file pins all dependency versions for reproducible installs across environments.

To update dependencies:

```bash
# Update all packages to latest versions (respecting constraints)
uv lock --upgrade

# Update specific packages
uv lock --upgrade-package numpy

# Commit uv.lock to track dependency state
git add uv.lock
git commit -m "Update dependencies"
```

### Troubleshooting UV

**Issue: UV command not found**

```bash
# Add UV to PATH or reinstall
pip install --upgrade uv
```

**Issue: Virtual environment conflicts**

```bash
# Remove and recreate
rm -rf .venv
uv sync
```

**Issue: Specific package version needed**

```bash
# Lock file will be updated automatically when you modify pyproject.toml
uv lock
```

## Resources

- [UV Documentation](https://docs.astral.sh/uv/)
- [ProperDocs Documentation](https://timcera.bitbucket.io/properdocs/full_description.html)
- [GitHub Pages Documentation](https://docs.github.com/en/pages)
- [Markdown Guide](https://www.markdownguide.org/)

## Support

For issues or questions about the documentation setup, please open an issue in the repository.
