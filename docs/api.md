# API Reference

## PropertyStratifiedSplitter

The main class for performing stratified splits based on molecular properties.

```python
from stratosampler import PropertyStratifiedSplitter
```

### Class: PropertyStratifiedSplitter

#### Constructor

```python
PropertyStratifiedSplitter(
    properties: List[str],
    n_bins: int = 5,
    test_size: float = 0.2,
    random_state: Optional[int] = None
)
```

**Parameters:**

- `properties` (List[str]): List of property names to stratify on
- `n_bins` (int): Number of bins for discretizing properties (default: 5)
- `test_size` (float): Fraction of data for test set (default: 0.2)
- `random_state` (Optional[int]): Random seed for reproducibility

#### Methods

##### split()

Perform stratified split on data.

```python
X_train, X_test, y_train, y_test = splitter.split(X, y=None)
```

**Parameters:**

- `X` (pd.DataFrame): Input features containing molecular data
- `y` (Optional[pd.Series]): Target variable

**Returns:**

- Tuple of (X_train, X_test, y_train, y_test)

## Metrics

Distribution metrics for validating splits.

### Function: calculate_distribution_metrics()

```python
from stratosampler.metrics import calculate_distribution_metrics

metrics = calculate_distribution_metrics(X_train, X_test, y_train, y_test)
```

Calculate statistical metrics comparing train and test distributions.

**Parameters:**

- `X_train` (pd.DataFrame): Training features
- `X_test` (pd.DataFrame): Test features
- `y_train` (pd.Series): Training targets
- `y_test` (pd.Series): Test targets

**Returns:**

- Dictionary containing:
  - `ks_statistic`: Kolmogorov-Smirnov test statistic
  - `ad_statistic`: Anderson-Darling test statistic
  - `mean_diff`: Difference in means between distributions

## Visualization

Plot distributions and comparisons.

### Function: plot_distribution_comparison()

```python
from stratosampler.visualisation import plot_distribution_comparison

plot_distribution_comparison(X_train, X_test)
```

Create comparison plots of train and test set distributions.

**Parameters:**

- `X_train` (pd.DataFrame): Training features
- `X_test` (pd.DataFrame): Test features

### Function: plot_stratification_results()

```python
from stratosampler.visualisation import plot_stratification_results

plot_stratification_results(X, splitter)
```

Visualize stratification results.

**Parameters:**

- `X` (pd.DataFrame): Full dataset
- `splitter` (PropertyStratifiedSplitter): Fitted splitter instance

## Molecular Loaders

Load molecules from various file formats using RDKit.

### Function: load_smiles()

```python
from stratosampler import load_smiles

mols, data = load_smiles(
    "compounds.smi",
    smiles_column=0,
    delimiter=None,
    id_column=None,
)
```

Load molecules from a SMILES file.

**Parameters:**

- `filepath` (str or Path): Path to SMILES file
- `smiles_column` (int): Column index containing SMILES (default: 0)
- `delimiter` (str, optional): Column delimiter (auto-detected if None)
- `id_column` (int, optional): Column index for compound IDs
- `keep_properties` (bool): Keep all file columns as metadata (default: False)
- `sanitize` (bool): Sanitize molecules (default: True)
- `add_hydrogens` (bool): Add explicit hydrogens (default: False)
- `compute_2d_coords` (bool): Compute 2D coordinates (default: False)
- `raise_on_invalid` (bool): Raise on invalid SMILES (default: False)

**Returns:**

- `mols` (list): RDKit molecule objects
- `data` (pd.DataFrame): DataFrame with SMILES and metadata

**Examples:**

```python
# Simple SMILES file (one per line)
mols, data = load_smiles("compounds.smi")

# CSV with custom columns
mols, data = load_smiles(
    "compounds.csv",
    delimiter=",",
    smiles_column=1,
    id_column=0,
)

# With 2D coordinates for visualization
mols, data = load_smiles(
    "compounds.smi",
    compute_2d_coords=True,
)
```

### Function: load_sdf()

```python
from stratosampler import load_sdf

mols, data = load_sdf(
    "compounds.sdf",
    include_properties=True,
    sanitize=True,
    remove_hs=True,
)
```

Load molecules from an SDF (Structure Data Format) file.

**Parameters:**

- `filepath` (str or Path): Path to SDF file
- `include_properties` (bool): Extract SDF properties (default: True)
- `sanitize` (bool): Sanitize molecules (default: True)
- `remove_hs` (bool): Remove explicit hydrogens (default: True)
- `compute_2d_coords` (bool): Compute 2D coordinates (default: False)
- `raise_on_invalid` (bool): Raise on invalid molecules (default: False)

**Returns:**

- `mols` (list): RDKit molecule objects
- `data` (pd.DataFrame): DataFrame with molecule properties

**Examples:**

```python
# Basic SDF loading
mols, data = load_sdf("compounds.sdf")

# Keep explicit hydrogens
mols, data = load_sdf(
    "compounds_with_h.sdf",
    remove_hs=False,
)

# Extract properties
mols, data = load_sdf("compounds.sdf", include_properties=True)
print(data.columns)
```

### Class: SmilesLoader

Advanced loader class for SMILES files with fine-grained control.

```python
from stratosampler import SmilesLoader

loader = SmilesLoader(raise_on_invalid=False)
mols, data = loader.load(
    "compounds.csv",
    delimiter=",",
    smiles_column=1,
    add_hydrogens=True,
)
```

### Class: SdfLoader

Advanced loader class for SDF files with fine-grained control.

```python
from stratosampler import SdfLoader

loader = SdfLoader(raise_on_invalid=False)
mols, data = loader.load(
    "compounds.sdf",
    include_properties=True,
    compute_2d_coords=True,
)
```
