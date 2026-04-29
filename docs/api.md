# API Reference

---

## Splitter

### `PropertyStratifiedSplitter`

```python
from stratosampler import PropertyStratifiedSplitter
```

**Constructor**

```python
PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    val_size=0.0,
    scaffold_aware=False,
    min_bin_size=1,
    random_state=None,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `properties` | `list[str]` | `["MolLogP","MolWt","TPSA"]` | Properties to stratify on |
| `n_bins` | `int` | `5` | Equal-frequency bins per property |
| `test_size` | `float` | `0.2` | Fraction of data for test set |
| `val_size` | `float` | `0.0` | Fraction for validation set (0 = no val) |
| `scaffold_aware` | `bool` | `False` | Keep Murcko scaffold groups together |
| `min_bin_size` | `int` | `1` | Rare strata below this are merged into a catch-all bin |
| `random_state` | `int \| None` | `None` | Random seed |

---

**`split(data, smiles_col=None, property_cols=None)`**

Split the dataset and return integer index arrays.

| Parameter | Description |
|-----------|-------------|
| `data` | `pd.DataFrame` |
| `smiles_col` | Column containing SMILES strings. Required unless `property_cols` is given. |
| `property_cols` | Pre-computed property columns already in `data`. |

Returns `(train_idx, test_idx)` or `(train_idx, val_idx, test_idx)` — NumPy integer arrays.

---

**`get_split_dataframes(data, smiles_col=None, property_cols=None)`**

Same as `split()` but returns DataFrames instead of index arrays.

---

### `compute_properties(smiles, properties)`

Compute RDKit molecular properties for a list of SMILES strings.

```python
from stratosampler import compute_properties

prop_df = compute_properties(df["SMILES"], ["MolLogP", "MolWt", "TPSA"])
```

Returns a `pd.DataFrame` with one row per molecule and one column per property. Invalid SMILES produce `NaN`.

---

### `BUILTIN_PROPERTIES`

Dictionary mapping built-in property names to their RDKit calculator functions. Inspect with:

```python
from stratosampler import BUILTIN_PROPERTIES
print(list(BUILTIN_PROPERTIES))
```

---

## Metrics

```python
from stratosampler import split_summary, distribution_report, coverage_score, ks_distance, js_divergence
```

---

### `split_summary(data, train_idx, test_idx, property_cols, val_idx=None, strata=None)`

Convenience wrapper returning a dict with all key quality metrics.

```python
summary = split_summary(df, train_idx, test_idx, ["MolLogP", "MolWt", "TPSA"])
```

**Returns** a dict with:

| Key | Description |
|-----|-------------|
| `n_total`, `n_train`, `n_test` | Set sizes |
| `train_frac`, `test_frac` | Fractions |
| `mean_ks_stat` | Mean KS statistic across properties (lower = better) |
| `mean_js_div` | Mean JS divergence across properties (lower = better) |
| `per_property` | `pd.DataFrame` with per-property metrics |
| `coverage_score` | Present if `strata` is supplied |

---

### `distribution_report(data, train_idx, test_idx, property_cols, val_idx=None)`

Per-property distribution metrics as a DataFrame.

**Returns** a `pd.DataFrame` with columns: `property`, `train_mean`, `test_mean`, `ks_stat`, `ks_pval`, `js_divergence` (plus `val_*` columns if `val_idx` is given).

---

### `coverage_score(train_idx, test_idx, strata)`

Fraction of unique strata that appear in **both** train and test sets.
A score of `1.0` means full applicability domain coverage.

```python
from stratosampler import coverage_score
score = coverage_score(train_idx, test_idx, strata_array)
```

---

### `ks_distance(train_vals, test_vals)`

Two-sample Kolmogorov-Smirnov test. Returns `(statistic, p_value)`.

---

### `js_divergence(train_vals, test_vals, n_bins=20)`

Jensen-Shannon divergence between train and test histograms. Returns a float in `[0, 1]`.

---

## Visualisation

```python
from stratosampler import plot_property_distributions, plot_split_comparison, plot_chemical_space
```

All functions return a `matplotlib.figure.Figure`.

---

### `plot_property_distributions(data, train_idx, test_idx, property_cols, val_idx=None, n_bins=30, figsize=None)`

Overlaid histograms of train/test (and optionally val) distributions, one subplot per property.

```python
fig = plot_property_distributions(df, train_idx, test_idx, ["MolLogP", "MolWt", "TPSA"])
fig.savefig("distributions.png", dpi=150, bbox_inches="tight")
```

---

### `plot_split_comparison(metrics_dict, property_cols, metric='ks_stat', figsize=None)`

Bar chart comparing a metric across multiple split strategies.

```python
results = {
    "random":     split_summary(df, rand_train, rand_test, props),
    "stratified": split_summary(df, train_idx, test_idx, props),
}
fig = plot_split_comparison(results, props, metric="ks_stat")
```

`metric` accepts `"ks_stat"` or `"js_divergence"`.

---

### `plot_chemical_space(data, train_idx, test_idx, x_col, y_col, val_idx=None, figsize=(7,6), alpha=0.5, s=18)`

2D scatter plot of two properties, coloured by split assignment.

```python
fig = plot_chemical_space(df, train_idx, test_idx, x_col="MolLogP", y_col="MolWt")
```

---

## IO (RDKit loaders)

```python
from stratosampler import load_smiles, load_sdf, SmilesLoader, SdfLoader
```

Requires `pip install "stratosampler[rdkit]"`.

---

### `load_smiles(filepath, smiles_column=0, delimiter=None, id_column=None, keep_properties=False, sanitize=True, add_hydrogens=False, compute_2d_coords=False, raise_on_invalid=False)`

Load molecules from a SMILES file (plain `.smi` or delimited `.csv`/`.tsv`).

Returns `(mols, data)` — a list of RDKit molecule objects and a `pd.DataFrame`.

---

### `load_sdf(filepath, include_properties=True, sanitize=True, remove_hs=True, compute_2d_coords=False, raise_on_invalid=False)`

Load molecules from an SDF file.

Returns `(mols, data)`.

---

### `SmilesLoader` / `SdfLoader`

Class-based loaders for fine-grained control:

```python
loader = SmilesLoader(raise_on_invalid=False)
mols, data = loader.load("compounds.csv", delimiter=",", smiles_column=1)
```
