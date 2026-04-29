# Examples

---

## 1. Basic stratified split

```python
import pandas as pd
from stratosampler import PropertyStratifiedSplitter, compute_properties, split_summary

df = pd.read_csv("qsar_data.csv")  # must have a "SMILES" column

splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42,
)
train_idx, test_idx = splitter.split(df, smiles_col="SMILES")

# Compute properties for metrics
props = ["MolLogP", "MolWt", "TPSA"]
prop_df = compute_properties(df["SMILES"], props)
df = pd.concat([df, prop_df], axis=1)

summary = split_summary(df, train_idx, test_idx, props)
print(f"Train: {summary['n_train']}  Test: {summary['n_test']}")
print(f"Mean KS:  {summary['mean_ks_stat']:.3f}")
print(f"Mean JSD: {summary['mean_js_div']:.3f}")
```

---

## 2. Compare random vs stratified vs scaffold-aware

```python
import numpy as np
from stratosampler import (
    PropertyStratifiedSplitter,
    split_summary,
    plot_split_comparison,
)

props = ["MolLogP", "MolWt", "TPSA"]

# Random split
rng = np.random.default_rng(42)
idx = np.arange(len(df))
rng.shuffle(idx)
n_test = int(0.2 * len(df))
random_results = split_summary(df, idx[n_test:], idx[:n_test], props)

# Stratified split
strat = PropertyStratifiedSplitter(test_size=0.2, random_state=42)
tr, te = strat.split(df, smiles_col="SMILES")
strat_results = split_summary(df, tr, te, props)

# Scaffold-aware stratified split
sc = PropertyStratifiedSplitter(test_size=0.2, scaffold_aware=True, random_state=42)
sc_tr, sc_te = sc.split(df, smiles_col="SMILES")
sc_results = split_summary(df, sc_tr, sc_te, props)

fig = plot_split_comparison(
    {"random": random_results, "stratified": strat_results, "scaffold+strat": sc_results},
    props,
    metric="ks_stat",
)
fig.savefig("strategy_comparison.png", dpi=150, bbox_inches="tight")
```

---

## 3. Visualise property distributions

```python
from stratosampler import plot_property_distributions

fig = plot_property_distributions(df, train_idx, test_idx, props)
fig.savefig("distributions.png", dpi=150, bbox_inches="tight")
```

---

## 4. Chemical space scatter

```python
from stratosampler import plot_chemical_space

fig = plot_chemical_space(df, train_idx, test_idx, x_col="MolLogP", y_col="MolWt")
fig.savefig("chemical_space.png", dpi=150, bbox_inches="tight")
```

---

## 5. Three-way split with validation set

```python
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    test_size=0.1,
    val_size=0.1,
    random_state=42,
)
train_idx, val_idx, test_idx = splitter.split(df, smiles_col="SMILES")

summary = split_summary(df, train_idx, test_idx, props, val_idx=val_idx)
print(f"Train: {summary['n_train']}  Val: {summary['n_val']}  Test: {summary['n_test']}")
```

---

## 6. Loading molecules from files

```python
from stratosampler import load_smiles, load_sdf

# From a plain SMILES file
mols, data = load_smiles("compounds.smi")

# From a CSV (SMILES in column 1, IDs in column 0)
mols, data = load_smiles("compounds.csv", delimiter=",", smiles_column=1, id_column=0)

# From an SDF with properties
mols, data = load_sdf("compounds.sdf", include_properties=True)

# data is a DataFrame — feed it straight into the splitter
train_idx, test_idx = splitter.split(data, smiles_col="SMILES")
```
