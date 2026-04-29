# Getting Started

## Installation

```bash
pip install stratosampler
# with RDKit (required to compute properties from SMILES):
pip install "stratosampler[rdkit]"
```

For development:

```bash
pip install -e ".[dev,rdkit]"
pytest tests/ -v
```

---

## Basic usage

### 1. Load your data

```python
import pandas as pd

df = pd.read_csv("molecules.csv")
# df must have a SMILES column, e.g. "SMILES"
```

### 2. Create a splitter

```python
from stratosampler import PropertyStratifiedSplitter

splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42,
)
```

### 3. Split

```python
# Returns integer index arrays into df
train_idx, test_idx = splitter.split(df, smiles_col="SMILES")

train_df = df.iloc[train_idx]
test_df  = df.iloc[test_idx]
```

Or get DataFrames directly:

```python
train_df, test_df = splitter.get_split_dataframes(df, smiles_col="SMILES")
```

### 4. Validate the split

```python
from stratosampler import split_summary, compute_properties

props = ["MolLogP", "MolWt", "TPSA"]

# Compute properties and add to DataFrame for metric functions
prop_df = compute_properties(df["SMILES"], props)
df = pd.concat([df, prop_df], axis=1)

summary = split_summary(df, train_idx, test_idx, props)
print(f"Mean KS statistic:  {summary['mean_ks_stat']:.3f}")   # lower = better
print(f"Mean JS divergence: {summary['mean_js_div']:.3f}")    # lower = better
```

---

## Split from pre-computed property columns

If your DataFrame already has property columns, skip SMILES entirely:

```python
train_idx, test_idx = splitter.split(df, property_cols=["logP", "MW", "TPSA"])
```

---

## Three-way split (train / val / test)

```python
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    test_size=0.1,
    val_size=0.1,
    random_state=42,
)
train_idx, val_idx, test_idx = splitter.split(df, smiles_col="SMILES")
```

---

## Scaffold-aware mode

Keeps molecules sharing a Murcko scaffold together in the same split,
preventing analogue leakage while still preserving property distributions.

```python
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    test_size=0.2,
    scaffold_aware=True,
    random_state=42,
)
train_idx, test_idx = splitter.split(df, smiles_col="SMILES")
```

---

## Built-in properties

These names are computed automatically from SMILES strings:

| Name | Description |
|------|-------------|
| `MolLogP` | Wildman-Crippen LogP |
| `MolWt` | Molecular weight |
| `TPSA` | Topological polar surface area |
| `NumHDonors` | H-bond donors |
| `NumHAcceptors` | H-bond acceptors |
| `NumRotBonds` | Rotatable bonds |
| `NumRings` | Total ring count |
| `NumAromaticRings` | Aromatic ring count |
| `FractionCSP3` | Fraction of sp3 carbons |
| `NumHeavyAtoms` | Heavy atom count |

Any valid `rdkit.Chem.Descriptors` attribute name also works.
