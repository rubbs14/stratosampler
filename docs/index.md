# stratosampler

**Stratified molecular dataset splitting for QSAR model development.**

---

## The problem with naive splits

A standard random `train_test_split` gives you a test set that *happens* to look like the training set — because it was drawn from the same distribution. This inflates your metrics.

```
Full dataset     →  logP range: -2 to 7,  MW range: 150–700
Random test set  →  logP range: -1.8 to 6.9,  MW range: 160–690  ← almost identical
```

**stratosampler** fixes this by stratifying across multiple molecular properties simultaneously, so train and test sets each represent the full property distribution of your dataset.

```
Stratified test set  →  logP range: -1.9 to 6.8,  MW range: 155–685
KS statistic (logP): 0.04 vs 0.18 for random split
```

---

## Install

```bash
pip install stratosampler
# with RDKit (required to compute properties from SMILES):
pip install "stratosampler[rdkit]"
```

---

## Quick start

```python
from stratosampler import PropertyStratifiedSplitter

splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42,
)

# Returns integer index arrays
train_idx, test_idx = splitter.split(df, smiles_col="SMILES")

# Or get DataFrames directly
train_df, test_df = splitter.get_split_dataframes(df, smiles_col="SMILES")
```

---

## Links

- **Repository**: [GitHub](https://github.com/rubbs14/stratosampler)
- **Issues**: [Report a bug](https://github.com/rubbs14/stratosampler/issues)
- **PyPI**: [stratosampler](https://pypi.org/project/stratosampler/)
- **License**: MIT
