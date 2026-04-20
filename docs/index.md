# stratosampler

**Stratified molecular dataset splitting for QSAR model development.**

## The Problem

When building QSAR models, the way you split your dataset directly determines the reliability of your evaluation.

A standard random `train_test_split` gives you a test set that *happens* to look like the training set — because it was drawn from the same distribution. This inflates your metrics and creates a false sense of generalization.

```
Full dataset     →  logP range: -2 to 7,  MW range: 150–700
Random test set  →  logP range: -1.8 to 6.9,  MW range: 160–690  ← almost identical
```

**stratosampler** fixes this by stratifying across multiple molecular properties simultaneously, ensuring that train and test sets each represent the full property distribution of your dataset.

```
Stratified test set  →  logP range: -1.9 to 6.8,  MW range: 155–685  ← intentionally representative
KS statistic (logP): 0.04 vs 0.18 for random split
```

## Key Features

- **Multi-property stratification**: Split across multiple molecular properties simultaneously
- **Property-aware sampling**: Ensures distributions match across train/test splits
- **RDKit integration**: Compute properties directly from SMILES strings
- **Comprehensive metrics**: Built-in statistical tests (KS, AD) to validate splits

## Installation

```bash
pip install stratosampler

# with RDKit support (required for computing properties from SMILES):
pip install "stratosampler[rdkit]"
```

## Quick Start

```python
from stratosampler import PropertyStratifiedSplitter

splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42
)

X_train, X_test, y_train, y_test = splitter.split(
    X, y=targets
)
```

## Links

- **Repository**: [GitHub](https://github.com/rubbs14/stratosampler)
- **Issues**: [Report a bug](https://github.com/rubbs14/stratosampler/issues)
- **License**: [MIT](https://github.com/rubbs14/stratosampler/blob/main/LICENSE)
