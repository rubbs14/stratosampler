# Getting Started

## Installation

### Basic Installation

```bash
pip install stratosampler
```

### With RDKit Support

To compute molecular properties directly from SMILES, install the RDKit extra:

```bash
pip install "stratosampler[rdkit]"
```

### Development Installation

For development and testing:

```bash
pip install -e ".[dev,rdkit]"
```

## Basic Usage

### 1. Create a Splitter

```python
from stratosampler import PropertyStratifiedSplitter

splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42
)
```

### 2. Prepare Your Data

```python
import pandas as pd

# Load your data
df = pd.read_csv("molecules.csv")
X = df[["SMILES", "activity"]]
y = df["bioactivity"]
```

### 3. Perform the Split

```python
X_train, X_test, y_train, y_test = splitter.split(X, y=y)
```

### 4. Validate the Split

```python
from stratosampler.metrics import calculate_distribution_metrics

metrics = calculate_distribution_metrics(X_train, X_test, y_train, y_test)
print(metrics)
```

## Configuration

### PropertyStratifiedSplitter Parameters

- **properties** (list): Properties to stratify on
- **n_bins** (int): Number of bins for stratification
- **test_size** (float): Fraction of data for test set (default: 0.2)
- **random_state** (int): Random seed for reproducibility

## Common Patterns

### Handling Missing Values

```python
# Drop rows with missing values
X = X.dropna()

# Or use fillna
X = X.fillna(X.mean())
```

### Cross-Validation

```python
from stratosampler import PropertyStratifiedSplitter

splitter = PropertyStratifiedSplitter(n_splits=5)

for train_idx, test_idx in splitter.split(X):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    # Train your model...
```

## Troubleshooting

### Property Computation Issues

If molecular properties can't be computed from SMILES:

1. Ensure RDKit is installed: `pip install rdkit`
2. Check SMILES validity
3. Verify property names match RDKit conventions

### Memory Issues with Large Datasets

For datasets with millions of molecules:

```python
# Process in batches
splitter = PropertyStratifiedSplitter(n_bins=10)
X_train, X_test = splitter.split(X)
```

## Next Steps

- Check the [API Reference](api.md) for detailed documentation
- Explore [Examples](examples.md) for more use cases
