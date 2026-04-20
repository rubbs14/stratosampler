# Examples

## Basic Molecular Dataset Split

### Example 1: Simple QSAR Dataset Split

```python
import pandas as pd
from stratosampler import PropertyStratifiedSplitter
from stratosampler.metrics import calculate_distribution_metrics

# Load data
df = pd.read_csv("qsar_data.csv")
X = df[["SMILES", "MolLogP", "MolWt", "TPSA"]]
y = df["activity"]

# Create splitter
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42
)

# Perform split
X_train, X_test, y_train, y_test = splitter.split(X, y=y)

# Validate split quality
metrics = calculate_distribution_metrics(X_train, X_test, y_train, y_test)
print(f"KS Statistic: {metrics['ks_statistic']:.4f}")
print(f"AD Statistic: {metrics['ad_statistic']:.4f}")
```

## Cross-Validation with Stratification

### Example 2: K-Fold Stratified Cross-Validation

```python
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from stratosampler import PropertyStratifiedSplitter

# Initialize model
model = RandomForestRegressor(random_state=42)

# Create stratified splitter
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt"],
    n_bins=5,
    test_size=0.2
)

# Perform stratified split
X_train, X_test, y_train, y_test = splitter.split(X, y=y)

# Train model
model.fit(X_train, y_train)

# Evaluate
train_score = model.score(X_train, y_train)
test_score = model.score(X_test, y_test)

print(f"Train R²: {train_score:.4f}")
print(f"Test R²: {test_score:.4f}")
```

## Visualization of Splits

### Example 3: Plotting Distribution Comparison

```python
import matplotlib.pyplot as plt
from stratosampler.visualisation import plot_distribution_comparison

# Create the split
X_train, X_test, y_train, y_test = splitter.split(X, y=y)

# Plot comparison
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
plot_distribution_comparison(X_train, X_test)
plt.tight_layout()
plt.show()
```

## Comparing Stratified vs Random Splits

### Example 4: Quality Comparison

```python
from sklearn.model_selection import train_test_split
from stratosampler.metrics import calculate_distribution_metrics

# Random split
X_train_random, X_test_random, y_train_r, y_test_r = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Stratified split
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42
)
X_train_strat, X_test_strat, y_train_s, y_test_s = splitter.split(X, y=y)

# Compare metrics
random_metrics = calculate_distribution_metrics(
    X_train_random, X_test_random, y_train_r, y_test_r
)
strat_metrics = calculate_distribution_metrics(
    X_train_strat, X_test_strat, y_train_s, y_test_s
)

print("Random Split KS Statistic:", random_metrics["ks_statistic"])
print("Stratified Split KS Statistic:", strat_metrics["ks_statistic"])
```

## Multi-Property Stratification

### Example 5: Advanced Configuration

```python
# Stratify on multiple chemical properties
properties = [
    "MolLogP",      # Lipophilicity
    "MolWt",        # Molecular weight
    "TPSA",         # Topological polar surface area
    "NumHBD",       # Hydrogen bond donors
    "NumHBA",       # Hydrogen bond acceptors
]

splitter = PropertyStratifiedSplitter(
    properties=properties,
    n_bins=4,  # Fewer bins for more properties
    test_size=0.2,
    random_state=42
)

X_train, X_test, y_train, y_test = splitter.split(X, y=y)
```

## Tips and Best Practices

### Property Selection

- Choose properties that reflect chemical diversity in your dataset
- Common choices: LogP, MW, TPSA, polar surface area
- Avoid properties with high correlation

### Bin Configuration

- **n_bins = 3-4**: Good for small datasets or many properties
- **n_bins = 5-6**: Default, good balance
- **n_bins > 6**: Use for large, diverse datasets

### Validation

Always validate your splits:

```python
# Check distributions match
metrics = calculate_distribution_metrics(X_train, X_test, y_train, y_test)

# Inspect property distributions
print(X_train[properties].describe())
print(X_test[properties].describe())
```
