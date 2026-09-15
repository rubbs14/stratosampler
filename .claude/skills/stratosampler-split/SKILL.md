---
name: stratosampler-split
description: Guide for running and evaluating property-stratified molecular dataset splits with stratosampler's PropertyStratifiedSplitter — loading a SMILES dataset, choosing which molecular properties (MolLogP, MolWt, TPSA, etc.) to stratify on, running plain vs scaffold-aware splits, computing and interpreting KS statistic / JS divergence / coverage score, and comparing random vs stratified vs scaffold-aware strategies via compare_strategies. Use this whenever the user wants to split a QSAR/cheminformatics dataset for train/test/val, mentions PropertyStratifiedSplitter, stratified splitting, scaffold_aware, split_dataset, compare_strategies, or asks whether a train/test split is "good" or representative — even if they just say "split this dataset" or "check my train/test split" without naming the library.
---

# stratosampler: property-stratified splitting workflow

## Why this exists

A random `train_test_split` on a molecular dataset tends to produce a test
set that mirrors the training set's property distribution purely by chance
— it inflates apparent model performance because the "test" is barely a
test. `PropertyStratifiedSplitter` (in
`stratosampler/splitters/property_stratified.py`) instead bins molecules by
one or more physicochemical properties and samples each bin proportionally
into train/test/val, so both sides represent the full property space. An
optional `scaffold_aware` mode additionally keeps molecules that share a
Murcko scaffold together in one split, preventing analogue leakage.

Follow the steps below in order. Steps 1–3 are required for any split;
step 4 is optional (only when leakage across close analogues matters);
steps 5–6 are required to say anything credible about split quality;
step 7 is for when the user wants a strategy comparison rather than a
single split.

## Step 1 — Load and inspect the dataset

```python
import pandas as pd

df = pd.read_csv("compounds.csv")
print(df.shape, df.columns.tolist())
# Identify the SMILES column (commonly "SMILES", "smiles", "Canonical_SMILES")
```

If an agent tool call is more convenient than raw pandas in this repo's
chat/agent context, `stratosampler.agent.tools.handle_tool_call("load_data",
{"path": ..., "smiles_col": ...})` returns shape/columns/dtypes/sample as
JSON — it's a thin wrapper, not a different code path.

## Step 2 — Choose properties to stratify on

Built-in properties (computed automatically from SMILES via RDKit, in
`BUILTIN_PROPERTIES`):

| Name | Meaning |
|---|---|
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

Any valid `rdkit.Chem.Descriptors` attribute name also works (e.g.
`"NumValenceElectrons"`), and if the dataset already has a target/activity
column (e.g. pIC50) it's often worth stratifying on that too, passed via
`property_cols` (see Step 3) rather than `properties`.

**Default properties everywhere in this codebase are
`["MolLogP", "MolWt", "TPSA"]`** — the splitter's `__init__`, the
`split_dataset` and `compare_strategies` agent tools, and the README all
default to this triple. Don't invent a different default unless the user's
dataset or question calls for it.

Gotcha: strata are formed by concatenating each property's bin label
(`n_bins` per property), so the number of joint strata grows as
`n_bins ** len(properties)`. Stratifying on more than ~3–4 properties with
the default `n_bins=5` on a small dataset (a few hundred rows) will make
most strata tiny; `min_bin_size` (default 1) merges anything below that
count into a catch-all `"other"` stratum, which quietly weakens the
stratification. Prefer 2–3 properties unless the dataset is large.

## Step 3 — Run PropertyStratifiedSplitter

Real constructor signature (`stratosampler/splitters/property_stratified.py`):

```python
PropertyStratifiedSplitter(
    properties: list[str] | None = None,   # default ["MolLogP", "MolWt", "TPSA"]
    n_bins: int = 5,
    test_size: float = 0.2,
    val_size: float = 0.0,                 # >0 enables a 3-way split
    scaffold_aware: bool = False,
    min_bin_size: int = 1,
    random_state: int | None = None,
)
```

`.split(data, smiles_col=None, property_cols=None)` returns
`(train_idx, test_idx)`, or `(train_idx, val_idx, test_idx)` when
`val_size > 0`. Exactly one of `smiles_col` / `property_cols` must be
usable — `property_cols` skips property computation if the columns already
exist in the DataFrame.

```python
from stratosampler import PropertyStratifiedSplitter

splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
    random_state=42,
)

# From SMILES (properties computed internally via RDKit)
train_idx, test_idx = splitter.split(df, smiles_col="SMILES")

# From property columns already present in df
train_idx, test_idx = splitter.split(df, property_cols=["logP", "MW", "TPSA"])

# Three-way split — set val_size > 0 in the constructor
splitter3 = PropertyStratifiedSplitter(test_size=0.1, val_size=0.1, random_state=42)
train_idx, val_idx, test_idx = splitter3.split(df, smiles_col="SMILES")

# Or get DataFrames directly instead of index arrays
train_df, test_df = splitter.get_split_dataframes(df, smiles_col="SMILES")
```

Via the agent tool layer instead (`stratosampler.agent.tools`), the
equivalent is `split_dataset` with params `path, smiles_col, properties,
n_bins=5, test_size=0.2, val_size=0.0, scaffold_aware=False,
random_state=42, output_dir`. Note the tool's default `random_state` is
`42` (reproducible by default), whereas the raw `PropertyStratifiedSplitter`
class defaults `random_state=None`. It writes `<stem>_train.csv` /
`<stem>_test.csv` (/ `<stem>_val.csv`) next to the input (or to
`output_dir`) and returns sizes, paths, and a `metrics` block from
`split_summary` (see Step 5) in one call — useful when the user just wants
files on disk plus a quality readout, not a Python session.

## Step 4 — scaffold_aware mode (optional)

```python
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    test_size=0.2,
    scaffold_aware=True,
    random_state=42,
)
train_idx, test_idx = splitter.split(df, smiles_col="SMILES")  # smiles_col is REQUIRED here
```

What it does: every molecule is assigned a generic Murcko scaffold; all
molecules sharing a scaffold are treated as one group and go to the same
split (never split across train/test). Each scaffold-group is then
assigned the *mode* (most common) property-stratum among its members, and
groups are stratified-split by that group-level stratum. This prevents
analogue leakage — a model can't get credit for "generalizing" to a
near-identical compound of one it trained on — while still roughly
respecting property coverage.

Trade-off worth telling the user: because whole scaffold groups move
together, per-molecule property stratification is less precise than the
non-scaffold mode (a large scaffold group can only sit entirely on one
side, dragging its stratum's balance). Expect slightly worse KS/JS numbers
than plain stratified splitting — that's the cost of avoiding leakage, not
a bug. Calling `.split(...)` with `scaffold_aware=True` and no
`smiles_col` raises `ValueError("scaffold_aware=True requires smiles_col.")`.

## Step 5 — Evaluate split quality

Functions live in `stratosampler/metrics/distribution.py`:

```python
from stratosampler.metrics.distribution import (
    distribution_report, split_summary, coverage_score, ks_distance, js_divergence,
)
from stratosampler import compute_properties

# Properties must exist as columns to evaluate — compute them if the split
# only used smiles_col:
prop_df = compute_properties(df["SMILES"].tolist(), ["MolLogP", "MolWt", "TPSA"])
df_eval = pd.concat([df.reset_index(drop=True), prop_df], axis=1)

report = distribution_report(df_eval, train_idx, test_idx, property_cols=["MolLogP", "MolWt", "TPSA"])
# -> DataFrame: property, train_mean, test_mean, ks_stat, ks_pval, js_divergence
#    (+ val_mean/val_ks_stat/val_js_divergence if val_idx is passed)

summary = split_summary(df_eval, train_idx, test_idx, property_cols=["MolLogP", "MolWt", "TPSA"])
# -> dict: n_total, n_train, n_test, train_frac, test_frac,
#          mean_ks_stat, mean_js_div, per_property (the report DataFrame)
```

**Gotcha — `coverage_score` needs a `strata` array that `.split()` never
returns.** `coverage_score(train_idx, test_idx, strata)` and
`split_summary(..., strata=strata)` both require the per-molecule stratum
labels, but `PropertyStratifiedSplitter.split()` only returns index
arrays — strata are computed internally by the private
`_assign_strata(prop_df)` method and discarded. In fact, the built-in
`split_dataset` agent tool handler (`_split_dataset` in
`stratosampler/agent/tools.py`) calls `split_summary(...)` **without** a
`strata` argument, so its JSON output never includes `coverage_score`
despite the README and system prompt listing it as a key metric. If the
user specifically wants coverage, either:

```python
strata = splitter._assign_strata(prop_df)   # private method, but the only way to recover strata
cov = coverage_score(train_idx, test_idx, strata)
```

or fold that into `split_summary(df_eval, train_idx, test_idx, property_cols, strata=strata)`
so it lands in the summary dict directly. Flag this to the user rather than
silently reporting KS/JS only when they explicitly asked about coverage.

## Step 6 — Interpret the numbers

From this repo's own documented thresholds (`stratosampler/agent/prompts.py`):

- **KS statistic** (train vs test, per property): 0 = identical distributions, 1 = maximally different.
  - `< 0.05` — excellent
  - `0.05–0.15` — good
  - `> 0.2` — poor (test set does not represent the property range of train)
- **JS divergence**: in `[0, 1]`, 0 = identical histograms. No fixed named bands are documented here — treat it the same direction as KS (lower is better) and compare relatively across strategies rather than reading it against an absolute cutoff.
- **Coverage score**: fraction of unique strata present in *both* train and test. `1.0` means every property-region of the dataset is represented on both sides — this is the number that speaks directly to QSAR applicability-domain coverage, not just central-tendency matching.

Always report numbers with the interpretation attached, e.g. "KS=0.04 —
excellent, train/test are nearly indistinguishable on MolWt" rather than a
bare number.

## Step 7 — Compare strategies (random vs stratified vs scaffold-aware)

The agent tool `compare_strategies` (`stratosampler/agent/tools.py`,
params: `path, smiles_col, properties, test_size=0.2, random_state=42`)
runs all three strategies on the same dataset and returns
`mean_ks_stat` / `mean_js_div` / `n_train` / `n_test` for each. It can be
invoked directly:

```python
from stratosampler.agent.tools import handle_tool_call
import json

result = json.loads(handle_tool_call("compare_strategies", {
    "path": "compounds.csv",
    "smiles_col": "SMILES",
    "properties": ["MolLogP", "MolWt", "TPSA"],  # optional, this is the default anyway
    "test_size": 0.2,
    "random_state": 42,
}))
print(result["strategies"])  # {"random": {...}, "stratified": {...}, "scaffold_stratified": {...}}
```

Or replicate it manually (this is exactly what the tool does internally,
useful when you need the actual index arrays rather than just summary
numbers):

```python
import numpy as np
from stratosampler import PropertyStratifiedSplitter, compute_properties, split_summary

props = ["MolLogP", "MolWt", "TPSA"]
prop_df = compute_properties(df["SMILES"].tolist(), props)
df_eval = pd.concat([df.reset_index(drop=True), prop_df], axis=1)

# Random baseline
rng = np.random.default_rng(42)
idx = np.arange(len(df)); rng.shuffle(idx)
n_test = int(0.2 * len(df))
random_summary = split_summary(df_eval, idx[n_test:], idx[:n_test], props)

# Property-stratified
strat = PropertyStratifiedSplitter(properties=props, test_size=0.2, random_state=42)
tr, te = strat.split(df, smiles_col="SMILES")
strat_summary = split_summary(df_eval, tr, te, props)

# Scaffold-aware stratified
sc = PropertyStratifiedSplitter(properties=props, test_size=0.2, scaffold_aware=True, random_state=42)
sc_tr, sc_te = sc.split(df, smiles_col="SMILES")
sc_summary = split_summary(df_eval, sc_tr, sc_te, props)
```

Lower `mean_ks_stat` / `mean_js_div` = better distribution match between
train and test for that strategy. Expect roughly: random split has the
*lowest* (misleadingly good-looking) KS/JS because it isn't a real
challenge; plain stratified split should show the best (lowest) KS/JS by
construction; scaffold-aware stratified typically sits a bit higher than
plain stratified (see Step 4's trade-off) but is still meant to beat random
on distribution match while also preventing analogue leakage. If
scaffold-aware comes out *worse* than random on KS/JS, that's a legitimate
finding worth surfacing, not just noise — say so.

## Step 8 — Visualize (optional)

```python
from stratosampler import plot_property_distributions, plot_split_comparison

fig = plot_property_distributions(df_eval, train_idx, test_idx, props)  # val_idx optional kwarg
fig.savefig("split_distributions.png", dpi=150, bbox_inches="tight")

fig2 = plot_split_comparison(
    {"random": random_summary, "stratified": strat_summary, "scaffold+strat": sc_summary},
    props,
    metric="ks_stat",  # or "js_divergence"
)
fig2.savefig("strategy_comparison.png", dpi=150, bbox_inches="tight")
```

The agent tool equivalent is `visualize_split` (`train_path, test_path,
smiles_col, properties, output_path`), which reads train/test CSVs back
from disk (e.g. the ones `split_dataset` just wrote) rather than taking
DataFrames directly.

## Quick-reference defaults

| Parameter | `PropertyStratifiedSplitter` default | Agent tool default |
|---|---|---|
| `properties` | `["MolLogP", "MolWt", "TPSA"]` | same |
| `n_bins` | `5` | `5` |
| `test_size` | `0.2` | `0.2` |
| `val_size` | `0.0` | `0.0` |
| `scaffold_aware` | `False` | `False` |
| `min_bin_size` | `1` | n/a |
| `random_state` | `None` | `42` (both `split_dataset` and `compare_strategies`) |

## Common mistakes to avoid

- Calling `.split(..., scaffold_aware=True)` without `smiles_col` — raises `ValueError`, no silent fallback.
- Passing a `properties` name that isn't in `BUILTIN_PROPERTIES` and isn't a valid `rdkit.Chem.Descriptors` attribute — `compute_properties` raises `ValueError` naming the available built-ins.
- Asking for `coverage_score` after only calling `split_summary` without `strata=` — it will simply be absent from the result dict (no error), which can look like it "should have shown up." Compute strata explicitly (Step 5) when coverage matters.
- Judging strategies on KS/JS alone — a random split's flattering-looking KS/JS is the *problem* the whole library exists to catch, not evidence it's the better split. Coverage-of-strata and resistance to analogue leakage matter just as much as raw distribution matching.
- Stratifying on too many properties for the dataset size (Step 2) — leads to sparse strata silently absorbed into `"other"`, diluting what the split actually guarantees.
