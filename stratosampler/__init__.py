"""
stratosampler
=============
Stratified molecular dataset splitting for QSAR model development.

Quick start
-----------
    from stratosampler import PropertyStratifiedSplitter

    splitter = PropertyStratifiedSplitter(
        properties=["MolLogP", "MolWt", "TPSA"],
        n_bins=5,
        test_size=0.2,
        random_state=42,
    )
    train_idx, test_idx = splitter.split(df, smiles_col="SMILES")
"""

from stratosampler.splitters.property_stratified import (
    PropertyStratifiedSplitter,
    compute_properties,
    BUILTIN_PROPERTIES,
)
from stratosampler.metrics.distribution import (
    distribution_report,
    split_summary,
    coverage_score,
    ks_distance,
    js_divergence,
)
from stratosampler.visualisation.plots import (
    plot_property_distributions,
    plot_split_comparison,
    plot_chemical_space,
)
from stratosampler.IO import (
    load_smiles,
    load_sdf,
    SmilesLoader,
    SdfLoader,
)

__version__ = "0.1.1"
__author__ = "Roberto Fino"
__all__ = [
    "PropertyStratifiedSplitter",
    "compute_properties",
    "BUILTIN_PROPERTIES",
    "distribution_report",
    "split_summary",
    "coverage_score",
    "ks_distance",
    "js_divergence",
    "plot_property_distributions",
    "plot_split_comparison",
    "plot_chemical_space",
    "load_smiles",
    "load_sdf",
    "SmilesLoader",
    "SdfLoader",
]
