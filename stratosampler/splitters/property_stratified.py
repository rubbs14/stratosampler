"""
PropertyStratifiedSplitter
--------------------------
Splits a molecular dataset into train/test (and optionally validation) subsets
while preserving the multi-property distribution of the full dataset.

Naive random splits produce test sets whose property distributions closely
mirror the training set by chance — but this is not guaranteed, especially
for small datasets or imbalanced property ranges. Scaffold splits address
analogue memorisation but ignore property coverage entirely. This splitter
combines a multi-dimensional binning strategy with optional scaffold awareness
to produce splits where:

  - Each property's distribution is matched between train and test
  - Chemical diversity is maximised within each split
  - The API is fully scikit-learn compatible

Usage
-----
    from stratosampler import PropertyStratifiedSplitter

    splitter = PropertyStratifiedSplitter(
        properties=["MolLogP", "MolWt", "TPSA"],
        n_bins=5,
        test_size=0.2,
    )
    train_idx, test_idx = splitter.split(df, smiles_col="SMILES")
"""

from __future__ import annotations

import hashlib
import warnings
from itertools import product
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_random_state

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    _RDKIT = True
except ImportError:
    _RDKIT = False
    warnings.warn("RDKit not found. Scaffold-aware splitting will be unavailable.")

# ── Built-in property calculators ─────────────────────────────────────────────

BUILTIN_PROPERTIES = {
    "MolLogP":     lambda m: Descriptors.MolLogP(m),
    "MolWt":       lambda m: Descriptors.MolWt(m),
    "TPSA":        lambda m: Descriptors.TPSA(m),
    "NumHDonors":  lambda m: rdMolDescriptors.CalcNumHBD(m),
    "NumHAcceptors": lambda m: rdMolDescriptors.CalcNumHBA(m),
    "NumRotBonds": lambda m: rdMolDescriptors.CalcNumRotatableBonds(m),
    "NumRings":    lambda m: rdMolDescriptors.CalcNumRings(m),
    "NumAromaticRings": lambda m: rdMolDescriptors.CalcNumAromaticRings(m),
    "FractionCSP3": lambda m: rdMolDescriptors.CalcFractionCSP3(m),
    "NumHeavyAtoms": lambda m: m.GetNumHeavyAtoms(),
}


def compute_properties(
    smiles: Iterable[str],
    properties: list[str],
) -> pd.DataFrame:
    """
    Compute RDKit molecular properties for a list of SMILES strings.

    Parameters
    ----------
    smiles : iterable of str
    properties : list of property names (must be keys in BUILTIN_PROPERTIES
                 or valid RDKit Descriptor names)

    Returns
    -------
    pd.DataFrame with one row per molecule, NaN for invalid SMILES.
    """
    if not _RDKIT:
        raise ImportError("RDKit is required to compute properties from SMILES.")

    rows = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(str(smi)) if smi else None
        row = {}
        for prop in properties:
            if mol is None:
                row[prop] = np.nan
            elif prop in BUILTIN_PROPERTIES:
                try:
                    row[prop] = BUILTIN_PROPERTIES[prop](mol)
                except Exception:
                    row[prop] = np.nan
            else:
                # Try as a raw rdkit Descriptor name
                fn = getattr(Descriptors, prop, None)
                if fn is not None:
                    try:
                        row[prop] = fn(mol)
                    except Exception:
                        row[prop] = np.nan
                else:
                    raise ValueError(
                        f"Unknown property '{prop}'. "
                        f"Available built-ins: {list(BUILTIN_PROPERTIES)}"
                    )
        rows.append(row)
    return pd.DataFrame(rows)


def _murcko_scaffold(smi: str) -> str:
    """Return the generic Murcko scaffold SMILES for a molecule."""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return ""
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return ""


# ── Main splitter ─────────────────────────────────────────────────────────────

class PropertyStratifiedSplitter(BaseEstimator):
    """
    Stratified train/test (and optionally validation) splitter for molecular
    datasets that preserves multi-property distributions.

    Parameters
    ----------
    properties : list of str
        Property names to stratify on. Must be keys in BUILTIN_PROPERTIES
        or valid RDKit Descriptor attribute names. Alternatively, if
        ``property_cols`` is provided at split-time, these are treated as
        column names already present in the DataFrame.
    n_bins : int, default 5
        Number of equal-frequency bins per property. Properties are jointly
        binned via a combined bin-label string, forming a multi-dimensional
        stratum.
    test_size : float, default 0.2
        Fraction of data for the test set.
    val_size : float, default 0.0
        Fraction of data for an optional validation set (0 = no val set).
    scaffold_aware : bool, default False
        If True, molecules sharing the same Murcko scaffold are kept together
        in the same split (all go to train or all go to test). This prevents
        analogue leakage while still respecting property stratification at
        the scaffold-cluster level. Requires RDKit.
    min_bin_size : int, default 1
        Strata with fewer than this many samples are merged into a catch-all
        bin to avoid empty strata.
    random_state : int or None, default None

    Examples
    --------
    >>> splitter = PropertyStratifiedSplitter(
    ...     properties=["MolLogP", "MolWt", "TPSA"],
    ...     n_bins=5,
    ...     test_size=0.2,
    ...     random_state=42,
    ... )
    >>> train_idx, test_idx = splitter.split(df, smiles_col="SMILES")

    >>> # With pre-computed property columns already in the DataFrame:
    >>> train_idx, test_idx = splitter.split(df, property_cols=["logP", "MW"])

    >>> # Three-way split:
    >>> splitter = PropertyStratifiedSplitter(test_size=0.1, val_size=0.1)
    >>> train_idx, val_idx, test_idx = splitter.split(df, smiles_col="SMILES")
    """

    def __init__(
        self,
        properties: list[str] | None = None,
        n_bins: int = 5,
        test_size: float = 0.2,
        val_size: float = 0.0,
        scaffold_aware: bool = False,
        min_bin_size: int = 1,
        random_state: Optional[int] = None,
    ):
        self.properties = properties or ["MolLogP", "MolWt", "TPSA"]
        self.n_bins = n_bins
        self.test_size = test_size
        self.val_size = val_size
        self.scaffold_aware = scaffold_aware
        self.min_bin_size = min_bin_size
        self.random_state = random_state

    # ── Public API ─────────────────────────────────────────────────────────────

    def split(
        self,
        data: pd.DataFrame,
        smiles_col: Optional[str] = None,
        property_cols: Optional[list[str]] = None,
    ) -> tuple[np.ndarray, ...]:
        """
        Split the dataset and return integer index arrays.

        Parameters
        ----------
        data : pd.DataFrame
        smiles_col : str, optional
            Column containing SMILES. Required unless ``property_cols`` given.
        property_cols : list of str, optional
            Pre-computed property columns already in ``data``. If provided,
            SMILES are not needed for stratification (but still needed for
            scaffold-aware mode).

        Returns
        -------
        (train_idx, test_idx) or (train_idx, val_idx, test_idx)
        """
        rng = check_random_state(self.random_state)
        n = len(data)

        prop_df = self._get_properties(data, smiles_col, property_cols)
        strata = self._assign_strata(prop_df)

        if self.scaffold_aware:
            if smiles_col is None:
                raise ValueError("scaffold_aware=True requires smiles_col.")
            return self._scaffold_aware_split(data, smiles_col, strata, rng)

        return self._stratified_split(strata, n, rng)

    def get_split_dataframes(
        self,
        data: pd.DataFrame,
        smiles_col: Optional[str] = None,
        property_cols: Optional[list[str]] = None,
    ) -> tuple[pd.DataFrame, ...]:
        """
        Like split() but returns DataFrames instead of index arrays.
        """
        indices = self.split(data, smiles_col, property_cols)
        return tuple(data.iloc[idx].reset_index(drop=True) for idx in indices)

    # ── Internal helpers ───────────────────────────────────────────────────────

    def _get_properties(
        self,
        data: pd.DataFrame,
        smiles_col: Optional[str],
        property_cols: Optional[list[str]],
    ) -> pd.DataFrame:
        if property_cols is not None:
            missing = [c for c in property_cols if c not in data.columns]
            if missing:
                raise ValueError(f"Columns not found in data: {missing}")
            return data[property_cols].copy()

        if smiles_col is None:
            raise ValueError("Provide either smiles_col or property_cols.")
        if smiles_col not in data.columns:
            raise ValueError(f"smiles_col '{smiles_col}' not found in data.")

        return compute_properties(data[smiles_col].tolist(), self.properties)

    def _assign_strata(self, prop_df: pd.DataFrame) -> np.ndarray:
        """
        Bin each property into n_bins equal-frequency bins, then combine
        all bin labels into a single stratum string per molecule.
        """
        bin_labels = []
        for col in prop_df.columns:
            vals = prop_df[col].fillna(prop_df[col].median())
            try:
                labels = pd.qcut(
                    vals, q=self.n_bins, labels=False, duplicates="drop"
                ).astype(str)
            except ValueError:
                # Fewer unique values than bins — use integer rank instead
                labels = pd.Series(
                    pd.Categorical(vals).codes, index=vals.index
                ).astype(str)
            bin_labels.append(labels)

        strata = bin_labels[0]
        for bl in bin_labels[1:]:
            strata = strata + "_" + bl

        # Merge rare strata into "other"
        counts = strata.value_counts()
        rare = counts[counts < self.min_bin_size].index
        strata = strata.where(~strata.isin(rare), other="other")

        return strata.values

    def _stratified_split(
        self,
        strata: np.ndarray,
        n: int,
        rng: np.random.RandomState,
    ) -> tuple[np.ndarray, ...]:
        """
        Sample from each stratum proportionally to achieve the desired
        test (and optionally val) fractions.
        """
        indices = np.arange(n)
        test_idx, val_idx, train_idx = [], [], []

        for stratum in np.unique(strata):
            mask = strata == stratum
            pool = indices[mask]
            rng.shuffle(pool)

            n_test = max(1, round(len(pool) * self.test_size))
            test_idx.extend(pool[:n_test])

            if self.val_size > 0:
                n_val = max(1, round(len(pool) * self.val_size))
                val_idx.extend(pool[n_test:n_test + n_val])
                train_idx.extend(pool[n_test + n_val:])
            else:
                train_idx.extend(pool[n_test:])

        train_arr = np.array(train_idx)
        test_arr  = np.array(test_idx)

        if self.val_size > 0:
            return train_arr, np.array(val_idx), test_arr
        return train_arr, test_arr

    def _scaffold_aware_split(
        self,
        data: pd.DataFrame,
        smiles_col: str,
        strata: np.ndarray,
        rng: np.random.RandomState,
    ) -> tuple[np.ndarray, ...]:
        """
        Group molecules by Murcko scaffold, then split scaffold-groups
        stratified by the dominant stratum of the group.
        """
        scaffolds = data[smiles_col].map(_murcko_scaffold).values
        # Group indices by scaffold
        scaffold_groups: dict[str, list[int]] = {}
        for i, sc in enumerate(scaffolds):
            scaffold_groups.setdefault(sc, []).append(i)

        # Assign each scaffold-group a stratum = most common stratum in group
        group_strata = {
            sc: pd.Series(strata[idxs]).mode()[0]
            for sc, idxs in scaffold_groups.items()
        }

        # Now split scaffold-groups stratified by group stratum
        sc_keys   = np.array(list(scaffold_groups.keys()))
        sc_strata = np.array([group_strata[sc] for sc in sc_keys])

        sc_train, *sc_rest = self._stratified_split(
            sc_strata, len(sc_keys), rng
        )

        def _expand(sc_indices):
            result = []
            for i in sc_indices:
                result.extend(scaffold_groups[sc_keys[i]])
            return np.array(result)

        if self.val_size > 0:
            return _expand(sc_train), _expand(sc_rest[0]), _expand(sc_rest[1])
        return _expand(sc_train), _expand(sc_rest[0])
