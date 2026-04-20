"""
Split quality metrics
---------------------
Functions to quantify how well a train/test split preserves the property
distribution of the full dataset. Use these to compare split strategies
or to validate a split before training a QSAR model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon


def ks_distance(
    train_vals: np.ndarray,
    test_vals: np.ndarray,
) -> tuple[float, float]:
    """
    Kolmogorov-Smirnov two-sample test between train and test property
    distributions.

    Returns
    -------
    statistic : float — KS statistic (0 = identical, 1 = maximally different)
    p_value   : float — p-value of the test
    """
    stat, pval = stats.ks_2samp(
        train_vals[~np.isnan(train_vals)],
        test_vals[~np.isnan(test_vals)],
    )
    return float(stat), float(pval)


def js_divergence(
    train_vals: np.ndarray,
    test_vals: np.ndarray,
    n_bins: int = 20,
) -> float:
    """
    Jensen-Shannon divergence between train and test property histograms.
    Returns a value in [0, 1] where 0 = identical distributions.
    """
    all_vals = np.concatenate([train_vals, test_vals])
    bins = np.linspace(np.nanmin(all_vals), np.nanmax(all_vals), n_bins + 1)

    p, _ = np.histogram(train_vals, bins=bins, density=True)
    q, _ = np.histogram(test_vals,  bins=bins, density=True)

    # Avoid zero-division in JSD
    p = p + 1e-10
    q = q + 1e-10
    p /= p.sum()
    q /= q.sum()

    return float(jensenshannon(p, q))


def distribution_report(
    data: pd.DataFrame,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    property_cols: list[str],
    val_idx: np.ndarray | None = None,
) -> pd.DataFrame:
    """
    Compute KS statistic and JS divergence for each property, comparing
    train vs test (and optionally train vs val).

    Returns
    -------
    pd.DataFrame with columns:
        property, ks_stat, ks_pval, js_div, [val_ks_stat, val_js_div]
    """
    rows = []
    for prop in property_cols:
        if prop not in data.columns:
            continue
        train_v = data.iloc[train_idx][prop].values.astype(float)
        test_v  = data.iloc[test_idx][prop].values.astype(float)

        ks_stat, ks_pval = ks_distance(train_v, test_v)
        js_div = js_divergence(train_v, test_v)

        row = dict(
            property=prop,
            train_mean=float(np.nanmean(train_v)),
            test_mean=float(np.nanmean(test_v)),
            ks_stat=ks_stat,
            ks_pval=ks_pval,
            js_divergence=js_div,
        )

        if val_idx is not None:
            val_v = data.iloc[val_idx][prop].values.astype(float)
            row["val_mean"] = float(np.nanmean(val_v))
            row["val_ks_stat"], _ = ks_distance(train_v, val_v)
            row["val_js_divergence"] = js_divergence(train_v, val_v)

        rows.append(row)

    return pd.DataFrame(rows)


def coverage_score(
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    strata: np.ndarray,
) -> float:
    """
    Fraction of unique strata in the full dataset that appear in BOTH
    train and test sets. A score of 1.0 means all strata are represented
    in both splits — ideal for QSAR applicability domain coverage.
    """
    all_strata   = set(np.unique(strata))
    train_strata = set(np.unique(strata[train_idx]))
    test_strata  = set(np.unique(strata[test_idx]))
    shared = train_strata & test_strata
    return len(shared) / len(all_strata) if all_strata else 0.0


def split_summary(
    data: pd.DataFrame,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    property_cols: list[str],
    val_idx: np.ndarray | None = None,
    strata: np.ndarray | None = None,
) -> dict:
    """
    Convenience wrapper returning a dict with all key quality metrics.
    Suitable for logging or display.
    """
    n = len(data)
    report = distribution_report(
        data, train_idx, test_idx, property_cols, val_idx
    )
    summary = {
        "n_total":       n,
        "n_train":       len(train_idx),
        "n_test":        len(test_idx),
        "train_frac":    len(train_idx) / n,
        "test_frac":     len(test_idx) / n,
        "mean_ks_stat":  float(report["ks_stat"].mean()),
        "mean_js_div":   float(report["js_divergence"].mean()),
        "per_property":  report,
    }
    if val_idx is not None:
        summary["n_val"]    = len(val_idx)
        summary["val_frac"] = len(val_idx) / n

    if strata is not None:
        summary["coverage_score"] = coverage_score(train_idx, test_idx, strata)

    return summary
