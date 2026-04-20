"""
Visualisation utilities for split quality inspection.
All functions return matplotlib Figure objects so they can be
saved, shown, or embedded in notebooks.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


_PALETTE = {
    "train": "#2E86AB",
    "test":  "#E84855",
    "val":   "#F9A03F",
    "full":  "#AAAAAA",
}

_SERIES_COLOURS = [
    "#E84855", "#2E86AB", "#F9A03F", "#3BB273", "#7B2D8B",
    "#FF6B35", "#1A936F", "#C84B31", "#4D9DE0", "#E15554",
]


def plot_series_sample(
    full_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    series_col: str = "series",
    pic50_col: str = "pIC50",
    pic50_bins: Optional[list] = None,
    jitter: float = 0.25,
    figsize: tuple = (12, 6),
    seed: int = 42,
) -> plt.Figure:
    """
    Scatter plot of pIC50 vs chemical series, with sampled points highlighted.

    Each series occupies a labelled column on the x-axis. The full population
    is shown as small, semi-transparent grey dots; sampled compounds are
    overlaid as larger, series-coloured filled circles with a white edge.
    Horizontal dashed lines mark pIC50 activity-range boundaries.

    Parameters
    ----------
    full_df    : DataFrame with all compounds (must contain series_col and pic50_col)
    sample_df  : DataFrame with the sampled subset (same column names)
    series_col : column name for the series / scaffold-cluster label
    pic50_col  : column name for pIC50 values
    pic50_bins : pIC50 boundary values for horizontal guide lines; defaults to
                 [5, 6, 7, 8, 9]
    jitter     : half-width of uniform horizontal jitter applied to each point
    figsize    : (width, height) in inches
    seed       : random seed for reproducible jitter

    Returns
    -------
    matplotlib.figure.Figure
    """
    rng = np.random.default_rng(seed)

    if pic50_bins is None:
        pic50_bins = [5, 6, 7, 8, 9]

    # Sort series: named series first (alphabetical), then "minor_series" last
    all_series = sorted(full_df[series_col].unique())
    minor = [s for s in all_series if "minor" in s.lower()]
    named = [s for s in all_series if s not in minor]
    series_order = sorted(named) + sorted(minor)

    x_pos = {s: i for i, s in enumerate(series_order)}
    n_series = len(series_order)
    colour_map = {s: _SERIES_COLOURS[i % len(_SERIES_COLOURS)]
                  for i, s in enumerate(series_order)}

    sampled_ids = set(sample_df["molecule_chembl_id"].values) \
        if "molecule_chembl_id" in sample_df.columns else set()

    fig, ax = plt.subplots(figsize=figsize)

    # --- background: full population (grey) -----------------------------------
    xs_bg, ys_bg = [], []
    for _, row in full_df.iterrows():
        if "molecule_chembl_id" in full_df.columns and row["molecule_chembl_id"] in sampled_ids:
            continue  # draw sampled points on top later
        s = row[series_col]
        xs_bg.append(x_pos[s] + rng.uniform(-jitter, jitter))
        ys_bg.append(row[pic50_col])

    ax.scatter(
        xs_bg, ys_bg,
        c="#CCCCCC", s=18, alpha=0.5, linewidths=0,
        label="unsampled", zorder=2,
    )

    # --- foreground: sampled points (series-coloured) -------------------------
    for series in series_order:
        mask = (sample_df[series_col] == series)
        sub = sample_df[mask]
        if sub.empty:
            continue
        xs = [x_pos[series] + rng.uniform(-jitter * 0.7, jitter * 0.7)
              for _ in range(len(sub))]
        ax.scatter(
            xs, sub[pic50_col].values,
            c=colour_map[series], s=55, alpha=0.92,
            edgecolors="white", linewidths=0.8,
            label=f"{series} (sampled)", zorder=4,
        )

    # --- pIC50 bin boundary lines ---------------------------------------------
    for boundary in pic50_bins:
        ax.axhline(boundary, color="#DDDDDD", linewidth=0.8,
                   linestyle="--", zorder=1)
        ax.text(
            n_series - 0.42, boundary + 0.05,
            f"pIC50={boundary}", fontsize=7, color="#999999", va="bottom",
        )

    # --- axes & labels --------------------------------------------------------
    ax.set_xticks(range(n_series))
    ax.set_xticklabels(series_order, rotation=30, ha="right", fontsize=9)
    ax.set_xlim(-0.6, n_series - 0.4)
    ax.set_ylabel("pIC50", fontsize=11)
    ax.set_xlabel("Chemical series (Murcko scaffold cluster)", fontsize=10)
    ax.set_title(
        "Stratified sample by series x pIC50 range  "
        f"(n={len(sample_df)} of {len(full_df)})",
        fontsize=12, fontweight="bold",
    )

    # count labels above each series column
    for series in series_order:
        n_pop    = (full_df[series_col] == series).sum()
        n_samp   = (sample_df[series_col] == series).sum()
        ax.text(
            x_pos[series], ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 11,
            f"n={n_samp}/{n_pop}",
            ha="center", va="bottom", fontsize=7.5, color="#555555",
        )

    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        loc="upper left", bbox_to_anchor=(1.01, 1),
        fontsize=8, frameon=False,
    )
    fig.tight_layout()
    return fig


def plot_property_distributions(
    data: pd.DataFrame,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    property_cols: list[str],
    val_idx: Optional[np.ndarray] = None,
    n_bins: int = 30,
    figsize: Optional[tuple] = None,
) -> plt.Figure:
    """
    Overlay train/test (and optionally val) property distributions as
    histograms, one subplot per property.

    Parameters
    ----------
    data         : DataFrame containing property columns
    train_idx    : integer index array for training set
    test_idx     : integer index array for test set
    property_cols: list of column names to plot
    val_idx      : optional integer index array for validation set
    n_bins       : number of histogram bins
    figsize      : (width, height); auto-computed if None

    Returns
    -------
    matplotlib.figure.Figure
    """
    n_props = len(property_cols)
    ncols = min(3, n_props)
    nrows = int(np.ceil(n_props / ncols))
    fw = figsize or (5 * ncols, 4 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=fw)
    axes = np.array(axes).flatten()

    splits = {"train": train_idx, "test": test_idx}
    if val_idx is not None:
        splits["val"] = val_idx

    for i, prop in enumerate(property_cols):
        ax = axes[i]
        if prop not in data.columns:
            ax.set_visible(False)
            continue

        all_vals = data[prop].dropna().values
        bins = np.linspace(np.nanmin(all_vals), np.nanmax(all_vals), n_bins + 1)

        for split_name, idx in splits.items():
            vals = data.iloc[idx][prop].dropna().values
            ax.hist(
                vals, bins=bins, alpha=0.55,
                color=_PALETTE[split_name], label=split_name,
                density=True, edgecolor="none",
            )

        ax.set_title(prop, fontsize=11, fontweight="bold")
        ax.set_xlabel(prop, fontsize=9)
        ax.set_ylabel("Density", fontsize=9)
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Property distributions by split", fontsize=13, y=1.01)
    fig.tight_layout()
    return fig


def plot_split_comparison(
    metrics_dict: dict[str, dict],
    property_cols: list[str],
    metric: str = "ks_stat",
    figsize: Optional[tuple] = None,
) -> plt.Figure:
    """
    Bar chart comparing a distribution metric across multiple split strategies.

    Parameters
    ----------
    metrics_dict  : {strategy_name: summary_dict} from metrics.split_summary()
    property_cols : properties to show on x-axis
    metric        : 'ks_stat' or 'js_divergence'
    figsize       : (width, height)

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> from stratosampler.metrics.distribution import split_summary
    >>> results = {
    ...     "random":     split_summary(data, *random_split, props),
    ...     "stratified": split_summary(data, *strat_split, props),
    ...     "scaffold":   split_summary(data, *scaffold_split, props),
    ... }
    >>> fig = plot_split_comparison(results, props, metric="ks_stat")
    """
    strategies = list(metrics_dict.keys())
    n_props = len(property_cols)
    x = np.arange(n_props)
    width = 0.8 / len(strategies)

    fw = figsize or (max(8, n_props * 1.5), 5)
    fig, ax = plt.subplots(figsize=fw)

    colours = plt.cm.tab10(np.linspace(0, 0.8, len(strategies)))

    for k, (name, summary) in enumerate(metrics_dict.items()):
        report = summary["per_property"].set_index("property")
        vals = [
            report.loc[p, metric] if p in report.index else 0.0
            for p in property_cols
        ]
        bars = ax.bar(
            x + k * width, vals, width,
            label=name, color=colours[k], alpha=0.85, edgecolor="white",
        )

    metric_label = {
        "ks_stat":      "KS statistic (↓ better)",
        "js_divergence":"JS divergence (↓ better)",
    }.get(metric, metric)

    ax.set_xticks(x + width * (len(strategies) - 1) / 2)
    ax.set_xticklabels(property_cols, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel(metric_label, fontsize=10)
    ax.set_title(
        f"Split strategy comparison — {metric_label}", fontsize=12
    )
    ax.legend(fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def plot_chemical_space(
    data: pd.DataFrame,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    x_col: str,
    y_col: str,
    val_idx: Optional[np.ndarray] = None,
    figsize: tuple = (7, 6),
    alpha: float = 0.5,
    s: float = 18,
) -> plt.Figure:
    """
    2D scatter plot of two properties, coloured by split assignment.
    Useful for a quick visual sanity check of coverage.

    Parameters
    ----------
    x_col, y_col : column names for the two axes
    """
    fig, ax = plt.subplots(figsize=figsize)

    splits = {"train": train_idx, "test": test_idx}
    if val_idx is not None:
        splits["val"] = val_idx

    zorders = {"train": 1, "val": 2, "test": 3}

    for split_name, idx in splits.items():
        sub = data.iloc[idx]
        ax.scatter(
            sub[x_col], sub[y_col],
            c=_PALETTE[split_name], label=split_name,
            alpha=alpha, s=s, edgecolors="none",
            zorder=zorders.get(split_name, 1),
        )

    ax.set_xlabel(x_col, fontsize=10)
    ax.set_ylabel(y_col, fontsize=10)
    ax.set_title(f"Chemical space coverage: {x_col} vs {y_col}", fontsize=11)
    ax.legend(fontsize=9, markerscale=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig
