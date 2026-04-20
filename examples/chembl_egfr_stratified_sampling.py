"""
ChEMBL EGFR Stratified Sampling Example
=========================================
Fetches IC50 activity data for EGFR (CHEMBL203) from the ChEMBL REST API,
assigns chemical series via Murcko scaffold clustering, then applies
stratified sampling across (series x pIC50 range) strata with a guaranteed
minimum of 5 representatives per series.

Usage
-----
    python examples/chembl_egfr_stratified_sampling.py

Dependencies
------------
    pip install requests rdkit pandas numpy
    # stratosampler must be on sys.path (run from repo root)
"""

from __future__ import annotations

import sys
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

# -- Allow running from repo root without install -------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from stratosampler.splitters.property_stratified import PropertyStratifiedSplitter
from stratosampler.visualisation.plots import plot_series_sample

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: RDKit not found - series assignment will fall back to ChEMBL ID prefix clusters.")


# -- Configuration --------------------------------------------------------------

TARGET_CHEMBL_ID  = "CHEMBL203"          # EGFR, Homo sapiens
TARGET_NAME       = "EGFR"
CHEMBL_BASE       = "https://www.ebi.ac.uk/chembl/api/data"
MAX_RECORDS       = 500                  # cap on how many activities to fetch
TOP_N_SERIES      = 5                    # named scaffold series (+ "minor_series")
SAMPLE_FRACTION   = 0.30                 # fraction to draw per stratum
MIN_PER_SERIES    = 5                    # guaranteed minimum representatives per series
PIC50_BINS        = [0, 5, 6, 7, 8, 9, 15]
PIC50_LABELS      = ["<5", "5-6", "6-7", "7-8", "8-9", ">=9"]
RANDOM_STATE      = 42


# -- 1. Fetch data from ChEMBL REST API ----------------------------------------

def fetch_chembl_activity(target_id: str, max_records: int = 500) -> pd.DataFrame:
    """
    Pull IC50 records with pChEMBL values and SMILES for a given ChEMBL target.
    Paginates automatically; deduplicates by molecule_chembl_id (keeps best IC50).
    """
    url = f"{CHEMBL_BASE}/activity.json"
    base_params = {
        "target_chembl_id": target_id,
        "standard_type": "IC50",
        "standard_units": "nM",
        "pchembl_value__isnull": False,
        "limit": 100,
    }

    records: list[dict] = []
    offset = 0

    print(f"Fetching IC50 data for {target_id} from ChEMBL ...")
    while len(records) < max_records:
        resp = requests.get(url, params={**base_params, "offset": offset}, timeout=30)
        resp.raise_for_status()
        payload = resp.json()

        activities = payload.get("activities", [])
        for act in activities:
            smiles    = act.get("canonical_smiles")
            chembl_id = act.get("molecule_chembl_id")
            std_val   = act.get("standard_value")
            pchembl   = act.get("pchembl_value")
            if smiles and chembl_id and std_val and pchembl:
                try:
                    records.append({
                        "molecule_chembl_id": chembl_id,
                        "smiles": smiles,
                        "ic50_nm": float(std_val),
                        "pIC50": float(pchembl),
                    })
                except (TypeError, ValueError):
                    continue

        meta = payload.get("page_meta", {})
        if meta.get("next") is None or len(records) >= max_records:
            break
        offset += base_params["limit"]

    df = (
        pd.DataFrame(records)
        .drop_duplicates("molecule_chembl_id")
        .sort_values("pIC50", ascending=False)
        .reset_index(drop=True)
    )
    print(f"  Retrieved {len(df)} unique compounds  "
          f"(pIC50 range: {df['pIC50'].min():.1f} - {df['pIC50'].max():.1f})")
    return df


# -- 2. Assign chemical series via Murcko scaffold clustering ------------------

def _murcko(smi: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return "__invalid__"
        sc = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(sc)
    except Exception:
        return "__invalid__"


def _fallback_series(chembl_id: str) -> str:
    """Rough series from ChEMBL ID numeric range (no RDKit)."""
    try:
        num = int(chembl_id.replace("CHEMBL", ""))
        bucket = (num // 50000) * 50000
        return f"id_group_{bucket}"
    except ValueError:
        return "misc"


def assign_series(df: pd.DataFrame, top_n: int = TOP_N_SERIES) -> pd.DataFrame:
    """
    Label each compound with a series name.

    With RDKit: Murcko scaffold → top N scaffolds each get their own label;
    the rest are collapsed into 'minor_series'.

    Without RDKit: ChEMBL ID numeric-range bucketing (rough proxy).
    """
    df = df.copy()

    if HAS_RDKIT:
        print("\nComputing Murcko scaffolds ...")
        df["scaffold"] = df["smiles"].map(_murcko)
        sc_counts = df["scaffold"].value_counts()
        top_scaffolds = sc_counts.head(top_n).index.tolist()

        series_map = {sc: f"series_{i+1}" for i, sc in enumerate(top_scaffolds)}

        def _label(sc):
            return series_map.get(sc, "minor_series")

        df["series"] = df["scaffold"].map(_label)
    else:
        df["series"] = df["molecule_chembl_id"].map(_fallback_series)

    counts = df["series"].value_counts()
    print("\nSeries composition (before sampling):")
    for series, n in counts.items():
        print(f"  {series:20s}  {n:4d} compounds")

    return df


# -- 3. Bin pIC50 ---------------------------------------------------------------

def assign_pic50_bins(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pIC50_bin"] = pd.cut(
        df["pIC50"],
        bins=PIC50_BINS,
        labels=PIC50_LABELS,
        include_lowest=True,
    ).astype(str)
    return df


# -- 4. Stratified sampling across (series x pIC50 bin) strata ----------------

def stratified_sample(
    df: pd.DataFrame,
    sample_fraction: float = SAMPLE_FRACTION,
    min_per_series: int = MIN_PER_SERIES,
    random_state: int = RANDOM_STATE,
) -> pd.DataFrame:
    """
    Draw a stratified sample from df using (series, pIC50_bin) as strata.

    Algorithm
    ---------
    1. Within each stratum take floor(n * sample_fraction) molecules (min 1).
    2. After the initial draw, check each series for the minimum guarantee.
    3. If a series is under-represented, top up by sampling from its
       unselected molecules (spreading top-ups across pIC50 bins if possible).

    Uses stratosampler's PropertyStratifiedSplitter internally for the
    pIC50-only binning validation step (reports distribution quality).
    """
    rng = np.random.default_rng(random_state)

    # -- 4a. Initial proportional draw per stratum ------------------------------
    selected: list[int] = []
    for (series, pic50_bin), grp in df.groupby(["series", "pIC50_bin"], observed=True):
        n_draw = max(1, math.floor(len(grp) * sample_fraction))
        chosen = rng.choice(grp.index.values, size=min(n_draw, len(grp)), replace=False)
        selected.extend(chosen.tolist())

    sampled_idx = set(selected)

    # -- 4b. Enforce minimum per series ----------------------------------------
    for series in df["series"].unique():
        series_mask  = df["series"] == series
        current_n    = sum(i in sampled_idx for i in df[series_mask].index)
        shortfall    = min_per_series - current_n

        if shortfall <= 0:
            continue

        # Pool: unselected molecules from this series, sorted by pIC50_bin
        # to spread top-ups across activity ranges
        unselected = (
            df[series_mask & ~df.index.isin(sampled_idx)]
            .sort_values("pIC50_bin")
        )
        top_up = rng.choice(
            unselected.index.values,
            size=min(shortfall, len(unselected)),
            replace=False,
        )
        sampled_idx.update(top_up.tolist())

    sample = df.loc[sorted(sampled_idx)].copy().reset_index(drop=True)
    return sample


# -- 5. Validate using stratosampler's built-in metrics ------------------------

def validate_with_stratosampler(full_df: pd.DataFrame, sample_df: pd.DataFrame) -> None:
    """
    Use PropertyStratifiedSplitter to produce an equivalent split on pIC50
    alone and compare it with our series-aware sample to quantify how well
    the pIC50 distribution is preserved.
    """
    try:
        from stratosampler.metrics.distribution import split_summary
    except ImportError:
        print("\n(stratosampler.metrics not available - skipping distribution report)")
        return

    # Build a reference split at the same size using PropertyStratifiedSplitter
    splitter = PropertyStratifiedSplitter(
        properties=["pIC50"],
        n_bins=len(PIC50_LABELS),
        test_size=len(sample_df) / len(full_df),
        random_state=RANDOM_STATE,
    )
    train_idx, test_idx = splitter.split(full_df, property_cols=["pIC50"])
    reference_sample = full_df.iloc[test_idx]

    print("\n-- Distribution quality (pIC50) --------------------------------------")
    print(f"  Full dataset   mean pIC50 = {full_df['pIC50'].mean():.3f}  "
          f"std = {full_df['pIC50'].std():.3f}")
    print(f"  Series sample  mean pIC50 = {sample_df['pIC50'].mean():.3f}  "
          f"std = {sample_df['pIC50'].std():.3f}")
    print(f"  pSS reference  mean pIC50 = {reference_sample['pIC50'].mean():.3f}  "
          f"std = {reference_sample['pIC50'].std():.3f}")

    # Use stratosampler's split_summary to score the sample vs the full set.
    # We treat "train = unsampled" and "test = sample" so the metrics describe
    # how faithfully the sample mirrors the full pIC50 distribution.
    sample_ids = set(sample_df["molecule_chembl_id"].values)
    in_sample = full_df["molecule_chembl_id"].isin(sample_ids).values
    test_positions  = np.where(in_sample)[0]
    train_positions = np.where(~in_sample)[0]

    try:
        summary = split_summary(
            full_df, train_positions, test_positions, property_cols=["pIC50"]
        )
        print("\n  stratosampler split_summary (series-aware sample vs remainder):")
        for k, v in summary.items():
            if k == "per_property":
                continue
            print(f"    {k}: {v:.4f}" if isinstance(v, float) else f"    {k}: {v}")
        if "per_property" in summary:
            print(f"    pIC50  KS={summary['per_property']['ks_stat'].iloc[0]:.4f}  "
                  f"JS={summary['per_property']['js_divergence'].iloc[0]:.4f}")
    except Exception as exc:
        print(f"  (split_summary error: {exc})")


# -- 6. Reporting ---------------------------------------------------------------

def report(full_df: pd.DataFrame, sample_df: pd.DataFrame) -> None:
    total = len(sample_df)
    print(f"\n{'='*60}")
    print(f"STRATIFIED SAMPLE  -  {TARGET_NAME} (ChEMBL {TARGET_CHEMBL_ID})")
    print(f"{'='*60}")
    print(f"Full dataset : {len(full_df):4d} unique compounds")
    print(f"Sample size  : {total:4d} ({100*total/len(full_df):.1f} % of full set)")

    print("\n-- Counts by series and pIC50 bin ------------------------------------")
    pivot = (
        sample_df.groupby(["series", "pIC50_bin"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=PIC50_LABELS, fill_value=0)
    )
    pivot["TOTAL"] = pivot.sum(axis=1)
    print(pivot.to_string())

    print("\n-- pIC50 statistics per series ---------------------------------------")
    stats = (
        sample_df.groupby("series")["pIC50"]
        .agg(n="count", mean="mean", std="std", min="min", max="max")
        .round(2)
    )
    print(stats.to_string())

    print("\n-- Top 5 most potent sampled compounds -------------------------------")
    top = sample_df.nlargest(5, "pIC50")[
        ["molecule_chembl_id", "series", "pIC50_bin", "pIC50", "smiles"]
    ]
    for _, row in top.iterrows():
        print(f"  {row['molecule_chembl_id']:14s}  {row['series']:20s}  "
              f"pIC50={row['pIC50']:.2f}  {row['smiles'][:50]}")


# -- Main -----------------------------------------------------------------------

def main() -> pd.DataFrame:
    # 1. Fetch
    df = fetch_chembl_activity(TARGET_CHEMBL_ID, max_records=MAX_RECORDS)

    # 2. Assign series
    df = assign_series(df, top_n=TOP_N_SERIES)

    # 3. Bin pIC50
    df = assign_pic50_bins(df)

    # 4. Stratified sample
    print(f"\nSampling (fraction={SAMPLE_FRACTION}, min_per_series={MIN_PER_SERIES}) ...")
    sample = stratified_sample(
        df,
        sample_fraction=SAMPLE_FRACTION,
        min_per_series=MIN_PER_SERIES,
        random_state=RANDOM_STATE,
    )

    # 5. Validate
    validate_with_stratosampler(df, sample)

    # 6. Report
    report(df, sample)

    # 7. Save
    out_path = Path(__file__).parent / "egfr_stratified_sample.csv"
    sample.to_csv(out_path, index=False)
    print(f"\nSample saved to: {out_path}")

    # 8. Plot
    fig = plot_series_sample(
        df, sample,
        pic50_bins=PIC50_BINS[1:-1],  # drop the 0 and 15 sentinels
        figsize=(13, 6),
    )
    plot_path = Path(__file__).parent / "egfr_stratified_sample_plot.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to:   {plot_path}")
    plt.show()

    return sample


if __name__ == "__main__":
    result = main()
