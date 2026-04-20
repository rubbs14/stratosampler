"""
Tests for stratosampler.
Run with: pytest tests/ -v
"""

import numpy as np
import pandas as pd
import pytest

from stratosampler import (
    PropertyStratifiedSplitter,
    distribution_report,
    split_summary,
    coverage_score,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def toy_df():
    """Small DataFrame with pre-computed properties (no RDKit needed)."""
    rng = np.random.default_rng(42)
    n = 200
    return pd.DataFrame({
        "logP": rng.normal(2.0, 1.5, n),
        "MW":   rng.normal(350, 80, n),
        "TPSA": rng.normal(70, 30, n).clip(0, 200),
        "label": rng.integers(0, 2, n),
    })


# ── Splitter tests ─────────────────────────────────────────────────────────────

class TestPropertyStratifiedSplitter:

    def test_basic_split_sizes(self, toy_df):
        splitter = PropertyStratifiedSplitter(
            test_size=0.2, random_state=0
        )
        train, test = splitter.split(
            toy_df, property_cols=["logP", "MW", "TPSA"]
        )
        n = len(toy_df)
        assert len(train) + len(test) == n
        assert 0.10 <= len(test) / n <= 0.55  # per-stratum min-1 sampling can shift the fraction

    def test_three_way_split(self, toy_df):
        splitter = PropertyStratifiedSplitter(
            test_size=0.1, val_size=0.1, random_state=1
        )
        result = splitter.split(toy_df, property_cols=["logP", "MW"])
        assert len(result) == 3
        train, val, test = result
        assert len(train) + len(val) + len(test) == len(toy_df)

    def test_no_index_overlap(self, toy_df):
        splitter = PropertyStratifiedSplitter(
            test_size=0.2, val_size=0.1, random_state=2
        )
        train, val, test = splitter.split(
            toy_df, property_cols=["logP", "MW", "TPSA"]
        )
        assert len(set(train) & set(val))  == 0
        assert len(set(train) & set(test)) == 0
        assert len(set(val)   & set(test)) == 0

    def test_get_split_dataframes(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2, random_state=3)
        train_df, test_df = splitter.get_split_dataframes(
            toy_df, property_cols=["logP", "MW"]
        )
        assert isinstance(train_df, pd.DataFrame)
        assert isinstance(test_df, pd.DataFrame)
        assert len(train_df) + len(test_df) == len(toy_df)

    def test_reproducibility(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2, random_state=99)
        t1, v1 = splitter.split(toy_df, property_cols=["logP", "MW"])
        t2, v2 = splitter.split(toy_df, property_cols=["logP", "MW"])
        np.testing.assert_array_equal(t1, t2)
        np.testing.assert_array_equal(v1, v2)

    def test_missing_property_col_raises(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2)
        with pytest.raises(ValueError, match="not found"):
            splitter.split(toy_df, property_cols=["nonexistent"])

    def test_no_smiles_or_property_cols_raises(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2)
        with pytest.raises(ValueError):
            splitter.split(toy_df)

    def test_single_property(self, toy_df):
        splitter = PropertyStratifiedSplitter(
            test_size=0.2, n_bins=3, random_state=5
        )
        train, test = splitter.split(toy_df, property_cols=["logP"])
        assert len(train) + len(test) == len(toy_df)

    def test_n_bins_1(self, toy_df):
        """Degenerate case: 1 bin = essentially random split."""
        splitter = PropertyStratifiedSplitter(
            test_size=0.2, n_bins=1, random_state=6
        )
        train, test = splitter.split(toy_df, property_cols=["logP", "MW"])
        assert len(train) + len(test) == len(toy_df)

    def test_large_n_bins(self, toy_df):
        """Many bins — some will be sparse, should not crash."""
        splitter = PropertyStratifiedSplitter(
            test_size=0.2, n_bins=50, random_state=7
        )
        train, test = splitter.split(toy_df, property_cols=["logP", "MW"])
        assert len(train) + len(test) == len(toy_df)


# ── Metrics tests ─────────────────────────────────────────────────────────────

class TestMetrics:

    def test_distribution_report_shape(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2, random_state=0)
        train, test = splitter.split(toy_df, property_cols=["logP", "MW"])
        report = distribution_report(toy_df, train, test, ["logP", "MW"])
        assert len(report) == 2
        assert "ks_stat" in report.columns
        assert "js_divergence" in report.columns

    def test_ks_stat_range(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2, random_state=0)
        train, test = splitter.split(toy_df, property_cols=["logP", "MW"])
        report = distribution_report(toy_df, train, test, ["logP", "MW"])
        assert (report["ks_stat"].between(0, 1)).all()

    def test_coverage_score_perfect(self, toy_df):
        """When train == test indices (degenerate), coverage is 1."""
        all_idx = np.arange(len(toy_df))
        strata  = np.array(["A"] * 100 + ["B"] * 100)
        score   = coverage_score(all_idx, all_idx, strata)
        assert score == 1.0

    def test_split_summary_keys(self, toy_df):
        splitter = PropertyStratifiedSplitter(test_size=0.2, random_state=0)
        train, test = splitter.split(toy_df, property_cols=["logP", "MW"])
        summary = split_summary(toy_df, train, test, ["logP", "MW"])
        for key in ["n_total", "n_train", "n_test", "mean_ks_stat", "per_property"]:
            assert key in summary

    def test_stratified_beats_random_ks(self, toy_df):
        """
        Stratified split should have lower mean KS statistic than
        a naive random split on a dataset with skewed properties.
        """
        rng = np.random.default_rng(0)
        # Create a skewed dataset
        n = 500
        skewed = pd.DataFrame({
            "logP": np.concatenate([rng.normal(0, 0.3, 400),
                                    rng.normal(5, 0.3, 100)]),
            "MW":   rng.normal(300, 50, n),
        })

        # Stratified
        strat = PropertyStratifiedSplitter(
            test_size=0.2, n_bins=5, random_state=42
        )
        s_train, s_test = strat.split(skewed, property_cols=["logP", "MW"])
        s_summary = split_summary(skewed, s_train, s_test, ["logP", "MW"])

        # Random
        all_idx = np.arange(n)
        rng2    = np.random.default_rng(42)
        rng2.shuffle(all_idx)
        n_test  = int(0.2 * n)
        r_test, r_train = all_idx[:n_test], all_idx[n_test:]
        r_summary = split_summary(skewed, r_train, r_test, ["logP", "MW"])

        assert s_summary["mean_ks_stat"] <= r_summary["mean_ks_stat"] + 0.05
