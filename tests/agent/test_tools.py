"""
Tests for stratosampler.agent.tools — handler functions and dispatcher.

No network access and no GROQ_API_KEY are required: `_fetch_pdb_structures`
mocks `urllib.request.urlopen` directly, and none of these tests touch
stratosampler.agent.agent (the Groq-backed StratoAgent).
"""

import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rdkit")

from stratosampler.agent.tools import (
    _load_data,
    _split_dataset,
    _compare_strategies,
    _compute_properties,
    _fetch_pdb_structures,
    _find_mcs,
    handle_tool_call,
)
from stratosampler.mcs.graph_mcs import PRESETS
from stratosampler.splitters.property_stratified import BUILTIN_PROPERTIES


# ── Shared fixtures ──────────────────────────────────────────────────────────

# A diverse pool of realistic, small drug-like/common SMILES: aliphatic
# alcohols/acids/amines, simple aromatics, heteroaromatics, and a few
# well-known drugs. Enough structural variety that property stratification
# and MCS comparisons behave meaningfully rather than degenerating.
DIVERSE_SMILES = [
    "CCO",                                  # ethanol
    "CCCO",                                 # propanol
    "CCCCO",                                # butanol
    "CCCCCO",                               # pentanol
    "CC(C)O",                               # isopropanol
    "CC(C)(C)O",                            # tert-butanol
    "CCOCC",                                # diethyl ether
    "CCC(=O)O",                             # propanoic acid
    "CCCC(=O)O",                            # butanoic acid
    "CC(=O)O",                              # acetic acid
    "c1ccccc1",                             # benzene
    "Cc1ccccc1",                            # toluene
    "Clc1ccccc1",                           # chlorobenzene
    "Cc1ccc(C)cc1",                         # p-xylene
    "c1ccncc1",                             # pyridine
    "c1ccoc1",                              # furan
    "c1ccsc1",                              # thiophene
    "CC(=O)Oc1ccccc1C(=O)O",                # aspirin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",         # caffeine
    "CC(=O)Nc1ccc(O)cc1",                   # paracetamol
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",           # ibuprofen
    "CCN(CC)CC",                            # triethylamine
    "CCN",                                  # ethylamine
    "CCCN",                                 # propylamine
    "NCCO",                                 # ethanolamine
    "OCCO",                                 # ethylene glycol
    "OCC(O)CO",                             # glycerol
    "CC#N",                                 # acetonitrile
    "CCOC(=O)C",                            # ethyl acetate
    "c1ccc2ccccc2c1",                       # naphthalene
]


@pytest.fixture
def diverse_csv(tmp_path):
    """A ~30-row molecular dataset CSV for split/compare/property tests."""
    df = pd.DataFrame({
        "SMILES": DIVERSE_SMILES,
        "activity": np.linspace(0.1, 9.5, len(DIVERSE_SMILES)),
    })
    path = tmp_path / "molecules.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def small_csv(tmp_path):
    """A tiny 3-row CSV for basic load_data tests."""
    df = pd.DataFrame({
        "SMILES": ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
        "activity": [0.1, 0.5, 0.9],
    })
    path = tmp_path / "small.csv"
    df.to_csv(path, index=False)
    return path


# ── _load_data ───────────────────────────────────────────────────────────────

class TestLoadData:
    def test_happy_path_no_smiles_col(self, small_csv):
        result = json.loads(_load_data(str(small_csv)))
        assert "error" not in result
        assert result["shape"] == [3, 2]
        assert result["columns"] == ["SMILES", "activity"]
        assert isinstance(result["dtypes"], dict)
        assert len(result["sample"]) == 3
        # No smiles_col requested -> these keys should be absent.
        assert "smiles_col" not in result
        assert "n_valid_smiles" not in result

    def test_happy_path_with_valid_smiles_col(self, small_csv):
        result = json.loads(_load_data(str(small_csv), smiles_col="SMILES"))
        assert "error" not in result
        assert result["smiles_col"] == "SMILES"
        assert result["n_valid_smiles"] == 3

    def test_smiles_col_not_in_columns_is_silently_ignored(self, small_csv):
        """smiles_col validity check: an unknown column name doesn't error,
        it just skips the smiles-specific fields (per the `if smiles_col and
        smiles_col in df.columns` guard)."""
        result = json.loads(_load_data(str(small_csv), smiles_col="not_a_column"))
        assert "error" not in result
        assert "smiles_col" not in result
        assert "n_valid_smiles" not in result

    def test_missing_file_returns_error(self, tmp_path):
        missing = tmp_path / "does_not_exist.csv"
        result = json.loads(_load_data(str(missing)))
        assert result == {"error": f"File not found: {missing}"}

    def test_sdf_path_loads_dataframe(self, tmp_path):
        """Regression test: `load_sdf` returns a (mols, DataFrame) tuple;
        `_load_data` must unpack it rather than treating the tuple as the
        DataFrame directly (previously raised AttributeError internally)."""
        from rdkit import Chem

        sdf_path = tmp_path / "mols.sdf"
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", "ethanol")
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        result = json.loads(_load_data(str(sdf_path)))
        assert "error" not in result
        assert result["shape"][0] == 1


# ── _split_dataset ───────────────────────────────────────────────────────────

class TestSplitDataset:
    def test_happy_path(self, diverse_csv):
        result = json.loads(_split_dataset(
            str(diverse_csv), smiles_col="SMILES",
            properties=["MolWt", "TPSA"], n_bins=3,
        ))
        assert "error" not in result, result

        n = len(DIVERSE_SMILES)
        assert result["train_size"] + result["test_size"] == n
        assert Path(result["train_path"]).exists()
        assert Path(result["test_path"]).exists()
        assert "val_size" not in result
        assert "val_path" not in result

        train_df = pd.read_csv(result["train_path"])
        test_df = pd.read_csv(result["test_path"])
        assert len(train_df) == result["train_size"]
        assert len(test_df) == result["test_size"]

        metrics = result["metrics"]
        for key in ("n_total", "n_train", "n_test", "train_frac", "test_frac",
                    "mean_ks_stat", "mean_js_div", "per_property", "coverage_score"):
            assert key in metrics
        assert metrics["n_total"] == n
        assert isinstance(metrics["mean_ks_stat"], float)
        assert isinstance(metrics["mean_js_div"], float)
        assert 0.0 <= metrics["coverage_score"] <= 1.0
        assert isinstance(metrics["per_property"], list)
        assert len(metrics["per_property"]) == 2
        for row in metrics["per_property"]:
            assert set(row) >= {"property", "train_mean", "test_mean",
                                 "ks_stat", "ks_pval", "js_divergence"}

    def test_scaffold_aware_branch(self, diverse_csv):
        result = json.loads(_split_dataset(
            str(diverse_csv), smiles_col="SMILES",
            properties=["MolWt"], n_bins=3, scaffold_aware=True,
        ))
        assert "error" not in result, result
        n = len(DIVERSE_SMILES)
        assert result["train_size"] + result["test_size"] == n
        assert Path(result["train_path"]).exists()
        assert Path(result["test_path"]).exists()

    def test_val_size_branch(self, diverse_csv):
        result = json.loads(_split_dataset(
            str(diverse_csv), smiles_col="SMILES",
            properties=["MolWt"], n_bins=3,
            test_size=0.2, val_size=0.2,
        ))
        assert "error" not in result, result
        n = len(DIVERSE_SMILES)
        assert "val_size" in result
        assert "val_path" in result
        assert Path(result["val_path"]).exists()
        assert result["train_size"] + result["val_size"] + result["test_size"] == n

        metrics = result["metrics"]
        assert "n_val" in metrics
        assert "val_frac" in metrics
        for row in metrics["per_property"]:
            assert "val_mean" in row
            assert "val_ks_stat" in row
            assert "val_js_divergence" in row

    def test_output_dir_override(self, diverse_csv, tmp_path):
        out_dir = tmp_path / "outputs"
        result = json.loads(_split_dataset(
            str(diverse_csv), smiles_col="SMILES",
            properties=["MolWt"], n_bins=3, output_dir=str(out_dir),
        ))
        assert "error" not in result, result
        assert result["train_path"].startswith(str(out_dir))
        assert result["test_path"].startswith(str(out_dir))
        assert Path(result["train_path"]).exists()

    def test_missing_smiles_col_in_csv_returns_error(self, diverse_csv):
        result = json.loads(_split_dataset(str(diverse_csv), smiles_col="nope"))
        assert "error" in result
        assert "traceback" in result



# ── _compare_strategies ──────────────────────────────────────────────────────

class TestCompareStrategies:
    def test_happy_path(self, diverse_csv):
        result = json.loads(_compare_strategies(
            str(diverse_csv), smiles_col="SMILES", properties=["MolWt"],
        ))
        assert "error" not in result, result
        assert result["properties_evaluated"] == ["MolWt"]
        assert "note" in result

        strategies = result["strategies"]
        n = len(DIVERSE_SMILES)
        for name in ("random", "stratified", "scaffold_stratified"):
            assert name in strategies
            strat = strategies[name]
            assert "mean_ks_stat" in strat
            assert "mean_js_div" in strat
            assert isinstance(strat["mean_ks_stat"], float)
            assert isinstance(strat["mean_js_div"], float)
            assert strat["n_train"] + strat["n_test"] == n

    def test_missing_smiles_col_returns_error(self, diverse_csv):
        result = json.loads(_compare_strategies(str(diverse_csv), smiles_col="nope"))
        assert "error" in result
        assert "traceback" in result


# ── _compute_properties ──────────────────────────────────────────────────────

class TestComputeProperties:
    def test_happy_path_default_properties(self, diverse_csv):
        result = json.loads(_compute_properties(str(diverse_csv), smiles_col="SMILES"))
        assert "error" not in result, result
        assert set(result.keys()) == set(BUILTIN_PROPERTIES.keys())

        n = len(DIVERSE_SMILES)
        for prop, stats in result.items():
            assert set(stats.keys()) == {"mean", "std", "min", "max", "n_valid", "n_missing"}
            assert stats["n_valid"] + stats["n_missing"] == n
            assert stats["n_missing"] == 0  # all SMILES in the fixture are valid
            assert stats["min"] <= stats["mean"] <= stats["max"]

    def test_explicit_property_subset(self, diverse_csv):
        result = json.loads(_compute_properties(
            str(diverse_csv), smiles_col="SMILES", properties=["MolWt", "TPSA"],
        ))
        assert "error" not in result, result
        assert set(result.keys()) == {"MolWt", "TPSA"}

    def test_invalid_smiles_counted_as_missing(self, tmp_path):
        df = pd.DataFrame({
            "SMILES": ["CCO", "c1ccccc1", "not_a_smiles", "also_not_valid", "CCN"],
        })
        path = tmp_path / "mixed.csv"
        df.to_csv(path, index=False)

        result = json.loads(_compute_properties(str(path), smiles_col="SMILES", properties=["MolWt"]))
        assert "error" not in result, result
        assert result["MolWt"]["n_valid"] == 3
        assert result["MolWt"]["n_missing"] == 2

    def test_missing_smiles_col_returns_error(self, diverse_csv):
        result = json.loads(_compute_properties(str(diverse_csv), smiles_col="nope"))
        assert "error" in result


# ── _find_mcs ────────────────────────────────────────────────────────────────

class TestFindMcs:
    def test_happy_path_default_configs_runs_all_presets(self):
        result = json.loads(_find_mcs(["CCO", "CCN"]))
        assert "error" not in result, result
        assert result["n_mols"] == 2
        assert result["n_invalid"] == 0
        assert len(result["results"]) == len(PRESETS)
        for row in result["results"]:
            assert set(row) >= {"config", "description", "mcs_size",
                                 "min_mol_size", "coverage", "jaccard", "approximate"}
            assert "mapping" not in row
            assert "mol_sizes" not in row

    def test_happy_path_explicit_config_names(self):
        result = json.loads(_find_mcs(
            ["c1ccccc1", "Cc1ccccc1"], config_names=["strict", "element_only"],
        ))
        assert "error" not in result, result
        assert result["n_mols"] == 2
        configs_returned = {row["config"] for row in result["results"]}
        assert configs_returned == {"strict", "element_only"}
        # Benzene ring is fully shared -> MCS should cover at least 6 atoms.
        for row in result["results"]:
            assert row["mcs_size"] >= 6

    def test_skips_invalid_smiles_but_counts_them(self):
        result = json.loads(_find_mcs(
            ["CCO", "not_a_smiles", "CCN"], config_names=["strict"],
        ))
        assert "error" not in result, result
        assert result["n_mols"] == 2
        assert result["n_invalid"] == 1

    def test_fewer_than_two_valid_smiles_returns_error(self):
        result = json.loads(_find_mcs(["not_a_smiles", "CCO"]))
        assert result == {"error": "Need at least 2 valid SMILES."}

    def test_all_invalid_smiles_returns_error(self):
        result = json.loads(_find_mcs(["nope", "still_nope"]))
        assert result == {"error": "Need at least 2 valid SMILES."}

    def test_unknown_config_name_returns_error(self):
        result = json.loads(_find_mcs(["CCO", "CCN"], config_names=["bogus_config"]))
        assert "error" in result
        assert "Unknown configs" in result["error"]
        assert "bogus_config" in result["error"]


# ── _fetch_pdb_structures (network mocked) ──────────────────────────────────

class _FakeResponse:
    """Minimal stand-in for the object returned by urllib.request.urlopen."""

    def __init__(self, payload: dict):
        self._payload = json.dumps(payload).encode()

    def read(self):
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def _metadata_payload(title="A structure", method="X-RAY DIFFRACTION",
                       resolution=1.8, deposit_date="2020-01-01T00:00:00Z"):
    return {
        "struct": {"title": title},
        "exptl": [{"method": method}],
        "rcsb_entry_info": {"resolution_combined": [resolution]},
        "rcsb_accession_info": {"deposit_date": deposit_date},
    }


class TestFetchPdbStructures:
    def test_explicit_pdb_ids_skips_search(self, monkeypatch):
        calls = []

        def fake_urlopen(url_or_req, timeout=None):
            calls.append(url_or_req)
            assert not isinstance(url_or_req, urllib.request.Request), \
                "explicit pdb_ids must skip the RCSB search endpoint"
            return _FakeResponse(_metadata_payload())

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        result = json.loads(_fetch_pdb_structures(
            target_name="EGFR", pdb_ids=["1m17", "2itp"],
        ))
        assert "error" not in result, result
        assert result["pdb_ids"] == ["1M17", "2ITP"]
        assert result["n_found"] == 2
        assert result["viewer_ready"] is True
        assert len(calls) == 2
        for s in result["structures"]:
            assert s["method"] == "X-RAY DIFFRACTION"
            assert s["resolution"] == 1.8
            assert s["deposit_date"] == "2020-01-01"

    def test_search_path_used_when_no_pdb_ids(self, monkeypatch):
        def fake_urlopen(url_or_req, timeout=None):
            if isinstance(url_or_req, urllib.request.Request):
                return _FakeResponse({
                    "result_set": [{"identifier": "5C9C"}, {"identifier": "3OG7"}],
                })
            return _FakeResponse(_metadata_payload(title="BRAF kinase domain"))

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        result = json.loads(_fetch_pdb_structures(target_name="BRAF", max_results=2))
        assert "error" not in result, result
        assert result["target"] == "BRAF"
        assert result["pdb_ids"] == ["5C9C", "3OG7"]
        assert result["n_found"] == 2
        assert all(s["title"] == "BRAF kinase domain" for s in result["structures"])

    def test_search_with_no_hits_returns_error(self, monkeypatch):
        def fake_urlopen(url_or_req, timeout=None):
            return _FakeResponse({"result_set": []})

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        result = json.loads(_fetch_pdb_structures(target_name="totally_unknown_target"))
        assert "error" in result

    def test_per_structure_metadata_failure_degrades_gracefully(self, monkeypatch):
        def fake_urlopen(url_or_req, timeout=None):
            if isinstance(url_or_req, urllib.request.Request):
                return _FakeResponse({
                    "result_set": [{"identifier": "1ABC"}, {"identifier": "2XYZ"}],
                })
            if "1ABC" in url_or_req:
                raise TimeoutError("simulated network failure")
            return _FakeResponse(_metadata_payload(title="Good structure", resolution=1.5))

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        result = json.loads(_fetch_pdb_structures(target_name="KRAS", max_results=2))
        assert "error" not in result, result
        assert result["n_found"] == 2

        by_id = {s["pdb_id"]: s for s in result["structures"]}
        assert by_id["1ABC"] == {
            "pdb_id": "1ABC", "title": "", "method": "", "resolution": None, "deposit_date": "",
        }
        assert by_id["2XYZ"]["title"] == "Good structure"
        assert by_id["2XYZ"]["resolution"] == 1.5


# ── handle_tool_call dispatcher ──────────────────────────────────────────────

class TestHandleToolCall:
    def test_unknown_tool_returns_error(self):
        result = json.loads(handle_tool_call("no_such_tool", {}))
        assert result == {"error": "Unknown tool: no_such_tool"}

    def test_missing_required_kwarg_is_caught_as_typeerror(self):
        # load_data requires `path`; calling with no kwargs raises TypeError
        # at the handler(**inputs) call site, which handle_tool_call catches.
        result = json.loads(handle_tool_call("load_data", {}))
        assert "error" in result
        assert "traceback" in result
        assert "path" in result["error"]

    def test_unexpected_kwarg_is_caught_as_typeerror(self):
        result = json.loads(handle_tool_call("find_mcs", {"smiles_list": ["CCO", "CCN"], "bogus_kwarg": 1}))
        assert "error" in result
        assert "traceback" in result

    def test_known_tool_dispatches_successfully(self, diverse_csv):
        result = json.loads(handle_tool_call(
            "compute_properties",
            {"path": str(diverse_csv), "smiles_col": "SMILES", "properties": ["MolWt"]},
        ))
        assert "error" not in result, result
        assert "MolWt" in result
