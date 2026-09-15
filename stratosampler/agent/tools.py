"""
Tool schemas (Claude API format) and handler functions for StratoAgent.
Each handler returns a JSON string. Errors are returned as {"error": "..."}.
"""

from __future__ import annotations

import json
import traceback
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Tool schemas
# ---------------------------------------------------------------------------

TOOL_SCHEMAS = [
    {
        "name": "load_data",
        "description": (
            "Load a molecular dataset from a CSV or SDF file and return metadata: "
            "shape, column names, dtypes, and a small sample. Call this first."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "Absolute or relative path to the CSV or SDF file.",
                },
                "smiles_col": {
                    "type": "string",
                    "description": "SMILES column name. Optional — provide if known.",
                },
            },
            "required": ["path"],
        },
    },
    {
        "name": "split_dataset",
        "description": (
            "Split a molecular dataset into train/test (and optionally validation) sets "
            "using property-stratified splitting. Saves outputs as CSV and returns sizes, "
            "output paths, and quality metrics."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the input CSV file."},
                "smiles_col": {"type": "string", "description": "SMILES column name."},
                "properties": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Properties to stratify on. Defaults to ['MolLogP', 'MolWt', 'TPSA'].",
                },
                "n_bins": {
                    "type": "integer",
                    "description": "Bins per property (default 5).",
                    "default": 5,
                },
                "test_size": {
                    "type": "number",
                    "description": "Test fraction (default 0.2).",
                    "default": 0.2,
                },
                "val_size": {
                    "type": "number",
                    "description": "Validation fraction (default 0 = no val set).",
                    "default": 0.0,
                },
                "scaffold_aware": {
                    "type": "boolean",
                    "description": "Keep Murcko scaffold families together (default false).",
                    "default": False,
                },
                "random_state": {
                    "type": "integer",
                    "description": "Random seed (default 42).",
                    "default": 42,
                },
                "output_dir": {
                    "type": "string",
                    "description": "Directory for output CSVs. Defaults to same dir as input.",
                },
            },
            "required": ["path", "smiles_col"],
        },
    },
    {
        "name": "compare_strategies",
        "description": (
            "Compare three splitting strategies on the same dataset: random, "
            "property-stratified, and scaffold-aware stratified. Returns KS and JS metrics per strategy."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the input CSV file."},
                "smiles_col": {"type": "string", "description": "SMILES column name."},
                "properties": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Properties to evaluate. Defaults to ['MolLogP', 'MolWt', 'TPSA'].",
                },
                "test_size": {
                    "type": "number",
                    "description": "Test fraction (default 0.2).",
                    "default": 0.2,
                },
                "random_state": {
                    "type": "integer",
                    "description": "Random seed (default 42).",
                    "default": 42,
                },
            },
            "required": ["path", "smiles_col"],
        },
    },
    {
        "name": "visualize_split",
        "description": (
            "Generate a property distribution plot comparing train and test sets. "
            "Saves the figure as PNG and returns the output path."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "train_path": {"type": "string", "description": "Path to the train CSV."},
                "test_path": {"type": "string", "description": "Path to the test CSV."},
                "smiles_col": {
                    "type": "string",
                    "description": "SMILES column name (needed to compute properties if not pre-computed).",
                },
                "properties": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Properties to plot. Defaults to ['MolLogP', 'MolWt', 'TPSA'].",
                },
                "output_path": {
                    "type": "string",
                    "description": "Where to save the PNG. Defaults to same dir as train CSV.",
                },
            },
            "required": ["train_path", "test_path"],
        },
    },
    {
        "name": "compute_properties",
        "description": (
            "Compute molecular properties from SMILES and return summary statistics "
            "(mean, std, min, max, n_valid) for each property."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to a CSV file."},
                "smiles_col": {"type": "string", "description": "SMILES column name."},
                "properties": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Properties to compute. Defaults to all built-in properties.",
                },
            },
            "required": ["path", "smiles_col"],
        },
    },
    {
        "name": "fetch_pdb_structures",
        "description": (
            "Search the RCSB Protein Data Bank for 3D crystal structures of a protein target "
            "(e.g. EGFR, BRAF, CDK2, KRAS). Returns PDB IDs, titles, resolution, and method. "
            "The webapp renders an interactive 3D viewer automatically."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_name": {
                    "type": "string",
                    "description": "Gene name or protein name (e.g. 'EGFR', 'BRAF', 'CDK2').",
                },
                "max_results": {
                    "type": "integer",
                    "description": "Max structures to return (default 5, max 10).",
                    "default": 5,
                },
                "pdb_ids": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Explicit PDB IDs to fetch directly (skips search).",
                },
            },
            "required": ["target_name"],
        },
    },
    {
        "name": "find_mcs",
        "description": (
            "Find the Maximum Common Substructure (MCS) across a list of SMILES strings "
            "using one or more matching configs. Returns MCS size, coverage, and Jaccard score."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles_list": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings (2 or more).",
                },
                "config_names": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "MCS config names to run. Omit to run all presets.",
                },
            },
            "required": ["smiles_list"],
        },
    },
]


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------

def _load_data(path: str, smiles_col: str | None = None) -> str:
    p = Path(path)
    if not p.exists():
        return json.dumps({"error": f"File not found: {path}"})
    try:
        if p.suffix.lower() == ".sdf":
            from stratosampler import load_sdf
            _, df = load_sdf(path)
        else:
            df = pd.read_csv(path)

        result: dict = {
            "shape": list(df.shape),
            "columns": df.columns.tolist(),
            "dtypes": df.dtypes.astype(str).to_dict(),
            "sample": df.head(3).fillna("").to_dict(orient="records"),
        }
        if smiles_col and smiles_col in df.columns:
            result["smiles_col"] = smiles_col
            result["n_valid_smiles"] = int(df[smiles_col].dropna().shape[0])
        return json.dumps(result)
    except Exception as exc:
        return json.dumps({"error": str(exc)})


def _split_dataset(
    path: str,
    smiles_col: str,
    properties: list[str] | None = None,
    n_bins: int = 5,
    test_size: float = 0.2,
    val_size: float = 0.0,
    scaffold_aware: bool = False,
    random_state: int = 42,
    output_dir: str | None = None,
) -> str:
    from stratosampler import PropertyStratifiedSplitter, compute_properties
    from stratosampler.metrics.distribution import split_summary

    try:
        df = pd.read_csv(path)
        properties = properties or ["MolLogP", "MolWt", "TPSA"]

        splitter = PropertyStratifiedSplitter(
            properties=properties,
            n_bins=n_bins,
            test_size=test_size,
            val_size=val_size,
            scaffold_aware=scaffold_aware,
            random_state=random_state,
        )
        indices = splitter.split(df, smiles_col=smiles_col)

        prop_df = compute_properties(df[smiles_col].tolist(), properties)
        strata = splitter._assign_strata(prop_df)
        df_eval = pd.concat([df.reset_index(drop=True), prop_df], axis=1)

        out_dir = Path(output_dir) if output_dir else Path(path).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = Path(path).stem

        val_idx = None
        if val_size > 0:
            train_idx, val_idx, test_idx = indices
            val_path = out_dir / f"{stem}_val.csv"
            df.iloc[val_idx].to_csv(val_path, index=False)
        else:
            train_idx, test_idx = indices
            val_path = None

        train_path = out_dir / f"{stem}_train.csv"
        test_path = out_dir / f"{stem}_test.csv"
        df.iloc[train_idx].to_csv(train_path, index=False)
        df.iloc[test_idx].to_csv(test_path, index=False)

        avail = [p for p in properties if p in df_eval.columns]
        metrics = split_summary(
            df_eval, train_idx, test_idx, property_cols=avail, val_idx=val_idx, strata=strata
        )

        metrics_out = {k: v for k, v in metrics.items() if k != "per_property"}
        metrics_out["per_property"] = metrics["per_property"].to_dict(orient="records")

        result: dict = {
            "train_size": len(train_idx),
            "test_size": len(test_idx),
            "train_path": str(train_path),
            "test_path": str(test_path),
            "metrics": metrics_out,
        }
        if val_path:
            result["val_size"] = len(val_idx)
            result["val_path"] = str(val_path)

        return json.dumps(result)
    except Exception as exc:
        return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})


def _compare_strategies(
    path: str,
    smiles_col: str,
    properties: list[str] | None = None,
    test_size: float = 0.2,
    random_state: int = 42,
) -> str:
    from stratosampler import PropertyStratifiedSplitter, compute_properties
    from stratosampler.metrics.distribution import split_summary

    try:
        df = pd.read_csv(path)
        properties = properties or ["MolLogP", "MolWt", "TPSA"]

        prop_df = compute_properties(df[smiles_col].tolist(), properties)
        df_eval = pd.concat([df.reset_index(drop=True), prop_df], axis=1)
        avail = [p for p in properties if p in df_eval.columns]

        n = len(df)
        rng = np.random.default_rng(random_state)

        results: dict = {}

        # Random
        idx = np.arange(n)
        rng.shuffle(idx)
        n_test = int(test_size * n)
        m = split_summary(df_eval, idx[n_test:], idx[:n_test], avail)
        results["random"] = {
            "mean_ks_stat": round(m["mean_ks_stat"], 4),
            "mean_js_div": round(m["mean_js_div"], 4),
            "n_train": m["n_train"],
            "n_test": m["n_test"],
        }

        # Stratified
        strat = PropertyStratifiedSplitter(
            properties=properties, test_size=test_size, random_state=random_state
        )
        tr, te = strat.split(df, smiles_col=smiles_col)
        m = split_summary(df_eval, tr, te, avail)
        results["stratified"] = {
            "mean_ks_stat": round(m["mean_ks_stat"], 4),
            "mean_js_div": round(m["mean_js_div"], 4),
            "n_train": m["n_train"],
            "n_test": m["n_test"],
        }

        # Scaffold-aware stratified
        sc = PropertyStratifiedSplitter(
            properties=properties, test_size=test_size, scaffold_aware=True, random_state=random_state
        )
        tr, te = sc.split(df, smiles_col=smiles_col)
        m = split_summary(df_eval, tr, te, avail)
        results["scaffold_stratified"] = {
            "mean_ks_stat": round(m["mean_ks_stat"], 4),
            "mean_js_div": round(m["mean_js_div"], 4),
            "n_train": m["n_train"],
            "n_test": m["n_test"],
        }

        return json.dumps({
            "properties_evaluated": avail,
            "strategies": results,
            "note": "Lower KS stat and JS divergence = better distribution match between train and test.",
        })
    except Exception as exc:
        return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})


def _visualize_split(
    train_path: str,
    test_path: str,
    smiles_col: str | None = None,
    properties: list[str] | None = None,
    output_path: str | None = None,
) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from stratosampler.visualisation.plots import plot_property_distributions

    try:
        train_df = pd.read_csv(train_path)
        test_df = pd.read_csv(test_path)
        properties = properties or ["MolLogP", "MolWt", "TPSA"]

        # Compute properties if not present and smiles_col is available
        if smiles_col:
            from stratosampler import compute_properties
            missing = [p for p in properties if p not in train_df.columns]
            if missing and smiles_col in train_df.columns:
                tp = compute_properties(train_df[smiles_col].tolist(), missing)
                ep = compute_properties(test_df[smiles_col].tolist(), missing)
                train_df = pd.concat([train_df.reset_index(drop=True), tp], axis=1)
                test_df = pd.concat([test_df.reset_index(drop=True), ep], axis=1)

        n_train = len(train_df)
        combined = pd.concat([train_df, test_df], ignore_index=True)
        train_idx = np.arange(n_train)
        test_idx = np.arange(n_train, len(combined))
        avail = [p for p in properties if p in combined.columns]

        fig = plot_property_distributions(combined, train_idx, test_idx, avail)

        out = output_path or str(Path(train_path).parent / "split_distributions.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)

        return json.dumps({"figure_path": out, "properties": avail})
    except Exception as exc:
        return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})


def _compute_properties(
    path: str,
    smiles_col: str,
    properties: list[str] | None = None,
) -> str:
    from stratosampler import compute_properties, BUILTIN_PROPERTIES

    try:
        df = pd.read_csv(path)
        properties = properties or list(BUILTIN_PROPERTIES.keys())

        prop_df = compute_properties(df[smiles_col].tolist(), properties)

        stats: dict = {}
        for col in prop_df.columns:
            vals = prop_df[col].dropna()
            stats[col] = {
                "mean": round(float(vals.mean()), 4),
                "std": round(float(vals.std()), 4),
                "min": round(float(vals.min()), 4),
                "max": round(float(vals.max()), 4),
                "n_valid": int(len(vals)),
                "n_missing": int(prop_df[col].isna().sum()),
            }
        return json.dumps(stats)
    except Exception as exc:
        return json.dumps({"error": str(exc)})


def _fetch_pdb_structures(
    target_name: str,
    max_results: int = 5,
    pdb_ids: list[str] | None = None,
) -> str:
    import urllib.request as _req

    try:
        max_results = min(max_results, 10)

        if pdb_ids:
            ids = [i.upper().strip() for i in pdb_ids[:max_results]]
        else:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": target_name},
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {"start": 0, "rows": max_results},
                    "sort": [{"sort_by": "score", "direction": "desc"}],
                },
            }
            data = json.dumps(query).encode()
            req = _req.Request(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                data=data,
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            with _req.urlopen(req, timeout=10) as resp:
                search_result = json.loads(resp.read())
            ids = [hit["identifier"] for hit in search_result.get("result_set", [])]

        if not ids:
            return json.dumps({"error": f"No PDB structures found for '{target_name}'."})

        structures = []
        for pdb_id in ids:
            try:
                url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                with _req.urlopen(url, timeout=5) as resp:
                    meta = json.loads(resp.read())
                exptl = meta.get("exptl", [{}])
                res_list = meta.get("rcsb_entry_info", {}).get("resolution_combined", [])
                structures.append({
                    "pdb_id": pdb_id,
                    "title": meta.get("struct", {}).get("title", ""),
                    "method": exptl[0].get("method", "") if exptl else "",
                    "resolution": round(res_list[0], 2) if res_list else None,
                    "deposit_date": meta.get("rcsb_accession_info", {}).get("deposit_date", "")[:10],
                })
            except Exception:
                structures.append({"pdb_id": pdb_id, "title": "", "method": "", "resolution": None, "deposit_date": ""})

        return json.dumps({
            "target": target_name,
            "n_found": len(structures),
            "pdb_ids": [s["pdb_id"] for s in structures],
            "structures": structures,
            "viewer_ready": True,
        })
    except Exception as exc:
        return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})


def _find_mcs(
    smiles_list: list[str],
    config_names: list[str] | None = None,
) -> str:
    from rdkit import Chem
    from stratosampler.mcs.graph_mcs import compare_configs, PRESETS

    try:
        mols = []
        invalid = []
        for i, smi in enumerate(smiles_list):
            m = Chem.MolFromSmiles(smi)
            if m is None:
                invalid.append(i)
            else:
                mols.append(m)

        if len(mols) < 2:
            return json.dumps({"error": "Need at least 2 valid SMILES."})

        if config_names:
            unknown = [n for n in config_names if n not in PRESETS]
            if unknown:
                return json.dumps({"error": f"Unknown configs: {unknown}. Available: {list(PRESETS)}"})
            configs = [PRESETS[n] for n in config_names]
        else:
            configs = list(PRESETS.values())

        df = compare_configs(mols, configs=configs)
        drop_cols = [c for c in ("mapping", "mol_sizes") if c in df.columns]
        result_df = df.drop(columns=drop_cols)

        return json.dumps({
            "n_mols": len(mols),
            "n_invalid": len(invalid),
            "results": result_df.to_dict(orient="records"),
        })
    except Exception as exc:
        return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})


# ---------------------------------------------------------------------------
# OpenAI-compatible schema conversion (used by Groq backend)
# ---------------------------------------------------------------------------

def _to_openai(schema: dict) -> dict:
    return {
        "type": "function",
        "function": {
            "name": schema["name"],
            "description": schema["description"],
            "parameters": schema["input_schema"],
        },
    }

OPENAI_TOOL_SCHEMAS = [_to_openai(s) for s in TOOL_SCHEMAS]


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

_HANDLERS = {
    "load_data": _load_data,
    "split_dataset": _split_dataset,
    "compare_strategies": _compare_strategies,
    "visualize_split": _visualize_split,
    "compute_properties": _compute_properties,
    "fetch_pdb_structures": _fetch_pdb_structures,
    "find_mcs": _find_mcs,
}


def handle_tool_call(name: str, inputs: dict) -> str:
    handler = _HANDLERS.get(name)
    if handler is None:
        return json.dumps({"error": f"Unknown tool: {name}"})
    try:
        return handler(**inputs)
    except Exception as exc:
        return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})
