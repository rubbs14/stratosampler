---
name: stratosampler-mcs
description: Guidance for Maximum Common Substructure (MCS) analysis in stratosampler — choosing the right preset config (strict, element_only, aromatic, full_atom, heavy_atom, ring_aware, murcko_strict/element/heavy, align_strict/element/murcko), calling compare_configs() directly in Python vs the find_mcs agent tool, and interpreting mcs_size/coverage/jaccard output. Use this whenever the user asks about MCS, maximum common substructure, shared scaffolds, SAR (structure-activity relationship) analysis across a molecule series, scaffold identification/clustering, R-group variation, or preparing molecules for 3D alignment/overlay — even if they just paste a list of SMILES and ask "what do these have in common" or "align these compounds." Also use it when someone is confused about which MCS preset to pick, why two configs give different mcs_size, or how the find_mcs tool's output differs from calling compare_configs directly.
---

# stratosampler MCS analysis

stratosampler's MCS engine (`stratosampler/mcs/graph_mcs.py`) finds the largest common
substructure across a set of molecules using either a correspondence-graph max-clique
algorithm or RDKit's C++ `rdFMCS` engine, depending on the preset. The behavior of "what
counts as a match" is fully pluggable via `MCSConfig` — atom compatibility, bond
compatibility, optional Murcko-scaffold pre-reduction, and whether the result must be a
single connected fragment.

There are two ways to run it:
- **Direct Python**: `compare_configs(mols, configs=...)` — full DataFrame, including the
  raw atom-index `mapping` needed for 3D alignment.
- **Agent tool**: `find_mcs` (in `stratosampler/agent/tools.py`) — takes SMILES strings
  instead of RDKit mols, and the JSON it returns has `mapping` and `mol_sizes` stripped
  out. Use this path only when you need agent-tool JSON output; use direct Python whenever
  you need the atom mapping itself (e.g. for alignment).

## Choosing a preset

All presets live in `PRESETS: Dict[str, MCSConfig]` in `graph_mcs.py`. Each is a
`MCSConfig(name, atom_compat, bond_compat, use_murcko, connected, use_rdkit_fmcs, description)`.
Pick based on what question you're actually asking:

| Preset | Atom match | Bond match | Murcko | Connected + fast (rdFMCS) | When to use |
|---|---|---|---|---|---|
| `strict` | exact element | exact bond order | no | no | Default/tightest. Confirms the shared substructure is chemically identical, not just same atoms wired differently. |
| `element_only` | exact element | any | no | no | Loosen bond order — use when comparing tautomers/resonance forms/different kekulization where atom connectivity is the same but bond-order perception differs. |
| `aromatic` | element + aromaticity | exact | no | no | Like `strict` but also requires aromaticity to agree, so a benzene ring won't match a cyclohexane ring of the same element pattern. |
| `full_atom` | element + aromaticity + formal charge | exact | no | no | Strictest atom typing. Use when charge state matters (e.g. carboxylate vs. protonated acid should NOT be treated as the same match). |
| `heavy_atom` | any heavy atom | any | no | no | Maximum permissiveness — ignores element and bond order, matches on topology alone. Use for scaffold-hopping-style comparisons across very different chemotypes, or to find the largest common molecular "skeleton" shape. |
| `ring_aware` | exact element | exact bond order + ring membership must agree | no | no | Like `strict` but a shared bond must also agree on ring-vs-chain context — prevents matching a ring bond to a chemically-equivalent open-chain bond. |
| `murcko_strict` | exact element | exact | **yes** | no | Strict matching, but run on Murcko scaffolds (ring systems + linkers, substituents stripped) — for scaffold identification ignoring R-group diversity. |
| `murcko_element` | exact element | any | **yes** | no | Same idea, looser bond order — scaffold-level tautomer tolerance. |
| `murcko_heavy` | any | any | **yes** | no | Maximum scaffold-topology overlap — "do these molecules share the same ring-system shape at all," ignoring atom identity and bond order. |
| `align_strict` | exact element | exact | no | **yes** | 3D-alignment core preset. Runs on the fast RDKit `rdFMCS` C++ backend and restricts the result to the single largest **connected** atom set — a disconnected mapping can't be used as one rigid-body atom correspondence for overlay. |
| `align_element` | exact element | any | no | **yes** | Same as `align_strict` but tolerant of bond-order differences (e.g. aligning across tautomers). |
| `align_murcko` | exact element | exact | **yes** | **yes** | Scaffold-level alignment — fast for large/macrocyclic molecules where you only need to superpose the ring system. |

Practical picking rules:
- **SAR analysis** (what varies across a congeneric series / R-group survey): start with
  `element_only` or `aromatic` to see the shared core, then compare against `strict` — a
  big gap between `strict` and `element_only` mcs_size tells you bond-order/tautomer
  perception is masking a real match. `ring_aware` is useful when you need to be sure a
  shared fragment sits in the same ring/chain context across the series.
- **Scaffold identification / clustering**: use the `murcko_*` family. `murcko_strict` for
  a tight scaffold match, `murcko_heavy` when you just want to know if the ring topology
  is shared at all regardless of heteroatom placement.
- **3D alignment prep**: always use one of the `align_*` presets. They're the only ones
  with `connected=True`, which is a hard requirement for alignment — you need one
  contiguous fragment to compute a rigid-body superposition, not several disconnected
  matched islands. They also route through `find_mcs_rdkit()` (RDKit's C++ engine) rather
  than the correspondence-graph clique search, so they're fast even on large molecules.

You can also build a custom `MCSConfig` from the exported building blocks
(`atom_exact`, `atom_any`, `atom_aromatic`, `atom_full`, `bond_any`, `bond_strict`,
`bond_ring_aware`) if none of the 12 presets fit — see `graph_mcs.py` for the signatures.

## Calling it directly in Python

```python
from rdkit import Chem
from stratosampler.mcs.graph_mcs import compare_configs, PRESETS

mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]  # caller must parse SMILES and drop None

# Run every preset:
df = compare_configs(mols)

# Or just the ones you care about:
df = compare_configs(mols, configs=[PRESETS["strict"], PRESETS["murcko_strict"], PRESETS["align_strict"]])

# Optional: label rows with molecule names instead of mol0, mol1, ...
df = compare_configs(mols, mol_names=["compound_A", "compound_B"], configs=[PRESETS["align_strict"]])
```

`compare_configs(mols, mol_names=None, configs=None)` requires already-parsed
`Chem.Mol` objects (invalid-SMILES filtering is the caller's job — `compare_configs`
itself does no validation). If `configs` is omitted it runs **all 12 presets**.

For 3D alignment, pull the `mapping` column — it's a list of per-molecule atom-index
lists, aligned positionally across molecules (`mapping[0]` are the matched atom indices
in `mols[0]`, `mapping[1]` in `mols[1]`, etc.), ready to feed into an atom-map-based
alignment call (e.g. RDKit's `rdMolAlign.AlignMol` / `GetO3A` with an `atomMap`):

```python
row = df[df["config"] == "align_strict"].iloc[0]
atom_map = list(zip(row["mapping"][0], row["mapping"][1]))  # (idx_in_mol0, idx_in_mol1) pairs
```

## Calling it via the agent tool (`find_mcs`)

Defined in `stratosampler/agent/tools.py`. Schema:

```json
{
  "name": "find_mcs",
  "input_schema": {
    "smiles_list": ["array of SMILES strings, 2 or more"],
    "config_names": ["optional array of preset names; omit to run all presets"]
  }
}
```

The `_find_mcs` handler: parses each SMILES with `Chem.MolFromSmiles`, collects indices of
anything that fails to parse as `invalid`, requires at least 2 valid molecules, validates
any given `config_names` against `PRESETS` (unknown names return an error listing valid
ones), calls `compare_configs(mols, configs=configs)`, then **drops the `mapping` and
`mol_sizes` columns** before returning JSON:

```json
{"n_mols": 2, "n_invalid": 0, "results": [{"config": "strict", "description": "...", "mcs_size": 12, "min_mol_size": 15, "coverage": 0.8, "jaccard": 0.6667, "approximate": false}, ...]}
```

Because `mapping` is stripped, **never use the `find_mcs` tool output as the source of an
atom correspondence for alignment** — call `compare_configs` directly in Python for that.
The agent tool is for reporting/summary use (SAR discussion, scaffold comparison) where
you only need the size/coverage/Jaccard numbers, not the underlying atom indices.

## Interpreting the output columns

`compare_configs` returns one row per config with these columns:

- **`config`** — preset name.
- **`description`** — the preset's human-readable description (from `MCSConfig.description`).
- **`mcs_size`** — number of atoms in the found MCS.
- **`mol_sizes`** — heavy-atom count per molecule, *after* Murcko reduction if the preset
  uses it (dropped by the `find_mcs` agent tool, present when calling `compare_configs`
  directly).
- **`min_mol_size`** — the smallest molecule's heavy-atom count; this is the denominator
  for `coverage`.
- **`coverage`** — `mcs_size / min_mol_size`. Fraction of the *smallest* molecule that's
  contained in the shared substructure. `1.0` means the smallest molecule is fully
  embedded in the others (e.g. it's a substructure of a larger analogue).
- **`jaccard`** — generalized similarity across *all* N molecules:
  `mcs_size / (sum(mol_sizes) - (n_mols - 1) * mcs_size)`. Unlike `coverage`, this
  penalizes size mismatches across the whole set, not just the smallest molecule — use it
  when comparing more than 2 molecules or when molecule sizes vary a lot.
  A jaccard of `1.0` means the molecules are essentially identical in the compared region;
  low jaccard with high coverage means one molecule is a small fragment embedded in much
  larger ones.
- **`approximate`** — `True` if the max-clique search fell back to a greedy heuristic
  (only relevant for non-`align_*`/non-rdFMCS presets: the correspondence graph exceeds
  120 nodes or 0.55 edge density, per the thresholds in `graph_mcs.py`). When `True`,
  treat `mcs_size` as a **lower bound**, not the exact optimum — this shows up most often
  with the permissive presets (`heavy_atom`, `element_only`, murcko variants) on larger or
  more numerous molecules, since permissive atom matching blows up the correspondence
  graph. The `align_*` presets avoid this entirely by using RDKit's C++ `rdFMCS` engine.
- **`mapping`** — raw per-molecule atom index lists (only present via direct
  `compare_configs` calls, not via the `find_mcs` agent tool). This is what you need for
  3D alignment.

## Common use cases at a glance

1. **SAR analysis** — run `element_only`/`aromatic` alongside `strict` on a congeneric
   series; compare `mcs_size`/`coverage` across configs to see how much of the divergence
   is real substituent variation vs. bond-order/tautomer artifacts. Low `jaccard` across
   the set flags a genuinely diverse series rather than simple R-group substitution.
2. **Scaffold identification** — run the `murcko_*` presets to compare ring systems in
   isolation; `murcko_heavy` for a permissive "same ring topology at all" check,
   `murcko_strict` for a tight scaffold match used for clustering compounds into series.
3. **3D alignment prep** — run an `align_*` preset, pull `mapping` from the DataFrame
   returned by `compare_configs` (not from the `find_mcs` tool, which strips it), and use
   it as the atom map for a rigid-body superposition. `align_strict` is the default;
   `align_element` for cross-tautomer alignment; `align_murcko` for large/macrocyclic
   molecules where scaffold-level alignment is sufficient and faster.
