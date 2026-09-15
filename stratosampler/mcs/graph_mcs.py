"""
Graph-based Maximum Common Substructure (MCS) via correspondence graph + max clique.

The algorithm:
1. Convert each molecule to a NetworkX graph (atoms=nodes, bonds=edges).
2. Build the correspondence graph of two molecular graphs:
   - Nodes: compatible atom pairs (u, v) — u ∈ G1, v ∈ G2
   - Edges: (u1,v1)-(u2,v2) iff bond u1-u2 ∈ G1 ↔ bond v1-v2 ∈ G2 (connectivity parity)
     + if both bonds exist, the bond_compat predicate must also pass
3. Maximum clique in the correspondence graph = MCS atom mapping.
4. For N molecules: iterate, maintaining a running atom-correspondence table.

Pluggable atom_compat / bond_compat callables let you swap matching semantics
(exact element, aromatic-aware, bond-order-free, Murcko, etc.) without touching
the core algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
NodeData = Dict
EdgeData = Dict
NodeCompatFn = Callable[[NodeData, NodeData], bool]
EdgeCompatFn = Callable[[EdgeData, EdgeData], bool]


# ---------------------------------------------------------------------------
# Molecule → graph
# ---------------------------------------------------------------------------

def mol_to_graph(mol: Chem.Mol) -> nx.Graph:
    """Convert an RDKit molecule to a NetworkX graph.

    Node attributes: symbol, atomic_num, is_aromatic, formal_charge
    Edge attributes: bond_type (float: 1.0 single, 2.0 double, 1.5 aromatic, 3.0 triple),
                     is_in_ring
    """
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            symbol=atom.GetSymbol(),
            atomic_num=atom.GetAtomicNum(),
            is_aromatic=atom.GetIsAromatic(),
            formal_charge=atom.GetFormalCharge(),
        )
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondTypeAsDouble(),
            is_in_ring=bond.IsInRing(),
        )
    return G


# ---------------------------------------------------------------------------
# Murcko scaffold helper
# ---------------------------------------------------------------------------

def to_murcko(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Return the Murcko scaffold of mol, or None if mol has no rings."""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return None
    return scaffold


# ---------------------------------------------------------------------------
# Atom compatibility functions
# ---------------------------------------------------------------------------

def atom_exact(a1: NodeData, a2: NodeData) -> bool:
    """Same element symbol (C==C, N==N, …)."""
    return a1["symbol"] == a2["symbol"]


def atom_any(a1: NodeData, a2: NodeData) -> bool:
    """Any heavy atom matches any other — maximises common topology."""
    return True


def atom_aromatic(a1: NodeData, a2: NodeData) -> bool:
    """Same element AND same aromaticity flag."""
    return a1["symbol"] == a2["symbol"] and a1["is_aromatic"] == a2["is_aromatic"]


def atom_full(a1: NodeData, a2: NodeData) -> bool:
    """Same element, aromaticity, and formal charge — strictest atom typing."""
    return (
        a1["symbol"] == a2["symbol"]
        and a1["is_aromatic"] == a2["is_aromatic"]
        and a1["formal_charge"] == a2["formal_charge"]
    )


# ---------------------------------------------------------------------------
# Bond compatibility functions
# ---------------------------------------------------------------------------

def bond_any(b1: EdgeData, b2: EdgeData) -> bool:
    """Bond presence only — ignore bond order entirely."""
    return True


def bond_strict(b1: EdgeData, b2: EdgeData) -> bool:
    """Exact bond type must match (1.0, 1.5, 2.0, 3.0)."""
    return b1["bond_type"] == b2["bond_type"]


def bond_ring_aware(b1: EdgeData, b2: EdgeData) -> bool:
    """Bond type match AND ring membership must agree."""
    return b1["bond_type"] == b2["bond_type"] and b1["is_in_ring"] == b2["is_in_ring"]


# ---------------------------------------------------------------------------
# MCSConfig — bundles atom+bond matchers + optional Murcko preprocessing
# ---------------------------------------------------------------------------

@dataclass
class MCSConfig:
    name: str
    atom_compat: NodeCompatFn = atom_exact
    bond_compat: EdgeCompatFn = bond_strict
    use_murcko: bool = False
    connected: bool = False      # restrict MCS to largest connected component (required for 3D alignment)
    use_rdkit_fmcs: bool = False  # use RDKit's C++ FindMCS instead of correspondence graph (fast for drug-like molecules)
    description: str = ""

    def __repr__(self) -> str:
        murcko_tag = " [Murcko]" if self.use_murcko else ""
        conn_tag   = " [connected]" if self.connected else ""
        rdkit_tag  = " [rdFMCS]" if self.use_rdkit_fmcs else ""
        return f"MCSConfig({self.name!r}{murcko_tag}{conn_tag}{rdkit_tag})"


PRESETS: Dict[str, MCSConfig] = {
    "strict": MCSConfig(
        "strict", atom_exact, bond_strict,
        description="Exact element + exact bond order",
    ),
    "element_only": MCSConfig(
        "element_only", atom_exact, bond_any,
        description="Exact element, ignore bond order",
    ),
    "aromatic": MCSConfig(
        "aromatic", atom_aromatic, bond_strict,
        description="Element + aromaticity flag, exact bond order",
    ),
    "full_atom": MCSConfig(
        "full_atom", atom_full, bond_strict,
        description="Element + aromaticity + charge, exact bond order",
    ),
    "heavy_atom": MCSConfig(
        "heavy_atom", atom_any, bond_any,
        description="Any heavy atom, ignore bond order — maximum topology overlap",
    ),
    "ring_aware": MCSConfig(
        "ring_aware", atom_exact, bond_ring_aware,
        description="Exact element + bond order + ring membership must agree",
    ),
    "murcko_strict": MCSConfig(
        "murcko_strict", atom_exact, bond_strict, use_murcko=True,
        description="Strict matching on Murcko scaffolds only",
    ),
    "murcko_element": MCSConfig(
        "murcko_element", atom_exact, bond_any, use_murcko=True,
        description="Element-only matching on Murcko scaffolds",
    ),
    "murcko_heavy": MCSConfig(
        "murcko_heavy", atom_any, bond_any, use_murcko=True,
        description="Maximum topology overlap on Murcko scaffolds",
    ),
    # --- alignment-oriented presets (rdFMCS C++ backend + connected) ---------
    "align_strict": MCSConfig(
        "align_strict", atom_exact, bond_strict,
        connected=True, use_rdkit_fmcs=True,
        description="Largest connected MCS, exact element + bond order — 3D alignment core",
    ),
    "align_element": MCSConfig(
        "align_element", atom_exact, bond_any,
        connected=True, use_rdkit_fmcs=True,
        description="Largest connected MCS, exact element, any bond order",
    ),
    "align_murcko": MCSConfig(
        "align_murcko", atom_exact, bond_strict,
        use_murcko=True, connected=True, use_rdkit_fmcs=True,
        description="Largest connected MCS on Murcko scaffolds — fast alignment for large molecules",
    ),
}


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def build_correspondence_graph(
    G1: nx.Graph,
    G2: nx.Graph,
    atom_compat: NodeCompatFn,
    bond_compat: EdgeCompatFn,
) -> nx.Graph:
    """Build the correspondence (modular product) graph for G1 and G2.

    A clique in the returned graph maps directly to a common induced subgraph.
    Edges encode connectivity parity: two node pairs are connected iff
    - both pairs have a bond (same bond type per bond_compat), OR
    - neither pair has a bond (non-edge consistency required for induced MCS).

    Edges are built with numpy adjacency-matrix operations to avoid an O(k²)
    pure-Python loop (k = number of compatible atom pairs, up to n1*n2 ≈ 1500
    for drug-like molecules under permissive atom matching).
    """
    nodes1 = list(G1.nodes())
    nodes2 = list(G2.nodes())
    n1, n2 = len(nodes1), len(nodes2)

    # --- compatible atom pairs -----------------------------------------------
    compatible: List[Tuple[int, int]] = [
        (u, v)
        for u in nodes1
        for v in nodes2
        if atom_compat(G1.nodes[u], G2.nodes[v])
    ]
    if not compatible:
        return nx.Graph()

    C = nx.Graph()
    C.add_nodes_from(compatible)

    # --- adjacency matrices for vectorised edge building ---------------------
    idx1 = {n: i for i, n in enumerate(nodes1)}
    idx2 = {n: i for i, n in enumerate(nodes2)}

    A1 = np.zeros((n1, n1), dtype=np.float32)
    A2 = np.zeros((n2, n2), dtype=np.float32)
    for u, v, d in G1.edges(data=True):
        A1[idx1[u], idx1[v]] = A1[idx1[v], idx1[u]] = d["bond_type"]
    for u, v, d in G2.edges(data=True):
        A2[idx2[u], idx2[v]] = A2[idx2[v], idx2[u]] = d["bond_type"]

    # Map compatible pairs to their row/col positions in the adjacency matrices
    cp = np.array(compatible, dtype=np.int32)          # shape (k, 2) — original atom indices
    ci1 = np.array([idx1[u] for u, _ in compatible], dtype=np.int32)  # positions in A1
    ci2 = np.array([idx2[v] for _, v in compatible], dtype=np.int32)  # positions in A2

    k = len(compatible)
    # Submatrices: C_e1[i,j] = A1[ci1[i], ci1[j]], etc.
    C_e1 = A1[np.ix_(ci1, ci1)]   # bond types between compatible G1 atoms
    C_e2 = A2[np.ix_(ci2, ci2)]   # bond types between compatible G2 atoms

    both_edge  = (C_e1 > 0) & (C_e2 > 0)   # bond in both
    no_edge    = (C_e1 == 0) & (C_e2 == 0)  # bond in neither

    # bond_compat check for "both edge" pairs.
    # For bond_any the type is irrelevant; for bond_strict types must match.
    # We approximate: call bond_compat once with representative dicts to detect
    # type-agnostic configs, otherwise apply element-wise type comparison.
    _dummy = {"bond_type": 1.0, "is_in_ring": False}
    _dummy2 = {"bond_type": 2.0, "is_in_ring": False}
    _type_matters = not bond_compat(_dummy, _dummy2)

    if _type_matters:
        types_match = (C_e1 == C_e2)
        valid = (both_edge & types_match) | no_edge
    else:
        valid = both_edge | no_edge

    # Mask out same-atom pairs (injective mapping constraint)
    same_u = (cp[:, 0:1] == cp[:, 0])   # shape (k, k)
    same_v = (cp[:, 1:2] == cp[:, 1])
    valid  = valid & ~same_u & ~same_v

    # Only upper triangle to avoid duplicate edges; bulk-insert for speed
    ii, jj = np.where(np.triu(valid, k=1))
    C.add_edges_from((compatible[i], compatible[j]) for i, j in zip(ii.tolist(), jj.tolist()))

    return C


_EXACT_NODE_LIMIT    = 120   # correspondence-graph nodes above which we approximate
_EXACT_DENSITY_LIMIT = 0.55  # edge density above which BK is exponentially slow
_GREEDY_RESTARTS     = 12    # number of degree-ordered restarts for the greedy heuristic


def _greedy_max_clique(C: nx.Graph) -> List:
    """Degree-ordered greedy max clique with multiple restarts.

    Much faster than networkx's ramsey_R2 for large sparse graphs:
    O(k * restarts) where k = clique size.  Typically within 1–2 atoms of
    optimal for drug-like correspondence graphs.
    """
    if not C:
        return []

    adj = {n: set(C.neighbors(n)) for n in C.nodes()}
    top_starts = sorted(adj, key=lambda n: len(adj[n]), reverse=True)
    top_starts = top_starts[:_GREEDY_RESTARTS]

    best: List = []
    for start in top_starts:
        clique = [start]
        candidates = set(adj[start])
        while candidates:
            # pick candidate with the highest sub-degree within current candidates
            v = max(candidates, key=lambda n: len(adj[n] & candidates))
            clique.append(v)
            candidates &= adj[v]
        if len(clique) > len(best):
            best = clique

    return best


def _max_clique(C: nx.Graph) -> List[Tuple[int, int]]:
    """Return the largest clique in C.

    Strategy by regime:
    - Small / sparse (n ≤ 120, density ≤ 0.55): exact Bron-Kerbosch.
    - Very dense (density > 0.75): exact BK on the *complement* graph
      (MCS ≈ full graph minus a small independent set — BK terminates fast).
    - Otherwise (large, medium density): degree-ordered greedy heuristic
      (lower bound, typically ≤2 atoms from exact for drug-like molecules).
    """
    if len(C) == 0:
        return []

    n = len(C)
    ne = len(C.edges())
    max_edges = n * (n - 1) / 2
    density = ne / max_edges if max_edges > 0 else 0.0

    if n <= _EXACT_NODE_LIMIT and density <= _EXACT_DENSITY_LIMIT:
        return max(nx.find_cliques(C), key=len, default=[])

    if density > 0.75:
        # Dense: exact BK converges quickly (complement is sparse)
        return max(nx.find_cliques(C), key=len, default=[])

    import warnings
    warnings.warn(
        f"Correspondence graph ({n} nodes, density={density:.2f}) falls back to "
        "greedy max-clique approximation — result is a lower bound.",
        stacklevel=4,
    )
    return _greedy_max_clique(C)


def find_mcs_pair(
    G1: nx.Graph,
    G2: nx.Graph,
    atom_compat: NodeCompatFn,
    bond_compat: EdgeCompatFn,
) -> List[Tuple[int, int]]:
    """Return the MCS atom mapping between two graphs as a list of (idx_G1, idx_G2) pairs."""
    C = build_correspondence_graph(G1, G2, atom_compat, bond_compat)
    return _max_clique(C)


def find_mcs_iterative(
    mols: List[Chem.Mol],
    atom_compat: NodeCompatFn,
    bond_compat: EdgeCompatFn,
) -> List[List[int]]:
    """Find the MCS across N ≥ 2 molecules by iterative intersection.

    Returns a list of per-molecule atom index lists that all map onto the MCS.
    E.g. result[0] are the MCS atom indices in mols[0], result[1] in mols[1], etc.
    Indices refer to atoms in the molecules passed in (post Murcko if applicable).

    Strategy:
    - Maintain a running atom-correspondence table: each row is a tuple of atom
      indices (one per molecule seen so far) that correspond to each other.
    - For each new molecule, reduce the table to the atoms that survive MCS with it.
    """
    if len(mols) < 2:
        raise ValueError("Need at least 2 molecules")

    graphs = [mol_to_graph(m) for m in mols]

    # Seed with first pair
    clique = find_mcs_pair(graphs[0], graphs[1], atom_compat, bond_compat)
    # mapping: list of tuples (idx_mol0, idx_mol1, …)
    mapping: List[tuple] = [tuple(pair) for pair in clique]

    for k in range(2, len(mols)):
        if not mapping:
            break

        mol0_indices = [row[0] for row in mapping]

        # Induced subgraph of mol0 on current MCS atoms, re-indexed to 0..m
        G0_sub = graphs[0].subgraph(mol0_indices)
        orig_to_new = {orig: new for new, orig in enumerate(mol0_indices)}
        G0_re = nx.relabel_nodes(G0_sub, orig_to_new)

        clique = find_mcs_pair(G0_re, graphs[k], atom_compat, bond_compat)

        surviving_new = {p[0] for p in clique}
        new_to_molk = {p[0]: p[1] for p in clique}

        new_mapping = []
        for row in mapping:
            new_idx = orig_to_new[row[0]]
            if new_idx in surviving_new:
                new_mapping.append(row + (new_to_molk[new_idx],))

        mapping = new_mapping

    if not mapping:
        return [[] for _ in mols]

    return [[row[i] for row in mapping] for i in range(len(mols))]


# ---------------------------------------------------------------------------
# RDKit FMCS backend (fast path for drug-like molecules)
# ---------------------------------------------------------------------------

_RDKIT_ATOM_COMPARE = {
    "exact":    rdFMCS.AtomCompare.CompareElements,
    "any":      rdFMCS.AtomCompare.CompareAny,
    "isotopes": rdFMCS.AtomCompare.CompareIsotopes,
}
_RDKIT_BOND_COMPARE = {
    "strict":   rdFMCS.BondCompare.CompareOrder,
    "any":      rdFMCS.BondCompare.CompareAny,
    "stereo":   rdFMCS.BondCompare.CompareOrderExact,
}


def find_mcs_rdkit(
    mols: List[Chem.Mol],
    atom_compare: str = "exact",
    bond_compare: str = "strict",
    ring_matches_ring_only: bool = True,
    complete_rings_only: bool = False,
    timeout: int = 10,
) -> List[List[int]]:
    """Find MCS using RDKit's C++ FindMCS engine — fast for drug-like molecules.

    Returns per-molecule atom index lists (same shape as find_mcs_iterative).
    Uses GetSubstructMatch to recover atom indices from the SMARTS result.
    """
    ac = _RDKIT_ATOM_COMPARE.get(atom_compare, rdFMCS.AtomCompare.CompareElements)
    bc = _RDKIT_BOND_COMPARE.get(bond_compare, rdFMCS.BondCompare.CompareOrder)

    result = rdFMCS.FindMCS(
        mols,
        atomCompare=ac,
        bondCompare=bc,
        ringMatchesRingOnly=ring_matches_ring_only,
        completeRingsOnly=complete_rings_only,
        timeout=timeout,
    )

    if result.numAtoms == 0:
        return [[] for _ in mols]

    query = Chem.MolFromSmarts(result.smartsString)
    if query is None:
        return [[] for _ in mols]

    matches = [m.GetSubstructMatch(query) for m in mols]
    if any(len(m) == 0 for m in matches):
        return [[] for _ in mols]

    return [list(m) for m in matches]


# ---------------------------------------------------------------------------
# Multi-config comparison
# ---------------------------------------------------------------------------

def compare_configs(
    mols: List[Chem.Mol],
    mol_names: Optional[List[str]] = None,
    configs: Optional[List[MCSConfig]] = None,
) -> pd.DataFrame:
    """Run MCS with each MCSConfig and return a comparison DataFrame.

    Columns
    -------
    config          : config name
    description     : human-readable description
    mcs_size        : number of atoms in the MCS
    mol_sizes       : list of heavy-atom counts per molecule (after Murcko if applicable)
    min_mol_size    : smallest molecule size (denominator for coverage)
    coverage        : mcs_size / min_mol_size
    jaccard         : mcs_size / (sum of mol_sizes - (n_mols-1)*mcs_size)  [generalised]
    approximate     : True if any clique step fell back to greedy approximation
    mapping         : raw per-molecule atom index lists
    """
    import warnings as _warnings

    if configs is None:
        configs = list(PRESETS.values())
    if mol_names is None:
        mol_names = [f"mol{i}" for i in range(len(mols))]

    rows = []
    for cfg in configs:
        working_mols: List[Chem.Mol] = []
        for m in mols:
            if cfg.use_murcko:
                s = to_murcko(m)
                working_mols.append(s if s is not None else m)
            else:
                working_mols.append(m)

        approximate = False
        try:
            if cfg.use_rdkit_fmcs:
                # Fast C++ path — uses rdFMCS; connectivity is handled by ring_matches_ring_only
                ac = "any" if cfg.atom_compat is atom_any else "exact"
                bc = "any" if cfg.bond_compat is bond_any else "strict"
                mapping = find_mcs_rdkit(working_mols, atom_compare=ac, bond_compare=bc)
            else:
                with _warnings.catch_warnings(record=True) as caught:
                    _warnings.simplefilter("always")
                    mapping = find_mcs_iterative(working_mols, cfg.atom_compat, cfg.bond_compat)
                approximate = any("greedy" in str(w.message) for w in caught)
        except Exception as exc:
            mapping = [[] for _ in mols]
            print(f"[{cfg.name}] failed: {exc}")

        # Restrict to largest connected component when requested (required for 3D alignment).
        # For rdFMCS with ring_matches_ring_only=True the result is typically already connected,
        # but we apply the filter uniformly for correctness.
        if cfg.connected and mapping and mapping[0]:
            G0   = mol_to_graph(working_mols[0])
            sub  = G0.subgraph(mapping[0])
            ccs  = sorted(nx.connected_components(sub), key=len, reverse=True)
            if len(ccs) > 1:
                keep      = ccs[0]
                keep_pos  = [i for i, a in enumerate(mapping[0]) if a in keep]
                mapping   = [[lst[i] for i in keep_pos] for lst in mapping]

        mcs_size = len(mapping[0]) if mapping and mapping[0] else 0
        mol_sizes = [m.GetNumHeavyAtoms() for m in working_mols]
        min_size = min(mol_sizes)
        coverage = mcs_size / min_size if min_size > 0 else 0.0

        # Generalised Jaccard: |MCS| / |union estimate|
        n = len(mols)
        union_est = sum(mol_sizes) - (n - 1) * mcs_size
        jaccard = mcs_size / union_est if union_est > 0 else 0.0

        rows.append(
            {
                "config": cfg.name,
                "description": cfg.description,
                "mcs_size": mcs_size,
                "mol_sizes": mol_sizes,
                "min_mol_size": min_size,
                "coverage": round(coverage, 4),
                "jaccard": round(jaccard, 4),
                "approximate": approximate,
                "mapping": mapping,
            }
        )

    return pd.DataFrame(rows)
