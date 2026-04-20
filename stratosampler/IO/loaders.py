"""
Molecular Loaders
-----------------
Convenient wrappers for loading molecules from various file formats using RDKit.

Supports:
  - SMILES files (CSV, TSV, plain text)
  - SDF (Structure Data Format) files
  - Custom delimiter and column specifications

Examples
--------
Load molecules from SMILES file::

    from stratosampler.IO import load_smiles

    # Simple SMILES file (one per line)
    mols, smiles = load_smiles("compounds.smi")

    # CSV with SMILES in column 1
    mols, smiles = load_smiles("compounds.csv", delimiter=",", smiles_column=1)

Load molecules from SDF file::

    from stratosampler.IO import load_sdf

    # SDF file with optional properties
    mols = load_sdf("compounds.sdf")
    mols_with_props = load_sdf("compounds.sdf", include_properties=True)
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional, Union

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    _RDKIT = True
except ImportError:
    _RDKIT = False
    warnings.warn(
        "RDKit not installed. Install with: pip install 'stratosampler[rdkit]'"
    )


class MoleculeLoader:
    """Base class for molecular file loaders."""

    def __init__(self, raise_on_invalid: bool = False):
        """
        Initialize the loader.

        Parameters
        ----------
        raise_on_invalid : bool, default=False
            If True, raise an exception when invalid molecules are encountered.
            If False, skip invalid molecules with a warning.
        """
        self.raise_on_invalid = raise_on_invalid
        self._validate_rdkit()

    @staticmethod
    def _validate_rdkit() -> None:
        """Validate that RDKit is available."""
        if not _RDKIT:
            raise ImportError(
                "RDKit is required. Install with: pip install 'stratosampler[rdkit]'"
            )

    @staticmethod
    def _validate_molecule(
        mol: Optional["Chem.Mol"],
        identifier: str,
        raise_on_invalid: bool = False,
    ) -> Optional["Chem.Mol"]:
        """
        Validate a molecule object.

        Parameters
        ----------
        mol : Chem.Mol or None
            RDKit molecule object
        identifier : str
            Identifier for error messages (e.g., SMILES string, compound ID)
        raise_on_invalid : bool
            Whether to raise on invalid molecules

        Returns
        -------
        Chem.Mol or None
            Valid molecule or None if invalid
        """
        if mol is None:
            msg = f"Invalid molecule: {identifier}"
            if raise_on_invalid:
                raise ValueError(msg)
            else:
                warnings.warn(msg, stacklevel=3)
            return None
        return mol


class SmilesLoader(MoleculeLoader):
    """Load molecules from SMILES files."""

    def load(
        self,
        filepath: Union[str, Path],
        smiles_column: int = 0,
        delimiter: Optional[str] = None,
        id_column: Optional[int] = None,
        keep_properties: bool = False,
        sanitize: bool = True,
        add_hydrogens: bool = False,
        compute_2d_coords: bool = False,
    ) -> tuple[list["Chem.Mol"], pd.DataFrame]:
        """
        Load molecules from a SMILES file.

        Parameters
        ----------
        filepath : str or Path
            Path to SMILES file
        smiles_column : int, default=0
            Column index containing SMILES strings
        delimiter : str, optional
            Column delimiter. If None, auto-detect (space, comma, tab)
        id_column : int, optional
            Column index for compound IDs. If None, use row index
        keep_properties : bool, default=False
            If True, keep all columns from file as metadata
        sanitize : bool, default=True
            Whether to sanitize molecules
        add_hydrogens : bool, default=False
            Whether to add explicit hydrogens
        compute_2d_coords : bool, default=False
            Whether to compute 2D coordinates

        Returns
        -------
        mols : list of Chem.Mol
            RDKit molecule objects
        data : pd.DataFrame
            DataFrame with SMILES, IDs, and optional properties

        Examples
        --------
        >>> loader = SmilesLoader()
        >>> mols, data = loader.load("compounds.smi")
        >>> print(len(mols), "molecules loaded")

        >>> # With column specification
        >>> mols, data = loader.load(
        ...     "compounds.csv",
        ...     delimiter=",",
        ...     smiles_column=1,
        ...     id_column=0,
        ... )
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        # Auto-detect delimiter if not specified
        if delimiter is None:
            delimiter = self._detect_delimiter(filepath)

        # Read file
        df = pd.read_csv(filepath, delimiter=delimiter, header=None)

        if smiles_column >= len(df.columns):
            raise ValueError(
                f"SMILES column {smiles_column} not found. "
                f"File has {len(df.columns)} columns."
            )

        # Extract SMILES
        smiles_list = df.iloc[:, smiles_column].astype(str).tolist()

        # Extract IDs
        if id_column is not None:
            if id_column >= len(df.columns):
                raise ValueError(f"ID column {id_column} not found.")
            ids = df.iloc[:, id_column].tolist()
        else:
            ids = [f"mol_{i}" for i in range(len(smiles_list))]

        # Load molecules
        molecules = []
        valid_rows = []

        for idx, (smi, mol_id) in enumerate(zip(smiles_list, ids)):
            mol = Chem.MolFromSmiles(smi, sanitize=sanitize)

            if mol is None:
                if self.raise_on_invalid:
                    raise ValueError(f"Invalid SMILES at row {idx}: {smi}")
                else:
                    warnings.warn(f"Skipping invalid SMILES at row {idx}: {smi}")
                    continue

            # Add hydrogens if requested
            if add_hydrogens:
                mol = Chem.AddHs(mol)

            # Compute 2D coordinates if requested
            if compute_2d_coords:
                AllChem.Compute2DCoords(mol)

            molecules.append(mol)
            valid_rows.append(idx)

        # Build result dataframe
        if keep_properties:
            result_df = df.iloc[valid_rows].copy()
        else:
            result_df = pd.DataFrame(
                {
                    "mol_id": [ids[i] for i in valid_rows],
                    "smiles": [smiles_list[i] for i in valid_rows],
                }
            )

        if id_column is None and not keep_properties:
            # Add SMILES and ID if not already included
            pass  # Already added above

        return molecules, result_df

    @staticmethod
    def _detect_delimiter(filepath: Path) -> str:
        """
        Auto-detect delimiter from file.

        Parameters
        ----------
        filepath : Path
            Path to file

        Returns
        -------
        str
            Detected delimiter (space, comma, or tab)
        """
        with open(filepath, "r") as f:
            first_line = f.readline().strip()

        # Count occurrences of each delimiter
        delimiters = {",": first_line.count(","), "\t": first_line.count("\t")}

        # Default to space if no clear delimiter found
        if delimiters[","] > 0:
            return ","
        elif delimiters["\t"] > 0:
            return "\t"
        else:
            return " "


class SdfLoader(MoleculeLoader):
    """Load molecules from SDF (Structure Data Format) files."""

    def load(
        self,
        filepath: Union[str, Path],
        include_properties: bool = True,
        sanitize: bool = True,
        remove_hs: bool = True,
        compute_2d_coords: bool = False,
    ) -> tuple[list["Chem.Mol"], pd.DataFrame]:
        """
        Load molecules from an SDF file.

        Parameters
        ----------
        filepath : str or Path
            Path to SDF file
        include_properties : bool, default=True
            If True, extract SDF properties into DataFrame
        sanitize : bool, default=True
            Whether to sanitize molecules
        remove_hs : bool, default=True
            Whether to remove explicit hydrogens
        compute_2d_coords : bool, default=False
            Whether to compute 2D coordinates

        Returns
        -------
        mols : list of Chem.Mol
            RDKit molecule objects
        data : pd.DataFrame
            DataFrame with molecule properties

        Examples
        --------
        >>> loader = SdfLoader()
        >>> mols, data = loader.load("compounds.sdf")
        >>> print(len(mols), "molecules loaded")
        >>> print(data.columns)
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        # Read SDF
        suppl = Chem.SDMolSupplier(str(filepath), sanitize=sanitize, removeHs=remove_hs)

        molecules = []
        properties_list = []
        valid_indices = []

        for idx, mol in enumerate(suppl):
            if mol is None:
                if self.raise_on_invalid:
                    raise ValueError(f"Invalid molecule at index {idx} in SDF")
                else:
                    warnings.warn(f"Skipping invalid molecule at index {idx} in SDF")
                    continue

            # Compute 2D coordinates if requested
            if compute_2d_coords:
                AllChem.Compute2DCoords(mol)

            molecules.append(mol)
            valid_indices.append(idx)

            # Extract properties if requested
            if include_properties:
                props = {"mol_id": f"mol_{idx}"}
                if mol.HasProp("_Name"):
                    props["name"] = mol.GetProp("_Name")

                # Get all SDF properties
                for key in mol.GetPropNames():
                    if not key.startswith("_"):
                        try:
                            props[key] = mol.GetProp(key)
                        except Exception:
                            props[key] = None

                properties_list.append(props)
            else:
                properties_list.append({"mol_id": f"mol_{idx}"})

        # Build result dataframe
        result_df = pd.DataFrame(properties_list)

        return molecules, result_df


def load_smiles(
    filepath: Union[str, Path],
    smiles_column: int = 0,
    delimiter: Optional[str] = None,
    id_column: Optional[int] = None,
    keep_properties: bool = False,
    sanitize: bool = True,
    add_hydrogens: bool = False,
    compute_2d_coords: bool = False,
    raise_on_invalid: bool = False,
) -> tuple[list["Chem.Mol"], pd.DataFrame]:
    """
    Load molecules from a SMILES file.

    Convenience function wrapping SmilesLoader.

    Parameters
    ----------
    filepath : str or Path
        Path to SMILES file
    smiles_column : int, default=0
        Column index containing SMILES strings
    delimiter : str, optional
        Column delimiter. If None, auto-detect
    id_column : int, optional
        Column index for compound IDs
    keep_properties : bool, default=False
        Keep all columns as metadata
    sanitize : bool, default=True
        Whether to sanitize molecules
    add_hydrogens : bool, default=False
        Whether to add explicit hydrogens
    compute_2d_coords : bool, default=False
        Whether to compute 2D coordinates
    raise_on_invalid : bool, default=False
        Whether to raise on invalid molecules

    Returns
    -------
    mols : list of Chem.Mol
        RDKit molecule objects
    data : pd.DataFrame
        DataFrame with SMILES and metadata

    Examples
    --------
    Load a simple SMILES file::

        from stratosampler.IO import load_smiles

        mols, data = load_smiles("compounds.smi")
        print(f"Loaded {len(mols)} molecules")

    Load a CSV with SMILES and compound IDs::

        mols, data = load_smiles(
            "compounds.csv",
            delimiter=",",
            smiles_column=1,
            id_column=0,
        )
    """
    loader = SmilesLoader(raise_on_invalid=raise_on_invalid)
    return loader.load(
        filepath=filepath,
        smiles_column=smiles_column,
        delimiter=delimiter,
        id_column=id_column,
        keep_properties=keep_properties,
        sanitize=sanitize,
        add_hydrogens=add_hydrogens,
        compute_2d_coords=compute_2d_coords,
    )


def load_sdf(
    filepath: Union[str, Path],
    include_properties: bool = True,
    sanitize: bool = True,
    remove_hs: bool = True,
    compute_2d_coords: bool = False,
    raise_on_invalid: bool = False,
) -> tuple[list["Chem.Mol"], pd.DataFrame]:
    """
    Load molecules from an SDF (Structure Data Format) file.

    Convenience function wrapping SdfLoader.

    Parameters
    ----------
    filepath : str or Path
        Path to SDF file
    include_properties : bool, default=True
        Extract SDF properties into DataFrame
    sanitize : bool, default=True
        Whether to sanitize molecules
    remove_hs : bool, default=True
        Whether to remove explicit hydrogens
    compute_2d_coords : bool, default=False
        Whether to compute 2D coordinates
    raise_on_invalid : bool, default=False
        Whether to raise on invalid molecules

    Returns
    -------
    mols : list of Chem.Mol
        RDKit molecule objects
    data : pd.DataFrame
        DataFrame with molecule properties

    Examples
    --------
    Load molecules from SDF::

        from stratosampler.IO import load_sdf

        mols, data = load_sdf("compounds.sdf")
        print(f"Loaded {len(mols)} molecules")

    Load with explicit hydrogen atoms::

        mols, data = load_sdf("compounds_with_h.sdf", remove_hs=False)
    """
    loader = SdfLoader(raise_on_invalid=raise_on_invalid)
    return loader.load(
        filepath=filepath,
        include_properties=include_properties,
        sanitize=sanitize,
        remove_hs=remove_hs,
        compute_2d_coords=compute_2d_coords,
    )


__all__ = [
    "SmilesLoader",
    "SdfLoader",
    "MoleculeLoader",
    "load_smiles",
    "load_sdf",
]
