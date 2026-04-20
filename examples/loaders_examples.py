"""
Examples of loading molecules using stratosampler IO loaders.

This notebook demonstrates various ways to load molecular data
from different file formats using the stratosampler.IO module.
"""

import pandas as pd
from stratosampler.IO import load_smiles, load_sdf, SmilesLoader, SdfLoader


def example_1_simple_smiles_loading():
    """Example 1: Load molecules from a simple SMILES file (one per line)."""
    print("=" * 60)
    print("Example 1: Simple SMILES Loading")
    print("=" * 60)

    # Assume we have a file "compounds.smi" with one SMILES per line:
    # CCO
    # CC(=O)O
    # c1ccccc1

    # Load molecules
    # mols, data = load_smiles("compounds.smi")
    # print(f"Loaded {len(mols)} molecules")
    # print(data)

    print("Usage: mols, data = load_smiles('compounds.smi')")
    print()


def example_2_smiles_with_ids():
    """Example 2: Load SMILES file with compound IDs."""
    print("=" * 60)
    print("Example 2: SMILES with Compound IDs")
    print("=" * 60)

    # File format: "compound_id SMILES"
    # mol_1 CCO
    # mol_2 CC(=O)O
    # mol_3 c1ccccc1

    # Load with space delimiter and custom ID column
    # mols, data = load_smiles(
    #     "compounds_with_ids.txt",
    #     smiles_column=1,
    #     id_column=0,
    # )
    # print(data)

    print(
        """Usage:
    mols, data = load_smiles(
        "compounds_with_ids.txt",
        smiles_column=1,
        id_column=0,
    )"""
    )
    print()


def example_3_csv_smiles():
    """Example 3: Load SMILES from CSV file with multiple properties."""
    print("=" * 60)
    print("Example 3: CSV with SMILES and Properties")
    print("=" * 60)

    # CSV file format:
    # mol_id,SMILES,activity,logp
    # mol_1,CCO,0.5,0.3
    # mol_2,CC(=O)O,0.8,0.1
    # mol_3,c1ccccc1,0.9,1.5

    # Load keeping all properties
    # mols, data = load_smiles(
    #     "compounds.csv",
    #     delimiter=",",
    #     smiles_column=1,
    #     id_column=0,
    #     keep_properties=True,
    # )
    # print(f"Loaded {len(mols)} molecules with {len(data.columns)} properties")
    # print(data.columns)

    print(
        """Usage:
    mols, data = load_smiles(
        "compounds.csv",
        delimiter=",",
        smiles_column=1,
        id_column=0,
        keep_properties=True,
    )
    
Result: DataFrame with molecule properties"""
    )
    print()


def example_4_sdf_loading():
    """Example 4: Load molecules from SDF file."""
    print("=" * 60)
    print("Example 4: Loading from SDF Files")
    print("=" * 60)

    # Load molecules from SDF (includes 3D coordinates and properties)
    # mols, data = load_sdf("compounds.sdf")
    # print(f"Loaded {len(mols)} molecules")
    # print(data.head())

    print("Usage: mols, data = load_sdf('compounds.sdf')")
    print("Result: List of RDKit Mol objects + DataFrame with properties")
    print()


def example_5_sdf_with_hydrogens():
    """Example 5: Load SDF with different hydrogen handling."""
    print("=" * 60)
    print("Example 5: SDF with Explicit Hydrogens")
    print("=" * 60)

    # Load SDF keeping explicit hydrogens
    # mols, data = load_sdf(
    #     "compounds_with_h.sdf",
    #     remove_hs=False,
    # )

    print(
        """Usage:
    mols, data = load_sdf(
        "compounds_with_h.sdf",
        remove_hs=False,  # Keep H atoms
    )"""
    )
    print()


def example_6_smiles_with_2d_coords():
    """Example 6: Load SMILES and compute 2D coordinates."""
    print("=" * 60)
    print("Example 6: SMILES with 2D Coordinates")
    print("=" * 60)

    # Load SMILES and compute 2D coordinates for visualization
    # mols, data = load_smiles(
    #     "compounds.smi",
    #     compute_2d_coords=True,
    # )
    #
    # # Now molecules have coordinates for drawing
    # from rdkit.Chem import Draw
    # img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200))
    # img.show()

    print(
        """Usage:
    mols, data = load_smiles(
        "compounds.smi",
        compute_2d_coords=True,
    )
    
# Visualize molecules
from rdkit.Chem import Draw
img = Draw.MolsToGridImage(mols, molsPerRow=4)
img.show()"""
    )
    print()


def example_7_using_loaders_directly():
    """Example 7: Using loader classes directly for advanced control."""
    print("=" * 60)
    print("Example 7: Direct Loader Usage")
    print("=" * 60)

    # For more control, use loader classes directly
    print(
        """from stratosampler.IO import SmilesLoader, SdfLoader

# SMILES Loader with custom options
smiles_loader = SmilesLoader(raise_on_invalid=True)
mols, data = smiles_loader.load(
    "compounds.csv",
    delimiter=",",
    smiles_column=1,
    sanitize=True,
    add_hydrogens=False,
)

# SDF Loader
sdf_loader = SdfLoader(raise_on_invalid=False)
mols, data = sdf_loader.load(
    "compounds.sdf",
    include_properties=True,
    remove_hs=True,
)"""
    )
    print()


def example_8_error_handling():
    """Example 8: Error handling with invalid molecules."""
    print("=" * 60)
    print("Example 8: Handling Invalid Molecules")
    print("=" * 60)

    # By default, invalid molecules are skipped with warnings
    # mols, data = load_smiles("compounds_with_errors.smi")  # Warns about invalid

    # Raise on invalid molecules
    # mols, data = load_smiles(
    #     "compounds.smi",
    #     raise_on_invalid=True,  # Will raise ValueError on invalid SMILES
    # )

    print(
        """# Skip invalid molecules (default)
mols, data = load_smiles("compounds.smi")

# Raise exception on invalid molecules
mols, data = load_smiles(
    "compounds.smi",
    raise_on_invalid=True,
)"""
    )
    print()


def example_9_integration_with_stratified_splitter():
    """Example 9: Integration with stratified splitter."""
    print("=" * 60)
    print("Example 9: Integration with Stratified Splitter")
    print("=" * 60)

    # Complete workflow: Load SMILES -> Compute properties -> Split
    print(
        """from stratosampler.IO import load_smiles
from stratosampler import PropertyStratifiedSplitter
from rdkit.Chem import AllChem

# 1. Load molecules from SMILES file
mols, data = load_smiles("compounds.smi")

# 2. Add SMILES back to dataframe for splitting
from rdkit.Chem import MolToSmiles
data["SMILES"] = [MolToSmiles(mol) for mol in mols]

# 3. Create stratified splitter
splitter = PropertyStratifiedSplitter(
    properties=["MolLogP", "MolWt", "TPSA"],
    n_bins=5,
    test_size=0.2,
)

# 4. Perform stratified split
train_idx, test_idx = splitter.split(data, smiles_col="SMILES")

# 5. Access training and test sets
train_mols = [mols[i] for i in train_idx]
test_mols = [mols[i] for i in test_idx]
"""
    )
    print()


def example_10_batch_processing():
    """Example 10: Processing multiple files."""
    print("=" * 60)
    print("Example 10: Batch Processing Multiple Files")
    print("=" * 60)

    # Process multiple SMILES files
    print(
        """from pathlib import Path
from stratosampler.IO import load_smiles

# Find all SMILES files
smiles_files = Path("data/").glob("*.smi")

all_molecules = []
all_data = []

for filepath in smiles_files:
    mols, data = load_smiles(filepath)
    all_molecules.extend(mols)
    all_data.append(data)

# Combine results
import pandas as pd
combined_data = pd.concat(all_data, ignore_index=True)
print(f"Loaded {len(all_molecules)} molecules from {len(list(smiles_files))} files")
"""
    )
    print()


def example_11_file_format_reference():
    """Example 11: Reference for different file formats."""
    print("=" * 60)
    print("Example 11: File Format Reference")
    print("=" * 60)

    print(
        """
SMILES File Format
------------------
Simple SMILES (one per line):
    CCO
    CC(=O)O
    c1ccccc1

SMILES with IDs (space-separated):
    mol_1 CCO
    mol_2 CC(=O)O
    mol_3 c1ccccc1

CSV Format (comma-separated):
    mol_id,SMILES,activity,logp
    mol_1,CCO,0.5,0.3
    mol_2,CC(=O)O,0.8,0.1

TSV Format (tab-separated):
    mol_id	SMILES	activity	logp
    mol_1	CCO	0.5	0.3
    mol_2	CC(=O)O	0.8	0.1

SDF Format
----------
3D molecular structures in SDF format with embedded properties.
Can be created from ChemDraw, Maestro, PyMOL, etc.
Example: Load with load_sdf("compounds.sdf")
"""
    )
    print()


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("STRATOSAMPLER MOLECULAR LOADERS - EXAMPLES")
    print("=" * 60 + "\n")

    example_1_simple_smiles_loading()
    example_2_smiles_with_ids()
    example_3_csv_smiles()
    example_4_sdf_loading()
    example_5_sdf_with_hydrogens()
    example_6_smiles_with_2d_coords()
    example_7_using_loaders_directly()
    example_8_error_handling()
    example_9_integration_with_stratified_splitter()
    example_10_batch_processing()
    example_11_file_format_reference()

    print("=" * 60)
    print("For more information, see the documentation:")
    print("  - API Reference: docs/api.md")
    print("  - Full Examples: examples/")
    print("=" * 60)


if __name__ == "__main__":
    main()
