"""
Tests for molecular loaders (SMILES and SDF).
"""

import tempfile
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem

from stratosampler.IO import load_smiles, load_sdf, SmilesLoader, SdfLoader


class TestSmilesLoader:
    """Test cases for SmilesLoader."""

    @pytest.fixture
    def simple_smiles_file(self):
        """Create a temporary SMILES file for testing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write("CCO ethanol\n")
            f.write("CC(=O)O acetic_acid\n")
            f.write("invalid_smiles bad_molecule\n")
            f.write("c1ccccc1 benzene\n")
            return Path(f.name)

    @pytest.fixture
    def csv_smiles_file(self):
        """Create a temporary CSV file with SMILES."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("mol_id,smiles,activity\n")
            f.write("mol_1,CCO,0.5\n")
            f.write("mol_2,CC(=O)O,0.8\n")
            f.write("mol_3,c1ccccc1,0.9\n")
            return Path(f.name)

    def test_load_simple_smiles(self, simple_smiles_file):
        """Test loading simple SMILES file with space delimiter."""
        loader = SmilesLoader(raise_on_invalid=False)
        mols, data = loader.load(simple_smiles_file)

        # Should skip invalid SMILES, so we get 3 molecules
        assert len(mols) == 3
        assert len(data) == 3

    def test_load_csv_with_columns(self, csv_smiles_file):
        """Test loading CSV with specific columns."""
        loader = SmilesLoader()
        mols, data = loader.load(
            csv_smiles_file,
            delimiter=",",
            smiles_column=1,
            id_column=0,
            keep_properties=True,
        )

        assert len(mols) == 3
        assert "smiles" in data.columns or "activity" in data.columns

    def test_invalid_smiles_handling(self):
        """Test handling of invalid SMILES strings."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi") as f:
            f.write("invalid_smiles\n")
            f.flush()

            # With raise_on_invalid=False (default), should skip and warn
            loader = SmilesLoader(raise_on_invalid=False)
            mols, data = loader.load(Path(f.name))
            assert len(mols) == 0

            # With raise_on_invalid=True, should raise
            loader = SmilesLoader(raise_on_invalid=True)
            with pytest.raises(ValueError):
                loader.load(Path(f.name))

    def test_nonexistent_file(self):
        """Test loading from nonexistent file."""
        loader = SmilesLoader()
        with pytest.raises(FileNotFoundError):
            loader.load(Path("nonexistent_file.smi"))

    def test_add_hydrogens(self, simple_smiles_file):
        """Test adding explicit hydrogens."""
        loader = SmilesLoader()
        mols, _ = loader.load(simple_smiles_file, add_hydrogens=True)

        # Molecules should have explicit hydrogens
        for mol in mols:
            assert mol.GetNumAtoms() > 0

    def test_compute_2d_coords(self, simple_smiles_file):
        """Test computing 2D coordinates."""
        loader = SmilesLoader()
        mols, _ = loader.load(simple_smiles_file, compute_2d_coords=True)

        # Molecules should have coordinates
        for mol in mols:
            conf = mol.GetConformer()
            assert conf is not None

    def test_delimiter_detection(self):
        """Test automatic delimiter detection."""
        loader = SmilesLoader()

        # Test space delimiter (default)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt") as f:
            f.write("CCO ethanol\n")
            f.flush()
            detected = loader._detect_delimiter(Path(f.name))
            assert detected == " "

        # Test comma delimiter
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv") as f:
            f.write("SMILES,ID\n")
            f.flush()
            detected = loader._detect_delimiter(Path(f.name))
            assert detected == ","


class TestSdfLoader:
    """Test cases for SdfLoader."""

    @pytest.fixture
    def simple_sdf_file(self):
        """Create a temporary SDF file for testing."""
        # Create molecules
        mols = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("CC(=O)O"),
            Chem.MolFromSmiles("c1ccccc1"),
        ]

        # Add names and properties
        for i, mol in enumerate(mols):
            mol.SetProp("_Name", f"compound_{i}")
            mol.SetProp("activity", str(0.5 + i * 0.2))

        # Write to SDF
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sdf", delete=False) as f:
            writer = Chem.SDWriter(f.name)
            for mol in mols:
                writer.write(mol)
            writer.close()
            return Path(f.name)

    def test_load_sdf_basic(self, simple_sdf_file):
        """Test basic SDF loading."""
        loader = SdfLoader()
        mols, data = loader.load(simple_sdf_file)

        assert len(mols) == 3
        assert len(data) == 3
        assert "mol_id" in data.columns

    def test_load_sdf_properties(self, simple_sdf_file):
        """Test loading SDF with properties."""
        loader = SdfLoader()
        mols, data = loader.load(simple_sdf_file, include_properties=True)

        assert "name" in data.columns or "activity" in data.columns

    def test_load_sdf_no_properties(self, simple_sdf_file):
        """Test loading SDF without extracting properties."""
        loader = SdfLoader()
        mols, data = loader.load(simple_sdf_file, include_properties=False)

        assert len(mols) == 3
        assert "mol_id" in data.columns

    def test_sdf_nonexistent_file(self):
        """Test loading from nonexistent SDF file."""
        loader = SdfLoader()
        with pytest.raises(FileNotFoundError):
            loader.load(Path("nonexistent_file.sdf"))

    def test_sdf_compute_2d_coords(self, simple_sdf_file):
        """Test computing 2D coordinates for SDF molecules."""
        loader = SdfLoader()
        mols, _ = loader.load(simple_sdf_file, compute_2d_coords=True)

        for mol in mols:
            conf = mol.GetConformer()
            assert conf is not None


class TestConvenienceFunctions:
    """Test convenience functions (load_smiles, load_sdf)."""

    def test_load_smiles_function(self):
        """Test load_smiles convenience function."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi") as f:
            f.write("CCO ethanol\n")
            f.write("CC(=O)O acetic_acid\n")
            f.flush()

            mols, data = load_smiles(Path(f.name))
            assert len(mols) == 2

    def test_load_sdf_function(self):
        """Test load_sdf convenience function."""
        # Create molecule
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", "ethanol")
        mol.SetProp("property", "test")

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
            writer = Chem.SDWriter(f.name)
            writer.write(mol)
            writer.close()

            mols, data = load_sdf(Path(f.name))
            assert len(mols) == 1


class TestIntegration:
    """Integration tests combining loaders with stratified splitter."""

    def test_load_and_split_workflow(self):
        """Test complete workflow: load SMILES -> compute properties -> split."""
        pytest.importorskip("stratosampler.splitters")

        from stratosampler import PropertyStratifiedSplitter

        # Create test SMILES
        smiles_list = [
            "CCO",  # ethanol
            "CC(=O)O",  # acetic acid
            "c1ccccc1",  # benzene
            "CCCc1ccccc1",  # propylbenzene
            "CC(C)O",  # isopropanol
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi") as f:
            for smi in smiles_list:
                f.write(f"{smi}\n")
            f.flush()

            # Load molecules
            mols, data = load_smiles(Path(f.name))
            assert len(mols) == 5

            # Prepare data for stratified splitting
            data["SMILES"] = [Chem.MolToSmiles(mol) for mol in mols]

            # Create splitter
            splitter = PropertyStratifiedSplitter(
                properties=["MolWt"],
                n_bins=2,
            )

            # Split
            train_idx, test_idx = splitter.split(data, smiles_col="SMILES")
            assert len(train_idx) + len(test_idx) == len(mols)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
