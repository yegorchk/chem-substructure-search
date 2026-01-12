import logging

import pytest

from chem_search.main import substructure_search


def test_substructure_search_basic_smiles_query():
    molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    sub = "c1ccccc1"  # benzene ring in SMILES
    assert substructure_search(molecules, sub) == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


def test_invalid_substructure_smiles_raises():
    molecules = ["CCO", "c1ccccc1"]
    with pytest.raises(ValueError, match="Invalid substructure SMILES"):
        substructure_search(molecules, "C1(CC")  # broken SMILES


def test_invalid_molecule_smiles_skipped_and_logged(caplog):
    molecules = ["CCO", "not-a-smiles", "CC(=O)N"]
    sub = "CCO"
    caplog.set_level(logging.WARNING)

    assert substructure_search(molecules, sub) == ["CCO"]
    assert any(
        "Invalid molecule SMILES skipped" in rec.message for rec in caplog.records
    )


def test_invalid_molecule_smiles_raise_mode():
    molecules = ["CCO", "not-a-smiles"]
    sub = "CCO"
    with pytest.raises(ValueError, match="Invalid molecule SMILES skipped"):
        substructure_search(molecules, sub, invalid_smiles="raise")
