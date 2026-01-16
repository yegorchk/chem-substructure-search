from __future__ import annotations

import logging
from collections.abc import Iterable
from typing import Literal

from rdkit import Chem

logger = logging.getLogger(__name__)


def configure_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
    )


def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    return Chem.MolFromSmiles(smiles)


def substructure_search(
    molecules: Iterable[str],
    substructure: str,
    *,
    invalid_smiles: Literal["skip", "raise"] = "skip",
) -> list[str]:
    if molecules is None:
        raise TypeError("molecules must be an iterable of SMILES strings, got None")

    query = _mol_from_smiles(substructure)
    if query is None:
        raise ValueError(f"Invalid substructure SMILES: {substructure!r}")

    result: list[str] = []
    for s in molecules:
        mol = _mol_from_smiles(s)
        if mol is None:
            msg = f"Invalid molecule SMILES skipped: {s!r}"
            if invalid_smiles == "raise":
                raise ValueError(msg)
            logger.warning(msg)
            continue

        if mol.HasSubstructMatch(query):
            result.append(s)

    return result


if __name__ == "__main__":
    configure_logging()

    demo_molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    demo_substructure = "c1ccccc1"  # SMILES
    print(substructure_search(demo_molecules, demo_substructure))
