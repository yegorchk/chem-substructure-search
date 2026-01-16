from __future__ import annotations

import logging

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.app.db.models import Molecule

logger = logging.getLogger(__name__)


def substructure_search_db(
    db: Session,
    *,
    substructure_smiles: str,
    limit: int,
    offset: int,
) -> tuple[list[Molecule], int]:
    sub = Chem.MolFromSmiles(substructure_smiles)
    if sub is None:
        raise ValueError("Invalid substructure SMILES")

    molecules = db.scalars(select(Molecule).order_by(Molecule.id.asc())).all()

    matches: list[Molecule] = []
    for m in molecules:
        mol = Chem.MolFromSmiles(m.smiles)
        if mol is None:
            # по идее в БД такого быть не должно (мы валидируем на входе), но safety ok
            continue
        if mol.HasSubstructMatch(sub):
            matches.append(m)

    total = len(matches)
    page = matches[offset : offset + limit]

    logger.info(
        "Substructure search sub=%s total=%s returned=%s",
        substructure_smiles,
        total,
        len(page),
    )
    return page, total


def similarity_search_db(
    db: Session,
    *,
    query_smiles: str,
    top_k: int,
    threshold: float,
) -> list[tuple[Molecule, float]]:
    qmol = Chem.MolFromSmiles(query_smiles)
    if qmol is None:
        raise ValueError("Invalid query SMILES")

    qfp = AllChem.GetMorganFingerprintAsBitVect(qmol, radius=2, nBits=2048)

    molecules = db.scalars(select(Molecule).order_by(Molecule.id.asc())).all()

    scored: list[tuple[Molecule, float]] = []
    for m in molecules:
        mol = Chem.MolFromSmiles(m.smiles)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        score = float(DataStructs.TanimotoSimilarity(qfp, fp))
        if score >= threshold:
            scored.append((m, score))

    scored.sort(key=lambda x: x[1], reverse=True)
    return scored[:top_k]
