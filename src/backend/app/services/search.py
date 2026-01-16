from __future__ import annotations

import logging

from rdkit import Chem
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
