from __future__ import annotations

import logging

from rdkit import Chem
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from backend.app.db.models import Molecule

logger = logging.getLogger(__name__)


def _validate_smiles(smiles: str) -> None:
    # фундамент: SMILES
    if Chem.MolFromSmiles(smiles) is None:
        raise ValueError("Invalid SMILES")


def create_molecule(db: Session, *, smiles: str, name: str | None) -> Molecule:
    _validate_smiles(smiles)

    m = Molecule(smiles=smiles, name=name)
    db.add(m)
    db.commit()
    db.refresh(m)

    logger.info("Created molecule id=%s", m.id)
    return m


def get_molecule(db: Session, molecule_id: int) -> Molecule | None:
    return db.get(Molecule, molecule_id)


def update_molecule(
    db: Session, *, molecule_id: int, smiles: str | None, name: str | None
) -> Molecule | None:
    m = db.get(Molecule, molecule_id)
    if m is None:
        return None

    if smiles is not None:
        _validate_smiles(smiles)
        m.smiles = smiles

    # name может быть None — это валидное значение "стереть имя"
    if name is not None or name is None:
        m.name = name

    db.add(m)
    db.commit()
    db.refresh(m)

    logger.info("Updated molecule id=%s", m.id)
    return m


def delete_molecule(db: Session, molecule_id: int) -> bool:
    m = db.get(Molecule, molecule_id)
    if m is None:
        return False

    db.delete(m)
    db.commit()
    logger.info("Deleted molecule id=%s", molecule_id)
    return True


def list_molecules(
    db: Session, *, limit: int, offset: int
) -> tuple[list[Molecule], int]:
    limit = max(1, min(limit, 200))
    offset = max(0, offset)

    total = db.scalar(select(func.count()).select_from(Molecule)) or 0
    items = db.scalars(
        select(Molecule).order_by(Molecule.id.asc()).limit(limit).offset(offset)
    ).all()
    return items, total
