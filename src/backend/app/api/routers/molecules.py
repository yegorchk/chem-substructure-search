from __future__ import annotations

import logging

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from backend.app.api.deps import get_db
from backend.app.schemas.molecules import (
    MoleculeCreate,
    MoleculeList,
    MoleculeOut,
    MoleculeUpdate,
)
from backend.app.services.molecules import (
    create_molecule,
    delete_molecule,
    get_molecule,
    list_molecules,
    update_molecule,
)

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/molecules")


@router.post("", response_model=MoleculeOut, status_code=status.HTTP_201_CREATED)
def create(payload: MoleculeCreate, db: Session = Depends(get_db)):
    try:
        m = create_molecule(db, smiles=payload.smiles, name=payload.name)
        return m
    except ValueError as e:
        # SMILES invalid
        raise HTTPException(status_code=422, detail=str(e))
    except IntegrityError:
        # UNIQUE smiles violation
        logger.exception("IntegrityError on create molecule")
        raise HTTPException(
            status_code=409, detail="Molecule with this SMILES already exists"
        )


@router.get("/{molecule_id}", response_model=MoleculeOut)
def get_one(molecule_id: int, db: Session = Depends(get_db)):
    m = get_molecule(db, molecule_id)
    if m is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return m


@router.put("/{molecule_id}", response_model=MoleculeOut)
def update_one(
    molecule_id: int, payload: MoleculeUpdate, db: Session = Depends(get_db)
):
    try:
        m = update_molecule(
            db, molecule_id=molecule_id, smiles=payload.smiles, name=payload.name
        )
        if m is None:
            raise HTTPException(status_code=404, detail="Molecule not found")
        return m
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except IntegrityError:
        logger.exception("IntegrityError on update molecule")
        raise HTTPException(
            status_code=409, detail="Molecule with this SMILES already exists"
        )


@router.delete("/{molecule_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_one(molecule_id: int, db: Session = Depends(get_db)):
    ok = delete_molecule(db, molecule_id)
    if not ok:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return None


@router.get("", response_model=MoleculeList)
def list_(
    limit: int = Query(default=50, ge=1, le=200),
    offset: int = Query(default=0, ge=0),
    db: Session = Depends(get_db),
):
    items, total = list_molecules(db, limit=limit, offset=offset)
    return MoleculeList(items=items, total=total, limit=limit, offset=offset)
