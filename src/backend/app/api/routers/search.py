from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from backend.app.api.deps import get_db
from backend.app.schemas.search import (
    SubstructureSearchRequest,
    SubstructureSearchResponse,
    SubstructureSearchItem,
)
from backend.app.services.search import substructure_search_db

router = APIRouter()


@router.post(
    "/search", response_model=SubstructureSearchResponse, status_code=status.HTTP_200_OK
)
def substructure_search(
    payload: SubstructureSearchRequest, db: Session = Depends(get_db)
):
    try:
        items, total = substructure_search_db(
            db,
            substructure_smiles=payload.substructure_smiles,
            limit=payload.limit,
            offset=payload.offset,
        )
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))

    return SubstructureSearchResponse(
        items=[
            SubstructureSearchItem(id=m.id, smiles=m.smiles, name=m.name) for m in items
        ],
        total=total,
        limit=payload.limit,
        offset=payload.offset,
    )
