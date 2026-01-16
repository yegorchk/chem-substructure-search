from __future__ import annotations

from pydantic import BaseModel, Field


class SubstructureSearchRequest(BaseModel):
    substructure_smiles: str = Field(min_length=1, max_length=2048)
    limit: int = Field(default=50, ge=1, le=200)
    offset: int = Field(default=0, ge=0)


class SubstructureSearchItem(BaseModel):
    id: int
    smiles: str
    name: str | None


class SubstructureSearchResponse(BaseModel):
    items: list[SubstructureSearchItem]
    total: int
    limit: int
    offset: int
