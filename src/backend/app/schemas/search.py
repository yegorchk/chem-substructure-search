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


class SimilaritySearchRequest(BaseModel):
    query_smiles: str = Field(min_length=1, max_length=2048)
    top_k: int = Field(default=10, ge=1, le=100)
    threshold: float = Field(default=0.0, ge=0.0, le=1.0)


class SimilaritySearchItem(BaseModel):
    id: int
    smiles: str
    name: str | None
    score: float


class SimilaritySearchResponse(BaseModel):
    items: list[SimilaritySearchItem]
    query_smiles: str
    top_k: int
    threshold: float
