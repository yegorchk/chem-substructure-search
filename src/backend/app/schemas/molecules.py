from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field


class MoleculeCreate(BaseModel):
    smiles: str = Field(min_length=1, max_length=2048)
    name: str | None = Field(default=None, max_length=255)


class MoleculeUpdate(BaseModel):
    smiles: str | None = Field(default=None, min_length=1, max_length=2048)
    name: str | None = Field(default=None, max_length=255)


class MoleculeOut(BaseModel):
    id: int
    smiles: str
    name: str | None

    model_config = ConfigDict(from_attributes=True)


class MoleculeList(BaseModel):
    items: list[MoleculeOut]
    total: int
    limit: int
    offset: int
