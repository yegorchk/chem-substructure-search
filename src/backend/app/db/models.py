from __future__ import annotations

from sqlalchemy import String, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column

from backend.app.db.base import Base


class Molecule(Base):
    __tablename__ = "molecules"
    __table_args__ = (UniqueConstraint("smiles", name="uq_molecules_smiles"),)

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    smiles: Mapped[str] = mapped_column(String(2048), nullable=False)
    name: Mapped[str | None] = mapped_column(String(255), nullable=True)
