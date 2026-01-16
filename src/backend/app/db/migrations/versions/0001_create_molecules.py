"""create molecules table

Revision ID: 0001_create_molecules
Revises:
Create Date: 2026-01-16
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision = "0001_create_molecules"
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "molecules",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("smiles", sa.String(length=2048), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=True),
    )
    op.create_unique_constraint("uq_molecules_smiles", "molecules", ["smiles"])


def downgrade() -> None:
    op.drop_constraint("uq_molecules_smiles", "molecules", type_="unique")
    op.drop_table("molecules")
