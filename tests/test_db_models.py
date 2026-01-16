import pytest
from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from backend.app.db.base import Base
from backend.app.db.models import Molecule


@pytest.fixture()
def db_session():
    # Для unit-тестов используем SQLite
    # для проверки моделей/constraints.
    engine = create_engine("sqlite+pysqlite:///:memory:", future=True)
    Base.metadata.create_all(engine)
    SessionLocal = sessionmaker(bind=engine, expire_on_commit=False, future=True)

    with SessionLocal() as session:
        yield session


def test_create_and_read_molecule(db_session: Session):
    m = Molecule(smiles="CCO", name="ethanol")
    db_session.add(m)
    db_session.commit()

    got = db_session.query(Molecule).filter_by(smiles="CCO").one()
    assert got.name == "ethanol"
    assert got.id is not None


def test_name_is_optional(db_session: Session):
    m = Molecule(smiles="c1ccccc1")
    db_session.add(m)
    db_session.commit()

    got = db_session.query(Molecule).filter_by(smiles="c1ccccc1").one()
    assert got.name is None


def test_smiles_unique_constraint(db_session: Session):
    db_session.add(Molecule(smiles="CCO"))
    db_session.commit()

    db_session.add(Molecule(smiles="CCO"))
    with pytest.raises(IntegrityError):
        db_session.commit()
