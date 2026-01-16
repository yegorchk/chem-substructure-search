import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.app.api.deps import get_db
from backend.app.db.base import Base
from backend.app.db import models  # noqa: F401
from backend.app.main import create_app


@pytest.fixture()
def client():
    engine = create_engine(
        "sqlite+pysqlite:///:memory:",
        future=True,
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(engine)
    TestingSessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)

    app = create_app()

    def override_get_db():
        db = TestingSessionLocal()
        try:
            yield db
        finally:
            db.close()

    app.dependency_overrides[get_db] = override_get_db
    return TestClient(app)


def test_search_benzene_ring(client: TestClient):
    # create molecules
    assert client.post("/molecules", json={"smiles": "c1ccccc1"}).status_code == 201
    assert (
        client.post("/molecules", json={"smiles": "CC(=O)Oc1ccccc1C(=O)O"}).status_code
        == 201
    )
    assert client.post("/molecules", json={"smiles": "CCO"}).status_code == 201

    r = client.post("/search", json={"substructure_smiles": "c1ccccc1"})
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["total"] == 2
    smiles = [x["smiles"] for x in data["items"]]
    assert "c1ccccc1" in smiles
    assert "CC(=O)Oc1ccccc1C(=O)O" in smiles


def test_search_invalid_smiles_422(client: TestClient):
    r = client.post("/search", json={"substructure_smiles": "not-a-smiles"})
    assert r.status_code == 422


def test_search_pagination(client: TestClient):
    client.post("/molecules", json={"smiles": "c1ccccc1"})
    client.post("/molecules", json={"smiles": "Cc1ccccc1"})
    client.post("/molecules", json={"smiles": "CCO"})

    r1 = client.post(
        "/search", json={"substructure_smiles": "c1ccccc1", "limit": 1, "offset": 0}
    )
    assert r1.status_code == 200
    d1 = r1.json()
    assert d1["total"] == 2
    assert len(d1["items"]) == 1

    r2 = client.post(
        "/search", json={"substructure_smiles": "c1ccccc1", "limit": 1, "offset": 1}
    )
    d2 = r2.json()
    assert d2["total"] == 2
    assert len(d2["items"]) == 1


def test_similarity_search_returns_best_match(client: TestClient):
    # ethanol and ethylamine are more similar to each other than benzene
    assert (
        client.post("/molecules", json={"smiles": "CCO", "name": "ethanol"}).status_code
        == 201
    )
    assert (
        client.post(
            "/molecules", json={"smiles": "CCN", "name": "ethylamine"}
        ).status_code
        == 201
    )
    assert (
        client.post(
            "/molecules", json={"smiles": "c1ccccc1", "name": "benzene"}
        ).status_code
        == 201
    )

    r = client.post(
        "/similarity-search", json={"query_smiles": "CCO", "top_k": 2, "threshold": 0.0}
    )
    assert r.status_code == 200, r.text
    data = r.json()
    assert len(data["items"]) == 2

    # best match should be the same molecule (Tanimoto=1.0) if present
    assert data["items"][0]["smiles"] == "CCO"
    assert data["items"][0]["score"] == 1.0
