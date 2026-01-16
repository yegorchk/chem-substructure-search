import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.app.api.deps import get_db
from backend.app.db.base import Base
from backend.app.main import create_app


@pytest.fixture()
def client():
    # SQLite для тестов
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


def test_create_and_get(client: TestClient):
    r = client.post("/molecules", json={"smiles": "CCO", "name": "ethanol"})
    assert r.status_code == 201, r.text
    created = r.json()
    assert "id" in created
    assert created["smiles"] == "CCO"

    mid = created["id"]
    r2 = client.get(f"/molecules/{mid}")
    assert r2.status_code == 200
    assert r2.json()["name"] == "ethanol"


def test_create_invalid_smiles(client: TestClient):
    r = client.post("/molecules", json={"smiles": "not-a-smiles"})
    assert r.status_code == 422


def test_list_pagination(client: TestClient):
    client.post("/molecules", json={"smiles": "CCO"})
    client.post("/molecules", json={"smiles": "c1ccccc1"})
    client.post("/molecules", json={"smiles": "CC(=O)O"})

    r = client.get("/molecules?limit=2&offset=0")
    assert r.status_code == 200
    payload = r.json()
    assert payload["total"] == 3
    assert payload["limit"] == 2
    assert payload["offset"] == 0
    assert len(payload["items"]) == 2

    r2 = client.get("/molecules?limit=2&offset=2")
    payload2 = r2.json()
    assert payload2["total"] == 3
    assert len(payload2["items"]) == 1


def test_unique_smiles_conflict(client: TestClient):
    r1 = client.post("/molecules", json={"smiles": "CCO"})
    assert r1.status_code == 201

    r2 = client.post("/molecules", json={"smiles": "CCO"})
    assert r2.status_code == 409
