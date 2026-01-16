from fastapi.testclient import TestClient

from backend.app.main import create_app


def test_health_ok():
    app = create_app()
    client = TestClient(app)

    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


def test_openapi_available():
    app = create_app()
    client = TestClient(app)

    resp = client.get("/openapi.json")
    assert resp.status_code == 200
    payload = resp.json()
    assert "paths" in payload
    assert "/health" in payload["paths"]
