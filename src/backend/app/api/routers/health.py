from fastapi import APIRouter, status

from backend.app.db.session import db_ping

router = APIRouter()


@router.get("/health")
def health():
    return {"status": "ok"}


@router.get("/health/db", status_code=status.HTTP_200_OK)
def health_db():
    if db_ping():
        return {"status": "ok", "db": "up"}
    return {"status": "degraded", "db": "down"}
