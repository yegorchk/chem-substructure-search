from fastapi import APIRouter

from backend.app.api.routers.health import router as health_router
from backend.app.api.routers.molecules import router as molecules_router
from backend.app.api.routers.search import router as search_router

api_router = APIRouter()
api_router.include_router(health_router, tags=["health"])
api_router.include_router(molecules_router, tags=["molecules"])
api_router.include_router(search_router, tags=["search"])
