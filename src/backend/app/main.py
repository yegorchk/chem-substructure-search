from __future__ import annotations

from fastapi import FastAPI

from backend.app.api.router import api_router
from backend.app.core.config import settings
from backend.app.core.logging import setup_logging


def create_app() -> FastAPI:
    setup_logging(settings.log_level)

    app = FastAPI(title=settings.name)
    app.include_router(api_router)
    return app


app = create_app()
