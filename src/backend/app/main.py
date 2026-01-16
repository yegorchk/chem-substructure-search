from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from backend.app.api.router import api_router
from backend.app.core.config import settings
from backend.app.core.logging import setup_logging


def create_app() -> FastAPI:
    setup_logging(settings.log_level)

    app = FastAPI(title=settings.name)
    app.include_router(api_router)

    static_dir = Path(__file__).resolve().parent / "static"

    if static_dir.exists():
        app.mount("/static", StaticFiles(directory=static_dir), name="static")

        @app.get("/", include_in_schema=False)
        def ui() -> FileResponse:
            return FileResponse(static_dir / "index.html")

    return app


app = create_app()
