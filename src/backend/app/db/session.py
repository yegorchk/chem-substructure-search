from __future__ import annotations

import logging
from collections.abc import Iterator
from contextlib import contextmanager

from sqlalchemy import create_engine, text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session, sessionmaker

from backend.app.core.config import settings

logger = logging.getLogger(__name__)

engine: Engine = create_engine(settings.database_url, pool_pre_ping=True)
SessionLocal = sessionmaker(
    bind=engine, autoflush=False, autocommit=False, expire_on_commit=False
)


@contextmanager
def session_scope() -> Iterator[Session]:
    """
    Контекстный менеджер для безопасной работы с транзакцией:
    - commit если всё ок
    - rollback если ошибка
    """
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        logger.exception("DB transaction failed; rolled back")
        raise
    finally:
        session.close()


def db_ping() -> bool:
    """Простейшая проверка соединения."""
    try:
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
        return True
    except Exception:
        logger.exception("DB ping failed")
        return False
