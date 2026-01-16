from __future__ import annotations

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="APP_",
        extra="ignore",
        env_file=".env",
        env_file_encoding="utf-8",
    )

    name: str = "chem-substructure-search"
    log_level: str = "INFO"

    # ВАЖНО: SQLAlchemy URL
    # пример: postgresql+psycopg://postgres:postgres@localhost:5432/chemdb
    database_url: str = "postgresql+psycopg://postgres:postgres@localhost:5432/chemdb"


settings = Settings()
