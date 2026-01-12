# chem-substructure-search (SMILES-first)

Stage 0: pure python module with RDKit-based substructure search using SMILES queries.

## Quickstart

### 1) Create venv
python3.12 -m venv .venv
source .venv/bin/activate

### 2) Install
pip install -U pip
pip install -e ".[dev]"

### 3) Test
pytest

### 4) Lint
ruff check .