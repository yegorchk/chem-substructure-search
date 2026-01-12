from __future__ import annotations

import logging
from typing import Iterable, List, Literal

from chem_search.main import substructure_search as _substructure_search

logger = logging.getLogger(__name__)


def configure_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
    )


def substructure_search(
    molecules: Iterable[str],
    substructure: str,
    *,
    invalid_smiles: Literal["skip", "raise"] = "skip",
) -> List[str]:
    """
    Stage 1 API per method sheet: /src/main.py, 2 args (SMILES list + SMILES substructure).  [oai_citation:3‡yadro_project.pdf](sediment://file_00000000532c71f485a40b2f26b940d9)
    """
    return _substructure_search(molecules, substructure, invalid_smiles=invalid_smiles)


if __name__ == "__main__":
    configure_logging()
    demo = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    print(substructure_search(demo, "c1ccccc1"))
