"""Shared fixtures and paths."""

from __future__ import annotations

from pathlib import Path

import pytest

from shaftsim.properties import MissingDataError, Species, load_species

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
TEST_DATA_DIR = ROOT / "tests" / "data"


@pytest.fixture(scope="session")
def real_species() -> dict[str, Species]:
    """Species from data/species.yaml; skips the test until the human has filled it in."""
    try:
        return load_species(DATA_DIR / "species.yaml")
    except MissingDataError as exc:
        pytest.skip(f"data/species.yaml not filled in yet ({str(exc)[:120]}...)")
