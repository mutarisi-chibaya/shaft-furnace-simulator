"""Provenance checks: every quoted parameter in data/kinetics.yaml, data/furnace.yaml and
data/uncertainty.yaml must appear in its source.

For each entry with a `quote`, the test checks that
1. the quote occurs (whitespace-normalised) in references/<ref_file>, the text layer of the
   source PDF produced by scripts/fetch_references.py (pdftotext);
2. `as_written` occurs in the quote;
3. value == number(as_written) * factor + offset (offset defaults to 0, e.g. 273.15 for degC).
Skipped when the references have not been downloaded.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
import yaml

from .conftest import DATA_DIR, ROOT

REF_DIR = ROOT / "references"


FILES = ("kinetics.yaml", "furnace.yaml", "uncertainty.yaml", "emissions.yaml", "sensitivity.yaml")


def _entries() -> list[tuple[str, dict[str, Any]]]:
    out = []
    for f in FILES:
        raw = yaml.safe_load((DATA_DIR / f).read_text(encoding="utf-8"))
        for section, items in raw.items():
            for key, e in items.items():
                if isinstance(e, dict) and "quote" in e:
                    out.append((f"{f}:{section}.{key}", e))
    return out


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text)


_cache: dict[str, str] = {}


def _source_text(ref_file: str) -> str:
    path = REF_DIR / ref_file
    if not path.exists():
        pytest.skip(f"{path.name} not found; run scripts/fetch_references.py")
    if ref_file not in _cache:
        _cache[ref_file] = _norm(Path(path).read_bytes().decode("latin-1"))
    return _cache[ref_file]


def parse_number(as_written: str) -> float:
    """'30, 584' -> 30584 (', ' thousands separator); '0,53' -> 0.53 (decimal comma)."""
    return float(as_written.replace(", ", "").replace(",", "."))


def test_parse_number() -> None:
    assert parse_number("30, 584") == 30584.0
    assert parse_number("0,53") == 0.53
    assert parse_number("0.01282") == 0.01282


@pytest.mark.parametrize(("name", "entry"), _entries(), ids=[n for n, _ in _entries()])
def test_quoted_value_matches_source(name: str, entry: dict[str, Any]) -> None:
    text = _source_text(entry["ref_file"])
    quote = _norm(entry["quote"])
    assert quote in text, f"{name}: quote not found in {entry['ref_file']}"
    assert entry["as_written"] in quote, f"{name}: as_written not in quote"
    expected = parse_number(entry["as_written"]) * float(entry["factor"]) + float(
        entry.get("offset", 0.0)
    )
    assert float(entry["value"]) == pytest.approx(expected, rel=1e-12), name


@pytest.mark.parametrize("fname", FILES)
def test_every_value_has_a_source(fname: str) -> None:
    raw = yaml.safe_load((DATA_DIR / fname).read_text(encoding="utf-8"))
    for section, items in raw.items():
        for key, e in items.items():
            assert e.get("source"), f"{fname}: {section}.{key} has no source"


def test_wakao_kaguei_quote() -> None:
    """data/transport.yaml: Nu = 2 + 1.1 Pr^(1/3) Re^0.6 as printed by Wagner (2008)."""
    text = _source_text("wagner2008_thesis.txt")
    assert "Nu = 2 + 1,1Pr1 3 Re0,6" in text
    raw = yaml.safe_load((DATA_DIR / "transport.yaml").read_text(encoding="utf-8"))
    wk = raw["wakao_kaguei"]
    assert (wk["a"], wk["b"], wk["m"]) == (2.0, 1.1, 0.6)
    assert wk["n"] == pytest.approx(1 / 3)
