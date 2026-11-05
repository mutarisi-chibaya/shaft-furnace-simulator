"""Energy demand of case A from the gas-loop model (Marsigny 2025 PhD thesis, Table 20, p. 125).

Usage:
    python scripts/extract_thesis_table20.py

Table 20 is intact in the text layer (references/thesis2025_h2rich.txt, made by
scripts/fetch_references.py). The script checks that each row occurs verbatim
(whitespace-normalised), that the parts add up to the total, and writes the case A column to
tests/data/reference_energy_marsigny2025_thesis_t20.csv.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TXT = ROOT / "references" / "thesis2025_h2rich.txt"
OUT = ROOT / "tests" / "data" / "reference_energy_marsigny2025_thesis_t20.csv"

# (key, row as printed (case A, case 2, case 1), description)
ROWS = (
    ("metallization_pct", "Metallization (%) 93.3 94.4 93.2", "metallisation (%)"),
    ("compression_GJ_per_t", "Compression energy (GJ/tDRI) 0.25 0.07 0.08", "gas compression"),
    (
        "heating_GJ_per_t",
        "Heating and reforming energy (GJ/tDRI) 2.28 2.81 2.78",
        "gas heating (no reforming in case A)",
    ),
    ("electrolysis_GJ_per_t", "Electrolysis energy (GJ/tDRI) 9.73 0.81 0.73", "water electrolysis"),
    ("total_GJ_per_t", "Global energy demand (GJ/tDRI) 12.26 3.69 3.59", "total"),
)


def main() -> None:
    if not TXT.exists():
        raise SystemExit(f"{TXT.name} not found; run scripts/fetch_references.py")
    text = re.sub(r"\s+", " ", TXT.read_bytes().decode("latin-1"))
    values = {}
    for key, row, _ in ROWS:
        if row not in text:
            raise SystemExit(f"row not found in the text layer: {row!r}")
        values[key] = float(row.split()[-3])  # case A column
    parts = (
        values["compression_GJ_per_t"]
        + values["heating_GJ_per_t"]
        + values["electrolysis_GJ_per_t"]
    )
    print(f"parts {parts:.2f} vs total {values['total_GJ_per_t']:.2f} GJ/t")
    if abs(parts - values["total_GJ_per_t"]) > 0.011:
        raise SystemExit("parts do not add up to the total")
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write(
            "# Energy demand of case A (100 % H2) from the gas-loop post-processing of REDUCTOR\n"
            "# (Marsigny 2025, PhD thesis, Table 20, p. 125, CC BY 4.0), extracted and checked by\n"
            "# scripts/extract_thesis_table20.py (do not edit). GJ per tonne of DRI.\n"
        )
        fh.write("key,value,description\n")
        for key, _, desc in ROWS:
            fh.write(f"{key},{values[key]:g},{desc}\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
