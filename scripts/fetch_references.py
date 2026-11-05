"""Download the open-access references used for the model parameters into references/.

The PDFs are not committed (their licences differ); this script restores them so that every
number in data/*.yaml can be checked against its source. If ``pdftotext`` (poppler) is on the
PATH, a ``.txt`` copy of each PDF is written next to it; tests/test_provenance.py uses these
to check that values typed into the YAML files appear verbatim in the source text.

Usage:
    python scripts/fetch_references.py
"""

from __future__ import annotations

import shutil
import subprocess
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "references"
USER_AGENT = "shaftsim-reference-fetch/0.1 (research script)"

# filename: (URL, note). All are open-access copies (HAL, NASA NTRS, NIST, GitHub).
REFERENCES: dict[str, tuple[str, str]] = {
    "marsigny2025_metals15_1034.pdf": (
        "https://hal.science/hal-05287190/document",
        "Marsigny, Letz, Mirgaux, Patisson (2025) Metals 15:1034, CC BY 4.0",
    ),
    "marsigny2025_metals15_862.pdf": (
        "https://hal.science/hal-05287188/document",
        "Marsigny et al. (2025) Metals 15:862, CC BY 4.0 (shaft furnace, Milestone 5)",
    ),
    "patisson2020_metals10_922.pdf": (
        "https://hal.univ-lorraine.fr/hal-02968422/document",
        "Patisson, Mirgaux (2020) Metals 10:922, CC BY 4.0",
    ),
    "thesis2025_h2rich.pdf": (
        "https://theses.hal.science/tel-05795515/document",
        "Marsigny (2025) PhD thesis, Univ. Lorraine, CC BY 4.0",
    ),
    "wagner2008_thesis.pdf": (
        "https://theses.hal.science/tel-00280689/document",
        "Wagner (2008) PhD thesis, INPL Nancy (HAL tel-00280689)",
    ),
    "ranzani2011_thesis.pdf": (
        "https://theses.hal.science/tel-01204934/document",
        "Ranzani da Costa (2011) PhD thesis, INPL Nancy (HAL tel-01204934)",
    ),
    "ranzani2013_jclepro.pdf": (
        "https://hal.science/hal-00943356/document",
        "Ranzani da Costa, Wagner, Patisson (2013) J. Clean. Prod. 46:27 (HAL author copy)",
    ),
    "svehla1962_nasa_tr_r132.pdf": (
        "https://ntrs.nasa.gov/api/citations/19630012982/downloads/19630012982.pdf",
        "Svehla (1962) NASA TR R-132, public domain",
    ),
    "marrero1972_jpcrd1_3.pdf": (
        "https://srd.nist.gov/JPCRD/jpcrd1.pdf",
        "Marrero, Mason (1972) J. Phys. Chem. Ref. Data 1:3 (NIST-hosted scan)",
    ),
    "mmcollisionint_cantera_v3.0.0.cpp": (
        "https://raw.githubusercontent.com/Cantera/cantera/v3.0.0/src/transport/MMCollisionInt.cpp",
        "Cantera v3.0.0 collision-integral tables (Monchick & Mason 1961), BSD-3",
    ),
}


def main() -> None:
    DEST.mkdir(exist_ok=True)
    pdftotext = shutil.which("pdftotext")
    for name, (url, note) in REFERENCES.items():
        path = DEST / name
        if not path.exists():
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as resp:
                path.write_bytes(resp.read())
            print(f"downloaded {name}  ({note})")
        if pdftotext and name.endswith(".pdf"):
            txt = path.with_suffix(".txt")
            if not txt.exists():
                subprocess.run([pdftotext, "-layout", str(path), str(txt)], check=False)
    if not pdftotext:
        print("pdftotext not found: provenance tests on PDF text will be skipped")


if __name__ == "__main__":
    main()
