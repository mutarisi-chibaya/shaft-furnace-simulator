"""Case A outputs of REDUCTOR (Marsigny, Mirgaux & Patisson 2025, Metals 15:862, Table 1).

Usage:
    python scripts/extract_reductor_table1.py

Table 1 is scrambled in the PDF text layer (columns interleaved), so the case A column was
READ FROM THE PAGE IMAGE (page 7 of the article, PDF page 8). To guard against reading errors
this script
1. checks that every printed value occurs as a token in the text layer of the PDF
   (references/marsigny2025_metals15_862.txt, made by scripts/fetch_references.py);
2. checks the internal consistency of the outputs: the mass fractions of the DRI add up to
   100 %, the metallisation follows from the DRI composition (Fe / total Fe, molar masses
   from data/species.yaml), the GUD from the top-gas composition, and the top-gas flow
   from the inlets (equimolar reduction);
and writes tests/data/reference_model_marsigny2025b_caseA.csv.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TXT = ROOT / "references" / "marsigny2025_metals15_862.txt"
OUT = ROOT / "tests" / "data" / "reference_model_marsigny2025b_caseA.csv"

# (key, unit, printed in Table 1 (case A), bold = model output, description)
ROWS: tuple[tuple[str, str, str, bool, str], ...] = (
    ("top_gas_flow", "mol/s", "2796", False, "top gas molar flow"),
    ("top_gas_T", "degC", "344", True, "top gas temperature"),
    ("top_gas_x_H2", "%", "64.1", True, "top gas H2 (molar)"),
    ("top_gas_x_H2O", "%", "30.9", True, "top gas H2O (molar)"),
    ("top_gas_x_N2", "%", "5", True, "top gas N2 (molar)"),
    ("dri_flow", "t/h", "117.19", False, "DRI mass flow"),
    ("dri_T", "degC", "261", True, "DRI temperature at the solid outlet"),
    ("dri_w_Fe2O3", "%", "0", True, "DRI Fe2O3 (mass)"),
    ("dri_w_Fe3O4", "%", "0.5", True, "DRI Fe3O4 (mass)"),
    ("dri_w_FeO", "%", "7.71", True, "DRI FeO (mass)"),
    ("dri_w_Fe", "%", "87.16", True, "DRI Fe (mass)"),
    ("dri_w_gangue", "%", "4.63", True, "DRI gangue (mass)"),
    ("metallization", "%", "93.3", True, "metallic Fe / total Fe"),
    ("gud", "-", "0.33", True, "gas utilisation degree x_H2O / (x_H2 + x_H2O) at the top"),
    ("energy_demand", "GJ/t", "12.20", True, "energy demand per tonne of DRI (gas loop)"),
    ("bustle_pressure", "Pa", "180,400", True, "pressure at the bustle"),
    ("cooling_pressure", "Pa", "181,600", True, "pressure at the cooling gas inlet"),
)


def number(printed: str) -> float:
    return float(printed.replace(",", ""))


def main() -> None:
    if not TXT.exists():
        raise SystemExit(f"{TXT.name} not found; run scripts/fetch_references.py")
    tokens = set(re.split(r"\s+", TXT.read_bytes().decode("latin-1")))
    # thousands are printed "180,400"; pdftotext keeps the token
    missing = [p for _, _, p, _, _ in ROWS if p not in tokens]
    if missing:
        raise SystemExit(f"values not found in the text layer: {missing}")
    v = {k: number(p) for k, _, p, _, _ in ROWS}

    sp = yaml.safe_load((ROOT / "data" / "species.yaml").read_text(encoding="utf-8"))["species"]
    m = {k: float(sp[k]["molar_mass"]) for k in ("Fe2O3", "Fe3O4", "FeO", "Fe")}
    fe_per = {"Fe2O3": 2.0, "Fe3O4": 3.0, "FeO": 1.0, "Fe": 1.0}
    total = sum(v[f"dri_w_{k}"] for k in ("Fe2O3", "Fe3O4", "FeO", "Fe", "gangue"))
    fe_tot = sum(v[f"dri_w_{k}"] / m[k] * fe_per[k] for k in m)
    metal = v["dri_w_Fe"] / m["Fe"] / fe_tot
    gud = v["top_gas_x_H2O"] / (v["top_gas_x_H2"] + v["top_gas_x_H2O"])
    checks = {
        "DRI mass fractions sum to 100 %": (total, 100.0, 0.05),
        "metallisation from DRI composition": (100 * metal, v["metallization"], 0.1),
        "GUD from top-gas composition": (gud, v["gud"], 0.005),
        "top-gas flow = 2300 + 500 mol/s": (v["top_gas_flow"], 2800.0, 5.0),
    }
    for name, (got, ref, tol) in checks.items():
        status = "ok" if abs(got - ref) <= tol else "FAILED"
        print(f"{name}: {got:.4g} vs {ref:.4g} (tol {tol}) {status}")
        if status != "ok":
            raise SystemExit("consistency check failed")

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write(
            "# REDUCTOR case A, outputs and reference values (Marsigny, Mirgaux & Patisson 2025,\n"
            "# Metals 15:862, Table 1, CC BY 4.0), read from the page image and checked against\n"
            "# the PDF text layer by scripts/extract_reductor_table1.py (do not edit).\n"
            "# output = 1: printed in bold (a model output); 0: plain (an input or a reported\n"
            "# flow). Units as printed.\n"
        )
        fh.write("key,unit,value,printed,output,description\n")
        for k, unit, printed, bold, desc in ROWS:
            fh.write(f'{k},{unit},{number(printed):g},"{printed}",{int(bold)},{desc}\n')
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
