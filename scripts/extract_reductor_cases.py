"""Inputs and outputs of REDUCTOR cases A-E (Marsigny, Mirgaux & Patisson 2025, Metals 15:862,
Table 1), for the case-to-case comparison of Milestone 7.

Usage:
    python scripts/extract_reductor_cases.py

As for case A (scripts/extract_reductor_table1.py), Table 1 is scrambled in the PDF text
layer, so the values were READ FROM THE PAGE IMAGE (page 7 of the article, PDF page 8). The
script
1. checks that every printed value occurs as a token in the text layer
   (references/marsigny2025_metals15_862.txt, made by scripts/fetch_references.py);
2. checks the internal consistency of each column: DRI mass fractions add up to 100 %,
   metallisation from the DRI composition, GUD from the top-gas composition, top-gas flow
   equal to the sum of the inlets (equimolar reduction). Deviations of the source itself are
   reported and listed in KNOWN_DEVIATIONS rather than corrected;
3. checks that the case A column equals tests/data/reference_model_marsigny2025b_caseA.csv;
and writes tests/data/reference_model_marsigny2025b_cases.csv (long format).
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TXT = ROOT / "references" / "marsigny2025_metals15_862.txt"
OUT = ROOT / "tests" / "data" / "reference_model_marsigny2025b_cases.csv"
CASE_A = ROOT / "tests" / "data" / "reference_model_marsigny2025b_caseA.csv"
CASES = ("A", "B", "C", "D", "E")

# (key, unit, bold = model output, description, printed values for cases A..E)
ROWS: tuple[tuple[str, str, bool, str, tuple[str, ...]], ...] = (
    (
        "bustle_flow",
        "mol/s",
        False,
        "bustle gas molar flow",
        ("2300", "2100", "2300", "2300", "2300"),
    ),
    ("bustle_T", "degC", False, "bustle gas temperature", ("890", "890", "950", "815", "890")),
    ("bustle_x_H2", "%", False, "bustle gas H2 (molar)", ("93", "93", "66", "93", "93")),
    ("bustle_x_H2O", "%", False, "bustle gas H2O (molar)", ("2", "2", "2", "2", "2")),
    ("bustle_x_N2", "%", False, "bustle gas N2 (molar)", ("5", "5", "32", "5", "5")),
    (
        "bustle_pressure",
        "Pa",
        True,
        "pressure at the bustle",
        ("180,400", "179,600", "212,400", "178,900", "176,700"),
    ),
    ("cooling_flow", "mol/s", False, "cooling gas molar flow", ("500", "700", "500", "500", "500")),
    ("cooling_T", "degC", False, "cooling gas temperature", ("70", "70", "70", "70", "70")),
    ("cooling_x_H2", "%", False, "cooling gas H2 (molar)", ("93", "93", "66", "93", "93")),
    ("cooling_x_H2O", "%", False, "cooling gas H2O (molar)", ("2", "2", "2", "2", "2")),
    ("cooling_x_N2", "%", False, "cooling gas N2 (molar)", ("5", "5", "32", "5", "5")),
    (
        "cooling_pressure",
        "Pa",
        True,
        "pressure at the cooling gas inlet",
        ("181,600", "181,900", "215,500", "180,800", "178,200"),
    ),
    ("solid_flow", "kg/s", False, "solid burden inlet mass flow", ("45", "45", "45", "45", "43")),
    ("solid_T", "degC", False, "solid burden inlet temperature", ("20", "20", "20", "1200", "20")),
    ("solid_w_Fe2O3", "%", False, "burden Fe2O3 (mass)", ("96.65",) * 5),
    ("solid_w_gangue", "%", False, "burden gangue (mass)", ("3.35",) * 5),
    ("pellet_diameter", "mm", False, "pellet diameter", ("15",) * 5),
    ("bed_porosity", "-", False, "bed porosity", ("0.5",) * 5),
    (
        "top_gas_flow",
        "mol/s",
        False,
        "top gas molar flow",
        ("2796", "2800", "2793", "2802", "2766"),
    ),
    ("top_gas_T", "degC", True, "top gas temperature", ("344", "341", "394", "806", "353")),
    ("top_gas_x_H2", "%", True, "top gas H2 (molar)", ("64.1", "65.1", "36.1", "62.9", "64.6")),
    ("top_gas_x_H2O", "%", True, "top gas H2O (molar)", ("30.9", "29.9", "31.9", "32.1", "30.4")),
    ("top_gas_x_N2", "%", True, "top gas N2 (molar)", ("5", "5", "32", "5", "5")),
    ("top_pressure", "Pa", False, "top gas pressure", ("150,000",) * 5),
    ("dri_flow", "t/h", False, "DRI mass flow", ("117.19", "117.22", "117.18", "116.86", "109.88")),
    (
        "dri_T",
        "degC",
        True,
        "DRI temperature at the solid outlet",
        ("261", "166", "402", "330", "403"),
    ),
    ("dri_w_Fe2O3", "%", True, "DRI Fe2O3 (mass)", ("0", "0", "0", "0", "0")),
    ("dri_w_Fe3O4", "%", True, "DRI Fe3O4 (mass)", ("0.5", "2.82", "0", "0", "0")),
    ("dri_w_FeO", "%", True, "DRI FeO (mass)", ("7.71", "5.04", "8.26", "7.11", "0")),
    ("dri_w_Fe", "%", True, "DRI Fe (mass)", ("87.16", "87.51", "87.11", "88.24", "95.15")),
    ("dri_w_gangue", "%", True, "DRI gangue (mass)", ("4.63", "4.63", "4.63", "4.64", "4.85")),
    ("metallization", "%", True, "metallic Fe / total Fe", ("93.3", "93.7", "93.2", "94.1", "100")),
    (
        "gud",
        "-",
        True,
        "gas utilisation degree x_H2O / (x_H2 + x_H2O) at the top",
        ("0.33", "0.33", "0.47", "0.34", "0.32"),
    ),
    (
        "energy_demand",
        "GJ/t",
        True,
        "energy demand per tonne of DRI (gas loop)",
        ("12.20", "11.99", "12.47", "11.68", "12.98"),
    ),
)

#: Checks that the printed values of the source do not pass at the default tolerance, with
#: the size of the deviation. Kept as found (not corrected); see the printed report.
#: - B: GUD printed 0.33, the printed top-gas composition gives 0.315 (29.9 / 95.0).
#: - C, E: top-gas flow 7 and 34 mol/s (0.25 %, 1.2 %) below the inlet total, although
#:   H2 -> H2O is equimolar.
KNOWN_DEVIATIONS: dict[tuple[str, str], float] = {
    ("B", "GUD from top-gas composition"): -0.0153,
    ("C", "top-gas flow = bustle + cooling"): -7.0,
    ("E", "top-gas flow = bustle + cooling"): -34.0,
}


def number(printed: str) -> float:
    return float(printed.replace(",", ""))


def consistency(v: dict[str, float], m: dict[str, float]) -> dict[str, tuple[float, float, float]]:
    """(computed, printed, tolerance) of each consistency check of one column."""
    fe_per = {"Fe2O3": 2.0, "Fe3O4": 3.0, "FeO": 1.0, "Fe": 1.0}
    total = sum(v[f"dri_w_{k}"] for k in ("Fe2O3", "Fe3O4", "FeO", "Fe", "gangue"))
    fe_tot = sum(v[f"dri_w_{k}"] / m[k] * fe_per[k] for k in m)
    return {
        "DRI mass fractions sum to 100 %": (total, 100.0, 0.05),
        "metallisation from DRI composition": (
            100 * v["dri_w_Fe"] / m["Fe"] / fe_tot,
            v["metallization"],
            0.1,
        ),
        "GUD from top-gas composition": (
            v["top_gas_x_H2O"] / (v["top_gas_x_H2"] + v["top_gas_x_H2O"]),
            v["gud"],
            0.005,
        ),
        "top-gas flow = bustle + cooling": (
            v["top_gas_flow"],
            v["bustle_flow"] + v["cooling_flow"],
            5.0,
        ),
    }


def main() -> None:
    if not TXT.exists():
        raise SystemExit(f"{TXT.name} not found; run scripts/fetch_references.py")
    tokens = set(re.split(r"\s+", TXT.read_bytes().decode("latin-1")))
    missing = sorted({p for row in ROWS for p in row[4] if p not in tokens})
    if missing:
        raise SystemExit(f"values not found in the text layer: {missing}")
    print("all printed values occur in the text layer")

    sp = yaml.safe_load((ROOT / "data" / "species.yaml").read_text(encoding="utf-8"))["species"]
    m = {k: float(sp[k]["molar_mass"]) for k in ("Fe2O3", "Fe3O4", "FeO", "Fe")}
    found: dict[tuple[str, str], float] = {}
    for j, case in enumerate(CASES):
        v = {k: number(vals[j]) for k, _, _, _, vals in ROWS}
        for name, (got, ref, tol) in consistency(v, m).items():
            ok = abs(got - ref) <= tol
            status = "ok" if ok else "DEVIATES"
            print(f"case {case}: {name}: {got:.4g} vs {ref:.4g} (tol {tol}) {status}")
            if not ok:
                found[(case, name)] = round(got - ref, 4)
    if found != KNOWN_DEVIATIONS:
        raise SystemExit(f"deviations differ from KNOWN_DEVIATIONS: {found}")

    with open(CASE_A, encoding="utf-8") as fh:
        a = {
            r["key"]: float(r["value"])
            for r in csv.DictReader(line for line in fh if not line.startswith("#"))
        }
    for k, _, _, _, vals in ROWS:
        if k in a and a[k] != number(vals[0]):
            raise SystemExit(f"case A {k}: {vals[0]} differs from {CASE_A.name} ({a[k]})")
    print(f"case A column equals {CASE_A.name}")

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write(
            "# REDUCTOR cases A-E, inputs and outputs (Marsigny, Mirgaux & Patisson 2025, Metals\n"
            "# 15:862, Table 1, CC BY 4.0), read from the page image and checked against the PDF\n"
            "# text layer by scripts/extract_reductor_cases.py (do not edit).\n"
            "# output = 1: printed in bold (a model output); 0: plain (an input or a reported\n"
            "# flow). Units as printed.\n"
        )
        fh.write("case,key,unit,value,printed,output,description\n")
        for case_j, case in enumerate(CASES):
            for k, unit, bold, desc, vals in ROWS:
                fh.write(
                    f'{case},{k},{unit},{number(vals[case_j]):g},"{vals[case_j]}",{int(bold)},{desc}\n'
                )
    print(f"wrote {OUT.relative_to(ROOT)}")
    write_inputs()


#: Model inputs per case: (yaml key, Table 1 key, factor, offset) -> SI value.
INPUTS = (
    ("bustle_flow", "bustle_flow", 1.0, 0.0, "mol/s"),
    ("bustle_temperature", "bustle_T", 1.0, 273.15, "K"),
    ("cooling_flow", "cooling_flow", 1.0, 0.0, "mol/s"),
    ("cooling_temperature", "cooling_T", 1.0, 273.15, "K"),
    ("x_H2", "bustle_x_H2", 0.01, 0.0, "-"),
    ("x_H2O", "bustle_x_H2O", 0.01, 0.0, "-"),
    ("x_inert", "bustle_x_N2", 0.01, 0.0, "-"),
    ("solid_flow", "solid_flow", 1.0, 0.0, "kg/s"),
    ("solid_temperature", "solid_T", 1.0, 273.15, "K"),
    ("pellet_diameter", "pellet_diameter", 1e-3, 0.0, "m"),
    ("top_pressure", "top_pressure", 1.0, 0.0, "Pa"),
    ("bustle_pressure", "bustle_pressure", 1.0, 0.0, "Pa"),
    ("cooling_pressure", "cooling_pressure", 1.0, 0.0, "Pa"),
)


def write_inputs() -> None:
    """data/furnace_cases.yaml: the inputs of cases A-E in SI units (generated)."""
    rows = {k: vals for k, _, _, _, vals in ROWS}
    for j, case in enumerate(CASES):  # same gas at both inlets in every case
        for sp in ("x_H2", "x_H2O", "x_N2"):
            if rows[f"bustle_{sp}"][j] != rows[f"cooling_{sp}"][j]:
                raise SystemExit(f"case {case}: bustle and cooling {sp} differ")
    path = ROOT / "data" / "furnace_cases.yaml"
    lines = [
        "# Operating conditions of REDUCTOR cases A-E (Marsigny, Mirgaux & Patisson 2025, Metals",
        "# 15:862, Table 1), GENERATED by scripts/extract_reductor_cases.py (do not edit): read",
        "# from the page image and checked against the PDF text layer. SI units. Gas composition",
        "# is the same at the bustle and cooling inlets in every case. The bustle and cooling",
        "# pressures are REDUCTOR outputs, used as boundary values (one mean pressure, see",
        "# shaftsim.operating). Loaded by shaftsim.operating.load_cases.",
        "cases:",
    ]
    for j, case in enumerate(CASES):
        lines.append(f"  {case}:")
        for key, src, factor, offset, unit in INPUTS:
            val = round(number(rows[src][j]) * factor + offset, 10)
            lines.append(
                f'    {key}: {{value: {val:g}, unit: "{unit}", source: Marsigny2025b-T1-image, '
                f'printed: "{rows[src][j]}"}}'
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
