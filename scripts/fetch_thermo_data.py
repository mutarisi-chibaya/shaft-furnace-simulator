"""Download NIST WebBook Shomate coefficients and NIST-JANAF tables, and write the data files.

Writes:
    data/species.yaml              Shomate coefficients + molar masses (NIST WebBook)
    tests/data/janaf_reference.csv validation rows (NIST-JANAF tables, Chase 1998)

Every number is copied programmatically from the downloaded pages; nothing is typed by hand.

Wüstite: real wüstite is non-stoichiometric. JANAF tabulates it as Fe0.947O (table Fe-001), but
NIST WebBook gives Shomate coefficients only for stoichiometric FeO. Shomate coefficients for
Fe0.947O are therefore FITTED here to the JANAF Fe-001 table, the same way Shomate fits are
built: A-E by least squares on Cp (relative residuals), then F and G fixed exactly by
H(298.15) - H298 = 0 and S(298.15) = S298(JANAF); H = dfH298(JANAF).
Raw downloads are cached in ``--cache`` (default: build/thermo_cache) so the extraction can be
audited and rerun offline.

Usage:
    python scripts/fetch_thermo_data.py [--cache DIR] [--offline]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import html
import re
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
USER_AGENT = "shaftsim-data-fetch/0.1 (research script)"

NIST_URL = "https://webbook.nist.gov/cgi/cbook.cgi?ID={id}&Units=SI&Mask={mask}"
JANAF_URL = "https://janaf.nist.gov/tables/{table}.txt"

# name: (NIST WebBook ID, mask 1=gas / 2=condensed, Shomate table label, JANAF table,
#        phase, elements)
SPECIES = {
    "H2": ("C1333740", 1, "Gas Phase Heat Capacity (Shomate Equation)", "H-050", "gas", {"H": 2}),
    "H2O": (
        "C7732185",
        1,
        "Gas Phase Heat Capacity (Shomate Equation)",
        "H-064",
        "gas",
        {"H": 2, "O": 1},
    ),
    # Inert carrier gases of the TGA experiments (He) and of industrial gas (N2).
    "He": ("C7440597", 1, "Gas Phase Heat Capacity (Shomate Equation)", "He-001", "gas", {"He": 1}),
    "N2": ("C7727379", 1, "Gas Phase Heat Capacity (Shomate Equation)", "N-023", "gas", {"N": 2}),
    "Fe": (
        "C7439896",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        "Fe-004",
        "solid",
        {"Fe": 1},
    ),
    "FeO": (
        "C1345251",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        "Fe-018",
        "solid",
        {"Fe": 1, "O": 1},
    ),
    "Fe3O4": (
        "C1309382",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        "Fe-032",
        "solid",
        {"Fe": 3, "O": 4},
    ),
    "Fe2O3": (
        "C1317608",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        "Fe-030",
        "solid",
        {"Fe": 2, "O": 3},
    ),
    # Gangue oxides (heat capacity of the burden, Milestone 4). No JANAF validation rows
    # (None): 3.35 wt-% of the pellets, so their effect on the energy balance is small.
    "CaO": (
        "C1305788",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        None,
        "solid",
        {"Ca": 1, "O": 1},
    ),
    "SiO2": (
        "C14808607",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        None,
        "solid",
        {"Si": 1, "O": 2},
    ),
    "MgO": (
        "C1309484",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        None,
        "solid",
        {"Mg": 1, "O": 1},
    ),
    "Al2O3": (
        "C1344281",
        2,
        "Solid Phase Heat Capacity (Shomate Equation)",
        None,
        "solid",
        {"Al": 2, "O": 3},
    ),
}

# NIST also lists, for Fe, a single "gamma phase" fit over 298-1809 K that overlaps the four
# alpha-delta ranges. It describes metastable gamma-iron, so only columns whose comment
# mentions the phase we want are kept.
PHASE_COMMENT_FILTER = {
    "Fe": chr(0x3B1) + "-" + chr(0x3B4),  # "alpha-delta" as NIST writes it
    "Al2O3": chr(0x3B1) + " phase",  # corundum; NIST also lists gamma, delta, kappa fits
}

# Validation temperatures [K]; rows outside a species' Shomate range are left out.
JANAF_TEMPERATURES = [298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]

COEFFS = "ABCDEFGH"

# Wüstite fit to JANAF Fe-001 (Fe0.947O, crystal). Upper limit kept below its melting point
# (~1650 K); above ~1650 K JANAF continues the crystal table as a metastable extrapolation.
WUSTITE: dict[str, Any] = {
    "name": "Fe0.947O",
    "x_fe": 0.947,
    "janaf_table": "Fe-001",
    "T_min": 298.15,
    "T_max": 1600.0,
}
TIMES = chr(0xD7)  # multiplication sign used by NIST in "1.2 x 10^-8"


def fetch(url: str, path: Path, offline: bool) -> str:
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as resp:
            path.write_bytes(resp.read())
    return path.read_text(encoding="utf-8", errors="replace")


def _cell(raw: str) -> str:
    """NIST cell -> number string: '-1.884667&times;10<sup>-8</sup>' -> '-1.884667e-8'."""
    c = html.unescape(re.sub(r"<sup>(.*?)</sup>", r"^\1", raw)).strip()
    c = re.sub(r"<[^>]+>", "", c)
    m = re.fullmatch(r"(-?[\d.]+)\s*" + TIMES + r"\s*10\^(-?\d+)", c)
    return f"{m.group(1)}e{m.group(2)}" if m else c


def parse_shomate(page: str, label: str, phase_filter: str | None) -> list[dict[str, str]]:
    m = re.search(rf'<table class="data" aria-label="{re.escape(label)}"(.*?)</table>', page, re.S)
    if not m:
        raise ValueError(f"table '{label}' not found")
    rows: dict[str, list[str]] = {}
    for r in re.finditer(r"<tr><th[^>]*>(.*?)</th>(.*?)</tr>", m.group(1), re.S):
        key = re.sub(r"<[^>]+>", "", r.group(1)).strip()
        rows[key] = [_cell(c) for c in re.findall(r"<td[^>]*>(.*?)</td>", r.group(2), re.S)]
    out = []
    for i, trange in enumerate(rows["Temperature (K)"]):
        comment = rows.get("Comment", [""] * (i + 1))[i]
        if phase_filter and phase_filter not in comment:
            continue
        lo, hi = (float(x) for x in re.findall(r"[\d.]+", trange))
        col = {"T_min": repr(lo), "T_max": repr(hi)}
        col |= {k: rows[k][i] for k in COEFFS}
        for k in COEFFS:
            float(col[k])  # fail loudly if a cell did not parse
        col["reference"] = rows["Reference"][i]
        col["comment"] = comment
        out.append(col)
    return sorted(out, key=lambda c: float(c["T_min"]))


def parse_molar_mass(page: str) -> float:
    m = re.search(r"Molecular weight</a>:</strong>\s*([\d.]+)", page)
    if not m:
        raise ValueError("molecular weight not found")
    return float(f"{float(m.group(1)) / 1000.0:.10g}")  # g/mol -> kg/mol, no float noise


def parse_janaf(text: str) -> dict[float, list[float]]:
    """JANAF .txt -> {T: [Cp, S, -(G-H)/T, H-H298, dfH, dfG, logKf]} (numeric rows only)."""
    rows: dict[float, list[float]] = {}
    for line in text.splitlines()[2:]:
        parts = line.split()
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue  # transition / annotation rows
        if len(vals) == 8 and vals[0] not in rows:
            rows[vals[0]] = vals[1:]
    return rows


def fit_shomate(rows: dict[float, list[float]], t_min: float, t_max: float) -> dict[str, float]:
    """Fit Shomate A-E to JANAF Cp (relative least squares); F, G, H from the 298.15 K values.

    Returns coefficients in NIST units (kJ/mol for F, H; J/(mol K) for G).
    """
    import numpy as np

    T = np.array([t for t in sorted(rows) if t_min <= t <= t_max])
    cp = np.array([rows[t][0] for t in T])
    t = T / 1000.0
    basis = np.column_stack([np.ones_like(t), t, t**2, t**3, 1 / t**2])
    coef, *_ = np.linalg.lstsq(basis / cp[:, None], np.ones_like(cp), rcond=None)
    a, b, c, d, e = (float(v) for v in coef)
    r = rows[298.15]
    t0 = 0.29815
    dfh298 = r[4]
    f = dfh298 - (a * t0 + b * t0**2 / 2 + c * t0**3 / 3 + d * t0**4 / 4 - e / t0)
    g = r[1] - (a * np.log(t0) + b * t0 + c * t0**2 / 2 + d * t0**3 / 3 - e / (2 * t0**2))
    return {"A": a, "B": b, "C": c, "D": d, "E": e, "F": f, "G": float(g), "H": dfh298}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "build" / "thermo_cache")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    today = dt.date.today().isoformat()

    yaml_lines = [
        "# GENERATED by scripts/fetch_thermo_data.py on " + today + " - do not edit by hand.",
        "# Source: NIST Chemistry WebBook, SRD 69 (Shomate coefficients, all citing Chase 1998,",
        "# NIST-JANAF Thermochemical Tables 4th ed.). See data/SOURCES.md for URLs.",
        "#",
        "# Conventions (NIST units; the code converts to SI): t = T/1000,",
        "#   Cp = A + B t + C t^2 + D t^3 + E/t^2                     [J/(mol K)]",
        "#   H - H298 = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F - H [kJ/mol]",
        "#   S = A ln t + B t + C t^2/2 + D t^3/3 - E/(2 t^2) + G       [J/(mol K)]",
        "# molar_mass in kg/mol (NIST 'Molecular weight' / 1000).",
        "species:",
    ]
    csv_rows: list[list[object]] = []
    for name, (nid, mask, label, jtable, phase, elements) in SPECIES.items():
        url = NIST_URL.format(id=nid, mask=mask)
        page = fetch(url, args.cache / f"nist_{name}.html", args.offline)
        ranges = parse_shomate(page, label, PHASE_COMMENT_FILTER.get(name))
        mm = parse_molar_mass(page)
        els = ", ".join(f"{k}: {v}" for k, v in elements.items())
        yaml_lines += [
            f"  {name}:",
            f"    phase: {phase}",
            f"    elements: {{{els}}}",
            f"    molar_mass: {mm!r}",
            f"    source: NIST-{name}",
            "    shomate:",
        ]
        for r in ranges:
            coeffs = ", ".join(f"{k}: {r[k]}" for k in COEFFS)
            yaml_lines.append(
                f"      - {{T_min: {r['T_min']}, T_max: {r['T_max']}, {coeffs}}}"
                f"  # {r['reference']}; {r['comment']}"
            )

        if jtable is None:
            continue
        janaf = parse_janaf(
            fetch(JANAF_URL.format(table=jtable), args.cache / f"janaf_{jtable}.txt", args.offline)
        )
        t_lo, t_hi = float(ranges[0]["T_min"]), float(ranges[-1]["T_max"])
        for T in JANAF_TEMPERATURES:
            if t_lo <= T <= t_hi and float(T) in janaf:
                cp, s, _, dh, _, dfg, _ = janaf[float(T)]
                csv_rows.append([name, T, cp, s, dh, dfg, f"JANAF-{jtable}"])

    # --- wüstite Fe0.947O, fitted to JANAF (see module docstring)
    wu = WUSTITE
    jt = str(wu["janaf_table"])
    wrows = parse_janaf(
        fetch(JANAF_URL.format(table=jt), args.cache / f"janaf_{jt}.txt", args.offline)
    )
    fit = fit_shomate(wrows, float(wu["T_min"]), float(wu["T_max"]))
    # Molar mass from the NIST values already read: M(O) = M(FeO) - M(Fe).
    m_fe = parse_molar_mass(fetch("", args.cache / "nist_Fe.html", True))
    m_feo = parse_molar_mass(fetch("", args.cache / "nist_FeO.html", True))
    x = float(wu["x_fe"])
    m_w = float(f"{x * m_fe + (m_feo - m_fe):.10g}")
    coeffs = ", ".join(f"{k}: {fit[k]:.10g}" for k in COEFFS)
    yaml_lines += [
        f"  {wu['name']}:",
        "    phase: solid",
        f"    elements: {{Fe: {x}, O: 1}}",
        f"    molar_mass: {m_w!r}  # {x} M(Fe) + [M(FeO) - M(Fe)], NIST molecular weights",
        f"    source: JANAF-{jt}-fit",
        "    shomate:",
        f"      - {{T_min: {wu['T_min']}, T_max: {wu['T_max']}, {coeffs}}}"
        f"  # fitted by scripts/fetch_thermo_data.py to JANAF {jt} (Chase, 1998)",
    ]
    for T in JANAF_TEMPERATURES:
        if float(wu["T_min"]) <= T <= float(wu["T_max"]) and float(T) in wrows:
            cp, s_, _, dh, _, dfg, _ = wrows[float(T)]
            csv_rows.append([wu["name"], T, cp, s_, dh, dfg, f"JANAF-{jt} (fitted)"])

    (ROOT / "data" / "species.yaml").write_text("\n".join(yaml_lines) + "\n", encoding="utf-8")
    with open(
        ROOT / "tests" / "data" / "janaf_reference.csv", "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "species",
                "T_K",
                "Cp_J_molK",
                "S_J_molK",
                "H_minus_H298_kJ_mol",
                "dfG_kJ_mol",
                "source",
            ]
        )
        w.writerows(csv_rows)
    print(f"wrote data/species.yaml and tests/data/janaf_reference.csv ({len(csv_rows)} rows)")


if __name__ == "__main__":
    main()
