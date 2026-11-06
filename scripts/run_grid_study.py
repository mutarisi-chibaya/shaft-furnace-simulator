"""Milestone 8: grid study of the whole-furnace case A, 50 to 800 cells.

Usage:
    python scripts/run_grid_study.py [--cells 50 100 200 400 800]

Each grid is solved by the coupled solver (shaftsim.coupled; default smoothing, effect below
1e-6), started from the solution on the previous grid (the first from the alternating
solver). Saves results/cases/grid_caseA.json (one record per grid) and prints Richardson
estimates from the last three grids (observed order from their differences).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
import warnings
from pathlib import Path
from typing import Any

from shaftsim import coupled as cp
from shaftsim import energy as en

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "cases"
sys.path.insert(0, str(ROOT / "scripts"))
import run_uq  # noqa: E402

KEYS = ("metallization", "dri_T_C", "top_gas_T_C", "gud", "h2_consumption_kg_per_t")


def run(cells: list[int]) -> list[dict[str, Any]]:
    """Solve case A on each grid in turn (warm start from the previous grid)."""
    out: list[dict[str, Any]] = []
    prev: en.ThermalSolution | None = None
    for n in cells:
        base, _ = run_uq.setup(n)
        prob = en.ThermalProblem(base.case, base.species, base.td, base.triplets, n)
        t0 = time.perf_counter()
        if prev is None:
            sol, st = cp.solve_thermal_coupled(prob)
        else:
            sol, st = cp.finish_coupled(prob, prev)
        s = sol.species
        rec = {
            "cells": n,
            "converged": bool(sol.converged),
            "iterations": st.iterations,
            "seconds": time.perf_counter() - t0,
            "metallization": s.metallization,
            "dri_T_C": sol.T_solid_out - 273.15,
            "top_gas_T_C": sol.T_top_gas - 273.15,
            "gud": s.gud,
            "h2_consumption_kg_per_t": s.h2_consumption,
            "energy_closure": sol.energy_closure,
        }
        print(json.dumps(rec), flush=True)
        out.append(rec)
        prev = sol
    return out


def richardson(res: list[dict[str, Any]], key: str) -> tuple[float, float]:
    """(extrapolated value, observed order) from the last three grids (refinement ratio 2).

    Returns NaN where the differences do not shrink monotonically (no asymptotic range).
    """
    a, b, c = (r[key] for r in res[-3:])
    d1, d2 = b - a, c - b
    if d1 == 0.0 or d2 == 0.0 or d2 / d1 <= 0.0 or abs(d2) >= abs(d1):
        return math.nan, math.nan
    p = math.log2(d1 / d2)
    return c + d2 / (2.0**p - 1.0), p


def main() -> None:
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, nargs="+", default=[50, 100, 200, 400, 800])
    a = ap.parse_args()
    res = run(a.cells)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "grid_caseA.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    if len(res) >= 3:
        for k in KEYS:
            ext, p = richardson(res, k)
            print(f"{k}: extrapolated {ext:.5g} (observed order {p:.2f})")


if __name__ == "__main__":
    main()
