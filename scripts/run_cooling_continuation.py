"""Milestone 7: the cooling-gas flow sweep of Fig. 10 by natural-parameter continuation.

Usage:
    python scripts/run_cooling_continuation.py [--cells 200] [--step 25] [--max-outer 150]

Case B operating conditions with the cooling gas flow stepped from 500 to 800 mol/s (bustle
flow = 2800 - cooling flow), once upward and once downward. Each step starts from the last
converged temperatures of the same direction. Results: results/cases/
cooling_continuation_<cells>cells.json (one record per step and direction; unconverged steps
are kept with converged = false and are not solutions).
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import numpy as np

from shaftsim import energy as en
from shaftsim import operating as op
from shaftsim import uq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "cases"
sys.path.insert(0, str(ROOT / "scripts"))
import run_uq  # noqa: E402


def sweep(cells: int, flows: np.ndarray, max_outer: int, direction: str) -> list[dict[str, Any]]:
    base, _ = run_uq.setup(cells)
    b = op.load_cases(ROOT / "data" / "furnace_cases.yaml", base.case)["B"]
    total = b.bustle_flow + b.cooling_flow
    T: tuple[np.ndarray, np.ndarray] | None = None
    out = []
    for cool in flows:
        oc = dataclasses.replace(b, cooling_flow=float(cool), bustle_flow=total - float(cool))
        prob = en.ThermalProblem(
            op.apply(base.case, oc), base.species, base.td, base.triplets, cells
        )
        t0 = time.perf_counter()
        sol = en.solve_thermal(
            prob,
            tol_T=float(uq.SOLVER["tol_T"]),
            species_tol=float(uq.SOLVER["species_tol"]),
            T_start=T,
            max_outer=max_outer,
        )
        rec = {
            "direction": direction,
            "cooling_flow": float(cool),
            "converged": bool(sol.converged),
            "outer_iterations": sol.outer_iterations,
            "metallization": sol.species.metallization,
            "dri_T_C": sol.T_solid_out - 273.15,
            "top_gas_T_C": sol.T_top_gas - 273.15,
            "seconds": time.perf_counter() - t0,
        }
        print(json.dumps(rec), flush=True)
        out.append(rec)
        if sol.converged:
            T = (sol.Ts, sol.Tg)
    return out


def main() -> None:
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=200)
    ap.add_argument("--step", type=float, default=25.0)
    ap.add_argument("--max-outer", type=int, default=150)
    a = ap.parse_args()
    flows = np.arange(500.0, 800.0 + 1e-9, a.step)
    res = sweep(a.cells, flows, a.max_outer, "up") + sweep(
        a.cells, flows[::-1], a.max_outer, "down"
    )
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"cooling_continuation_{a.cells}cells.json"
    path.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
