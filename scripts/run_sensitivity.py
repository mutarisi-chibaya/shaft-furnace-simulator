"""Milestone 7: one-at-a-time sensitivity of the case A furnace to bustle gas temperature,
bustle gas flow and pellet diameter, with energy demand and CO2 per tonne of DRI.

Usage:
    python scripts/run_sensitivity.py run   [--cells 50 200] [--points 6] [--workers 6]
    python scripts/run_sensitivity.py table

run:   sweeps each input over its range in data/sensitivity.yaml (the others at case A) and
       saves results/sensitivity/sensitivity_<cells>cells.json.
table: prints the change of each output over each range, per grid.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from shaftsim import emissions as em
from shaftsim import energy as en
from shaftsim import operating as op
from shaftsim import uq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "sensitivity"
sys.path.insert(0, str(ROOT / "scripts"))
import run_uq  # noqa: E402

#: input key in data/sensitivity.yaml -> OperatingConditions field
INPUTS = {
    "bustle_temperature": "bustle_T",
    "bustle_flow": "bustle_flow",
    "pellet_diameter": "pellet_diameter",
}


def ranges() -> dict[str, tuple[float, float]]:
    raw = yaml.safe_load((DATA / "sensitivity.yaml").read_text(encoding="utf-8"))
    return {k: (float(raw[k]["low"]["value"]), float(raw[k]["high"]["value"])) for k in INPUTS}


def solve(args: tuple[str, float, int]) -> dict[str, Any]:
    """One run: input `key` set to `value`, everything else at case A."""
    warnings.simplefilter("ignore")
    key, value, cells = args
    base, _ = run_uq.setup(cells)
    change: dict[str, Any] = {INPUTS[key]: value}
    oc = dataclasses.replace(op.conditions_of(base.case), **change)
    case = op.apply(base.case, oc)
    prob = en.ThermalProblem(case, base.species, base.td, base.triplets, cells)
    sol = en.solve_thermal(prob, **uq.SOLVER)  # type: ignore[arg-type]
    if not sol.converged:
        sol = en.solve_thermal(
            prob, tol_T=uq.TOL_T_RETRY, species_tol=1e-8, T_start=(sol.Ts, sol.Tg)
        )
    data = em.load_emissions(DATA / "emissions.yaml", DATA / "grid_intensity.yaml")
    d = em.furnace_energy_demand(case, sol, base.species, data)
    s = sol.species
    return {
        "input": key,
        "value": value,
        "converged": bool(sol.converged),
        "metallization": s.metallization,
        "h2_consumption_kg_per_t": s.h2_consumption,
        "top_gas_T_C": sol.T_top_gas - 273.15,
        "dri_T_C": sol.T_solid_out - 273.15,
        "electrolysis_GJ_per_t": d.electrolysis / 1e6,
        "heating_GJ_per_t": d.heating / 1e6,
        "co2_south_africa_t_per_t": em.co2_per_dri(d, data.grids["south_africa"].value),
    }


def run(cells_list: list[int], points: int, workers: int) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs0 = [(k, float(v)) for k, (lo, hi) in ranges().items() for v in np.linspace(lo, hi, points)]
    for cells in cells_list:
        print(f"{len(jobs0)} runs on {cells} cells", flush=True)
        with ProcessPoolExecutor(workers) as ex:
            res = list(ex.map(solve, [(k, v, cells) for k, v in jobs0]))
        bad = [r for r in res if not r["converged"]]
        print(f"  {len(res) - len(bad)}/{len(res)} converged", flush=True)
        (OUT / f"sensitivity_{cells}cells.json").write_text(
            json.dumps(res, indent=1), encoding="utf-8"
        )


def load(cells: int) -> list[dict[str, Any]] | None:
    p = OUT / f"sensitivity_{cells}cells.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


SHOW = (
    "metallization",
    "h2_consumption_kg_per_t",
    "electrolysis_GJ_per_t",
    "heating_GJ_per_t",
    "co2_south_africa_t_per_t",
    "dri_T_C",
)


def table() -> None:
    for cells in (50, 200, 400):
        res = load(cells)
        if res is None:
            continue
        print(f"\n{cells} cells: output at the low / high end of each range")
        for key in INPUTS:
            rows = sorted((r for r in res if r["input"] == key), key=lambda r: r["value"])
            lo, hi = rows[0], rows[-1]
            print(f"  {key} {lo['value']:g} -> {hi['value']:g}")
            for o in SHOW:
                print(f"    {o:26s} {lo[o]:9.4g} -> {hi[o]:9.4g}  (change {hi[o] - lo[o]:+.4g})")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=("run", "table"))
    ap.add_argument("--cells", type=int, nargs="+", default=[50, 200])
    ap.add_argument("--points", type=int, default=6)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.action == "run":
        run(a.cells, a.points, a.workers)
    else:
        table()


if __name__ == "__main__":
    main()
