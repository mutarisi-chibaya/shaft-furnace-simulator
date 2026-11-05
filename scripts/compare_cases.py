"""Milestone 7: the 1D model on REDUCTOR cases A-E and on the cooling-flow sweep of Fig. 10.

Usage:
    python scripts/compare_cases.py run   [--cells 50 200 400] [--workers 6]
    python scripts/compare_cases.py table

run:   solves every case (data/furnace_cases.yaml) and the Fig. 10 sweep (case B with the
       cooling gas flow 500-800 mol/s at a constant total of 2800 mol/s) on each grid, in
       parallel, and saves results/cases/cases_<cells>cells.npz.
table: prints the comparison with Table 1 and Fig. 10 (tests/data), as changes from case A.

Nothing is adjusted: every run uses the parameters of case A except the operating conditions
of Table 1.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from shaftsim import energy as en
from shaftsim import operating as op
from shaftsim import uq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "cases"
TESTS = ROOT / "tests" / "data"
sys.path.insert(0, str(ROOT / "scripts"))
import run_uq  # noqa: E402

SWEEP = (500.0, 600.0, 700.0, 800.0)
#: Table 1 key -> (model output key, factor model -> Table 1 units)
OUTPUTS = {
    "metallization": ("metallization", 100.0),
    "top_gas_T": ("top_gas_T_C", 1.0),
    "dri_T": ("dri_T_C", 1.0),
    "gud": ("gud", 1.0),
}


def runs() -> dict[str, op.OperatingConditions]:
    """Named operating points: cases A-E and the Fig. 10 sweep."""
    base, _ = run_uq.setup(10)
    cases = op.load_cases(DATA / "furnace_cases.yaml", base.case)
    out = dict(cases)
    b = cases["B"]
    for flow in SWEEP:
        out[f"B_cool{flow:g}"] = dataclasses.replace(
            b, cooling_flow=flow, bustle_flow=b.bustle_flow + b.cooling_flow - flow
        )
    return out


def solve(args: tuple[str, op.OperatingConditions, int]) -> tuple[str, uq.RunResult, float]:
    warnings.simplefilter("ignore")
    name, oc, cells = args
    base, _ = run_uq.setup(cells)
    case = op.apply(base.case, oc)
    prob = en.ThermalProblem(case, base.species, base.td, base.triplets, cells)
    r = uq.run_one(prob)
    return name, r, r.energy_closure


def run(cells_list: list[int], workers: int) -> None:
    pts = runs()
    OUT.mkdir(parents=True, exist_ok=True)
    for cells in cells_list:
        print(f"{len(pts)} runs on {cells} cells", flush=True)
        jobs = [(n, oc, cells) for n, oc in pts.items()]
        with ProcessPoolExecutor(min(workers, len(jobs))) as ex:
            res = list(ex.map(solve, jobs))
        for name, r, closure in res:
            print(
                f"  {name:12s} converged {r.converged} (retried {r.retried}) {r.seconds:6.0f} s, "
                f"metallisation {r.outputs['metallization']:.4f}, energy closure {closure:.1e}",
                flush=True,
            )
        np.savez_compressed(
            OUT / f"cases_{cells}cells.npz",
            names=np.array([n for n, _, _ in res]),
            outputs=uq.output_matrix([r for _, r, _ in res]),
            output_names=np.array(uq.OUTPUTS),
            Ts=np.array([r.Ts for _, r, _ in res]),
            Tg=np.array([r.Tg for _, r, _ in res]),
            x_H2=np.array([r.x_H2 for _, r, _ in res]),
            X3=np.array([r.X3 for _, r, _ in res]),
            converged=np.array([r.converged for _, r, _ in res]),
            energy_closure=np.array([c for _, _, c in res]),
            seconds=np.array([r.seconds for _, r, _ in res]),
            conditions=np.array(json.dumps({n: dataclasses.asdict(oc) for n, oc in pts.items()})),
            z_depth=run_uq.setup(cells)[0].grid.z,
        )


def load_results(cells: int) -> dict[str, dict[str, float]] | None:
    """Model outputs per run name (Table 1 units), or None if not computed. Unconverged runs
    are NaN (they are not solutions; see docs/validation.md, Milestone 7)."""
    p = OUT / f"cases_{cells}cells.npz"
    if not p.exists():
        return None
    with np.load(p) as d:
        names = [str(n) for n in d["names"]]
        onames = [str(n) for n in d["output_names"]]
        conv = d["converged"]
        out = {}
        for i, n in enumerate(names):
            row = dict(zip(onames, d["outputs"][i], strict=True))
            out[n] = {
                k: float(row[m]) * f if conv[i] else float("nan") for k, (m, f) in OUTPUTS.items()
            }
        return out


def reference_cases() -> dict[str, dict[str, float]]:
    """Table 1 outputs per case (tests/data/reference_model_marsigny2025b_cases.csv)."""
    with open(TESTS / "reference_model_marsigny2025b_cases.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    out: dict[str, dict[str, float]] = {}
    for r in rows:
        if r["key"] in OUTPUTS:
            out.setdefault(r["case"], {})[r["key"]] = float(r["value"])
    return out


def reference_fig10() -> dict[float, float]:
    with open(TESTS / "reference_model_marsigny2025b_fig10.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    return {float(r["cooling_flow_mol_s"]): float(r["metallization_pct"]) for r in rows}


def table() -> None:
    ref = reference_cases()
    fig10 = reference_fig10()
    for cells in (50, 100, 200, 400):
        res = load_results(cells)
        if res is None:
            continue
        print(f"\n{cells} cells: change from case A (model / REDUCTOR)")
        print(f"{'case':6s}" + "".join(f"{k:>26s}" for k in OUTPUTS))
        for c in "BCDE":
            cells_txt = []
            for k in OUTPUTS:
                dm = res[c][k] - res["A"][k]
                dr = ref[c][k] - ref["A"][k]
                cells_txt.append(f"{dm:+10.3g} / {dr:+8.3g}")
            print(f"{c:6s}" + "".join(f"{t:>26s}" for t in cells_txt))
        print(
            f"{'abs A':6s}"
            + "".join(f"{res['A'][k]:>15.4g} / {ref['A'][k]:<8.4g}" for k in OUTPUTS)
        )
        print("Fig. 10 sweep, metallisation (%): cooling flow, model, REDUCTOR")
        for flow in SWEEP:
            print(
                f"  {flow:5g}  {res[f'B_cool{flow:g}']['metallization']:7.2f}  {fig10[flow]:7.2f}"
            )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=("run", "table"))
    ap.add_argument("--cells", type=int, nargs="+", default=[50, 200, 400])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.action == "run":
        run(a.cells, a.workers)
    else:
        table()


if __name__ == "__main__":
    main()
