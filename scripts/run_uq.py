"""Monte Carlo, Sobol and grid-check runs for Milestone 6 (whole furnace, case A).

Usage:
    python scripts/run_uq.py mc     [--n 512] [--cells 50] [--seed 2026] [--workers 6]
    python scripts/run_uq.py sobol  [--n 256] [--cells 50] [--seed 2026] [--workers 6]
    python scripts/run_uq.py grid   [--n 24]  [--cells 200] [--seed 2026] [--workers 6]
    python scripts/run_uq.py retry  --file results/uq/mc_50cells.npz [--workers 6]

mc:    Latin hypercube sample of the 8 uncertain inputs (shaftsim.uq).
sobol: Saltelli/Sobol design for the three groups (kinetics, tortuosity, transfer
       correlation), first-order and total indices with bootstrap confidence intervals.
grid:  the first n points of the mc design on a finer grid, to check that the spreads (not
       the offsets) are grid-independent.
retry: re-solve only the runs of an existing file that are flagged unconverged (with the
       fallback of shaftsim.uq.run_one), replace them in place and, for a Sobol file,
       recompute the indices.

Results go to results/uq/<kind>_<cells>cells.npz (inputs, outputs, profiles, run times).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import numpy as np
from SALib.analyze import sobol as sobol_analyze
from SALib.sample import sobol as sobol_sample

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr
from shaftsim import uq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "uq"
PROFILE_KEYS = ("Ts", "Tg", "x_H2", "X3")


def setup(cells: int) -> tuple[en.ThermalProblem, uq.UncertaintySpec]:
    warnings.simplefilter("ignore")
    sp = pr.load_species(DATA / "species.yaml")
    td = tr.load_transport(
        DATA / "transport.yaml",
        DATA / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    trip = pl.load_refit_triplets(DATA / "kinetics_refit.yaml")
    case = en.load_furnace_case(DATA / "furnace.yaml", DATA / "kinetics.yaml")
    spec = uq.load_uncertainty(
        DATA / "uncertainty.yaml", DATA / "kinetics_refit.yaml", case.shaft.pellet.tortuosity
    )
    return en.ThermalProblem(case, sp, td, trip, cells), spec


def salib_problem() -> dict[str, object]:
    return {
        "num_vars": len(uq.NAMES),
        "names": list(uq.NAMES),
        "bounds": [[0.0, 1.0]] * len(uq.NAMES),
        "groups": list(uq.GROUPS),
    }


def sobol_indices(outputs: np.ndarray, seed: int) -> dict[str, dict[str, list[float]]]:
    """Grouped first-order and total Sobol indices of every output (95 % bootstrap CIs)."""
    indices = {}
    for k, name in enumerate(uq.OUTPUTS):
        si = sobol_analyze.analyze(
            salib_problem(),
            outputs[:, k],
            calc_second_order=False,
            num_resamples=1000,
            conf_level=0.95,
            seed=seed,
        )
        indices[name] = {
            key: np.asarray(si[key]).tolist() for key in ("S1", "S1_conf", "ST", "ST_conf")
        }
    return indices


def progress_printer() -> tuple[float, object]:
    t0 = time.perf_counter()
    last = [t0]

    def progress(k: int, n: int) -> None:
        now = time.perf_counter()
        if now - last[0] > 30 or k == n:
            print(f"  {k}/{n} runs, {60 * k / (now - t0):.1f} runs/min", flush=True)
            last[0] = now

    return t0, progress


def run(
    base: en.ThermalProblem, spec: uq.UncertaintySpec, pts: np.ndarray, workers: int
) -> dict[str, np.ndarray]:
    t0, progress = progress_printer()
    res = uq.run_points(base, spec, pts, workers=workers, progress=progress)  # type: ignore[arg-type]
    wall = time.perf_counter() - t0
    params = [uq.to_parameters(spec, u) for u in pts]
    return {
        "u": pts,
        "theta": np.array([p.theta for p in params]),
        "tortuosity": np.array([p.tortuosity for p in params]),
        "transfer_mix": np.array([p.transfer_mix for p in params]),
        "outputs": uq.output_matrix(res),
        "output_names": np.array(uq.OUTPUTS),
        "Ts": np.array([r.Ts for r in res]),
        "Tg": np.array([r.Tg for r in res]),
        "x_H2": np.array([r.x_H2 for r in res]),
        "X3": np.array([r.X3 for r in res]),
        "z_depth": base.grid.z,
        "converged": np.array([r.converged for r in res]),
        "retried": np.array([r.retried for r in res]),
        "seconds": np.array([r.seconds for r in res]),
        "wall_seconds": np.array(wall),
        "workers": np.array(workers),
        "runs_per_minute": np.array(60 * len(pts) / wall),
    }


def retry(path: Path, workers: int) -> None:
    """Re-solve the flagged runs of an existing results file and update it in place."""
    with np.load(path) as d:
        data = {k: d[k] for k in d.files}
    meta = json.loads(str(data["meta"]))
    base, spec = setup(int(meta["cells"]))
    bad = np.where(~data["converged"])[0]
    print(f"retry: {bad.size} flagged runs in {path.name}: {bad.tolist()}", flush=True)
    if bad.size == 0:
        return
    _, progress = progress_printer()
    res = uq.run_points(base, spec, data["u"][bad], workers=workers, progress=progress)  # type: ignore[arg-type]
    if "retried" not in data:
        data["retried"] = np.zeros(data["converged"].size, dtype=bool)
    new = uq.output_matrix(res)
    for j, (i, r) in enumerate(zip(bad, res, strict=True)):
        print(
            f"  run {i}: converged {r.converged}, retried {r.retried}, "
            f"metallisation {data['outputs'][i, 0]:.8f} -> {new[j, 0]:.8f}"
        )
        data["outputs"][i] = new[j]
        for key in PROFILE_KEYS:
            data[key][i] = getattr(r, key)
        data["converged"][i] = r.converged
        data["retried"][i] = True
        data["seconds"][i] += r.seconds  # time of the flagged first solve plus this one
    if meta["kind"] == "sobol":
        meta["indices"] = sobol_indices(data["outputs"], int(meta["seed"]))
        print(json.dumps(meta["indices"], indent=1))
    meta["retried_runs"] = bad.tolist()
    data["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(path, **data)
    ok = data["converged"]
    print(f"updated {path.name}: {ok.sum()}/{ok.size} converged")
    if not ok.all():
        sys.exit("some runs still did not converge (kept in the file, flagged)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=("mc", "sobol", "grid", "retry"))
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--cells", type=int, default=None)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--file", type=Path, default=None)
    a = ap.parse_args()
    if a.kind == "retry":
        if a.file is None:
            ap.error("retry needs --file")
        retry(a.file, a.workers)
        return
    n = a.n or {"mc": 512, "sobol": 256, "grid": 24}[a.kind]
    cells = a.cells or {"mc": 50, "sobol": 50, "grid": 200}[a.kind]
    base, spec = setup(cells)
    if a.kind in ("mc", "grid"):
        pts = uq.latin_hypercube(512 if a.kind == "grid" else n, seed=a.seed)[:n]
    else:
        pts = sobol_sample.sample(salib_problem(), n, calc_second_order=False, seed=a.seed)
    print(f"{a.kind}: {len(pts)} runs on {cells} cells, {a.workers} workers", flush=True)
    data = run(base, spec, pts, a.workers)
    meta: dict[str, object] = {"kind": a.kind, "n": n, "cells": cells, "seed": a.seed}
    if a.kind == "sobol":
        meta["groups"] = sorted(set(uq.GROUPS), key=uq.GROUPS.index)
        meta["indices"] = sobol_indices(data["outputs"], a.seed)
        print(json.dumps(meta["indices"], indent=1))
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{a.kind}_{cells}cells.npz"
    arrays: dict[str, Any] = {"meta": np.array(json.dumps(meta)), **data}
    np.savez_compressed(path, **arrays)
    ok = data["converged"]
    print(
        f"wrote {path.relative_to(ROOT)}: {ok.sum()}/{ok.size} converged "
        f"({int(data['retried'].sum())} after the fallback tolerance), "
        f"{float(data['runs_per_minute']):.1f} runs/min, "
        f"median {np.median(data['seconds']):.1f} s/run"
    )
    if not ok.all():
        sys.exit("some runs did not converge (kept in the file, flagged)")


if __name__ == "__main__":
    main()
