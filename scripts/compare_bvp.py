"""Compare the finite-volume shaft solver with scipy.integrate.solve_bvp (PLAN.md, Milestone 3).

The continuous problem, z downward from the top:
    dX_j/dz = r_j(X, f) / u_s,   df/dz = c dXbar/dz = c sum_j w_j r_j / u_s
    X_j(0) = 0,                  f(H) = f_in
is solved with solve_bvp (collocation, adaptive mesh) from
  (a) a naive initial guess (no reduction, inlet gas everywhere), and
  (b) the converged finite-volume profile (an expert guess),
and compared with the finite-volume result extrapolated to zero cell size (Richardson,
first order). Results are printed and written to build/compare_bvp.txt for docs/validation.md.

Usage:
    python scripts/compare_bvp.py
"""

from __future__ import annotations

import dataclasses
import functools
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import solver as sv
from shaftsim import transport as tr

ROOT = Path(__file__).resolve().parents[1]
W = np.array(pl.OXYGEN_WEIGHTS)


def bvp_rhs(z: np.ndarray, y: np.ndarray, prob: sh.ShaftProblem) -> np.ndarray:
    """dX_j/dz = r_j/u_s, df/dz = c sum_j w_j r_j / u_s (continuous form of the FV model)."""
    X = np.clip(y[:3], 0.0, 1.0)
    f = np.clip(y[3], 0.0, prob.f_in)
    r = pl.step_rates_array(prob.coeffs, X, f, prob.g_in - f, prob.rate_multiplier)
    r = r / prob.case.solid_velocity
    return np.vstack([r, prob.c_gas * (W @ r)])


def bvp_bc(ya: np.ndarray, yb: np.ndarray, prob: sh.ShaftProblem) -> np.ndarray:
    """X_j(0) = 0 (fresh solid at the top), f(H) = f_in (gas inlet at the bottom)."""
    return np.array([ya[0], ya[1], ya[2], yb[3] - prob.f_in])


def main() -> None:
    warnings.simplefilter("ignore")
    sp = pr.load_species(ROOT / "data" / "species.yaml")
    td = tr.load_transport(
        ROOT / "data" / "transport.yaml",
        ROOT / "data" / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    trip = pl.load_refit_triplets(ROOT / "data" / "kinetics_refit.yaml")
    lines = []
    cases = (
        ("700 degC, 3 m", 973.15, 3.0, 2800.0, 1.0),
        ("890 degC, 10.1 m", 1163.15, 10.1, 2800.0, 1.0),
        ("800 degC, gas-limited 1157 mol/s", 1073.15, 10.1, 1157.0, 1.0),
        ("800 degC, gas-limited 1639 mol/s", 1073.15, 10.1, 1639.0, 1.0),
        ("800 degC, 1157 mol/s, rates x100", 1073.15, 10.1, 1157.0, 100.0),
    )
    for label, T, H, F, mult in cases:
        case = dataclasses.replace(
            sh.load_case(ROOT / "data" / "furnace.yaml", ROOT / "data" / "kinetics.yaml", T=T),
            height=H,
            gas_flow=F,
        )
        co = sh.shaft_coefficients(case, sp, td, trip)
        ref = sh.ShaftProblem(case, co, 1, mult)
        rhs = functools.partial(bvp_rhs, prob=ref)
        bc = functools.partial(bvp_bc, prob=ref)

        # finite volumes on 200 / 400 / 800 cells, Richardson extrapolation (first order)
        fv = {}
        t_fv = 0.0
        y_prev = None
        for n in (200, 400, 800):
            # fast kinetics: continuation from the physical-rate solution on the same grid
            if mult != 1.0:
                y_prev, _ = sv.solve(sh.ShaftProblem(case, co, n))
                for m_step in (10.0, 100.0, 1000.0):
                    if m_step >= mult:
                        break
                    y_prev, _ = sv.solve(sh.ShaftProblem(case, co, n, m_step), y0=y_prev)
            prob = sh.ShaftProblem(case, co, n, mult)
            t0 = time.time()
            y, st = sv.solve(prob, y0=y_prev)
            t_fv = time.time() - t0
            if not st.converged:
                lines.append(f"  FV {n} cells did not converge ({st.method})")
            fv[n] = (prob, y, sh.make_solution(prob, y, td, st.residuals[-1], st.iterations))
        m = {n: s.metallization for n, (_, _, s) in fv.items()}
        m_extrap = 2 * m[800] - m[400]
        lines.append(
            f"{label}: FV metallisation 200/400/800 cells = "
            f"{m[200]:.5f} / {m[400]:.5f} / {m[800]:.5f}, Richardson {m_extrap:.5f} "
            f"(800 cells: {t_fv:.1f} s)"
        )

        prob, y, s = fv[400]
        z_fv = np.concatenate([[0.0], s.z + 0.5 * prob.dz])
        guesses = {
            "naive guess": np.vstack([np.zeros((3, 50)), np.full(50, ref.f_in)]),
            "FV solution as guess": None,
        }
        for gname, g in guesses.items():
            if g is None:
                Xg = np.concatenate([np.zeros((3, 1)), s.X_steps], axis=1)
                fg = np.concatenate([[s.x_H2[0]], s.x_H2])
                z0, y0 = z_fv, np.vstack([Xg, fg])
            else:
                z0, y0 = np.linspace(0.0, H, g.shape[1]), g
            t0 = time.time()
            res = solve_bvp(rhs, bc, z0, y0, tol=1e-6, max_nodes=200000)
            dt = time.time() - t0
            if res.success:
                metal = float(np.clip(res.sol(H)[2], 0, 1))
                lines.append(
                    f"  solve_bvp, {gname}: success, {res.x.size} nodes, {dt:.1f} s, "
                    f"metallisation {metal:.5f} (FV Richardson {m_extrap:.5f})"
                )
            else:
                lines.append(f"  solve_bvp, {gname}: FAILED ({res.message}), {dt:.1f} s")
    text = "\n".join(lines)
    print(text)
    out = ROOT / "build" / "compare_bvp.txt"
    out.parent.mkdir(exist_ok=True)
    out.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
