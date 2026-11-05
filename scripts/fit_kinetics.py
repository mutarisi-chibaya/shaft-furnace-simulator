"""Re-fit the pellet-model kinetic triplets (k0_j, Ea_j) to the Marsigny et al. (2025) TGA data.

Why: the published triplets (Marsigny 2025, Table 2) were fitted with a model whose structural
and solid-state-diffusion parameters are not published. With the documented structure used in
shaftsim (data/kinetics.yaml), the published triplets reproduce 900 degC reasonably but are
far too slow at 600-700 degC (see docs/validation.md). We therefore repeat the authors'
procedure (Sec. 3.1: fit the triplets to the TG96 curves at 600-900 degC, X <= 0.99) with our
model. Nothing else is adjusted.

Data: tests/data/pellet_reduction_marsigny2025_fig6.csv (experimental curves, digitised).
Objective: least squares on X_model - X_exp, each curve resampled at 40 conversion levels
(X = 0.02 ... 0.98) so that every temperature has equal weight.
Parameters: per step, ln k_j(T_ref) and Ea_j, T_ref = 1023.15 K (decorrelates k0 and Ea).

Also runs a leave-one-temperature-out (LOTO) cross-validation: each temperature is predicted
by a fit to the other six. The seven LOTO parameter vectors also give a jackknife covariance
of the parameters (curve-level resampling; used for uncertainty quantification, Milestone 6,
because the Jacobian covariance treats the 40 points of a curve as independent and so
underestimates the spread). Results -> data/kinetics_refit.yaml (generated).

Usage:
    python scripts/fit_kinetics.py [--no-loto]
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr

ROOT = Path(__file__).resolve().parents[1]
Curves = dict[int, tuple[np.ndarray, np.ndarray]]
T_REF = 1023.15
X_LEVELS = np.linspace(0.02, 0.98, 40)
T_END = 150 * 60.0


def load() -> tuple[pl.PelletParameters, pl.TGAConditions, dict[str, pr.Species], tr.TransportData]:
    sp = pr.load_species(ROOT / "data" / "species.yaml")
    td = tr.load_transport(
        ROOT / "data" / "transport.yaml",
        ROOT / "data" / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    pp, cond = pl.load_kinetics(ROOT / "data" / "kinetics.yaml")
    return pp, cond, sp, td


def experimental_curves() -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """{T_degC: (X_levels, t_exp at those levels [s])} from the digitised curves."""
    d = pd.read_csv(ROOT / "tests" / "data" / "pellet_reduction_marsigny2025_fig6.csv")
    out = {}
    for T, g in d.groupby("temperature_C"):
        g = g.sort_values("time_min")
        x = np.maximum.accumulate(g.conversion.to_numpy())
        t = g.time_min.to_numpy() * 60.0
        x_pts = np.concatenate([[0.0], x])
        t_pts = np.concatenate([[0.0], t])
        keep = X_LEVELS <= x.max()
        out[int(str(T))] = (X_LEVELS[keep], np.interp(X_LEVELS[keep], x_pts, t_pts))
    return out


def triplets_from(theta: np.ndarray) -> pl.Triplets:
    trip = []
    for j in range(3):
        ln_kref, ea = theta[2 * j], theta[2 * j + 1] * 1e4
        trip.append(pl.KineticTriplet(math.exp(ln_kref) * math.exp(ea / (pr.R_GAS * T_REF)), ea))
    return (trip[0], trip[1], trip[2])


def theta_from(trip: pl.Triplets) -> np.ndarray:
    th = []
    for t in trip:
        th += [math.log(t.k(T_REF)), t.Ea / 1e4]
    return np.array(th)


def residuals(
    theta: np.ndarray,
    temps: list[int],
    curves: Curves,
    pp: pl.PelletParameters,
    cond: pl.TGAConditions,
    sp: dict[str, pr.Species],
    td: tr.TransportData,
) -> np.ndarray:
    trip = triplets_from(theta)
    res = []
    for T in temps:
        x_lv, t_exp = curves[T]
        r = pl.simulate_tga(pp, cond, T + 273.15, T_END, sp, td, triplets=trip, n_out=600)
        x_model = np.interp(t_exp, r.t, r.X)
        res.append(x_model - x_lv)
    return np.concatenate(res)


def fit(
    temps: list[int],
    curves: Curves,
    theta0: np.ndarray,
    pp: pl.PelletParameters,
    cond: pl.TGAConditions,
    sp: dict[str, pr.Species],
    td: tr.TransportData,
) -> Any:
    lo = np.array([-30, 0.5, -30, 0.5, -30, 0.5])
    hi = np.array([5, 30.0, 5, 30.0, 5, 30.0])
    return least_squares(
        residuals,
        theta0,
        args=(temps, curves, pp, cond, sp, td),
        bounds=(lo, hi),
        x_scale=np.array([1, 1, 1, 1, 1, 1.0]),
        diff_step=1e-3,
        xtol=1e-6,
        ftol=1e-8,
    )


def per_curve_rmse(
    theta: np.ndarray,
    temps: list[int],
    curves: Curves,
    pp: pl.PelletParameters,
    cond: pl.TGAConditions,
    sp: dict[str, pr.Species],
    td: tr.TransportData,
) -> dict[int, float]:
    res = residuals(theta, temps, curves, pp, cond, sp, td)
    out, i = {}, 0
    for T in temps:
        n = len(curves[T][0])
        out[T] = float(np.sqrt(np.mean(res[i : i + n] ** 2)))
        i += n
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-loto", action="store_true")
    args = ap.parse_args()
    pp, cond, sp, td = load()
    curves = experimental_curves()
    temps = sorted(curves)
    theta_pub = theta_from(pp.triplets)
    rmse_pub = per_curve_rmse(theta_pub, temps, curves, pp, cond, sp, td)
    print("published triplets, RMSE per T:", {k: round(v, 4) for k, v in rmse_pub.items()})

    t0 = time.time()
    sol = fit(temps, curves, theta_pub, pp, cond, sp, td)
    print(f"full fit: cost {sol.cost:.4g}, {sol.nfev} evaluations, {time.time() - t0:.0f} s")
    trip = triplets_from(sol.x)
    rmse_fit = per_curve_rmse(sol.x, temps, curves, pp, cond, sp, td)
    dof = max(len(sol.fun) - len(sol.x), 1)
    s2 = float(np.sum(sol.fun**2) / dof)
    cov = np.linalg.pinv(sol.jac.T @ sol.jac) * s2
    se = np.sqrt(np.diag(cov))

    loto: dict[int, float] = {}
    loto_theta: dict[int, np.ndarray] = {}
    if not args.no_loto:
        for T in temps:
            others = [u for u in temps if u != T]
            s = fit(others, curves, sol.x, pp, cond, sp, td)
            loto[T] = per_curve_rmse(s.x, [T], curves, pp, cond, sp, td)[T]
            loto_theta[T] = s.x
            print(f"LOTO {T} degC: predictive RMSE {loto[T]:.4f}")

    lines = [
        "# GENERATED by scripts/fit_kinetics.py on "
        + dt.date.today().isoformat()
        + " - do not edit by hand.",
        "# Kinetic triplets re-fitted to the Marsigny et al. (2025) TG96 curves (Fig. 6,",
        "# 600-900 degC, digitised: tests/data/pellet_reduction_marsigny2025_fig6.csv) with the",
        "# model and structure of shaftsim.pellet / data/kinetics.yaml. See docs/validation.md.",
        f"# Parameterisation: ln k_j(T_ref), Ea_j, T_ref = {T_REF} K.",
        "kinetic_triplets:",
    ]
    for j, t in enumerate(trip, start=1):
        lines += [
            f"  step_{j}:",
            f"    k0: {t.k0:.6e}  # m/s",
            f"    Ea: {t.Ea:.1f}  # J/mol",
            f"    ln_k_Tref: {sol.x[2 * j - 2]:.6f}",
            f"    se_ln_k_Tref: {se[2 * j - 2]:.6f}",
            f"    se_Ea: {se[2 * j - 1] * 1e4:.1f}  # J/mol, from the Jacobian",
        ]
    lines.append("correlation_matrix:  # order ln k1, Ea1, ln k2, Ea2, ln k3, Ea3")
    corr = cov / np.outer(se, se)
    for row in corr:
        lines.append("  - [" + ", ".join(f"{v:.3f}" for v in row) + "]")
    lines.append("rmse_fit:  # per temperature (degC), conversion units")
    lines += [f"  {T}: {v:.4f}" for T, v in rmse_fit.items()]
    lines.append("rmse_published_triplets:")
    lines += [f"  {T}: {v:.4f}" for T, v in rmse_pub.items()]
    if loto:
        lines.append("rmse_leave_one_temperature_out:")
        lines += [f"  {T}: {v:.4f}" for T, v in loto.items()]
        th = np.array([loto_theta[T] for T in temps])  # (n, 6), Ea in 1e4 J/mol
        n_j = th.shape[0]
        dev = th - th.mean(axis=0)
        cov_jk = (n_j - 1) / n_j * dev.T @ dev
        scale = np.array([1.0, 1e4, 1.0, 1e4, 1.0, 1e4])
        lines.append("loto_parameters:  # ln k_j(T_ref), Ea_j [J/mol] of each LOTO fit")
        for T in temps:
            v = loto_theta[T] * scale
            lines.append(f"  {T}: [" + ", ".join(f"{x:.6g}" for x in v) + "]")
        lines.append(
            "jackknife_covariance:  # (n-1)/n sum (theta_i - mean)(theta_i - mean)^T over the"
        )
        lines.append("  # LOTO fits; order ln k1, Ea1, ln k2, Ea2, ln k3, Ea3; Ea in J/mol")
        for row in cov_jk * np.outer(scale, scale):
            lines.append("  - [" + ", ".join(f"{v:.6g}" for v in row) + "]")
        lines.append("jacobian_covariance:  # same order and units, from the Jacobian at the fit")
        for row in cov * np.outer(scale, scale):
            lines.append("  - [" + ", ".join(f"{v:.6g}" for v in row) + "]")
    out = ROOT / "data" / "kinetics_refit.yaml"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}")
    for j, t in enumerate(trip, 1):
        print(
            f"step {j}: k0 {t.k0:.4g} m/s, Ea {t.Ea / 1e3:.1f} kJ/mol "
            f"(published {pp.triplets[j - 1].k0:.4g}, {pp.triplets[j - 1].Ea / 1e3:.1f})"
        )


if __name__ == "__main__":
    main()
