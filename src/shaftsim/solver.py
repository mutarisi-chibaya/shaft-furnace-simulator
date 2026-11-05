"""Steady-state solution of the discrete shaft model (shaftsim.shaft.ShaftProblem).

Strategy (three methods, tried in this order; the result says which one converged)
-----------------------------------------------------------------------------------
1. **Discrete shooting on the top gas** (primary). If the H2 fraction f_0 of the gas leaving
   the top cell is known, the discrete equations can be marched from the top: cell 0 receives
   fresh solid and contains gas f_0, so its solid equations give X_0; the gas balance then gives
   the gas entering it from below, f_1 = f_0 + c (Xbar_0 - 0); and so on down to f_n, which must
   equal the inlet. f_n(f_0) is increasing (more H2 everywhere means more reduction and more
   consumption), so a bracketing root finder on f_0 (Brent) always converges. Unlike
   shooting on the continuous ODEs (which PLAN.md warns about), every step solves the
   discrete cell equations exactly, so there is no error amplification; for physical
   kinetics the residual reaches ~1e-11.
   Limit: as kinetics become very fast, X(f) in a cell tends to a step function of f and
   f_n(f_0) to a jump; bisection then cannot locate f_0 in double precision.
2. **Symmetric block Gauss-Seidel sweeps with coupled cells** (polish / continuation): a
   downward then an upward pass, each cell solved exactly for its own solid conversions and
   outlet gas (:func:`coupled_cell`) given its upstream neighbours. Solving each phase for a
   frozen other phase instead lets a cell with fast kinetics consume more H2 than it receives,
   and the sweeps then oscillate. Robust, but sharp fronts move only about one cell per few
   sweeps, so it is used from a good starting point. Same idea as the cell-by-cell
   Gauss-Seidel of REDUCTOR (Marsigny et al. 2025).
3. **Pseudo-transient continuation** (last resort, Kelley & Keyes 1998): implicit Euler steps
   (I/dtau + J) dy = -R with a growing pseudo-time step (switched evolution relaxation), so
   the iteration turns into Newton's method. The rate law is not smooth (min() of the oxide
   availability constraint, unbounded derivative of (1-X)^(2/3) at X -> 1), so Newton stalls
   around 1e-7 near full conversion.

scripts/compare_bvp.py compares scipy.integrate.solve_bvp on the same problem.

With side inlets (bustle) and variable cross-section (shaftsim.shaft.Grid), the march uses
the per-cell residence time kappa_i, the per-cell c_i = N_O / F_i and the mixed inflow
f_mix,i = (F_{i+1} f_{i+1} + S_i f_s,i) / F_i; the monotonicity argument is unchanged.

The solid equations of a cell, g_j(X_j) = X_j - X_up,j - kappa_i r_j(X_j), are strictly
increasing in X_j (the rate decreases with conversion), so each has a unique root in
[X_up,j, 1] (Brent's method), solved in the order j = 1, 2, 3 because of the availability
constraint.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray
from scipy.optimize import brentq

from shaftsim import pellet as pl
from shaftsim.shaft import ShaftProblem

FloatArray = NDArray[np.float64]
_BAND = 12  # unknowns spanned by cells i-1..i+1 (4 unknowns per cell)


@dataclass(frozen=True)
class SolverStats:
    """Convergence record.

    Attributes:
        residuals: Max-norm of the scaled residual after each iteration.
        sweeps: Number of Gauss-Seidel sweeps.
        newton_steps: Number of accepted pseudo-transient (Newton) steps.
        converged: Whether the final residual is below the tolerance.
        method: Method that produced the result: "shooting", "sweeps",
            "grid continuation + sweeps" or "pseudo-transient".
    """

    residuals: tuple[float, ...]
    sweeps: int
    newton_steps: int
    converged: bool
    method: str = ""

    @property
    def iterations(self) -> int:
        """Total iterations (sweeps + Newton steps)."""
        return self.sweeps + self.newton_steps


def initial_guess(problem: ShaftProblem) -> FloatArray:
    """No reduction, inlet gas everywhere."""
    X = np.zeros((3, problem.n_cells))
    f = np.minimum(problem.f_in, problem.grid.g)
    return problem.pack(X, f)


def _solid_given_gas(
    problem: ShaftProblem, x_up: tuple[float, float, float], x_h2: float, x_h2o: float, i: int
) -> tuple[float, float, float]:
    """Solve the three solid equations of cell i for a given gas composition."""
    kappa = float(problem.grid.kappa[i])
    co, mult = problem.coeffs_of(i), problem.rate_multiplier
    X = list(x_up)
    for j in range(3):

        def g(v: float, j: int = j) -> float:
            trial = list(X)
            trial[j] = v
            r = pl.step_rates_scalar(co, trial[0], trial[1], trial[2], x_h2, x_h2o, mult)
            return v - x_up[j] - kappa * r[j]

        lo = x_up[j]
        if g(1.0) <= 0.0:
            X[j] = 1.0
        elif g(lo) < 0.0:
            X[j] = float(brentq(g, lo, 1.0, xtol=1e-15, rtol=1e-15))
    return X[0], X[1], X[2]


# --------------------------------------------------------------------------- 1. shooting


def march_from_top(problem: ShaftProblem, f0: float) -> tuple[FloatArray, FloatArray, float]:
    """March the discrete equations from the top for a given top-cell H2 fraction f0.

    Returns:
        (X shape (3, n), f shape (n,), f_n): conversions, cell H2 fractions, and the H2
        fraction implied at the bottom gas inlet face (to be compared with the true inlet
        value).
    """
    n = problem.n_cells
    gr = problem.grid
    c = problem.c_cell
    F_below = np.concatenate([gr.F[1:], [problem.case.gas_flow]])
    w1, w2, w3 = pl.OXYGEN_WEIGHTS
    X = np.zeros((3, n))
    f = np.zeros(n)
    x_up = (0.0, 0.0, 0.0)
    xbar_up = 0.0
    fi = f0
    for i in range(n):
        gi = float(gr.g[i])
        fi = min(max(fi, 0.0), gi)
        Xi = _solid_given_gas(problem, x_up, fi, gi - fi, i)
        X[:, i] = Xi
        f[i] = fi
        xbar = w1 * Xi[0] + w2 * Xi[1] + w3 * Xi[2]
        f_mix = fi + float(c[i]) * (xbar - xbar_up)
        fi = (float(gr.F[i]) * f_mix - float(gr.S[i] * gr.f_side[i])) / float(F_below[i])
        xbar_up, x_up = xbar, Xi
    return X, f, fi


def shoot(problem: ShaftProblem) -> FloatArray:
    """Discrete shooting: root of the monotone miss f_n(f0) - f_in on the top-cell H2 fraction
    f0 (Brent's method within the bracket [0, f_max]; see module docstring)."""
    lo, hi = 0.0, problem.f_max

    def miss(f0: float) -> float:
        return float(march_from_top(problem, f0)[2] - problem.f_in)

    m_lo, m_hi = miss(lo), miss(hi)
    if m_lo < 0.0 < m_hi:
        # Brent's method on the monotone miss function: ~10-20 marches instead of 60
        # bisections; it falls back to bisection steps where the function is step-like
        # (very fast kinetics), so the bracket guarantee is kept.
        f0 = float(brentq(miss, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=200))
    else:
        f0 = lo if abs(m_lo) <= abs(m_hi) else hi
    X, f, _ = march_from_top(problem, f0)
    return problem.pack(X, f)


# --------------------------------------------------------------------------- 2. sweeps


def coupled_cell(
    problem: ShaftProblem, x_up: tuple[float, float, float], f_down: float, i: int
) -> tuple[tuple[float, float, float], float]:
    """Solve one cell for its solid conversions and outlet H2 fraction together.

    Given the solid entering from above (x_up) and the H2 fraction of the gas entering from
    below (f_down, already mixed with any side gas: ShaftProblem.mixed_inflow), find f with
    h(f) = f - f_down + c_i (Xbar(X(f)) - Xbar_up) = 0, where X(f) solves the solid equations
    for gas (f, g_i - f). More H2 gives more reduction, so h is increasing; h(0) = -f_down
    <= 0 and h(f_down) >= 0 bracket the unique root, so a cell can never consume more H2
    than it receives.
    """
    w1, w2, w3 = pl.OXYGEN_WEIGHTS
    xbar_up = w1 * x_up[0] + w2 * x_up[1] + w3 * x_up[2]
    c, g_in = float(problem.c_cell[i]), float(problem.grid.g[i])

    def solid(f: float) -> tuple[float, float, float]:
        return _solid_given_gas(problem, x_up, f, g_in - f, i)

    def h(f: float) -> float:
        X = solid(f)
        return f - f_down + c * (w1 * X[0] + w2 * X[1] + w3 * X[2] - xbar_up)

    if f_down <= 0.0 or h(f_down) <= 0.0:
        return solid(max(f_down, 0.0)), max(f_down, 0.0)
    f = float(brentq(h, 0.0, f_down, xtol=1e-15, rtol=1e-15))
    return solid(f), f


def gauss_seidel_sweeps(problem: ShaftProblem, y: FloatArray, n_sweeps: int) -> FloatArray:
    """Symmetric block Gauss-Seidel with coupled cells (downward pass, then upward pass)."""
    X, f = problem.unpack(y)
    X, f = X.copy(), f.copy()
    n = problem.n_cells
    gr = problem.grid
    F_below = np.concatenate([gr.F[1:], [problem.case.gas_flow]])
    side = gr.S * gr.f_side
    for _ in range(n_sweeps):
        for order in (range(n), range(n - 1, -1, -1)):
            for i in order:
                x_up = (0.0, 0.0, 0.0) if i == 0 else (X[0, i - 1], X[1, i - 1], X[2, i - 1])
                f_below = problem.f_in if i == n - 1 else float(f[i + 1])
                f_down = (float(F_below[i]) * f_below + float(side[i])) / float(gr.F[i])
                X[:, i], f[i] = coupled_cell(problem, x_up, f_down, i)
    return problem.pack(X, f)


# --------------------------------------------------------------------------- 3. Newton


def banded_jacobian(problem: ShaftProblem, y: FloatArray, r0: FloatArray) -> sp.csc_matrix:
    """Finite-difference Jacobian dR/dy with grouped column perturbations (sparse)."""
    n = y.size
    rows, cols, vals = [], [], []
    idx = np.arange(n)
    upper = np.where(idx % 4 == 3, problem.f_max, 1.0)
    for g in range(_BAND):
        group = idx[g::_BAND]
        h = 1e-7 * np.maximum(1.0, np.abs(y[group]))
        sign = np.where(y[group] + h > upper[group], -1.0, 1.0)  # stay inside the bounds
        yp = y.copy()
        yp[group] += sign * h
        dr = problem.residual(yp) - r0
        for c, hc, s in zip(group, h, sign, strict=True):
            cell = c // 4
            lo, hi = max(0, (cell - 1) * 4), min(n, (cell + 2) * 4)
            rows.append(np.arange(lo, hi))
            cols.append(np.full(hi - lo, c))
            vals.append(dr[lo:hi] / (s * hc))
    return sp.csc_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
    )


def _project(problem: ShaftProblem, y: FloatArray) -> FloatArray:
    Y = y.reshape(problem.n_cells, 4).copy()
    Y[:, :3] = np.clip(Y[:, :3], 0.0, 1.0)
    Y[:, 3] = np.clip(Y[:, 3], 0.0, problem.f_max)
    return Y.ravel()


def pseudo_transient(
    problem: ShaftProblem,
    y: FloatArray,
    tol: float,
    max_iter: int,
    dtau0: float = 1.0,
    dtau_max: float = 1e12,
) -> tuple[FloatArray, list[float]]:
    """Pseudo-transient continuation with backtracking; returns the iterate and residuals."""
    r = problem.residual(y)
    norm = float(np.max(np.abs(r)))
    history: list[float] = []
    dtau = dtau0
    eye = sp.identity(y.size, format="csc")
    for _ in range(max_iter):
        if norm < tol:
            break
        J = banded_jacobian(problem, y, r)
        dy = spla.spsolve((eye / dtau + J).tocsc(), -r)
        lam, accepted = 1.0, False
        for _ls in range(8):
            y_new = _project(problem, y + lam * dy)
            r_new = problem.residual(y_new)
            if np.all(np.isfinite(r_new)) and np.linalg.norm(r_new) < np.linalg.norm(r):
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            dtau /= 10.0
            if dtau < 1e-6:
                break
            continue
        norm_new = float(np.max(np.abs(r_new)))
        dtau = min(dtau * norm / max(norm_new, 1e-300), dtau_max)
        y, r, norm = y_new, r_new, norm_new
        history.append(norm)
    return y, history


# --------------------------------------------------------------------------- driver


def coarsen(problem: ShaftProblem) -> ShaftProblem:
    """The same problem on half the cells (for grid continuation).

    Per-cell coefficients (non-isothermal problems) are taken from the fine cell nearest to
    each coarse cell centre: the coarse solution only serves as a starting point.
    """
    coarse = replace(problem, n_cells=problem.n_cells // 2, cell_coeffs=None)
    if problem.cell_coeffs is None:
        return coarse
    z_fine = problem.grid.z
    idx = [int(np.argmin(np.abs(z_fine - zc))) for zc in coarse.grid.z]
    cc = tuple(problem.cell_coeffs[i] for i in idx)
    return replace(coarse, coeffs=cc[0], cell_coeffs=cc)


def prolong(coarse: ShaftProblem, y_coarse: FloatArray, fine: ShaftProblem) -> FloatArray:
    """Interpolate a coarse-grid solution onto a finer grid (linear, at cell centres)."""
    Xc, fc = coarse.unpack(y_coarse)
    zc = coarse.grid.z
    zf = fine.grid.z
    Xf = np.vstack([np.interp(zf, zc, Xc[j]) for j in range(3)])
    ff = np.interp(zf, zc, fc)
    return fine.pack(Xf, ff)


def close_gas_balance(problem: ShaftProblem, y: FloatArray) -> FloatArray:
    """Recompute the gas profile from the solid profile by the exact upward march.

    After this, the discrete gas equations hold to round-off, so the O and H element balances
    close to machine precision whatever the remaining solid residual.
    """
    X, _ = problem.unpack(y)
    xbar = np.asarray(pl.overall_conversion(X))
    dxbar = np.diff(np.concatenate([[0.0], xbar]))
    gr = problem.grid
    n_o = problem.case.removable_oxygen_flow
    f = np.empty(problem.n_cells)
    h2_below = problem.case.gas_flow * problem.f_in  # H2 flow entering from below [mol/s]
    for i in range(problem.n_cells - 1, -1, -1):
        h2 = h2_below + gr.S[i] * gr.f_side[i] - n_o * dxbar[i]
        f[i] = h2 / gr.F[i]
        h2_below = h2
    return problem.pack(X, f)


def _finish(
    problem: ShaftProblem,
    y: FloatArray,
    history: list[float],
    sweeps: int,
    newton: int,
    tol: float,
    method: str,
) -> tuple[FloatArray, SolverStats]:
    y = close_gas_balance(problem, y)
    norm = float(np.max(np.abs(problem.residual(y))))
    history.append(norm)
    return y, SolverStats(tuple(history), sweeps, newton, norm < tol, method)


def solve(
    problem: ShaftProblem,
    y0: FloatArray | None = None,
    tol: float = 1e-9,
    max_sweeps: int = 300,
    max_newton: int = 100,
) -> tuple[FloatArray, SolverStats]:
    """Solve the discrete steady state R(y) = 0.

    Args:
        problem: Discrete shaft problem.
        y0: Starting point. If None, discrete shooting is tried first; if it fails on a fine
            grid, the problem is solved on half the cells (recursively) and interpolated. If
            given (e.g. the solution of a neighbouring case, for continuation), the sweeps
            start from it.
        tol: Threshold on the max-norm of the scaled residual [-] (conversion units).
        max_sweeps: Maximum Gauss-Seidel sweeps.
        max_newton: Maximum pseudo-transient steps if the sweeps stagnate.

    Returns:
        The solution vector and convergence statistics. The caller must check
        ``stats.converged``; an unconverged result is returned, not hidden.
    """
    history: list[float] = []
    method = "sweeps"
    if y0 is None:
        y = shoot(problem)
        norm = float(np.max(np.abs(problem.residual(y))))
        history.append(norm)
        if norm < tol:
            return _finish(problem, y, history, 0, 0, tol, "shooting")
        if norm >= 0.05 and problem.n_cells >= 50:
            # grid continuation: solve on half the cells, interpolate, then polish
            coarse = coarsen(problem)
            y_c, st_c = solve(coarse, None, tol, max_sweeps, max_newton)
            if st_c.converged:
                y, method = prolong(coarse, y_c, problem), "grid continuation + sweeps"
            else:
                y = initial_guess(problem)
        elif norm >= 0.05:
            y = initial_guess(problem)
    else:
        y = y0.copy()
    best, since_best, sweeps = np.inf, 0, 0
    for sweeps in range(1, max_sweeps + 1):
        y = gauss_seidel_sweeps(problem, y, 1)
        norm = float(np.max(np.abs(problem.residual(y))))
        history.append(norm)
        if norm < 0.1 * tol:  # margin: closing the gas balance perturbs the solid residual
            return _finish(problem, y, history, sweeps, 0, tol, method)
        if norm < 0.9 * best:
            best, since_best = norm, 0
        else:
            since_best += 1
            if since_best >= 10:  # no 10 % improvement in 10 sweeps: stagnated
                break
    y, newton_hist = pseudo_transient(problem, y, 0.1 * tol, max_newton)
    history += newton_hist
    return _finish(problem, y, history, sweeps, len(newton_hist), tol, "pseudo-transient")
