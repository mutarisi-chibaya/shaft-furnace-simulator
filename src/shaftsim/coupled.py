"""Fully coupled steady-state solver: species and energy in one Newton system (Milestone 8).

The alternating scheme of :func:`shaftsim.energy.solve_thermal` (species at frozen T, energy
at frozen species, Anderson-accelerated) cycles where a cell flips between "wustite left"
and "wustite exhausted" as its temperature moves (docs/validation.md, Milestone 7). Here all
unknowns of a cell, (X1, X2, X3, f, T_s, T_g), are solved together.

Discrete equations: exactly those of the alternating scheme, i.e. the species residual of
:class:`shaftsim.shaft.ShaftProblem` with the per-cell coefficients at the current
temperatures, and the energy residual of :class:`shaftsim.energy.ThermalProblem`. Scaling:
the gas-species row is negated (so every row has a positive diagonal) and the energy rows
are divided by the heat-capacity flow of the phase in the cell, so they read in kelvin.

Numerical choices:
- Unknowns: the conversions are solved for as unreacted-core radii u_j = (1 - X_j)^(1/3)
  (the shrinking-core variable of the grainy pellet model). In X the rate falls like
  (1 - X)^(2/3) towards complete conversion, so its derivative is infinite at X = 1 and
  Newton's method overshoots and stalls where a step completes; in u the rate behaves like
  u^2 and is smooth. The equations themselves are unchanged.
- Jacobian by finite differences with a colouring: a cell's residual depends only on its own
  unknowns and those of the cells above (solid upwind) and below (gas upwind), so perturbing
  one unknown in every third cell at once gives the block-tridiagonal Jacobian in 18 residual
  evaluations.
- Pseudo-transient continuation: solve (J + D / dtau) dz = -R with D = diag J (floored in
  magnitude), so each equation relaxes on its own scale and in its own direction; dtau
  grows as the residual falls (switched evolution relaxation, at least doubling after each
  decrease) and Newton's method is recovered as dtau -> infinity. A step that increases the
  scaled residual tenfold or gives non-finite values is retried with a smaller dtau.
- After each step the unknowns are projected onto their physical bounds: u >= U_MIN
  (complete conversion), 0 <= f <= g (no negative H2O), and the temperature range of the
  property tables.
- Active set (projected Newton): where a step is complete, u sits at U_MIN and its equation
  depends on u only through u^3 (sensitivity ~1e-16), which makes the Newton matrix
  singular. A core radius at its bound whose equation is satisfied or pushes it further down
  is frozen for that step (its row and column are replaced by the identity); it is released
  as soon as its equation asks it to grow. All rows, frozen or not, enter the convergence
  test.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

from . import energy as en

FloatArray = NDArray[np.float64]

N_VAR = 6  #: unknowns per cell: X1, X2, X3, f, T_s, T_g
SPECIES_ROWS = slice(0, 4)
ENERGY_ROWS = slice(4, 6)


def pack(X: FloatArray, f: FloatArray, Ts: FloatArray, Tg: FloatArray) -> FloatArray:
    """(X shape (3, n), f, T_s, T_g) -> unknown vector, ordered per cell."""
    return np.column_stack([X.T, f, Ts, Tg]).ravel()


def unpack(y: FloatArray, n: int) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Inverse of :func:`pack`."""
    Y = y.reshape(n, N_VAR)
    return Y[:, :3].T, Y[:, 3], Y[:, 4], Y[:, 5]


def to_internal(y: FloatArray, n: int) -> FloatArray:
    """Unknowns -> solver variables: conversions replaced by u = (1 - X)^(1/3)."""
    Z = y.reshape(n, N_VAR).copy()
    Z[:, :3] = np.cbrt(np.clip(1.0 - Z[:, :3], 0.0, 1.0))
    return Z.ravel()


def from_internal(z: FloatArray, n: int) -> FloatArray:
    """Inverse of :func:`to_internal`: X = 1 - u^3."""
    Y = z.reshape(n, N_VAR).copy()
    Y[:, :3] = 1.0 - Y[:, :3] ** 3
    return Y.ravel()


def energy_scales(prob: en.ThermalProblem, y: FloatArray) -> FloatArray:
    """Heat-capacity flows [W/K] of solid and gas per cell, shape (n, 2), used to express the
    energy residual in kelvin. Fixed for a solve (taken at its starting point)."""
    X, f, Ts, Tg = unpack(y, prob.n_cells)
    gr = prob.grid
    _, Cs = prob.solid_enthalpy(X, Ts)
    _, Cg = prob.gas_enthalpy(f, Tg, gr.F, gr.g)
    return np.column_stack([np.maximum(Cs, 1e-12), np.maximum(Cg, 1e-12)])


def residual(
    prob: en.ThermalProblem,
    y: FloatArray,
    scales: FloatArray,
    one_minus: FloatArray | None = None,
) -> FloatArray:
    """Scaled coupled residual, ordered per cell like the unknowns.

    Species rows: dimensionless (conversion units). Energy rows: kelvin. ``one_minus``
    optionally gives 1 - X at full precision (shape (3, n)).
    """
    n = prob.n_cells
    X, f, Ts, Tg = unpack(y, n)
    sprob = prob.species_problem(Ts, Tg)
    rs = sprob.residual(sprob.pack(X, f), one_minus=one_minus).reshape(n, 4).copy()
    rs[:, 3] *= -1.0  # gas row: positive diagonal
    UA = prob.heat_transfer_UA(Tg)
    re = prob.energy_residual(X, f, Ts, Tg, UA).reshape(n, 2) / scales
    return np.column_stack([rs, re]).ravel()


def residual_internal(prob: en.ThermalProblem, z: FloatArray, scales: FloatArray) -> FloatArray:
    """:func:`residual` as a function of the solver variables (u instead of X); the rates
    use 1 - X = u^3 directly, so no precision is lost near complete conversion."""
    n = prob.n_cells
    u = z.reshape(n, N_VAR)[:, :3].T
    return residual(prob, from_internal(z, n), scales, one_minus=u**3)


#: Smallest core radius: X <= 1 - 1e-24 (1 in floating point; the rate there is ~1e-16 of
#: its value at the start of the step). Keeps the u-equations non-singular where a step is
#: complete upstream.
U_MIN = 1e-8

#: Finite-difference steps per solver variable (core radius / mole fraction units, kelvin).
#: The core-radius step must be small against the smoothing width of the availability
#: constraint (near X = 0 a step h in u moves X by 3h): with 1e-7 the Jacobian was ~10 %
#: wrong inside a 1e-5 transition and Newton converged only linearly, one cell per
#: iteration; 1e-9 restores quadratic convergence (rounding error ~1e-7 relative).
FD_STEP = np.array([1e-9, 1e-9, 1e-9, 1e-8, 1e-4, 1e-4])


def jacobian(
    prob: en.ThermalProblem, z: FloatArray, scales: FloatArray, r0: FloatArray
) -> sp.csc_matrix:
    """Banded Jacobian of :func:`residual_internal` by coloured finite differences.

    A cell's residual depends on the cells within ``bw`` of it: bw = 1 for the first-order
    scheme (block-tridiagonal, 3 colours), bw = 2 for the second-order scheme (the face
    reconstruction reaches one cell further; 5 colours).
    """
    n = prob.n_cells
    bw = 1 if prob.scheme == 1 else 2
    ncol = 2 * bw + 1
    rows: list[NDArray[np.int_]] = []
    cols: list[NDArray[np.int_]] = []
    vals: list[FloatArray] = []
    for colour in range(ncol):
        cells = np.arange(colour, n, ncol)
        for v in range(N_VAR):
            h = FD_STEP[v]
            # one-sided step pointing into the admissible range
            idx = cells * N_VAR + v
            zp = z.copy()
            sign = np.where((v < 3) & (z[idx] + h > 1.0), -1.0, 1.0)
            zp[idx] += sign * h
            dr = (residual_internal(prob, zp, scales) - r0).reshape(n, N_VAR)
            for c, s in zip(cells, sign, strict=True):
                for nb in range(c - bw, c + bw + 1):
                    if 0 <= nb < n:
                        rows.append(np.arange(nb * N_VAR, (nb + 1) * N_VAR))
                        cols.append(np.full(N_VAR, c * N_VAR + v))
                        vals.append(dr[nb] / (s * h))
    return sp.csc_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n * N_VAR, n * N_VAR),
    )


def project(prob: en.ThermalProblem, z: FloatArray) -> FloatArray:
    """Clip the solver variables to their physical bounds (see the module docstring)."""
    n = prob.n_cells
    Y = z.reshape(n, N_VAR).copy()
    Y[:, :3] = np.clip(Y[:, :3], U_MIN, 1.0)
    Y[:, 3] = np.clip(Y[:, 3], 0.0, prob.grid.g)
    lo, hi = prob.T_TABLE
    Y[:, 4:] = np.clip(Y[:, 4:], lo, hi)
    return Y.ravel()


@dataclass(frozen=True)
class CoupledStats:
    """Outcome of :func:`solve_coupled`.

    Attributes:
        converged: Whether both residual blocks are below their tolerances.
        iterations: Accepted pseudo-time steps.
        rejected: Rejected (retried) steps.
        species_residual: Final max-norm of the species rows [-].
        energy_residual: Final max-norm of the energy rows [K].
        history: Scaled max-norm after each accepted step (max of block norm / tolerance).
        seconds: Wall time [s].
    """

    converged: bool
    iterations: int
    rejected: int
    species_residual: float
    energy_residual: float
    history: list[float] = field(default_factory=list)
    seconds: float = 0.0


def _norms(r: FloatArray, n: int) -> tuple[float, float]:
    R = r.reshape(n, N_VAR)
    return float(np.max(np.abs(R[:, SPECIES_ROWS]))), float(np.max(np.abs(R[:, ENERGY_ROWS])))


def _freeze_bound(
    J: sp.csc_matrix, z: FloatArray, r: FloatArray, n: int, tol: float
) -> tuple[sp.csc_matrix, FloatArray]:
    """Remove core radii at their lower bound from the Newton system (see module docstring).

    Frozen: u <= 2 U_MIN and residual <= tol. For a u-row the residual is X - X_up - kappa r;
    a positive value means X is too large, i.e. u must grow, so such a variable stays free.
    """
    Z = z.reshape(n, N_VAR)
    R = r.reshape(n, N_VAR)
    m = np.zeros((n, N_VAR), dtype=bool)
    m[:, :3] = (Z[:, :3] <= 2.0 * U_MIN) & (R[:, :3] <= tol)
    mask = m.ravel()
    if not mask.any():
        return J, r
    keep = sp.diags((~mask).astype(float))
    J = (keep @ J @ keep + sp.diags(mask.astype(float))).tocsc()
    return J, np.where(mask, 0.0, r)


def _line_search(
    prob: en.ThermalProblem,
    z: FloatArray,
    dz: FloatArray,
    scales: FloatArray,
    norm: float,
    scaled: Callable[[FloatArray], float],
    halvings: int,
) -> tuple[FloatArray, FloatArray, float]:
    """Backtracking along dz: the full step if it lowers the scaled residual, otherwise the
    best of up to ``halvings`` halved steps (projected onto the bounds). Returns
    (z_new, r_new, scaled norm); the caller decides whether to accept."""
    best: tuple[FloatArray, FloatArray, float] | None = None
    alpha = 1.0
    for _ in range(halvings + 1):
        z_new = project(prob, z + alpha * dz)
        r_new = residual_internal(prob, z_new, scales)
        norm_new = scaled(r_new) if np.all(np.isfinite(r_new)) else np.inf
        if best is None or norm_new < best[2]:
            best = (z_new, r_new, norm_new)
        if norm_new < norm:
            break
        alpha *= 0.5
    assert best is not None
    return best


def _release_bound(z: FloatArray, r: FloatArray, n: int, tol: float) -> tuple[FloatArray, bool]:
    """Lift core radii off their lower bound where their equation asks them to grow.

    A u-row reads R = (1 - u^3) - X_up - kappa r(u) with r ~ u^2, so dR/du ~ -3u^2 - 2 kappa k u
    vanishes at u = U_MIN: once a Newton step has overshot to complete conversion, the
    linearisation cannot bring u back even though R > 0 says it must grow. Such a u is reset
    to cbrt(u^3 + R), the root of the equation without its reaction term, which bounds the
    true root from above and lies where the slope is finite. Returns (z, changed)."""
    Z = z.reshape(n, N_VAR)[:, :3]
    R = r.reshape(n, N_VAR)[:, :3]
    m = (Z <= 2.0 * U_MIN) & (R > tol)
    if not m.any():
        return z, False
    z = z.copy()
    Zn = z.reshape(n, N_VAR)
    Zn[:, :3] = np.where(m, np.cbrt(np.clip(Z**3 + R, 0.0, 1.0)), Z)
    return z, True


def solve_coupled(
    prob: en.ThermalProblem,
    y0: FloatArray,
    species_tol: float = 1e-10,
    energy_tol: float = 1e-6,
    max_iter: int = 300,
    dtau0: float = 1.0,
    smoothing: tuple[float, float] = (0.0, 0.0),
    line_search: int = 0,
    max_retries: int = 12,
    callback: Callable[[int, FloatArray, FloatArray, float], None] | None = None,
) -> tuple[FloatArray, CoupledStats]:
    """Solve the coupled steady state from the starting point y0 (see the module docstring).

    Args:
        prob: Thermal problem.
        y0: Starting unknowns (:func:`pack`), e.g. from the alternating solver or a coarser
            grid.
        species_tol: Max-norm tolerance of the species rows [-].
        energy_tol: Max-norm tolerance of the energy rows [K].
        max_iter: Maximum accepted pseudo-time steps.
        dtau0: Initial pseudo-time step (relative to each equation's own scale).
        smoothing: (availability_width, min_width) of the availability constraint (see
            :func:`shaftsim.pellet.step_rates_array`); (0, 0) = the original discontinuous
            rule. Overrides the widths of ``prob``.
        line_search: Number of step halvings tried along each pseudo-time Newton direction
            before the pseudo-time step itself is reduced (0 = full steps only).
        max_retries: Pseudo-time step reductions per iteration before giving up.
        callback: Called after each accepted step as callback(iteration, z, r, dtau), with
            the solver variables z, the scaled residual r and the next pseudo-time step
            (diagnostics only).

    Returns:
        (y, stats); ``stats.converged`` must be checked by the caller.
    """
    t0 = time.perf_counter()
    prob = prob.with_parameters(
        availability_width=smoothing[0],
        min_width=smoothing[1],
    )
    n = prob.n_cells
    z = project(prob, to_internal(y0, n))
    scales = energy_scales(prob, from_internal(z, n))
    r = residual_internal(prob, z, scales)

    def scaled(rr: FloatArray) -> float:
        s, e = _norms(rr, n)
        return max(s / species_tol, e / energy_tol)

    norm = scaled(r)
    best = (norm, z, r)
    dtau = dtau0
    history = [norm]
    rejected = 0
    it = 0
    while norm > 1.0 and it < max_iter:
        z_rel, changed = _release_bound(z, r, n, species_tol)
        if changed:
            z = project(prob, z_rel)
            r = residual_internal(prob, z, scales)
            norm = scaled(r)
        J = jacobian(prob, z, scales, r)
        J, r_solve = _freeze_bound(J, z, r, n, species_tol)
        d = J.diagonal()
        # pseudo-time term along each equation's own diagonal (u-rows have d < 0)
        D = np.where(d < 0.0, -1.0, 1.0) * np.maximum(np.abs(d), 1e-10)
        for _ in range(max_retries):
            A = (J + sp.diags(D / dtau)).tocsc()
            dz = spla.spsolve(A, -r_solve)
            if not np.all(np.isfinite(dz)):
                rejected += 1
                dtau *= 0.25
                continue
            z_new, r_new, norm_new = _line_search(prob, z, dz, scales, norm, scaled, line_search)
            if norm_new < 10.0 * norm:
                break
            rejected += 1
            dtau *= 0.25
        else:
            break  # no acceptable step even with a tiny pseudo-time step
        if norm_new < best[0]:
            best = (norm_new, z_new, r_new)
        # SER with a floor on the growth: at least double after a decrease (so dtau recovers
        # after rejections even when the decrease is small), halve after an increase
        ratio = norm / max(norm_new, 1e-300)
        dtau = float(min(dtau * (max(2.0, min(ratio, 10.0)) if ratio >= 1.0 else 0.5), 1e12))
        z, r, norm = z_new, r_new, norm_new
        history.append(norm)
        it += 1
        if callback is not None:
            callback(it, z, r, dtau)
    if best[0] < norm:  # return the best state visited, not the last one
        norm, z, r = best
    s_res, e_res = _norms(r, n)
    return from_internal(z, n), CoupledStats(
        converged=norm <= 1.0,
        iterations=it,
        rejected=rejected,
        species_residual=s_res,
        energy_residual=e_res,
        history=history,
        seconds=time.perf_counter() - t0,
    )


#: Default smoothing (availability_width, min_width) of the coupled solver. Its effect on the
#: outputs is below 1e-7 in metallisation and 1e-4 K in temperature (case A, case B sweep;
#: docs/validation.md, Milestone 8); without smoothing the discrete equations have no exact
#: root where a cell sits on the availability switch (case B, 625-750 mol/s cooling flow).
DEFAULT_SMOOTHING = (1e-5, 1e-5)


def thermal_solution(
    prob: en.ThermalProblem, y: FloatArray, stats: CoupledStats
) -> en.ThermalSolution:
    """Wrap a coupled solution as a :class:`shaftsim.energy.ThermalSolution` (profiles,
    indicators and balance closures, as returned by :func:`shaftsim.energy.solve_thermal`).
    ``outer_iterations`` holds the coupled (pseudo-time) iterations."""
    from . import shaft as sh

    n = prob.n_cells
    X, f, Ts, Tg = unpack(y, n)
    sprob = prob.species_problem(Ts, Tg)
    sol = sh.make_solution(
        sprob, sprob.pack(X, f), prob.td, stats.species_residual, stats.iterations
    )
    closure, closure_s = en.energy_closure(prob, X, f, Ts, Tg)
    return en.ThermalSolution(
        species=sol,
        Ts=Ts.copy(),
        Tg=Tg.copy(),
        T_top_gas=float(Tg[0]),
        T_solid_out=float(Ts[-1]),
        energy_closure=closure,
        energy_closure_sensible=closure_s,
        outer_iterations=stats.iterations,
        converged=stats.converged,
    )


def solve_thermal_coupled(
    prob: en.ThermalProblem,
    smoothing: tuple[float, float] = DEFAULT_SMOOTHING,
    species_tol: float = 1e-10,
    energy_tol: float = 1e-6,
    coarse_cells: int = 25,
    max_iter: int = 300,
) -> tuple[en.ThermalSolution, CoupledStats]:
    """Steady state by the coupled solver, started from the alternating solver.

    Start: :func:`shaftsim.energy.solve_thermal` with loose tolerances on half the cells
    (on the grid itself below ``coarse_cells``), interpolated onto the cells; it need not
    have converged (in the band where it cycles it still gives a good start). Then
    :func:`solve_coupled` with backtracking (6 halvings) and up to 30 pseudo-time
    reductions per step.

    Args:
        prob: Thermal problem (its own smoothing widths are overridden by ``smoothing``).
        smoothing: Widths, see :func:`solve_coupled`.
        species_tol: Species tolerance [-].
        energy_tol: Energy tolerance [K].
        coarse_cells: Grid size below which the start is computed on the grid itself.
        max_iter: Maximum coupled iterations.

    Returns:
        (solution, stats); ``solution.converged`` must be checked by the caller.
    """
    from dataclasses import replace

    n = prob.n_cells
    n_c = n // 2 if n > coarse_cells else n
    cprob = replace(prob, n_cells=n_c, scheme=1)  # alternating solver: first order only
    cprob.__dict__["_tables"] = prob._tables  # tables do not depend on the grid
    start = en.solve_thermal(cprob, tol_T=1e-4, species_tol=1e-8)
    return finish_coupled(prob, start, smoothing, species_tol, energy_tol, max_iter)


def finish_coupled(
    prob: en.ThermalProblem,
    start: en.ThermalSolution,
    smoothing: tuple[float, float] = DEFAULT_SMOOTHING,
    species_tol: float = 1e-10,
    energy_tol: float = 1e-6,
    max_iter: int = 300,
) -> tuple[en.ThermalSolution, CoupledStats]:
    """Coupled solve started from a (possibly unconverged) alternating solution ``start``,
    on the same or another grid (profiles are interpolated onto the cells of ``prob``).

    Used where :func:`shaftsim.energy.solve_thermal` does not converge (Milestone 7 band):
    its end state is a good start even when it cycles. Arguments as for
    :func:`solve_thermal_coupled`.
    """
    sprob = prob.with_parameters(
        availability_width=smoothing[0],
        min_width=smoothing[1],
    )
    zc = start.species.z
    zf = prob.grid.z
    X0 = np.vstack([np.interp(zf, zc, start.species.X_steps[j]) for j in range(3)])
    y0 = pack(
        X0,
        np.interp(zf, zc, start.species.x_H2),
        np.interp(zf, zc, start.Ts),
        np.interp(zf, zc, start.Tg),
    )
    y, stats = solve_coupled(
        sprob,
        y0,
        species_tol,
        energy_tol,
        max_iter,
        smoothing=smoothing,
        line_search=6,
        max_retries=30,
    )
    return thermal_solution(sprob, y, stats), stats


def rescue(
    prob: en.ThermalProblem,
    end_state: en.ThermalSolution,
    neighbour: en.ThermalSolution | None = None,
    smoothing: tuple[float, float] = DEFAULT_SMOOTHING,
) -> tuple[en.ThermalSolution, CoupledStats, str]:
    """Coupled solve where the alternating solver failed, from the most reliable start.

    The coupled solver converges from good starts only (it is a finisher). The end state of
    a cycling alternating iteration can be far from the solution (e.g. fully reduced and
    cold), so it is tried last. Starts, in order: ``neighbour`` (a converged solution of a
    nearby operating point, e.g. in a continuation), the alternating solution on half the
    cells (:func:`solve_thermal_coupled`; the cooling front is smeared there and the
    alternating solver converges), then ``end_state``.

    Returns:
        (solution, stats, name of the start that was used: "neighbour", "coarse" or
        "end_state"); if none converges, the result of the last attempt.
    """
    tries: list[str] = (["neighbour"] if neighbour is not None else []) + ["coarse", "end_state"]
    out: tuple[en.ThermalSolution, CoupledStats] | None = None
    for name in tries:
        if name == "neighbour":
            assert neighbour is not None
            out = finish_coupled(prob, neighbour, smoothing)
        elif name == "coarse":
            out = solve_thermal_coupled(prob, smoothing)
        else:
            out = finish_coupled(prob, end_state, smoothing)
        if out[0].converged:
            return out[0], out[1], name
    assert out is not None
    return out[0], out[1], tries[-1]
