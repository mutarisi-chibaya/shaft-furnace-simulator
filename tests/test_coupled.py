"""Tests for shaftsim.coupled (fully coupled Newton solver, Milestone 8) and the smoothed
availability constraint of shaftsim.pellet.step_rates_array."""

from __future__ import annotations

import numpy as np
import pytest

from shaftsim import coupled as cp
from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import transport as tr

from .conftest import DATA_DIR

N = 40


@pytest.fixture(scope="module")
def prob() -> en.ThermalProblem:
    species = pr.load_species(DATA_DIR / "species.yaml")
    mm = {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None}
    td = tr.load_transport(DATA_DIR / "transport.yaml", DATA_DIR / "collision_integrals.csv", mm)
    case = en.load_thermal_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")
    trip = pl.load_refit_triplets(DATA_DIR / "kinetics_refit.yaml")
    return en.ThermalProblem(case, species, td, trip, N)


@pytest.fixture(scope="module")
def alternating(prob: en.ThermalProblem) -> en.ThermalSolution:
    sol = en.solve_thermal(prob, tol_T=1e-7, species_tol=1e-10)
    assert sol.converged
    return sol


def _y(sol: en.ThermalSolution) -> np.ndarray:
    return cp.pack(sol.species.X_steps, sol.species.x_H2, sol.Ts, sol.Tg)


@pytest.fixture(scope="module")
def smoothed(
    prob: en.ThermalProblem, alternating: en.ThermalSolution
) -> tuple[en.ThermalProblem, np.ndarray, cp.CoupledStats]:
    y, st = cp.solve_coupled(
        prob, _y(alternating), smoothing=cp.DEFAULT_SMOOTHING, line_search=6, max_retries=30
    )
    p = prob.with_parameters(
        availability_width=cp.DEFAULT_SMOOTHING[0], min_width=cp.DEFAULT_SMOOTHING[1]
    )
    return p, y, st


# ----------------------------------------------------------------------------- variables


def test_pack_and_internal_round_trip() -> None:
    rng = np.random.default_rng(1)
    n = 7
    X = rng.random((3, n))
    f, Ts, Tg = rng.random(n), 300 + 900 * rng.random(n), 300 + 900 * rng.random(n)
    y = cp.pack(X, f, Ts, Tg)
    for a, b in zip(cp.unpack(y, n), (X, f, Ts, Tg), strict=True):
        np.testing.assert_array_equal(a, b)
    np.testing.assert_allclose(cp.from_internal(cp.to_internal(y, n), n), y, rtol=0, atol=1e-14)


def test_jacobian_matches_directional_derivative(
    prob: en.ThermalProblem, alternating: en.ThermalSolution
) -> None:
    """Coloured FD Jacobian vs a central directional difference along a random direction,
    at a perturbed reacting state (all three steps active somewhere)."""
    p = prob.with_parameters(availability_width=1e-3, min_width=1e-3)
    n = prob.n_cells
    z = cp.to_internal(_y(alternating), n)
    rng = np.random.default_rng(7)
    z = cp.project(p, z + 1e-3 * rng.standard_normal(z.size) * np.tile([1, 1, 1, 0.1, 100, 100], n))
    scales = cp.energy_scales(p, cp.from_internal(z, n))
    r0 = cp.residual_internal(p, z, scales)
    J = cp.jacobian(p, z, scales, r0)
    v = rng.standard_normal(z.size) * np.tile([1e-6, 1e-6, 1e-6, 1e-7, 1e-3, 1e-3], n)
    # keep the direction admissible (away from the core-radius bounds)
    Z = z.reshape(n, cp.N_VAR)
    V = v.reshape(n, cp.N_VAR)
    V[:, :3] = np.where((Z[:, :3] < 1e-4) | (Z[:, :3] > 1 - 1e-4), 0.0, V[:, :3])
    d = (cp.residual_internal(p, z + v, scales) - cp.residual_internal(p, z - v, scales)) / 2
    np.testing.assert_allclose(J @ v, d, rtol=0, atol=2e-3 * np.max(np.abs(d)))


# ----------------------------------------------------------------------------- solver


def test_coupled_equals_alternating_without_smoothing(
    prob: en.ThermalProblem, alternating: en.ThermalSolution
) -> None:
    """PLAN.md, Milestone 8: same discrete equations, so the same solution where both
    solvers converge (1e-8)."""
    y, st = cp.solve_coupled(prob, _y(alternating))
    assert st.converged
    X, _f, Ts, Tg = cp.unpack(y, prob.n_cells)
    assert abs(X[2, -1] - alternating.species.metallization) < 1e-8
    np.testing.assert_allclose(X, alternating.species.X_steps, rtol=0, atol=1e-8)
    np.testing.assert_allclose(Ts, alternating.Ts, rtol=0, atol=1e-6)
    np.testing.assert_allclose(Tg, alternating.Tg, rtol=0, atol=1e-6)


def test_smoothed_solution_converges_and_balances_close(
    smoothed: tuple[en.ThermalProblem, np.ndarray, cp.CoupledStats],
) -> None:
    p, y, st = smoothed
    assert st.converged
    sol = cp.thermal_solution(p, y, st)
    assert abs(sol.energy_closure) < 1e-10
    assert abs(sol.energy_closure_sensible) < 1e-10
    sprob = p.species_problem(sol.Ts, sol.Tg)
    for err in sh.element_balances(sprob, sol.species).values():
        assert abs(err) < 1e-12


def test_smoothing_effect_below_tolerance(
    alternating: en.ThermalSolution,
    smoothed: tuple[en.ThermalProblem, np.ndarray, cp.CoupledStats],
) -> None:
    """PLAN.md, Milestone 8: effect of the default smoothing on the outputs below 1e-6
    (metallisation, absolute) and 1e-6 relative for temperatures."""
    p, y, st = smoothed
    sol = cp.thermal_solution(p, y, st)
    assert abs(sol.species.metallization - alternating.species.metallization) < 1e-6
    assert abs(sol.T_solid_out - alternating.T_solid_out) < 1e-6 * alternating.T_solid_out
    assert abs(sol.T_top_gas - alternating.T_top_gas) < 1e-6 * alternating.T_top_gas


def test_release_lifts_core_radius_off_bound() -> None:
    """A core radius at U_MIN whose equation asks it to grow is reset to cbrt(u^3 + R); one
    whose equation is satisfied stays."""
    n = 2
    z = np.zeros(n * cp.N_VAR)
    Z = z.reshape(n, cp.N_VAR)
    Z[:, :3] = [[cp.U_MIN, cp.U_MIN, 0.5], [cp.U_MIN, 0.4, 0.4]]
    r = np.zeros_like(z)
    R = r.reshape(n, cp.N_VAR)
    R[0, 1] = 7e-5  # X too large by 7e-5
    R[1, 0] = -1e-12  # pushes further down: stays
    z2, changed = cp._release_bound(z, r, n, 1e-10)
    assert changed
    Z2 = z2.reshape(n, cp.N_VAR)
    assert Z2[0, 1] == pytest.approx(np.cbrt(cp.U_MIN**3 + 7e-5), rel=1e-12)
    assert Z2[0, 0] == cp.U_MIN and Z2[1, 0] == cp.U_MIN
    np.testing.assert_array_equal(np.delete(z2, 1), np.delete(z, 1))


# ----------------------------------------------------------------------------- rate law


def _coeffs(
    n: int, a3: float = 2900.0
) -> tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients]:
    one = np.ones(n)
    return tuple(  # type: ignore[return-value]
        pl.StepCoefficients(a_chem=a * one, a_diff=50.0 * one, a_ext=2.0 * one, ratio=k * one)
        for a, k in ((400.0, 6e4), (1000.0, 0.35), (a3, 0.37))
    )


def test_smoothed_rates_reduce_to_original() -> None:
    rng = np.random.default_rng(5)
    n = 50
    X = np.sort(rng.random((3, n)), axis=0)[::-1]
    X[2, :10] = X[1, :10]  # some cells exactly on the availability switch
    xh2 = 0.5 + 0.4 * rng.random(n)
    xh2o = 0.95 - xh2
    c = _coeffs(n)
    r0 = pl.step_rates_array(c, X, xh2, xh2o)
    r1 = pl.step_rates_array(c, X, xh2, xh2o, one_minus=1.0 - X)
    np.testing.assert_allclose(r1, r0, rtol=1e-12, atol=0)
    r2 = pl.step_rates_array(c, X, xh2, xh2o, availability_width=1e-12, min_width=1e-12)
    np.testing.assert_allclose(r2, r0, rtol=1e-9, atol=0)


def test_availability_smoothing_is_continuous_across_switch() -> None:
    """Original rule: step 3 jumps from r3 to min(r3, r2) when X3 reaches X2. Smoothed: the
    rate changes continuously with the gap X2 - X3."""
    c = _coeffs(1, a3=100.0)  # step 3 intrinsically faster than step 2
    xh2, xh2o = np.array([0.8]), np.array([0.15])
    gaps = np.array([-1e-9, 0.0, 1e-9])

    def r3(gap: float, **kw: float) -> float:
        X = np.array([[0.9], [0.3], [0.3 - gap]])
        return float(pl.step_rates_array(c, X, xh2, xh2o, **kw)[2, 0])

    raw = [r3(g) for g in gaps]
    assert abs(raw[2] - raw[1]) > 0.1 * abs(raw[2])  # the discontinuity
    sm = [r3(g, availability_width=1e-5, min_width=1e-5) for g in gaps]
    assert abs(sm[2] - sm[0]) < 1e-3 * abs(sm[1])


def test_smooth_min_keeps_step_at_rest() -> None:
    """A step whose predecessor is at equilibrium (rate 0) must not move: the smooth
    minimum is exactly 0 if either rate is 0, and never exceeds the exact minimum."""
    c = _coeffs(3)
    X = np.zeros((3, 3))
    x_eq2 = 0.95 / 1.35  # step-2 equilibrium H2 fraction: r2 = 0 at and below it
    xh2 = np.array([x_eq2 - 0.01, x_eq2, 0.9])
    xh2o = 0.95 - xh2
    r = pl.step_rates_array(c, X, xh2, xh2o, availability_width=1e-5, min_width=1e-5)
    exact = pl.step_rates_array(c, X, xh2, xh2o)
    assert r[1, 0] == 0.0 and r[2, 0] == 0.0
    assert np.all(r[2] <= exact[2] * (1 + 1e-15))
    assert np.all(r >= 0.0)


# ----------------------------------------------------------------------------- scheme 2


def test_outflow_faces_exact_for_linear_profiles_on_nonuniform_grid() -> None:
    dz = np.array([0.1, 0.1, 0.3, 0.3, 0.3, 0.05, 0.05, 0.2])
    zf = np.concatenate([[0.0], np.cumsum(dz)])
    z = 0.5 * (zf[:-1] + zf[1:])
    phi = 2.0 + 3.0 * z
    down = sh.outflow_faces(phi, z, dz, downward=True)
    up = sh.outflow_faces(phi, z, dz, downward=False)
    np.testing.assert_allclose(down[1:-1], 2.0 + 3.0 * zf[2:-1], rtol=0, atol=1e-13)
    np.testing.assert_allclose(up[1:-1], 2.0 + 3.0 * zf[1:-2], rtol=0, atol=1e-13)
    # end cells: first order (face = cell value)
    assert down[0] == phi[0] and down[-1] == phi[-1] and up[0] == phi[0] and up[-1] == phi[-1]


def test_outflow_faces_stay_between_neighbours_for_monotone_profiles() -> None:
    rng = np.random.default_rng(11)
    n = 60
    z = np.arange(n) + 0.5
    dz = np.ones(n)
    phi = np.sort(rng.random(n)) ** 3  # monotone, strongly varying slopes
    phi[-10:] = 1.0  # flat at complete conversion
    down = sh.outflow_faces(phi, z, dz, downward=True)
    assert np.all(down[:-1] >= phi[:-1] - 1e-15)
    assert np.all(down[:-1] <= phi[1:] + 1e-15)
    assert np.all(down <= 1.0)


def test_scheme2_fluxes_telescope(prob: en.ThermalProblem, alternating: en.ThermalSolution) -> None:
    """Second-order faces are conservative for any state, converged or not: the sum of the
    cell energy residuals equals the boundary enthalpy flows, and the sum of the gas
    residuals (times F) equals H2 in - out - oxygen removed (in - out of the solid)."""
    p2 = prob.with_parameters(scheme=2)
    rng = np.random.default_rng(4)
    s = alternating
    n = prob.n_cells
    X = np.clip(s.species.X_steps + 0.01 * rng.standard_normal((3, n)), 0.0, 1.0)
    f = s.species.x_H2 * (1 + 0.01 * rng.standard_normal(n))
    Ts, Tg = s.Ts + rng.standard_normal(n), s.Tg + rng.standard_normal(n)
    gr = p2.grid
    r = p2.energy_residual(X, f, Ts, Tg, p2.heat_transfer_UA(Tg)).reshape(n, 2)
    Hs_in = p2.solid_enthalpy(np.zeros((3, 1)), np.array([p2.case.T_solid_in]))[0][0]
    Hg_in = p2.gas_enthalpy(np.array([p2.case.shaft.x_in["H2"]]), np.array([p2.case.T_gas_in]))[0][
        0
    ]
    Hs_out = p2.solid_enthalpy(X[:, -1:], Ts[-1:])[0][0]
    Hg_out = p2.gas_enthalpy(f[:1], Tg[:1], gr.F[:1], gr.g[:1])[0][0]
    expected = Hs_out - Hs_in + Hg_out - Hg_in - float(np.sum(p2.side_enthalpy))
    assert float(np.sum(r)) == pytest.approx(expected, rel=1e-11, abs=1e-6 * abs(Hs_in))
    sp2 = p2.species_problem(Ts, Tg)
    rg = sp2.residual(sp2.pack(X, f)).reshape(n, 4)[:, 3] * gr.F
    case = p2.case.shaft
    h2_in = case.gas_flow * case.x_in["H2"] + float(np.sum(gr.S * gr.f_side))
    o_out = case.removable_oxygen_flow * float(pl.overall_conversion(X[:, -1]))
    assert float(np.sum(rg)) == pytest.approx(h2_in - gr.F[0] * f[0] - o_out, rel=1e-11)


def test_scheme2_second_order_without_reactions(prob: en.ThermalProblem) -> None:
    """Reaction-free heat exchange (h x 3e-3, as in the shooting-reference test of
    test_energy): grid self-convergence ratio of the top gas temperature on 25/50/100 cells
    is about 4 for scheme 2 (2 for scheme 1). The reacting solution has kinks (availability,
    equilibrium switch-off), where any scheme is first order (docs/validation.md, M8)."""
    tops = []
    for n in (25, 50, 100):
        p = prob.with_parameters(n_cells=n, rate_multiplier=0.0, h_multiplier=3e-3)
        s1 = en.solve_thermal(p, tol_T=1e-9, species_tol=1e-11)
        y, st = cp.solve_coupled(p.with_parameters(scheme=2), _y(s1), energy_tol=1e-9)
        assert st.converged
        tops.append(cp.unpack(y, n)[3][0])
    ratio = (tops[0] - tops[1]) / (tops[1] - tops[2])
    assert 3.6 < ratio < 4.4


def test_scheme_guards(prob: en.ThermalProblem, alternating: en.ThermalSolution) -> None:
    """The marching species solver and the energy Newton refuse the second-order scheme
    (they assume donor-cell fluxes)."""
    from shaftsim import solver as sv

    p2 = prob.with_parameters(scheme=2)
    s = alternating
    with pytest.raises(ValueError):
        sv.solve(p2.species_problem(s.Ts, s.Tg))
    with pytest.raises(ValueError):
        p2.solve_energy(s.species.X_steps, s.species.x_H2, s.Ts, s.Tg)
