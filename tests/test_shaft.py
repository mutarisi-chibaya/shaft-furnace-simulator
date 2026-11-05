"""Tests for shaftsim.shaft and shaftsim.solver (isothermal 1D counter-current shaft)."""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import solver as sv
from shaftsim import transport as tr

from .conftest import DATA_DIR

T800 = 1073.15


@pytest.fixture(scope="module")
def species() -> dict[str, pr.Species]:
    return pr.load_species(DATA_DIR / "species.yaml")


@pytest.fixture(scope="module")
def td(species: dict[str, pr.Species]) -> tr.TransportData:
    mm = {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None}
    return tr.load_transport(DATA_DIR / "transport.yaml", DATA_DIR / "collision_integrals.csv", mm)


@pytest.fixture(scope="module")
def trip() -> pl.Triplets:
    return pl.load_refit_triplets(DATA_DIR / "kinetics_refit.yaml")


def _case(T: float, **changes: float) -> sh.ShaftCase:
    base = sh.load_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml", T=T)
    return dataclasses.replace(base, **changes)  # type: ignore[arg-type]


def _solve(
    case: sh.ShaftCase,
    species: dict[str, pr.Species],
    td: tr.TransportData,
    trip: pl.Triplets,
    n: int = 100,
    mult: float = 1.0,
    y0: np.ndarray | None = None,
) -> tuple[sh.ShaftProblem, np.ndarray, sh.ShaftSolution]:
    prob = sh.ShaftProblem(case, sh.shaft_coefficients(case, species, td, trip), n, mult)
    y, st = sv.solve(prob, y0=y0)
    assert st.converged, (st.method, st.residuals[-1])
    return prob, y, sh.make_solution(prob, y, td, st.residuals[-1], st.iterations)


def _pinch_metallization(case: sh.ShaftCase, species: dict[str, pr.Species]) -> float:
    """Thermodynamic limit with instantaneous reactions (counter-current pinch).

    The gas leaving the iron zone can carry at most the H2O fraction in equilibrium with
    Fe/FexO: x_H2O <= g r3/(1 + r3), g = x_H2 + x_H2O. So the oxygen removable from wustite is
    F (g r3/(1+r3) - x_H2O,in); metallisation = min(1, that / O_wustite).
    """
    r3 = float(pr.equilibrium_h2o_h2_ratio(pr.WUSTITE_TO_IRON, species, case.T))
    g = case.x_in["H2"] + case.x_in["H2O"]
    o3 = case.gas_flow * (g * r3 / (1 + r3) - case.x_in["H2O"])
    return min(1.0, o3 / (pl.OXYGEN_WEIGHTS[2] * case.removable_oxygen_flow))


# ----------------------------------------------------------------------------- case data


def test_case_a_geometry_and_flows() -> None:
    c = _case(1163.15)
    assert c.height == pytest.approx(10.1)
    assert c.area == pytest.approx(math.pi * 2.6**2)
    assert c.gas_flow == 2800.0
    assert c.x_in["H2"] + c.x_in["H2O"] + c.x_in["N2"] == pytest.approx(1.0)
    # solid velocity from the volume balance of the burden
    v_solid = c.solid_mass_flow / c.pellet.density
    assert c.solid_velocity * (1 - c.bed_porosity) * c.area == pytest.approx(v_solid)
    # axial solid mass flux ~2.1 kg m-2 s-1, consistent with Fig. 2 of Marsigny et al. 2025b
    assert c.solid_mass_flow / c.area == pytest.approx(2.12, abs=0.01)


def test_below_wustite_limit_solves_and_reduces_less(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """Below the wustite stability limit (893 K) the model still converges (metastable
    wustite) and the burden is less reduced than at higher temperature."""
    m = []
    for T in (823.15, 873.15, 973.15):
        m.append(_solve(_case(T, height=10.1), species, td, trip)[2].metallization)
    assert m[0] < m[1] < m[2]


# ----------------------------------------------------------------------------- balances


@pytest.mark.parametrize(("T", "flow"), [(T800, 2800.0), (T800, 1639.0), (1163.15, 2800.0)])
def test_element_balances_close(species, td, trip, T: float, flow: float) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: Fe, O, H balances close to 1e-6. By construction they close to round-off."""
    prob, _, sol = _solve(_case(T, gas_flow=flow), species, td, trip)
    for el, err in sh.element_balances(prob, sol).items():
        assert abs(err) < 1e-12, el


def test_top_gas_from_mass_balance(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """With full metallisation the top-gas H2O follows from stoichiometry alone."""
    case = _case(1163.15)
    _, _, sol = _solve(case, species, td, trip)
    assert sol.metallization == pytest.approx(1.0, abs=1e-9)
    x_h2o_top = case.x_in["H2O"] + case.removable_oxygen_flow / case.gas_flow
    assert sol.x_H2O[0] == pytest.approx(x_h2o_top, rel=1e-9)
    assert sol.gud == pytest.approx(x_h2o_top / (case.x_in["H2"] + case.x_in["H2O"]), rel=1e-9)


# ----------------------------------------------------------------------------- limits


def test_zero_kinetics_gives_no_reduction(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    case = _case(T800)
    _, _, sol = _solve(case, species, td, trip, n=50, mult=0.0)
    assert np.all(sol.X_steps == 0.0)
    np.testing.assert_allclose(sol.x_H2, case.x_in["H2"], rtol=0, atol=1e-15)


@pytest.mark.parametrize("frac", [0.6, 0.85])
def test_gas_limited_shaft_reaches_equilibrium_pinch(species, td, trip, frac: float) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: fast kinetics reproduce the equilibrium limit of Milestone 1.

    Below the minimum gas flow for full metallisation, the physical kinetics at 800 degC are
    already fast enough for the outlet to sit on the thermodynamic pinch.
    """
    case0 = _case(T800)
    g = case0.x_in["H2"] + case0.x_in["H2O"]
    r3 = float(pr.equilibrium_h2o_h2_ratio(pr.WUSTITE_TO_IRON, species, T800))
    f_min = (
        pl.OXYGEN_WEIGHTS[2] * case0.removable_oxygen_flow / (g * r3 / (1 + r3) - case0.x_in["H2O"])
    )
    case = dataclasses.replace(case0, gas_flow=frac * f_min)
    _, _, sol = _solve(case, species, td, trip)
    assert sol.metallization == pytest.approx(_pinch_metallization(case, species), abs=1e-5)
    assert sol.metallization == pytest.approx(frac, abs=1e-3)


def test_faster_kinetics_do_not_beat_thermodynamics(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """Rates x10, x100, x1000 give the same pinch value as the physical kinetics.

    Solved by continuation: each case starts from the previous solution (the discrete
    shooting cannot resolve the near-discontinuity of very fast kinetics, see solver.py).
    """
    case = dataclasses.replace(_case(T800), gas_flow=1157.0)
    pinch = _pinch_metallization(case, species)
    y = None
    for mult in (1.0, 10.0, 100.0, 1000.0):
        _, y, sol = _solve(case, species, td, trip, mult=mult, y0=y)
        assert sol.metallization == pytest.approx(pinch, abs=1e-5), mult


# ----------------------------------------------------------------------------- numerics


def test_grid_convergence_first_order(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: results converge with grid refinement. First-order upwind: errors halve."""
    case = _case(973.15, height=3.0)  # 700 degC, short zone: outlet not saturated
    m = [_solve(case, species, td, trip, n=n)[2].metallization for n in (50, 100, 200, 400)]
    d = np.diff(m)
    ratios = d[:-1] / d[1:]
    assert np.all((ratios > 1.8) & (ratios < 2.2)), ratios
    richardson_error_400 = abs(d[-1])  # first order: error(n) ~ m(n) - m(n/2)
    assert richardson_error_400 < 0.003


def test_solver_methods_agree(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """Shooting and coupled Gauss-Seidel sweeps converge to the same discrete solution."""
    case = _case(973.15, height=3.0)
    prob = sh.ShaftProblem(case, sh.shaft_coefficients(case, species, td, trip), 60)
    y_shoot, st1 = sv.solve(prob)
    assert st1.method == "shooting"
    y_sweep, st2 = sv.solve(prob, y0=sv.initial_guess(prob))
    assert st2.method == "sweeps" and st2.converged
    np.testing.assert_allclose(y_shoot, y_sweep, atol=1e-8)


def test_profiles_monotone(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    _, _, sol = _solve(_case(973.15, height=3.0), species, td, trip)
    assert np.all(np.diff(sol.X) >= -1e-12)  # conversion increases downward
    assert np.all(np.diff(sol.x_H2) >= -1e-12)  # H2 decreases upward
    assert np.all(sol.X_steps[1] <= sol.X_steps[0] + 1e-12)
    assert np.all(sol.X_steps[2] <= sol.X_steps[1] + 1e-12)


def test_residual_and_jacobian_consistent(species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """The banded finite-difference Jacobian matches a dense one (fallback solver)."""
    case = _case(T800, height=2.0)
    prob = sh.ShaftProblem(case, sh.shaft_coefficients(case, species, td, trip), 8)
    y = sv.shoot(prob) * 0.97
    r0 = prob.residual(y)
    J = sv.banded_jacobian(prob, y, r0).toarray()
    dense = np.empty_like(J)
    for k in range(y.size):
        h = 1e-7 * max(1.0, abs(y[k]))
        yp = y.copy()
        yp[k] += h
        dense[:, k] = (prob.residual(yp) - r0) / h
    np.testing.assert_allclose(J, dense, atol=1e-5)
