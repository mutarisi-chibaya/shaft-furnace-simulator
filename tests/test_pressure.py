"""Tests for shaftsim.pressure (Ergun pressure field, Milestone 8) and the local pressure in
the pellet kinetics (shaftsim.energy.ThermalProblem.pressure_cells)."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import pressure as ps
from shaftsim import properties as pr
from shaftsim import transport as tr

from .conftest import DATA_DIR

P_TOP = 150000.0


@pytest.fixture(scope="module")
def parts() -> tuple[dict[str, pr.Species], tr.TransportData, pl.Triplets]:
    species = pr.load_species(DATA_DIR / "species.yaml")
    mm = {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None}
    td = tr.load_transport(DATA_DIR / "transport.yaml", DATA_DIR / "collision_integrals.csv", mm)
    return species, td, pl.load_refit_triplets(DATA_DIR / "kinetics_refit.yaml")


@pytest.fixture(scope="module")
def furnace(parts) -> en.ThermalProblem:  # type: ignore[no-untyped-def]
    species, td, trip = parts
    case = en.load_furnace_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")
    return en.ThermalProblem(case, species, td, trip, 30)


def test_ergun_constants_loaded(parts) -> None:  # type: ignore[no-untyped-def]
    er = parts[1].ergun
    assert er is not None and (er.viscous, er.inertial) == (150.0, 1.75)


def test_profile_matches_direct_ergun_integration(furnace: en.ThermalProblem) -> None:
    """Exact p^2 integration per cell vs RK45 on the original form dp/dz = K_v mu u +
    K_i rho u^2, with u = F R T/(p A) and rho = p M/(R T) evaluated from p itself (piecewise
    constant flow, area, temperature and composition, as in the cells)."""
    prob = furnace
    gr = prob.grid
    n = prob.n_cells
    f = np.linspace(0.6, 0.9, n)
    Tg = np.linspace(600.0, 1150.0, n)
    prof = ps.pressure_profile(prob, f, Tg, P_TOP)
    case = prob.case.shaft
    er = prob.td.ergun
    assert er is not None
    eps, d = case.bed_porosity, case.pellet.diameter
    mm = prob.td.molar_mass
    p = P_TOP
    for i in range(n):
        x = {"H2": f[i], "H2O": gr.g[i] - f[i], case.inert: 1.0 - gr.g[i]}
        mu = tr.mixture_viscosity(prob.td, x, float(Tg[i]))
        m_mix = sum(v * mm[k] for k, v in x.items())

        def rhs(
            _z: float, y: np.ndarray, i: int = i, mu: float = mu, m_mix: float = m_mix
        ) -> list[float]:
            u = gr.F[i] * pr.R_GAS * Tg[i] / (y[0] * gr.area[i])
            rho = y[0] * m_mix / (pr.R_GAS * Tg[i])
            return [
                er.viscous * (1 - eps) ** 2 / (eps**3 * d**2) * mu * u
                + er.inertial * (1 - eps) / (eps**3 * d) * rho * u * u
            ]

        p = float(solve_ivp(rhs, (0.0, gr.dz[i]), [p], rtol=1e-12, atol=1e-6).y[0, -1])
        assert prof.p_faces[i + 1] == pytest.approx(p, rel=1e-9)
    assert prof.p_faces[0] == P_TOP
    assert np.all(np.diff(prof.p_faces) > 0.0)  # pressure rises downward
    assert np.all((prof.p_cells > prof.p_faces[:-1]) & (prof.p_cells < prof.p_faces[1:]))
    assert prof.cooling == prof.p_faces[-1]
    assert prof.bustle == prof.p_faces[gr.side_cells[0] + 1]


def test_local_pressure_rescaling_is_exact(furnace: en.ThermalProblem) -> None:
    """Rescaled table coefficients at 1.3 p equal a direct evaluation at 1.3 p (at a table
    node, so without interpolation error); a uniform pressure equal to the table pressure
    reproduces the original coefficients exactly."""
    prob = furnace
    n = prob.n_cells
    T = np.full(n, 1100.0)  # a table node
    p0 = prob.case.shaft.p
    ref = prob.cell_coefficients(T, T)
    same = prob.with_parameters(pressure_cells=tuple([p0] * n)).cell_coefficients(T, T)
    assert same == ref
    p1 = 1.3 * p0
    got = prob.with_parameters(pressure_cells=tuple([p1] * n)).cell_coefficients(T, T)[0]
    case = prob.case.shaft
    kg = float(prob.mass_transfer_kg(T)[0]) * p0 / p1  # k_g ~ 1/p
    d_mix = tr.mixture_diffusivity(prob.td, "H2", case.x_in, 1100.0, p1, "H2O")
    direct = pl.step_coefficients(
        case.pellet, 1100.0, p1, kg, d_mix, prob.species, prob.td, prob.triplets
    )
    for g, dr in zip(got, direct, strict=True):
        assert g.a_chem == pytest.approx(dr.a_chem, rel=1e-9)
        assert g.a_diff == pytest.approx(dr.a_diff, rel=1e-9)
        assert g.a_ext == pytest.approx(dr.a_ext, rel=1e-9)
        assert g.ratio == pytest.approx(dr.ratio, rel=1e-12)


def test_molecular_diffusivity_scales_as_inverse_pressure(parts) -> None:  # type: ignore[no-untyped-def]
    td = parts[1]
    x = {"H2": 0.93, "H2O": 0.02, "N2": 0.05}
    d1 = tr.mixture_diffusivity(td, "H2", x, 1000.0, 1.5e5, "H2O")
    d2 = tr.mixture_diffusivity(td, "H2", x, 1000.0, 3.0e5, "H2O")
    assert d1 == pytest.approx(2.0 * d2, rel=1e-12)
    assert math.isfinite(d1)


def test_solve_with_pressure_converges(parts) -> None:  # type: ignore[no-untyped-def]
    """Fixed point pressure <-> solution (reduction zone alone, 20 cells, for speed)."""
    species, td, trip = parts
    case = en.load_thermal_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")
    zone = en.ThermalProblem(case, species, td, trip, 20)

    def solve(p: en.ThermalProblem, prev: en.ThermalSolution | None) -> en.ThermalSolution:
        t0 = None if prev is None else (prev.Ts, prev.Tg)
        return en.solve_thermal(p, tol_T=1e-6, species_tol=1e-10, T_start=t0)

    sol, prof, it = ps.solve_with_pressure(zone, P_TOP, solve, tol=1.0)
    assert sol.converged and it < 10
    again = ps.pressure_profile(zone, sol.species.x_H2, sol.Tg, P_TOP)
    assert np.max(np.abs(again.p_cells - prof.p_cells)) == 0.0
    assert abs(sol.energy_closure) < 1e-10
