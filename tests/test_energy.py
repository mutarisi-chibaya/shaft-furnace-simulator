"""Tests for shaftsim.energy (non-isothermal 1D shaft)."""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import transport as tr

from .conftest import DATA_DIR, TEST_DATA_DIR


@pytest.fixture(scope="module")
def species() -> dict[str, pr.Species]:
    return pr.load_species(DATA_DIR / "species.yaml")


@pytest.fixture(scope="module")
def td(species: dict[str, pr.Species]) -> tr.TransportData:
    mm = {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None}
    return tr.load_transport(DATA_DIR / "transport.yaml", DATA_DIR / "collision_integrals.csv", mm)


@pytest.fixture(scope="module")
def case() -> en.ThermalCase:
    return en.load_thermal_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")


@pytest.fixture(scope="module")
def trip() -> pl.Triplets:
    return pl.load_refit_triplets(DATA_DIR / "kinetics_refit.yaml")


@pytest.fixture(scope="module")
def no_reaction(case, species, td, trip) -> en.ThermalSolution:  # type: ignore[no-untyped-def]
    prob = en.ThermalProblem(case, species, td, trip, 40, rate_multiplier=0.0)
    return en.solve_thermal(prob)


@pytest.fixture(scope="module")
def reacting(case, species, td, trip) -> en.ThermalSolution:  # type: ignore[no-untyped-def]
    return en.solve_thermal(en.ThermalProblem(case, species, td, trip, 40))


# ----------------------------------------------------------------------------- properties


def test_enthalpy_extension_continuous(species: dict[str, pr.Species]) -> None:
    h2o = species["H2O"]
    lo = h2o.T_min
    below = en.enthalpy_ext(h2o, lo - 1e-6)
    at = float(pr.enthalpy(h2o, lo))
    assert float(below) == pytest.approx(at, abs=1e-3)
    # linear below the range with slope Cp(T_min)
    slope = (en.enthalpy_ext(h2o, 400.0) - en.enthalpy_ext(h2o, 300.0)) / 100.0
    assert float(slope) == pytest.approx(float(pr.heat_capacity(h2o, lo)), rel=1e-12)


def test_enthalpy_bridged_across_shomate_ranges(species: dict[str, pr.Species]) -> None:
    """H is continuous and increasing across every NIST range boundary, cp_ext is its derivative,
    and outside the +/-10 K windows it is NIST's H unchanged."""
    w = en.TRANSITION_HALF_WIDTH
    for sp_ in species.values():
        for r in sp_.ranges[1:]:
            Tb = r.T_min
            jump = en.enthalpy_ext(sp_, Tb + 1e-9) - en.enthalpy_ext(sp_, Tb - 1e-9)
            assert abs(float(jump)) < 1e-3, sp_.name
            T = np.linspace(Tb - 2 * w, Tb + 2 * w, 401)
            assert np.all(np.diff(en.enthalpy_ext(sp_, T)) > 0), sp_.name
            mid = np.array([Tb - 0.4 * w, Tb + 0.3 * w])
            dh = (en.enthalpy_ext(sp_, mid + 1e-4) - en.enthalpy_ext(sp_, mid - 1e-4)) / 2e-4
            np.testing.assert_allclose(dh, en.cp_ext(sp_, mid), rtol=1e-6)
            out = np.array([Tb - 1.5 * w, Tb + 1.5 * w])
            out = out[(out > sp_.T_min) & (out < sp_.T_max)]
            np.testing.assert_allclose(en.enthalpy_ext(sp_, out), pr.enthalpy(sp_, out), rtol=1e-14)


def test_property_tables_accurate(case, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """Tabulated coefficients (1 K grid, log-linear) match the direct calculation to 2e-5."""
    prob = en.ThermalProblem(case, species, td, trip, 10)
    for T in (655.3, 893.7, 1101.2):
        c = prob.cell_coefficients(np.array([T]), np.array([T]))[0]
        cs = dataclasses.replace(case.shaft, T=T)
        kg = sh.bed_mass_transfer_coefficient(cs, td)
        d = tr.mixture_diffusivity(td, "H2", case.shaft.x_in, T, case.shaft.p, "H2O")
        direct = pl.step_coefficients(case.shaft.pellet, T, case.shaft.p, kg, d, species, td, trip)
        for j in range(3):
            for k in ("a_chem", "a_diff", "a_ext", "ratio"):
                assert getattr(c[j], k) == pytest.approx(getattr(direct[j], k), rel=2e-5)
        ua = prob.heat_transfer_UA(np.array([T]))[0]
        assert ua == pytest.approx(prob._ua_per_dz_direct(T) * prob.dz, rel=1e-6)


def test_svehla_h2_viscosity_and_conductivity(species, td) -> None:  # type: ignore[no-untyped-def]
    """Implementation check against Svehla (1962) Table III (same theory and constants):
    viscosity and modified-Eucken conductivity of H2 within 1 % (the collision integrals
    here come from Monchick-Mason tables, Svehla's from Hirschfelder et al.)."""
    d = pd.read_csv(TEST_DATA_DIR / "svehla1962_table3_H2.csv", comment="#")
    for T, cp_r, eta_up, lam in zip(
        d.T_K.to_numpy(float),
        d.Cp_over_R.to_numpy(float),
        d.eta_uP.to_numpy(float),
        d.lambda_1e6.to_numpy(float),
        strict=True,
    ):
        assert tr.viscosity(td, "H2", T) == pytest.approx(eta_up * 1e-7, rel=0.01)
        lam_si = lam * 1e-6 * 418.4  # cal/(cm s K) -> W/(m K)
        assert tr.thermal_conductivity(td, "H2", T, cp_r * pr.R_GAS) == pytest.approx(
            lam_si, rel=0.01
        )


def test_energy_jacobian_matches_finite_differences(case, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    prob = en.ThermalProblem(case, species, td, trip, 6)
    rng = np.random.default_rng(3)
    X = np.sort(rng.random((3, 6)), axis=0)[::-1] * 0.5
    X = np.maximum.accumulate(X, axis=1)
    f = np.linspace(0.6, 0.9, 6)
    # away from NIST Shomate segment boundaries (e.g. Fe at 700 K, where H jumps 141 J/mol)
    Ts = np.array([705.0, 785.0, 865.0, 925.0, 985.0, 1065.0])
    Tg = Ts + 40.0
    UA = prob.heat_transfer_UA(Tg)
    r0 = prob.energy_residual(X, f, Ts, Tg, UA)
    _, Cs = prob.solid_enthalpy(X, Ts)
    _, Cg = prob.gas_enthalpy(f, Tg)
    f_down = np.concatenate([f[1:], [case.shaft.x_in["H2"]]])
    n_o = case.shaft.gas_flow * (f_down - f)
    J = en.energy_jacobian(
        Cs, Cg, UA, n_o, en.cp_ext(species["H2"], Tg), en.cp_ext(species["H2O"], Ts)
    ).toarray()
    for k in range(12):
        h = 1e-3
        Ts2, Tg2 = Ts.copy(), Tg.copy()
        (Ts2 if k % 2 == 0 else Tg2)[k // 2] += h
        dr = (prob.energy_residual(X, f, Ts2, Tg2, UA) - r0) / h
        np.testing.assert_allclose(J[:, k], dr, rtol=1e-4, atol=1e-3 * np.abs(J).max() * 1e-6)


# ----------------------------------------------------------------------------- balances


def test_energy_balance_closes(no_reaction, reacting) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: energy balance closes to 1e-6 of the inlet enthalpy flow (here: round-off)."""
    for s in (no_reaction, reacting):
        assert s.converged
        assert abs(s.energy_closure) < 1e-10
        assert abs(s.energy_closure_sensible) < 1e-10


def test_element_balances_still_close(case, species, td, trip, reacting) -> None:  # type: ignore[no-untyped-def]
    prob = en.ThermalProblem(case, species, td, trip, 40)
    sprob = prob.species_problem(reacting.Ts, reacting.Tg)
    for err in sh.element_balances(sprob, reacting.species).values():
        assert abs(err) < 1e-12


# ----------------------------------------------------------------------------- heat exchange


def test_no_transfer_keeps_inlet_temperatures(case, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    prob = en.ThermalProblem(case, species, td, trip, 20, rate_multiplier=0.0, h_multiplier=0.0)
    s = en.solve_thermal(prob)
    np.testing.assert_allclose(s.Ts, case.T_solid_in, atol=1e-6)
    np.testing.assert_allclose(s.Tg, case.T_gas_in, atol=1e-6)


def test_counter_current_exchange_profiles(no_reaction, case) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: with reactions off, gas and solid temperatures approach each other as in a
    counter-current exchanger: gas hotter than solid everywhere, both increasing downward, and
    the stream with the smaller heat-capacity flow (the burden) leaves near the inlet
    temperature of the other stream (pinch at the bottom)."""
    s = no_reaction
    assert np.all(s.Tg > s.Ts)
    assert np.all(np.diff(s.Ts) > 0) and np.all(np.diff(s.Tg) > 0)
    assert s.T_solid_out == pytest.approx(case.T_gas_in, abs=1.0)
    gap = s.Tg - s.Ts
    assert gap[-1] < gap[0]  # the streams approach each other toward the bottom


def test_large_transfer_reaches_pinch(case, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """h x 1e3: the burden leaves at the gas inlet temperature; the top gas temperature then
    follows from the energy balance alone."""
    prob = en.ThermalProblem(case, species, td, trip, 40, rate_multiplier=0.0, h_multiplier=1e3)
    s = en.solve_thermal(prob)
    assert s.T_solid_out == pytest.approx(case.T_gas_in, abs=1e-3)
    X0 = np.zeros((3, 1))
    f_in = np.array([case.shaft.x_in["H2"]])
    hs_in = prob.solid_enthalpy(X0, np.array([case.T_solid_in]))[0][0]
    hs_out = prob.solid_enthalpy(X0, np.array([case.T_gas_in]))[0][0]
    hg_in = prob.gas_enthalpy(f_in, np.array([case.T_gas_in]))[0][0]
    # top gas enthalpy = hg_in - (hs_out - hs_in); invert by bisection
    target = hg_in - (hs_out - hs_in)
    lo, hi = 300.0, case.T_gas_in
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if prob.gas_enthalpy(f_in, np.array([mid]))[0][0] > target:
            hi = mid
        else:
            lo = mid
    assert s.T_top_gas == pytest.approx(0.5 * (lo + hi), abs=1e-3)


def test_heat_exchange_matches_shooting_reference(case, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """Reactions off, h x 3e-3 (at the real h the bed pinches and the top gas follows from the
    overall balance alone, which would make the comparison trivial): finite volumes (50/100
    cells, Richardson) agree with an independent shooting solution of
    dTs/dz = q/Cs, dTg/dz = q/Cg, q = UA'(Tg - Ts) (adaptive RK45, rtol 1e-10; brentq on the
    top gas temperature). Cs, Cg, UA' are tabulated on a 0.5 K grid refined at the edges of
    the enthalpy bridging windows of enthalpy_ext, where Cp jumps. The shooting solution
    itself agrees with scipy's solve_bvp to 1e-6 K (checked once; too slow for the suite)."""
    hm = 3e-3
    fv = {}
    for n in (50, 100):
        prob = en.ThermalProblem(case, species, td, trip, n, rate_multiplier=0.0, h_multiplier=hm)
        s = en.solve_thermal(prob)
        assert s.converged
        fv[n] = (s.T_top_gas, s.T_solid_out)
    top_fv, sout_fv = (2 * fv[100][k] - fv[50][k] for k in range(2))

    prob = en.ThermalProblem(case, species, td, trip, 1, rate_multiplier=0.0)
    f_in = case.shaft.x_in["H2"]
    # Cp jumps at the edges of the enthalpy bridging windows: put nodes just either side of
    # every edge so linear interpolation does not smear the jumps (~0.05 K each otherwise)
    w = en.TRANSITION_HALF_WIDTH
    edges = np.array([r.T_min + d for v in species.values() for r in v.ranges[1:] for d in (-w, w)])
    grid = np.unique(np.concatenate([np.arange(250.0, 1300.25, 0.5), edges - 1e-6, edges + 1e-6]))
    cs = prob.solid_enthalpy(np.zeros((3, grid.size)), grid)[1]
    cg = prob.gas_enthalpy(np.full_like(grid, f_in), grid)[1]
    ua = prob.ua_per_length_ref(grid) * hm

    def rhs(_z: float, y: np.ndarray) -> list[float]:
        Ts, Tg = y
        q = np.interp(Tg, grid, ua) * (Tg - Ts)
        return [q / np.interp(Ts, grid, cs), q / np.interp(Tg, grid, cg)]

    def shoot(t_top: float) -> tuple[float, float]:
        r = solve_ivp(
            rhs, (0.0, case.shaft.height), [case.T_solid_in, t_top], rtol=1e-10, atol=1e-8
        )
        return float(r.y[1, -1]), float(r.y[0, -1])

    top = brentq(lambda t: shoot(t)[0] - case.T_gas_in, case.T_solid_in, case.T_gas_in, xtol=1e-9)
    assert top_fv == pytest.approx(top, abs=0.05)
    assert sout_fv == pytest.approx(shoot(top)[1], abs=0.08)


# ----------------------------------------------------------------------------- reactions


def test_endothermic_reduction_cools_top_gas(no_reaction, reacting) -> None:  # type: ignore[no-untyped-def]
    """H2 reduction is endothermic overall: with reactions the top gas leaves colder."""
    assert reacting.T_top_gas < no_reaction.T_top_gas - 50.0
    assert 0.0 < reacting.species.metallization < 1.0
