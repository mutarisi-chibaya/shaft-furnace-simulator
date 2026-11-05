"""Tests for shaftsim.pellet (grainy pellet model, law of additive reaction times)."""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pandas as pd
import pytest

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr

from .conftest import DATA_DIR, TEST_DATA_DIR

T900 = 1173.15


@pytest.fixture(scope="module")
def species() -> dict[str, pr.Species]:
    return pr.load_species(DATA_DIR / "species.yaml")


@pytest.fixture(scope="module")
def td(species: dict[str, pr.Species]) -> tr.TransportData:
    mm = {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None}
    return tr.load_transport(DATA_DIR / "transport.yaml", DATA_DIR / "collision_integrals.csv", mm)


@pytest.fixture(scope="module")
def base() -> tuple[pl.PelletParameters, pl.TGAConditions]:
    return pl.load_kinetics(DATA_DIR / "kinetics.yaml")


@pytest.fixture(scope="module")
def refit() -> pl.Triplets:
    return pl.load_refit_triplets(DATA_DIR / "kinetics_refit.yaml")


@pytest.fixture(scope="module")
def run900(base, species, td, refit) -> pl.PelletResult:  # type: ignore[no-untyped-def]
    pp, cond = base
    return pl.simulate_tga(pp, cond, T900, 3600.0, species, td, triplets=refit, n_out=1200)


# ----------------------------------------------------------------------------- structure


def test_pellet_mass_matches_paper(base) -> None:  # type: ignore[no-untyped-def]
    """Marsigny 2025 Sec. 2.1: a 13 mm pellet is 'approximately 4 g'."""
    pp, _ = base
    mass = pp.density * math.pi * pp.diameter**3 / 6
    assert mass == pytest.approx(4.0e-3, rel=0.05)


def test_layer_porosities_increase_and_bounded(base) -> None:  # type: ignore[no-untyped-def]
    pp, _ = base
    eps = (pp.porosity, *pl.layer_porosities(pp))
    assert all(0 < e < 1 for e in eps)
    assert list(eps) == sorted(eps)


def test_grain_diameter_from_specific_area(base) -> None:  # type: ignore[no-untyped-def]
    pp, _ = base
    d = pl.grain_diameter(pp)
    sphere_area_per_mass = 6 / (pp.true_density["Fe2O3"] * d)
    assert sphere_area_per_mass == pytest.approx(pp.specific_area)


def test_oxygen_weights_sum_to_one() -> None:
    assert sum(pl.OXYGEN_WEIGHTS) == pytest.approx(1.0)


def test_solid_moles_conserve_iron_and_oxygen(base) -> None:  # type: ignore[no-untyped-def]
    pp, _ = base
    n0 = pl.initial_moles_fe2o3(pp)
    for X in ([0, 0, 0], [0.7, 0.4, 0.1], [1, 1, 1]):
        n = pl.solid_moles(pp, np.array(X, float))
        fe = 2 * n["Fe2O3"] + 3 * n["Fe3O4"] + pl.X_FE_WUSTITE_KINETIC * n["Fe0.95O"] + n["Fe"]
        assert fe == pytest.approx(2 * n0)
        o = 3 * n["Fe2O3"] + 4 * n["Fe3O4"] + n["Fe0.95O"]
        removed = pl.O_PER_FE2O3 * n0 * float(pl.overall_conversion(np.array(X, float)))
        assert o + removed == pytest.approx(3 * n0)


# ----------------------------------------------------------------------------- rate law


def test_rate_zero_at_and_beyond_equilibrium() -> None:
    sc = pl.StepCoefficients(a_chem=100.0, a_diff=10.0, a_ext=1.0, ratio=0.5)
    x_eq = 1.0 / (1.0 + 0.5)  # H2 fraction at equilibrium in pure H2/H2O
    for x_h2 in (x_eq, x_eq - 0.05):
        times = sc.times(x_h2, 1.0 - x_h2)
        assert pl.lart_rate(times, 0.3) == 0.0
    assert pl.lart_rate(sc.times(x_eq + 0.05, 1.0 - x_eq - 0.05), 0.3) > 0.0


def test_rate_zero_at_full_conversion() -> None:
    assert pl.lart_rate(pl.StepTimes(10.0, 1.0, 1.0), 1.0) == 0.0


def test_chemical_regime_analytic_solution() -> None:
    """With only chemical control, X(t) = 1 - (1 - t/tau)^3 (shrinking sphere)."""
    tau = 500.0
    times = pl.StepTimes(tau, 0.0, 0.0)
    t = np.linspace(0, tau, 2001)
    X = np.zeros_like(t)
    for i in range(1, len(t)):  # midpoint integration
        dt = t[i] - t[i - 1]
        xm = X[i - 1] + 0.5 * dt * pl.lart_rate(times, X[i - 1])
        X[i] = min(X[i - 1] + dt * pl.lart_rate(times, xm), 1.0)
    exact = 1 - (1 - t / tau) ** 3
    assert np.max(np.abs(X - exact)) < 2e-3


def test_resistances_add() -> None:
    """Adding a resistance can only slow the step."""
    r_chem = pl.lart_rate(pl.StepTimes(100.0, 0.0, 0.0), 0.5)
    r_all = pl.lart_rate(pl.StepTimes(100.0, 50.0, 20.0), 0.5)
    assert r_all < r_chem


def test_availability_constraint() -> None:
    times = (pl.StepTimes(1e4, 0, 0), pl.StepTimes(10.0, 0, 0), pl.StepTimes(10.0, 0, 0))
    r = pl.step_rates(times, np.array([0.2, 0.2, 0.2]))
    assert r[1] <= r[0] + 1e-15 and r[2] <= r[1] + 1e-15


# ----------------------------------------------------------------------------- simulation


def test_simulation_bounds_and_ordering(run900: pl.PelletResult) -> None:
    X = run900.X_steps
    assert np.all((X >= 0) & (X <= 1))
    assert np.all(X[1] <= X[0] + 1e-6) and np.all(X[2] <= X[1] + 1e-6)
    assert np.all(np.diff(run900.X) >= -1e-9)
    assert run900.X[-1] > 0.99


def test_oxygen_mass_balance(base, run900: pl.PelletResult) -> None:  # type: ignore[no-untyped-def]
    """H2 consumed (= H2O produced) integrated over time equals the oxygen removed."""
    pp, _ = base
    n0 = pl.initial_moles_fe2o3(pp)
    h2_consumed = np.trapezoid(run900.o_rate, run900.t)
    assert h2_consumed == pytest.approx(pl.O_PER_FE2O3 * n0 * run900.X[-1], rel=2e-3)


def test_gas_balance_identity(base, td, run900: pl.PelletResult) -> None:  # type: ignore[no-untyped-def]
    """'mean' mode: x_H2 seen = x_in - 0.5 * consumption / total molar flow."""
    _, cond = base
    f_tot = td.normal_pressure * cond.gas_flow_normal / (pr.R_GAS * td.normal_temperature)
    np.testing.assert_allclose(
        run900.x_H2_bulk, cond.x_H2 - 0.5 * run900.o_rate / f_tot, rtol=1e-6, atol=1e-9
    )


def test_depletion_modes_ordered(base, species, td, refit) -> None:  # type: ignore[no-untyped-def]
    pp, cond = base
    t50 = {}
    for mode in ("static", "mean", "outlet"):
        r = pl.simulate_tga(pp, cond, T900, 1800.0, species, td, gas_mode=mode, triplets=refit)
        t50[mode] = float(np.interp(0.5, r.X, r.t))
    assert t50["static"] < t50["mean"] < t50["outlet"]


def test_no_reduction_without_hydrogen(base, species, td, refit) -> None:  # type: ignore[no-untyped-def]
    pp, cond = base
    r = pl.simulate_tga(
        pp, dataclasses.replace(cond, x_H2=0.0), T900, 600.0, species, td, triplets=refit
    )
    assert np.all(r.X == 0.0)


# ----------------------------------------------------------------------------- validation

FIG6 = TEST_DATA_DIR / "pellet_reduction_marsigny2025_fig6.csv"
FIG9 = TEST_DATA_DIR / "pellet_reduction_marsigny2025_fig9.csv"
REPEAT = TEST_DATA_DIR / "pellet_t50_t90_marsigny2025.csv"


def _rmse(res: pl.PelletResult, t_min: np.ndarray, x: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.interp(t_min * 60, res.t, res.X) - x) ** 2)))


def test_repeatability_times_within_two_sigma(run900: pl.PelletResult) -> None:
    """HELD-OUT: t50 and t90 at 900 degC (5 runs, mean +- std) were not used in the fit."""
    data = pd.read_csv(REPEAT, comment="#")
    cols = zip(
        data.conversion.to_numpy(float),
        data.time_min.to_numpy(float),
        data.std_min.to_numpy(float),
        strict=True,
    )
    for x, t_exp, std in cols:
        t_model = float(np.interp(x, run900.X, run900.t)) / 60
        assert abs(t_model - t_exp) <= 2 * std, (x, t_model)


@pytest.mark.parametrize("T_c", [600, 650, 700, 750, 800, 850, 900])
def test_fit_quality_fig6(base, species, td, refit, T_c: int) -> None:  # type: ignore[no-untyped-def]
    """FIT QUALITY (these curves were used to fit the triplets): RMSE <= 0.03 in X.

    0.03 is about the scatter between repeated runs: at 900 degC the run in Fig. 6 and the
    nominally identical run in Fig. 9 differ by up to 0.04 in X. Exception, 0.04 at 700 and
    750 degC: Marsigny 2025 (Sec. 3.1) report a kinetic slowdown above X = 0.6 there that
    their own model also over-predicts; it is not part of either model.
    """
    pp, cond = base
    d = pd.read_csv(FIG6)
    g = d[d.temperature_C == T_c]
    r = pl.simulate_tga(pp, cond, T_c + 273.15, 150 * 60.0, species, td, triplets=refit, n_out=900)
    tol = 0.04 if T_c in (700, 750) else 0.03
    assert _rmse(r, g.time_min.to_numpy(), g.conversion.to_numpy()) <= tol


@pytest.mark.parametrize("x_h2", [0.20, 0.55, 0.80, 1.00])
def test_prediction_fig9_dilution(base, species, td, refit, x_h2: float) -> None:  # type: ignore[no-untyped-def]
    """HELD-OUT: H2-He dilution series at 900 degC (Fig. 9), not used in the fit.

    Tolerance RMSE <= 0.06 in X: twice the run-to-run scatter above. The measured RMSE values
    (0.02-0.05, docs/validation.md) show the model is somewhat too fast at every H2 fraction.
    """
    pp, cond = base
    d = pd.read_csv(FIG9)
    g = d[np.isclose(d.x_H2_in, x_h2)]
    c = dataclasses.replace(cond, x_H2=x_h2)
    r = pl.simulate_tga(pp, c, T900, 40 * 60.0, species, td, triplets=refit, n_out=800)
    assert _rmse(r, g.time_min.to_numpy(), g.conversion.to_numpy()) <= 0.06


def test_published_triplets_known_deviation(base, species, td) -> None:  # type: ignore[no-untyped-def]
    """Guards docs/validation.md: with this model structure the published triplets
    (Marsigny 2025, Table 2) are far too slow at 600 degC."""
    pp, cond = base
    d = pd.read_csv(FIG6)
    g = d[d.temperature_C == 600]
    r = pl.simulate_tga(pp, cond, 873.15, 150 * 60.0, species, td, n_out=900)
    assert _rmse(r, g.time_min.to_numpy(), g.conversion.to_numpy()) > 0.2
