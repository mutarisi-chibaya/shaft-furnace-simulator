"""Tests for shaftsim.uq and the parameter factors of ThermalProblem (Milestone 6)."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr
from shaftsim import uq

from .conftest import DATA_DIR


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


@pytest.fixture(scope="module")
def furnace() -> en.ThermalCase:
    return en.load_furnace_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")


@pytest.fixture(scope="module")
def spec(furnace: en.ThermalCase) -> uq.UncertaintySpec:
    return uq.load_uncertainty(
        DATA_DIR / "uncertainty.yaml",
        DATA_DIR / "kinetics_refit.yaml",
        furnace.shaft.pellet.tortuosity,
    )


@pytest.fixture(scope="module")
def base(furnace, species, td, trip) -> en.ThermalProblem:  # type: ignore[no-untyped-def]
    return en.ThermalProblem(furnace, species, td, trip, 20)


# ----------------------------------------------------------------------------- factors


def _max_rel(a, b) -> float:  # type: ignore[no-untyped-def]
    return float(
        max(
            abs(getattr(x[j], k) / getattr(y[j], k) - 1)
            for x, y in zip(a, b, strict=True)
            for j in range(3)
            for k in ("a_chem", "a_diff", "a_ext", "ratio")
        )
    )


def test_kinetic_factor_matches_rebuilt_tables(base, furnace, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    new = tuple(
        dataclasses.replace(t, k0=t.k0 * f, Ea=t.Ea + d)
        for t, f, d in zip(trip, (1.3, 0.8, 1.1), (3000.0, -2000.0, 1500.0), strict=True)
    )
    Ts = np.linspace(400.0, 1150.0, 20)
    a = base.with_parameters(kinetics=new).cell_coefficients(Ts, Ts + 15)
    b = en.ThermalProblem(furnace, species, td, new, 20).cell_coefficients(Ts, Ts + 15)
    assert _max_rel(a, b) < 1e-5  # log-linear table interpolation


def test_tortuosity_factor_matches_rebuilt_tables(base, furnace, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    pp = dataclasses.replace(furnace.shaft.pellet, tortuosity=furnace.shaft.pellet.tortuosity * 1.4)
    case2 = dataclasses.replace(furnace, shaft=dataclasses.replace(furnace.shaft, pellet=pp))
    Ts = np.linspace(400.0, 1150.0, 20)
    a = base.with_parameters(tortuosity_factor=1.4).cell_coefficients(Ts, Ts + 15)
    b = en.ThermalProblem(case2, species, td, trip, 20).cell_coefficients(Ts, Ts + 15)
    assert _max_rel(a, b) < 1e-12


def test_transfer_mix_endpoints(base, furnace, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """lambda = 1 equals tables built with Ranz-Marshall in place of Wakao-Kaguei."""
    tdr = dataclasses.replace(td, wakao_kaguei=td.ranz_marshall)
    ref = en.ThermalProblem(furnace, species, tdr, trip, 20)
    Tg = np.linspace(350.0, 1160.0, 20)
    one = base.with_parameters(transfer_mix=1.0)
    np.testing.assert_allclose(one.mass_transfer_kg(Tg), ref.mass_transfer_kg(Tg), rtol=1e-12)
    np.testing.assert_allclose(one.heat_transfer_UA(Tg), ref.heat_transfer_UA(Tg), rtol=1e-12)
    half = base.with_parameters(transfer_mix=0.5)
    np.testing.assert_allclose(
        half.heat_transfer_UA(Tg),
        0.5 * (base.heat_transfer_UA(Tg) + ref.heat_transfer_UA(Tg)),
        rtol=1e-12,
    )


def test_with_parameters_rejects_table_fields(base: en.ThermalProblem) -> None:
    with pytest.raises(ValueError, match="property tables"):
        base.with_parameters(triplets=None)
    assert base.with_parameters(tortuosity_factor=2.0)._tables is base._tables


# ----------------------------------------------------------------------------- sampling


def test_theta_triplet_round_trip(spec: uq.UncertaintySpec, trip: pl.Triplets) -> None:
    back = uq.triplets_from_theta(spec.theta_mean)
    for a, b in zip(back, trip, strict=True):
        assert a.k0 == pytest.approx(b.k0, rel=1e-12)
        assert a.Ea == pytest.approx(b.Ea, rel=1e-9)


def test_latin_hypercube_stratified_and_reproducible() -> None:
    u = uq.latin_hypercube(50, seed=7)
    assert u.shape == (50, len(uq.NAMES))
    for k in range(u.shape[1]):  # one point per stratum in every dimension
        assert sorted(np.floor(u[:, k] * 50).astype(int)) == list(range(50))
    np.testing.assert_array_equal(u, uq.latin_hypercube(50, seed=7))
    assert not np.array_equal(u, uq.latin_hypercube(50, seed=8))


def test_samples_stay_in_their_ranges(spec: uq.UncertaintySpec) -> None:
    """PLAN.md: sampled parameters stay inside their ranges (kinetics within the truncation
    along the principal axes)."""
    s, v = spec.principal_axes
    for u in uq.latin_hypercube(400, seed=1):
        p = uq.to_parameters(spec, u)
        assert spec.tortuosity_range[0] <= p.tortuosity <= spec.tortuosity_range[1]
        assert spec.transfer_range[0] <= p.transfer_mix <= spec.transfer_range[1]
        z = (v.T @ (p.theta - spec.theta_mean))[s > 0] / s[s > 0]
        assert np.all(np.abs(z) <= spec.truncation + 1e-9)
    with pytest.raises(ValueError):
        uq.to_parameters(spec, np.full(len(uq.NAMES), 1.5))


def test_centre_of_cube_is_nominal_kinetics(spec: uq.UncertaintySpec) -> None:
    p = uq.to_parameters(spec, np.full(len(uq.NAMES), 0.5))
    np.testing.assert_allclose(p.theta, spec.theta_mean, rtol=1e-12)


# ----------------------------------------------------------------------------- runs


def test_zero_uncertainty_collapses_to_deterministic(base, spec) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: with zero uncertainty every Monte Carlo run equals the deterministic result."""
    det = uq.run_one(base)
    zero = spec.zero()
    res = uq.run_points(base, zero, uq.latin_hypercube(3, seed=0), workers=1)
    for r in res:
        assert r.converged
        for k in uq.OUTPUTS:
            assert r.outputs[k] == pytest.approx(det.outputs[k], rel=1e-10, abs=1e-10)


def test_runs_reproducible_with_seed_in_parallel(base, spec) -> None:  # type: ignore[no-untyped-def]
    """PLAN.md: results are reproducible with a fixed seed (also across processes)."""
    pts = uq.latin_hypercube(4, seed=3)
    a = uq.output_matrix(uq.run_points(base, spec, pts, workers=2))
    b = uq.output_matrix(uq.run_points(base, spec, uq.latin_hypercube(4, seed=3), workers=1))
    np.testing.assert_array_equal(a, b)
    assert np.ptp(a[:, 0]) > 0  # the samples do differ


# ----------------------------------------------------------------------------- statistics


def test_summary_and_bands() -> None:
    v = np.arange(101, dtype=float)
    s = uq.summary(v)
    assert s["mean"] == pytest.approx(50.0)
    assert (s["p05"], s["p50"], s["p95"]) == pytest.approx((5.0, 50.0, 95.0))
    assert s["std"] == pytest.approx(np.std(v, ddof=1))
    prof = np.vstack([v, 2 * v]).T  # 101 runs, 2 cells
    b = uq.quantile_bands(prof)
    np.testing.assert_allclose(b[:, 0], [5.0, 50.0, 95.0])
    np.testing.assert_allclose(b[:, 1], [10.0, 100.0, 190.0])


def test_running_statistics_converge_to_full_sample() -> None:
    rng = np.random.default_rng(0)
    v = rng.normal(size=1000)
    r = uq.running_statistics(v, [10, 100, 1000])
    assert r["n"].tolist() == [10.0, 100.0, 1000.0]
    assert r["mean"][-1] == pytest.approx(uq.summary(v)["mean"])
    assert r["p05"][-1] < r["mean"][-1] < r["p95"][-1]
