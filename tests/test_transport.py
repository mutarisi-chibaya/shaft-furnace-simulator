"""Tests for shaftsim.transport."""

from __future__ import annotations

import math

import pytest

from shaftsim import properties as pr
from shaftsim import transport as tr

from .conftest import DATA_DIR


@pytest.fixture(scope="module")
def td() -> tr.TransportData:
    sp = pr.load_species(DATA_DIR / "species.yaml")
    mm = {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None}
    return tr.load_transport(DATA_DIR / "transport.yaml", DATA_DIR / "collision_integrals.csv", mm)


P_ATM = 101325.0


# ----------------------------------------------------------------------------- collision integrals


def test_collision_table_nodes_reproduced(td: tr.TransportData) -> None:
    ci = td.collision
    for t, o11, o22 in zip(ci.tstar, ci.omega11, ci.omega22, strict=True):
        assert ci.omega11_at(t) == pytest.approx(o11, rel=1e-9)
        assert ci.omega22_at(t) == pytest.approx(o22, rel=1e-9)


def test_collision_extrapolation_continuous_and_bounded(td: tr.TransportData) -> None:
    ci = td.collision
    assert ci.omega22_at(100.0 * (1 + 1e-9)) == pytest.approx(ci.omega22_at(100.0), rel=1e-6)
    assert ci.omega22_at(300.0) < ci.omega22_at(100.0)
    with pytest.raises(ValueError):
        ci.omega22_at(1000.0)
    with pytest.raises(ValueError):
        ci.omega11_at(0.05)


# ----------------------------------------------------------------------------- diffusivities


def test_chapman_enskog_symmetric_and_inverse_pressure(td: tr.TransportData) -> None:
    d1 = tr.chapman_enskog_diffusivity(td, "H2", "H2O", 1000.0, P_ATM)
    assert tr.chapman_enskog_diffusivity(td, "H2O", "H2", 1000.0, P_ATM) == pytest.approx(d1)
    assert tr.chapman_enskog_diffusivity(td, "H2", "H2O", 1000.0, 2 * P_ATM) == pytest.approx(
        d1 / 2
    )


def test_chapman_enskog_matches_classic_formula(td: tr.TransportData) -> None:
    """Same as D = 0.0018583 T^1.5 sqrt(1/M1+1/M2) / (p sigma^2 Omega) [cm2/s, atm, A, g/mol]
    (independent check of the constants in the kinetic-theory expression, to 0.1 %)."""
    T = 1000.0
    a, b = td.lj["H2"], td.lj["H2O"]
    sigma_a = 0.5 * (a.sigma + b.sigma) * 1e10
    om = td.collision.omega11_at(T / math.sqrt(a.eps_over_k * b.eps_over_k))
    ma, mb = td.molar_mass["H2"] * 1e3, td.molar_mass["H2O"] * 1e3
    d_cgs = 0.0018583 * T**1.5 * math.sqrt(1 / ma + 1 / mb) / (sigma_a**2 * om)
    assert tr.chapman_enskog_diffusivity(td, "H2", "H2O", T, P_ATM) * 1e4 == pytest.approx(
        d_cgs, rel=1e-3
    )


def test_marrero_mason_equation(td: tr.TransportData) -> None:
    """Eq. (4.3-1) evaluated by hand for He-H2 at 300 K."""
    T = 300.0
    pd = 27.0e-3 * T**1.510 / math.log(0.0534e8 / T) ** 2  # atm cm2/s
    assert tr.binary_diffusivity(td, "He", "H2", T, P_ATM) == pytest.approx(pd * 1e-4, rel=1e-12)


@pytest.mark.parametrize(
    ("pair", "T", "max_dev"),
    # Measured deviation of Chapman-Enskog (Svehla force constants) from the experimental
    # Marrero-Mason correlations; frozen here so that docs/validation.md stays true.
    [
        (("He", "H2"), 300.0, 0.05),
        (("He", "H2"), 1000.0, 0.06),
        (("H2", "N2"), 300.0, 0.05),
        (("H2", "N2"), 1000.0, 0.12),
        (("H2", "N2"), 1500.0, 0.15),
    ],
)
def test_chapman_enskog_vs_marrero_mason(
    td: tr.TransportData, pair: tuple[str, str], T: float, max_dev: float
) -> None:
    ce = tr.chapman_enskog_diffusivity(td, *pair, T, P_ATM)
    mm = tr.binary_diffusivity(td, *pair, T, P_ATM)
    assert abs(ce / mm - 1) < max_dev


def test_correlation_range_enforced(td: tr.TransportData) -> None:
    with pytest.raises(ValueError):
        tr.binary_diffusivity(td, "He", "H2", 20000.0, P_ATM)


def test_blanc_limits(td: tr.TransportData) -> None:
    T = 1173.15
    binary = tr.binary_diffusivity(td, "H2", "He", T, P_ATM)
    assert tr.mixture_diffusivity(td, "H2", {"H2": 0.3, "He": 0.7}, T, P_ATM, "H2O") == (
        pytest.approx(binary)
    )
    pure = tr.mixture_diffusivity(td, "H2", {"H2": 1.0}, T, P_ATM, "H2O")
    assert pure == pytest.approx(tr.binary_diffusivity(td, "H2", "H2O", T, P_ATM))
    mix = tr.mixture_diffusivity(td, "H2", {"H2": 0.5, "He": 0.25, "H2O": 0.25}, T, P_ATM, "H2O")
    lo, hi = sorted([binary, tr.binary_diffusivity(td, "H2", "H2O", T, P_ATM)])
    assert lo <= mix <= hi


# ----------------------------------------------------------------------------- viscosity


def test_wilke_pure_limit(td: tr.TransportData) -> None:
    assert tr.mixture_viscosity(td, {"N2": 1.0}, 900.0) == pytest.approx(
        tr.viscosity(td, "N2", 900.0)
    )


def test_wilke_mixture_between_pure_values(td: tr.TransportData) -> None:
    T = 1173.15
    mix = tr.mixture_viscosity(td, {"H2": 0.55, "He": 0.45}, T)
    lo, hi = sorted([tr.viscosity(td, "H2", T), tr.viscosity(td, "He", T)])
    assert lo < mix < hi


def test_viscosity_increases_with_temperature(td: tr.TransportData) -> None:
    for s in ("H2", "H2O", "He", "N2"):
        assert tr.viscosity(td, s, 1200.0) > tr.viscosity(td, s, 600.0)


# ----------------------------------------------------------------------------- pores, external


def test_knudsen_formula() -> None:
    d, T, M = 1e-7, 1000.0, 2.016e-3
    assert tr.knudsen_diffusivity(d, T, M) == pytest.approx(
        d / 3 * math.sqrt(8 * pr.R_GAS * T / (math.pi * M))
    )


def test_bosanquet_limits() -> None:
    assert tr.effective_diffusivity(0.4, 2.0, 1e-4, 1e9) == pytest.approx(0.2 * 1e-4, rel=1e-9)
    assert tr.effective_diffusivity(0.4, 2.0, 1e9, 1e-6) == pytest.approx(0.2 * 1e-6, rel=1e-9)


def test_ranz_marshall_stagnant_limit(td: tr.TransportData) -> None:
    """At zero velocity Sh = 2, so k_g = 2 D / d."""
    kg = tr.ranz_marshall_kg(td.ranz_marshall, 0.012, 0.0, 0.03, 4e-5, 1e-3)
    assert kg == pytest.approx(2 * 1e-3 / 0.012)
