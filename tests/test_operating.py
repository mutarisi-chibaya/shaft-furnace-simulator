"""Tests for shaftsim.operating (Milestone 7)."""

from __future__ import annotations

import dataclasses

import pytest

from shaftsim import energy as en
from shaftsim import operating as op

from .conftest import DATA_DIR


@pytest.fixture(scope="module")
def furnace() -> en.ThermalCase:
    return en.load_furnace_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")


@pytest.fixture(scope="module")
def cases(furnace: en.ThermalCase) -> dict[str, op.OperatingConditions]:
    return op.load_cases(DATA_DIR / "furnace_cases.yaml", furnace)


def test_mean_pressure_of_constant_profile_and_case_a(furnace: en.ThermalCase) -> None:
    assert op.mean_pressure(10.1, 25.0, 1e5, 1e5, 1e5) == pytest.approx(1e5)
    # case A: same mean as energy.load_furnace_case (174.6 kPa)
    assert op.mean_pressure(10.1, 25.0, 150e3, 180.4e3, 181.6e3) == pytest.approx(
        furnace.shaft.p, rel=1e-12
    )
    with pytest.raises(ValueError):
        op.mean_pressure(30.0, 25.0, 1e5, 1e5, 1e5)


def test_case_a_from_cases_file_equals_furnace_yaml(
    furnace: en.ThermalCase, cases: dict[str, op.OperatingConditions]
) -> None:
    a = op.conditions_of(furnace)
    b = cases["A"]
    for f in dataclasses.fields(op.OperatingConditions):
        x, y = getattr(a, f.name), getattr(b, f.name)
        if isinstance(x, dict):
            assert x == pytest.approx(y, abs=1e-12)
        else:
            assert x == pytest.approx(y, rel=1e-12), f.name


def test_apply_round_trip_and_changes(furnace: en.ThermalCase) -> None:
    oc = op.conditions_of(furnace)
    assert op.apply(furnace, oc) == furnace
    new = dataclasses.replace(
        oc, bustle_flow=2100.0, cooling_flow=700.0, solid_T=1473.15, pellet_diameter=0.012
    )
    c = op.apply(furnace, new)
    assert op.conditions_of(c) == new
    assert c.shaft.side_inlets[0].flow == 2100.0 and c.shaft.gas_flow == 700.0
    assert c.T_solid_in == 1473.15 and c.shaft.pellet.diameter == 0.012
    assert c.shaft.segments == furnace.shaft.segments  # geometry unchanged
    assert c.shaft.pellet.porosity == furnace.shaft.pellet.porosity


def test_reductor_cases_differ_as_in_table_1(cases: dict[str, op.OperatingConditions]) -> None:
    """Spot checks of the generated data/furnace_cases.yaml against Table 1."""
    a = cases["A"]
    assert (cases["B"].bustle_flow, cases["B"].cooling_flow) == (2100.0, 700.0)
    assert cases["C"].bustle_T == pytest.approx(950 + 273.15)
    assert cases["C"].x_gas["N2"] == pytest.approx(0.32)
    assert cases["D"].solid_T == pytest.approx(1200 + 273.15)
    assert cases["E"].solid_flow == 43.0
    for c in cases.values():
        assert c.bustle_flow + c.cooling_flow == 2800.0
        assert c.pellet_diameter == a.pellet_diameter


def test_invalid_conditions_rejected(furnace: en.ThermalCase) -> None:
    oc = op.conditions_of(furnace)
    with pytest.raises(ValueError, match="add up"):
        dataclasses.replace(oc, x_gas={"H2": 0.9, "H2O": 0.02, "N2": 0.05})
    with pytest.raises(ValueError, match="positive"):
        dataclasses.replace(oc, solid_flow=0.0)
    with pytest.raises(ValueError, match="species"):
        op.apply(furnace, dataclasses.replace(oc, x_gas={"H2": 0.95, "H2O": 0.05}))
