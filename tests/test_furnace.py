"""Tests for the whole-furnace model (segments, cone, bustle side inlet) and the comparison
with REDUCTOR case A (Milestone 5)."""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pandas as pd
import pytest

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import solver as sv
from shaftsim import transport as tr
from shaftsim import validation as va

from .conftest import DATA_DIR, TEST_DATA_DIR

REF_OUT = TEST_DATA_DIR / "reference_model_marsigny2025b_caseA.csv"
REF_PROF = TEST_DATA_DIR / "reference_model_marsigny2025b_caseA_profiles.csv"


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
def zone() -> en.ThermalCase:
    return en.load_thermal_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")


@pytest.fixture(scope="module")
def furnace() -> en.ThermalCase:
    return en.load_furnace_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")


@pytest.fixture(scope="module")
def furnace_solution(furnace, species, td, trip):  # type: ignore[no-untyped-def]
    prob = en.ThermalProblem(furnace, species, td, trip, 50)
    return prob, en.solve_thermal(prob)


def _split_bottom_inlet(zone: en.ThermalCase, bottom_flow: float) -> en.ThermalCase:
    """The same gas, split into a bottom inlet and a side inlet at the bottom face."""
    s = zone.shaft
    side = sh.SideInlet(s.height, s.gas_flow - bottom_flow, dict(s.x_in), zone.T_gas_in)
    return dataclasses.replace(
        zone, shaft=dataclasses.replace(s, gas_flow=bottom_flow, side_inlets=(side,))
    )


# ----------------------------------------------------------------------------- grid


def test_grid_geometry_and_flows(furnace: en.ThermalCase) -> None:
    case = furnace.shaft
    gr = sh.make_grid(case, 100)
    segs = case.segments
    assert segs is not None
    # segment boundaries (and the bustle) on cell faces
    for depth in np.cumsum([sg.length for sg in segs]):
        assert np.min(np.abs(gr.z_faces - depth)) < 1e-9
    # bed volume: cylinders + frustum, exact
    v_cyl = math.pi * segs[0].r_top ** 2 * (segs[0].length + segs[1].length)
    c = segs[2]
    v_cone = math.pi * c.length * (c.r_top**2 + c.r_top * c.r_bottom + c.r_bottom**2) / 3.0
    assert float(np.sum(gr.area * gr.dz)) == pytest.approx(v_cyl + v_cone, rel=1e-12)
    t_res = (1 - case.bed_porosity) * (v_cyl + v_cone) / case.solid_volume_flow
    assert float(np.sum(gr.kappa)) == pytest.approx(t_res, rel=1e-12)
    # flows: bottom 500 mol/s below the bustle, 2800 above; same composition everywhere
    bustle = gr.side_cells[0]
    assert gr.z_faces[bustle + 1] == pytest.approx(case.side_inlets[0].depth)
    np.testing.assert_allclose(gr.F[bustle + 1 :], case.gas_flow)
    np.testing.assert_allclose(gr.F[: bustle + 1], case.total_gas_flow)
    np.testing.assert_allclose(gr.g, case.x_in["H2"] + case.x_in["H2O"], rtol=1e-14)


def test_grid_rejects_side_inlet_off_boundary(furnace: en.ThermalCase) -> None:
    s = furnace.shaft
    bad = dataclasses.replace(s, side_inlets=(dataclasses.replace(s.side_inlets[0], depth=3.0),))
    with pytest.raises(ValueError, match="segment boundary"):
        sh.make_grid(bad, 50)


def test_single_segment_grid_is_uniform(zone: en.ThermalCase) -> None:
    gr = sh.make_grid(zone.shaft, 40)
    np.testing.assert_allclose(gr.dz, zone.shaft.height / 40, rtol=1e-12)
    np.testing.assert_allclose(gr.F, zone.shaft.gas_flow)
    np.testing.assert_allclose(gr.kappa, (zone.shaft.height / 40) / zone.shaft.solid_velocity)


# ----------------------------------------------------------------------------- equivalences


def test_side_inlet_at_bottom_equals_single_inlet(zone, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """2300 mol/s injected at the bottom face as a side inlet + 500 mol/s bottom inlet gives
    the same discrete solution as 2800 mol/s at the bottom (species and energy)."""
    a = en.solve_thermal(en.ThermalProblem(zone, species, td, trip, 30))
    split = _split_bottom_inlet(zone, 500.0)
    b = en.solve_thermal(en.ThermalProblem(split, species, td, trip, 30))
    assert a.converged and b.converged
    np.testing.assert_allclose(b.Ts, a.Ts, atol=1e-6)
    np.testing.assert_allclose(b.Tg, a.Tg, atol=1e-6)
    np.testing.assert_allclose(b.species.X_steps, a.species.X_steps, atol=1e-8)
    np.testing.assert_allclose(b.species.x_H2, a.species.x_H2, atol=1e-8)


def test_split_cylinder_equals_single_cylinder(zone) -> None:  # type: ignore[no-untyped-def]
    """Two cylinder segments of the same radius: same grid, same isothermal solution."""
    s = dataclasses.replace(zone.shaft, T=1100.0)
    r = math.sqrt(s.area / math.pi)
    two = dataclasses.replace(
        s, segments=(sh.Segment(s.height / 2, r, r), sh.Segment(s.height / 2, r, r))
    )
    g1, g2 = sh.make_grid(s, 40), sh.make_grid(two, 40)
    for k in ("dz", "area", "kappa", "F", "g"):
        np.testing.assert_allclose(getattr(g2, k), getattr(g1, k), rtol=1e-12)


def test_isothermal_side_inlet_species(zone, species, td, trip) -> None:  # type: ignore[no-untyped-def]
    """Isothermal species problem with a side inlet half-way: element balances close and
    the gas below the inlet carries only the bottom flow."""
    s = dataclasses.replace(zone.shaft, T=1100.0)
    r = math.sqrt(s.area / math.pi)
    side = sh.SideInlet(s.height / 2, 2300.0, dict(s.x_in), 1100.0)
    s2 = dataclasses.replace(
        s,
        gas_flow=500.0,
        segments=(sh.Segment(s.height / 2, r, r), sh.Segment(s.height / 2, r, r)),
        side_inlets=(side,),
    )
    co = sh.shaft_coefficients(s, species, td, trip)
    prob = sh.ShaftProblem(s2, co, 40)
    y, st = sv.solve(prob)
    assert st.converged
    sol = sh.make_solution(prob, y, td, st.residuals[-1], st.iterations)
    for err in sh.element_balances(prob, sol).values():
        assert abs(err) < 1e-12
    assert 0.0 < sol.metallization <= 1.0
    # H2 consumed equals oxygen removed from the solid
    o_removed = s2.removable_oxygen_flow * float(sol.X[-1])
    assert float(np.sum(prob.oxygen_removed(prob.unpack(y)[1]))) == pytest.approx(o_removed)


# ----------------------------------------------------------------------------- whole furnace


def test_whole_furnace_balances(furnace_solution) -> None:  # type: ignore[no-untyped-def]
    prob, sol = furnace_solution
    assert sol.converged
    assert abs(sol.energy_closure) < 1e-10
    assert abs(sol.energy_closure_sensible) < 1e-10
    sprob = prob.species_problem(sol.Ts, sol.Tg)
    for err in sh.element_balances(sprob, sol.species).values():
        assert abs(err) < 1e-12


def test_whole_furnace_physics(furnace_solution, furnace) -> None:  # type: ignore[no-untyped-def]
    """Qualitative checks: the cold cooling gas cools the DRI well below the bustle gas
    temperature, the burden is heated above the inlet, and metallisation is substantial."""
    _, sol = furnace_solution
    assert furnace.T_gas_in < sol.T_solid_out < 0.5 * (furnace.T_gas_in + 1163.15)
    assert sol.T_top_gas > furnace.T_solid_in
    assert np.max(sol.Ts) < furnace.shaft.side_inlets[0].T
    assert sol.species.metallization > 0.8


# ----------------------------------------------------------------------------- mass fractions


def test_solid_mass_fractions(zone: en.ThermalCase) -> None:
    pp = zone.shaft.pellet
    X = np.array([[0.0, 1.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]])
    w = pl.solid_mass_fractions(pp, X)
    np.testing.assert_allclose(sum(w.values()), 1.0, rtol=1e-14)
    assert w["Fe2O3"][0] == pytest.approx(pp.w_Fe2O3)
    assert w["Fe3O4"][1] + w["gangue"][1] == pytest.approx(1.0)
    assert w["Fe0.95O"][2] + w["gangue"][2] == pytest.approx(1.0)
    # full metallisation: Fe / (Fe + gangue) from the hematite content
    fe = 2 * pp.w_Fe2O3 / pp.molar_mass["Fe2O3"] * pp.molar_mass["Fe"]
    assert w["Fe"][3] == pytest.approx(fe / (fe + 1 - pp.w_Fe2O3))


# ----------------------------------------------------------------------------- references


def test_reference_table1_consistent(species: dict[str, pr.Species]) -> None:
    """Table 1 (as extracted): DRI composition adds up and gives the reported metallisation."""
    ref = va.load_reference_outputs(REF_OUT)
    w = {k: ref[f"dri_w_{k}"] for k in ("Fe2O3", "Fe3O4", "FeO", "Fe", "gangue")}
    assert sum(w.values()) == pytest.approx(100.0, abs=0.05)
    m = {k: float(species[k].molar_mass or np.nan) for k in ("Fe2O3", "Fe3O4", "FeO", "Fe")}
    fe_per = {"Fe2O3": 2.0, "Fe3O4": 3.0, "FeO": 1.0, "Fe": 1.0}
    fe_tot = sum(w[k] / m[k] * fe_per[k] for k in m)
    assert 100 * w["Fe"] / m["Fe"] / fe_tot == pytest.approx(ref["metallization"], abs=0.1)


def test_reference_profiles_sane(furnace: en.ThermalCase) -> None:
    """Digitised maps: detected edge radius follows the geometry, most pixels classified."""
    d = va.load_reference_profiles(REF_PROF)
    assert set(d.quantity) == {"T_solid", "x_H2", "w_Fe2O3", "w_Fe3O4", "w_FeO", "w_Fe"}
    segs = furnace.shaft.segments
    assert segs is not None
    cone = segs[2]
    z = d.z.to_numpy()
    r_geo = np.where(z >= 0, cone.r_top, cone.r_top + (cone.r_bottom - cone.r_top) * (-z / 10))
    assert np.max(np.abs(d.r_edge.to_numpy() - r_geo)) < 0.12
    assert d.classified.min() > 0.8
    assert np.all(d.mean_lo <= d["mean"]) and np.all(d["mean"] <= d.mean_hi)


def test_compare_profile_band_logic() -> None:
    ref = pd.DataFrame(
        {"quantity": "q", "z": [0.0, 1.0, 2.0], "mean": 1.0, "mean_lo": 0.5, "mean_hi": 1.5}
    )
    model = pd.DataFrame({"z": [0.0, 1.0, 2.0], "q": [1.0, 2.0, 0.0]})
    c = va.compare_profile(ref, model, "q", (0.0, 2.0))
    assert c.n_points == 3
    assert c.within_band == pytest.approx(1 / 3)
    assert c.mean_excess == pytest.approx((0.0 + 0.5 + 0.5) / 3)
    assert c.max_excess == pytest.approx(0.5)
    assert c.rms_midpoint == pytest.approx(math.sqrt(2 / 3))


def test_model_outputs_consistent(furnace_solution) -> None:  # type: ignore[no-untyped-def]
    """Model outputs in Table 1 units: mass fractions add to 100 %, metallisation follows."""
    prob, sol = furnace_solution
    out = va.model_outputs(prob, sol)
    total = sum(out[f"dri_w_{k}"] for k in ("Fe2O3", "Fe3O4", "FeO", "Fe", "gangue"))
    assert total == pytest.approx(100.0, rel=1e-12)
    assert out["top_gas_flow"] == pytest.approx(2800.0)
    assert out["top_gas_x_H2"] + out["top_gas_x_H2O"] + out["top_gas_x_N2"] == pytest.approx(100)
    assert out["metallization"] == pytest.approx(100 * sol.species.metallization)


def test_regression_case_a_50_cells(furnace_solution) -> None:  # type: ignore[no-untyped-def]
    """Pins the 50-cell row of the grid study in docs/validation.md (Milestone 5): if this
    fails, the model changed and the documented comparison must be regenerated."""
    prob, sol = furnace_solution
    out = va.model_outputs(prob, sol)
    assert out["metallization"] == pytest.approx(94.56, abs=0.05)
    assert out["top_gas_T"] == pytest.approx(363.0, abs=0.3)
    assert out["dri_T"] == pytest.approx(268.3, abs=0.3)
    assert out["gud"] == pytest.approx(0.3165, abs=5e-4)
