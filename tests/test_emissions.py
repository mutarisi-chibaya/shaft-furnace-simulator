"""Tests for shaftsim.emissions (Milestone 7)."""

from __future__ import annotations

import csv

import numpy as np
import pytest
import yaml
from scipy.integrate import quad

from shaftsim import emissions as em
from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr

from .conftest import DATA_DIR, ROOT

TEST_DATA = ROOT / "tests" / "data"


@pytest.fixture(scope="module")
def species() -> dict[str, pr.Species]:
    return pr.load_species(DATA_DIR / "species.yaml")


@pytest.fixture(scope="module")
def data() -> em.EmissionsData:
    return em.load_emissions(DATA_DIR / "emissions.yaml", DATA_DIR / "grid_intensity.yaml")


def test_normal_molar_volume_and_electrolysis_energy(data: em.EmissionsData) -> None:
    assert em.normal_molar_volume(273.15, 101325.0) == pytest.approx(22.414e-3, rel=1e-4)
    # 18 MJ/Nm3 -> 403 kJ/mol -> 55.6 kWh per kg H2
    assert data.electrolysis_per_mol == pytest.approx(18e6 * 22.414e-3, rel=1e-4)
    assert data.electrolysis_per_mol / 2.01588e-3 / 3.6e6 == pytest.approx(55.6, abs=0.05)


def test_grid_intensities_match_generated_file(data: em.EmissionsData) -> None:
    raw = yaml.safe_load((DATA_DIR / "grid_intensity.yaml").read_text(encoding="utf-8"))["areas"]
    for key, g in data.grids.items():
        assert g.g_per_kWh == pytest.approx(raw[key]["latest"]["g_per_kWh"], rel=1e-6)
        assert g.year == raw[key]["latest"]["year"]
    assert data.grids["south_africa"].g_per_kWh > data.grids["world"].g_per_kWh


def test_gas_heating_equals_integral_of_cp(species: dict[str, pr.Species]) -> None:
    x = {"H2": 0.93, "H2O": 0.02, "N2": 0.05}
    q = em.gas_heating_power(species, x, 2.0, 600.0, 1163.15)
    ref = 2.0 * sum(
        xi * quad(lambda T, s=s: float(en.cp_ext(species[s], T)), 600.0, 1163.15, limit=200)[0]
        for s, xi in x.items()
    )
    assert q == pytest.approx(ref, rel=1e-4)
    assert em.gas_heating_power(species, x, 2.0, 900.0, 900.0) == 0.0


def test_energy_demand_and_co2_arithmetic(data: em.EmissionsData) -> None:
    d = em.energy_demand(h2_flow=2.0, dri_flow=4.0, heating_power=7.0e6, data=data)
    assert d.electrolysis == pytest.approx(2.0 * data.electrolysis_per_mol / 4.0)
    assert d.heating == pytest.approx(7.0e6 / data.heater_efficiency / 4.0)
    assert d.total == pytest.approx(d.electrolysis + d.heating)
    assert em.co2_per_dri(d, 0.0) == 0.0
    assert em.co2_per_dri(d, 2e-7) == pytest.approx(2 * em.co2_per_dri(d, 1e-7))
    with pytest.raises(ValueError):
        em.energy_demand(1.0, 0.0, 1.0, data)


def test_case_a_energy_demand_vs_gas_loop_model(species, data) -> None:  # type: ignore[no-untyped-def]
    """Case A (50 cells) against Table 20 of Marsigny's thesis (tests/data): electrolysis
    within 2 % (it follows from the H2 consumption); electrolysis + heating within 5 % of the
    reference total minus compression (heating here has no heat recovery, an upper bound)."""
    with open(TEST_DATA / "reference_energy_marsigny2025_thesis_t20.csv", encoding="utf-8") as fh:
        ref = {
            r["key"]: float(r["value"])
            for r in csv.DictReader(line for line in fh if not line.startswith("#"))
        }
    td = tr.load_transport(
        DATA_DIR / "transport.yaml",
        DATA_DIR / "collision_integrals.csv",
        {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None},
    )
    case = en.load_furnace_case(DATA_DIR / "furnace.yaml", DATA_DIR / "kinetics.yaml")
    trip = pl.load_refit_triplets(DATA_DIR / "kinetics_refit.yaml")
    prob = en.ThermalProblem(case, species, td, trip, 50)
    sol = en.solve_thermal(prob, tol_T=1e-4, species_tol=1e-8)
    assert sol.converged
    d = em.furnace_energy_demand(case, sol, species, data)
    assert d.electrolysis / 1e6 == pytest.approx(ref["electrolysis_GJ_per_t"], rel=0.02)
    assert d.total / 1e6 == pytest.approx(
        ref["total_GJ_per_t"] - ref["compression_GJ_per_t"], rel=0.05
    )
    assert np.isfinite(em.co2_per_dri(d, data.grids["south_africa"].value))
