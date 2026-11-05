"""Electricity demand and CO2 per tonne of DRI (Milestone 7).

Scope, as in the gas-loop model of Marsigny et al. (2025b, Sec. 2.2), whose parameters are
used (data/emissions.yaml):

- **Water electrolysis** supplies the fresh H2, equal to the H2 consumed in the furnace
  (the rest of the gas is recycled). Electricity: 18 MJ per normal cubic metre of H2.
- **Gas heating**: the bustle gas is heated from the condenser exit temperature (70 degC) to
  the bustle temperature by a heater of efficiency 0.7, assumed electric. Heat recovery
  from the top gas (the heat exchanger of the paper's gas loop) is NOT counted, so this term
  is an upper bound. The cooling gas enters at the condenser temperature and is not heated.
- **Not included**: gas compression (about 2 % of the total in the paper), pellet making, the
  EAF, and the emissions of building the plant and the electrolysers.

CO2 = electricity per kg DRI x grid intensity (life-cycle, data/grid_intensity.yaml).
Everything is per kg of DRI in SI units (J/kg, kg CO2/kg). For display: 1 GJ/t = 1 MJ/kg =
1e6 J/kg, and 1 kg CO2/kg = 1 t CO2/t.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import yaml

from . import energy as en
from . import properties as pr


@dataclass(frozen=True)
class GridIntensity:
    """Life-cycle CO2 intensity of a power grid.

    Attributes:
        name: Area name.
        year: Year of the value.
        value: Intensity [kg CO2e per J of electricity].
    """

    name: str
    year: int
    value: float

    @property
    def g_per_kWh(self) -> float:
        """The same in g CO2e/kWh (display unit)."""
        return self.value * 3.6e9


@dataclass(frozen=True)
class EmissionsData:
    """Parameters of the energy and CO2 accounting (see the module docstring).

    Attributes:
        electrolysis_per_mol: Electricity per mol of H2 produced [J/mol].
        heater_efficiency: Efficiency of the bustle gas heater [-].
        recycle_T: Temperature of the recycled gas before heating [K].
        grids: Grid intensities by key (south_africa, world, france).
        benchmarks: Benchmark intensities by key [kg CO2/kg] (different system boundaries,
            see data/emissions.yaml).
    """

    electrolysis_per_mol: float
    heater_efficiency: float
    recycle_T: float
    grids: Mapping[str, GridIntensity]
    benchmarks: Mapping[str, float]


def normal_molar_volume(T_n: float, p_n: float) -> float:
    """Ideal-gas molar volume at the normal state [m3/mol] (0.0224 at 273.15 K, 101325 Pa)."""
    return pr.R_GAS * T_n / p_n


def electrolysis_energy_per_mol(per_normal_m3: float, T_n: float, p_n: float) -> float:
    """Electricity per mol of H2 [J/mol] from a value per normal cubic metre [J/Nm3]."""
    return per_normal_m3 * normal_molar_volume(T_n, p_n)


def load_emissions(emissions_yaml: str | Path, grid_yaml: str | Path) -> EmissionsData:
    """Read data/emissions.yaml and data/grid_intensity.yaml."""
    with open(emissions_yaml, encoding="utf-8") as fh:
        e = yaml.safe_load(fh)
    with open(grid_yaml, encoding="utf-8") as fh:
        g = yaml.safe_load(fh)["areas"]

    def v(section: str, key: str) -> float:
        return float(e[section][key]["value"])

    return EmissionsData(
        electrolysis_per_mol=electrolysis_energy_per_mol(
            v("electrolysis", "specific_energy"),
            v("normal_conditions", "temperature"),
            v("normal_conditions", "pressure"),
        ),
        heater_efficiency=v("gas_heating", "heater_efficiency"),
        recycle_T=v("gas_heating", "recycle_temperature"),
        grids={
            k: GridIntensity(str(a["name"]), int(a["latest"]["year"]), float(a["latest"]["value"]))
            for k, a in g.items()
        },
        benchmarks={k: float(b["value"]) for k, b in e["benchmarks"].items()},
    )


def gas_heating_power(
    species: Mapping[str, pr.Species],
    x: Mapping[str, float],
    flow: float,
    T_from: float,
    T_to: float,
) -> float:
    """Heat needed to bring a gas stream from T_from to T_to [W].

    Args:
        species: Thermodynamic data (energy.enthalpy_ext; extended linearly below a species'
            data range, e.g. H2O below 500 K).
        x: Mole fractions [-].
        flow: Molar flow [mol/s].
        T_from, T_to: Temperatures [K].
    """
    dh = sum(
        xi * float(en.enthalpy_ext(species[s], T_to) - en.enthalpy_ext(species[s], T_from))
        for s, xi in x.items()
    )
    return flow * dh


@dataclass(frozen=True)
class EnergyDemand:
    """Electricity per kg of DRI [J/kg] (1 MJ/kg = 1 GJ/t).

    Attributes:
        electrolysis: Water electrolysis for the fresh H2 [J/kg].
        heating: Bustle gas heating (electricity, no heat recovery) [J/kg].
    """

    electrolysis: float
    heating: float

    @property
    def total(self) -> float:
        """Electrolysis + heating [J/kg]."""
        return self.electrolysis + self.heating


def energy_demand(
    h2_flow: float, dri_flow: float, heating_power: float, data: EmissionsData
) -> EnergyDemand:
    """Electricity per kg of DRI.

    Args:
        h2_flow: Fresh H2 needed = H2 consumed in the furnace [mol/s].
        dri_flow: DRI mass flow [kg/s].
        heating_power: Heat taken up by the bustle gas [W].
        data: Emissions parameters.
    """
    if dri_flow <= 0.0:
        raise ValueError("DRI flow must be positive")
    return EnergyDemand(
        electrolysis=h2_flow * data.electrolysis_per_mol / dri_flow,
        heating=heating_power / data.heater_efficiency / dri_flow,
    )


def furnace_energy_demand(
    case: en.ThermalCase,
    sol: en.ThermalSolution,
    species: Mapping[str, pr.Species],
    data: EmissionsData,
) -> EnergyDemand:
    """Electricity per kg of DRI for a solved whole-furnace case (one bustle inlet)."""
    s = sol.species
    m_h2 = species["H2"].molar_mass
    if m_h2 is None:
        raise ValueError("H2 molar mass missing")
    h2_flow = s.h2_consumption * (s.product_flow / 1000.0) / m_h2
    bustle = case.shaft.side_inlets[0]
    q = gas_heating_power(species, bustle.x, bustle.flow, data.recycle_T, bustle.T)
    return energy_demand(h2_flow, s.product_flow, q, data)


def co2_per_dri(demand: EnergyDemand, grid_intensity: float) -> float:
    """CO2 per kg of DRI [kg/kg = t/t] for a grid intensity [kg CO2e/J]."""
    return demand.total * grid_intensity
