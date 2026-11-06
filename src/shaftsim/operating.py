"""Operating conditions of the whole furnace (Milestone 7): read, change and apply them.

The whole-furnace case (energy.load_furnace_case) has two gas inlets of the same composition
(bustle side inlet, cooling gas at the bottom) and one solid inlet. ``OperatingConditions``
collects the quantities that the REDUCTOR cases A-E of Marsigny et al. (2025b, Table 1) and
the sensitivity studies change. Geometry, pellet structure and kinetics stay as loaded.

Pressure: the model uses one pressure for the whole furnace (energy.load_furnace_case). It is
the length-weighted mean of the piecewise-linear profile through the top, bustle and cooling
pressures (``mean_pressure``); the bustle and cooling pressures are REDUCTOR outputs, so a
changed operating point keeps the pressures it is given (by default those of the base case).
Pressure enters only the transport properties (the equilibria are equimolar).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path

import yaml

from . import energy as en


def mean_pressure(
    zone: float, length: float, p_top: float, p_bustle: float, p_cooling: float
) -> float:
    """Length-weighted mean pressure [Pa] of the piecewise-linear profile top -> bustle ->
    cooling inlet.

    Args:
        zone: Length of the zone above the bustle [m].
        length: Total length of the furnace [m].
        p_top, p_bustle, p_cooling: Pressures at the top, the bustle and the bottom [Pa].
    """
    if not 0.0 < zone < length:
        raise ValueError("the bustle must lie inside the furnace")
    return (
        zone * 0.5 * (p_top + p_bustle) + (length - zone) * 0.5 * (p_bustle + p_cooling)
    ) / length


@dataclass(frozen=True)
class OperatingConditions:
    """Operating point of the whole furnace (SI units).

    Attributes:
        bustle_flow: Bustle gas molar flow [mol/s].
        bustle_T: Bustle gas temperature [K].
        cooling_flow: Cooling gas molar flow at the bottom [mol/s].
        cooling_T: Cooling gas temperature [K].
        x_gas: Mole fractions of both inlet gases (keys H2, H2O and the inert) [-].
        solid_flow: Burden mass flow at the top [kg/s].
        solid_T: Burden temperature at the top [K].
        pellet_diameter: Pellet diameter [m].
        p: Mean pressure of the furnace [Pa] (see the module docstring).
        pressures: (top, bustle, cooling) pressures [Pa] for the compression energy; None =
            keep those of the case the conditions are applied to.
    """

    bustle_flow: float
    bustle_T: float
    cooling_flow: float
    cooling_T: float
    x_gas: Mapping[str, float] = field(hash=False)
    solid_flow: float
    solid_T: float
    pellet_diameter: float
    p: float
    pressures: tuple[float, float, float] | None = None

    def __post_init__(self) -> None:
        if abs(sum(self.x_gas.values()) - 1.0) > 1e-9:
            raise ValueError("inlet gas mole fractions must add up to 1")
        positive = (
            self.bustle_flow,
            self.bustle_T,
            self.cooling_flow,
            self.cooling_T,
            self.solid_flow,
            self.solid_T,
            self.pellet_diameter,
            self.p,
        )
        if min(positive) <= 0.0:
            raise ValueError("flows, temperatures, diameter and pressure must be positive")


def conditions_of(case: en.ThermalCase) -> OperatingConditions:
    """The operating conditions of a whole-furnace case (one bustle side inlet)."""
    s = case.shaft
    if len(s.side_inlets) != 1:
        raise ValueError("expected a whole-furnace case with one bustle inlet")
    b = s.side_inlets[0]
    if dict(b.x) != dict(s.x_in):
        raise ValueError("bustle and cooling gas compositions differ")
    return OperatingConditions(
        bustle_flow=b.flow,
        bustle_T=b.T,
        cooling_flow=s.gas_flow,
        cooling_T=case.T_gas_in,
        x_gas=dict(s.x_in),
        solid_flow=s.solid_mass_flow,
        solid_T=case.T_solid_in,
        pellet_diameter=s.pellet.diameter,
        p=s.p,
        pressures=case.pressures,
    )


def apply(case: en.ThermalCase, oc: OperatingConditions) -> en.ThermalCase:
    """A copy of a whole-furnace case with the given operating conditions.

    The result needs a new ``ThermalProblem`` (its property tables depend on the composition,
    pressure, pellet diameter and reference gas flow).
    """
    s = case.shaft
    conditions_of(case)  # validates the structure
    if set(oc.x_gas) != set(s.x_in):
        raise ValueError(f"gas species must be {sorted(s.x_in)}")
    x = dict(oc.x_gas)
    bustle = replace(s.side_inlets[0], flow=oc.bustle_flow, x=x, T=oc.bustle_T)
    shaft = replace(
        s,
        gas_flow=oc.cooling_flow,
        x_in=x,
        solid_mass_flow=oc.solid_flow,
        pellet=replace(s.pellet, diameter=oc.pellet_diameter),
        p=oc.p,
        side_inlets=(bustle,),
    )
    return replace(
        case,
        shaft=shaft,
        T_solid_in=oc.solid_T,
        T_gas_in=oc.cooling_T,
        pressures=case.pressures if oc.pressures is None else oc.pressures,
    )


def load_cases(cases_yaml: str | Path, base: en.ThermalCase) -> dict[str, OperatingConditions]:
    """Operating conditions of the REDUCTOR cases (data/furnace_cases.yaml).

    The pellet diameter is taken from each case; the mean pressure from its top, bustle and
    cooling pressures with the geometry of ``base``.
    """
    with open(cases_yaml, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)["cases"]
    segs = base.shaft.segments
    if segs is None:
        raise ValueError("base must be a whole-furnace case")
    zone = base.shaft.side_inlets[0].depth
    length = sum(sg.length for sg in segs)
    inert = base.shaft.inert
    out = {}
    for name, c in raw.items():

        def v(key: str, c: Mapping[str, Mapping[str, float]] = c) -> float:
            return float(c[key]["value"])

        x = {"H2": v("x_H2"), "H2O": v("x_H2O"), inert: v("x_inert")}
        out[name] = OperatingConditions(
            bustle_flow=v("bustle_flow"),
            bustle_T=v("bustle_temperature"),
            cooling_flow=v("cooling_flow"),
            cooling_T=v("cooling_temperature"),
            x_gas=x,
            solid_flow=v("solid_flow"),
            solid_T=v("solid_temperature"),
            pellet_diameter=v("pellet_diameter"),
            p=mean_pressure(
                zone, length, v("top_pressure"), v("bustle_pressure"), v("cooling_pressure")
            ),
            pressures=(v("top_pressure"), v("bustle_pressure"), v("cooling_pressure")),
        )
    return out
