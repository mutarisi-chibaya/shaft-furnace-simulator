"""One furnace scenario from operating conditions to KPIs and profiles (used by the app).

``load_base`` reads every data file once; ``run_scenario`` applies operating conditions to the
case A furnace, solves it (Milestone 4 solver with the Milestone 6 sampling tolerances and
fallback) and collects outlet KPIs, profiles, energy demand and balance closures.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from . import coupled as cp
from . import emissions as em
from . import energy as en
from . import operating as op
from . import pellet as pl
from . import properties as pr
from . import shaft as sh
from . import transport as tr
from . import uq

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class BaseModel:
    """Everything loaded from data/ that a scenario needs.

    Attributes:
        case: Case A whole furnace.
        species: Thermodynamic data.
        td: Transport data.
        triplets: Refitted kinetic constants.
        emissions: Energy and CO2 parameters.
        reference_cases: Operating conditions of REDUCTOR cases A-E.
    """

    case: en.ThermalCase
    species: dict[str, pr.Species]
    td: tr.TransportData
    triplets: pl.Triplets
    emissions: em.EmissionsData
    reference_cases: dict[str, op.OperatingConditions]


def load_base(data_dir: str | Path) -> BaseModel:
    """Read the data files of the case A furnace (a few hundred ms)."""
    d = Path(data_dir)
    species = pr.load_species(d / "species.yaml")
    td = tr.load_transport(
        d / "transport.yaml",
        d / "collision_integrals.csv",
        {k: v.molar_mass for k, v in species.items() if v.molar_mass is not None},
    )
    case = en.load_furnace_case(d / "furnace.yaml", d / "kinetics.yaml")
    return BaseModel(
        case=case,
        species=species,
        td=td,
        triplets=pl.load_refit_triplets(d / "kinetics_refit.yaml"),
        emissions=em.load_emissions(d / "emissions.yaml", d / "grid_intensity.yaml"),
        reference_cases=op.load_cases(d / "furnace_cases.yaml", case),
    )


@dataclass(frozen=True)
class ScenarioResult:
    """Solution of one scenario (SI units).

    Attributes:
        height: Cell-centre height above the bottom of the upper cylinder [m] (bustle at
            4.9 m, cone below 0).
        Ts, Tg: Solid and gas temperature [K].
        x_H2, x_H2O: Gas mole fractions [-].
        X_steps: Conversions of the three reduction steps, shape (3, n) [-].
        metallization: Metallic Fe / total Fe at the outlet [-].
        h2_consumption: H2 consumed per tonne of DRI [kg/t].
        top_gas_T, dri_T: Top gas and DRI temperatures [K].
        gud: Gas utilisation degree [-].
        dri_flow: DRI mass flow [kg/s].
        energy: Electricity demand per kg of DRI.
        energy_closure: Relative energy balance closure [-].
        element_closure: Largest relative Fe/O/H balance closure [-].
        converged: Whether the solver converged.
        seconds: Wall time [s].
    """

    height: FloatArray
    Ts: FloatArray
    Tg: FloatArray
    x_H2: FloatArray
    x_H2O: FloatArray
    X_steps: FloatArray
    metallization: float
    h2_consumption: float
    top_gas_T: float
    dri_T: float
    gud: float
    dri_flow: float
    energy: em.EnergyDemand
    energy_closure: float
    element_closure: float
    converged: bool
    seconds: float


def run_scenario(base: BaseModel, oc: op.OperatingConditions, n_cells: int = 50) -> ScenarioResult:
    """Solve the furnace at the given operating conditions.

    Uses the sampling tolerances of Milestone 6 (uq.SOLVER) and its fallback
    (uq.TOL_T_RETRY). If the alternating scheme still does not converge (it cycles where the
    cone's cooling front is balanced, Milestone 7), the coupled solver finishes from its end
    state (:func:`shaftsim.coupled.finish_coupled`, smoothing effect below 1e-6). On 50 cells
    a run takes about 10-15 s (property tables included); the result is not grid-converged
    (docs/validation.md, Milestones 5-6).
    """
    t0 = time.perf_counter()
    case = op.apply(base.case, oc)
    prob = en.ThermalProblem(case, base.species, base.td, base.triplets, n_cells)
    sol = en.solve_thermal(prob, **uq.SOLVER)  # type: ignore[arg-type]
    if not sol.converged:
        sol = en.solve_thermal(
            prob,
            tol_T=uq.TOL_T_RETRY,
            species_tol=float(uq.SOLVER["species_tol"]),
            T_start=(sol.Ts, sol.Tg),
        )
    if not sol.converged:  # Milestone 7 band: finish with the coupled solver (Milestone 8)
        sol, _ = cp.finish_coupled(prob, sol)
    s = sol.species
    sprob = prob.species_problem(sol.Ts, sol.Tg)
    closure = sh.element_balances(sprob, s)
    segs = case.shaft.segments or ()
    upper = sum(sg.length for sg in segs[:2]) if len(segs) == 3 else case.shaft.height
    return ScenarioResult(
        height=upper - prob.grid.z,
        Ts=sol.Ts.copy(),
        Tg=sol.Tg.copy(),
        x_H2=s.x_H2.copy(),
        x_H2O=s.x_H2O.copy(),
        X_steps=s.X_steps.copy(),
        metallization=s.metallization,
        h2_consumption=s.h2_consumption,
        top_gas_T=sol.T_top_gas,
        dri_T=sol.T_solid_out,
        gud=s.gud,
        dri_flow=s.product_flow,
        energy=em.furnace_energy_demand(case, sol, base.species, base.emissions),
        energy_closure=float(sol.energy_closure),
        element_closure=float(max(abs(v) for v in closure.values())),
        converged=bool(sol.converged),
        seconds=time.perf_counter() - t0,
    )
