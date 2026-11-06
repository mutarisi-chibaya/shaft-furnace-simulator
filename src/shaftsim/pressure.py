"""Gas pressure field of the shaft from the Ergun equation (Milestone 8).

Ergun (1952), as given by Wagner (2008, Sec. III.5.4; constants in data/transport.yaml):

    -dP/dh = [K_v (1 - eps)^2 / (eps^3 d^2) mu + K_i (1 - eps) / (eps^3 d) rho u] u,

with u the superficial gas velocity, eps the bed porosity and d the pellet diameter. With z the
depth (gas flowing upward, towards smaller z) the pressure rises downward. In a cell with
molar gas flow F, area A, temperature T and mean molar mass M, the ideal gas gives
u = F R T / (p A) and rho u = F M / A, so

    dp/dz = C / p,   C = R T F / A * [K_v (1 - eps)^2 / (eps^3 d^2) mu
                                       + K_i (1 - eps) / (eps^3 d) F M / A],

and p^2 is integrated exactly cell by cell (C constant in a cell): p_(k+1)^2 = p_k^2 + 2 C dz,
from the top face (top pressure, a boundary value) downward. Cell pressures are taken at the
cell centre. The bustle and cooling-gas pressures (REDUCTOR outputs in Table 1 of Marsigny et
al. 2025) are the values at the bustle face and at the bottom face.

Coupling: the pressure depends on the gas temperature and composition, the solution on the
pressure only through the pellet kinetics (shaftsim.energy.ThermalProblem.pressure_cells).
:func:`solve_with_pressure` iterates the two (fixed point; the kinetics depend weakly on p).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from . import energy as en
from . import properties as pr
from . import transport as tr

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class PressureProfile:
    """Pressure along the shaft.

    Attributes:
        z_faces: Depth of the cell faces below the top [m].
        p_faces: Pressure at the faces [Pa] (p_faces[0] = top pressure).
        p_cells: Pressure at the cell centres [Pa].
        bustle: Pressure at the bustle face [Pa] (NaN without side inlet).
        cooling: Pressure at the bottom face, where the cooling gas enters [Pa].
    """

    z_faces: FloatArray
    p_faces: FloatArray
    p_cells: FloatArray
    bustle: float
    cooling: float


def ergun_coefficient(prob: en.ThermalProblem, f: FloatArray, Tg: FloatArray) -> FloatArray:
    """C per cell [Pa^2/m] such that dp/dz = C / p (see the module docstring).

    Args:
        prob: Thermal problem (geometry, flows, pellet diameter, porosity, inert gas).
        f: H2 mole fraction per cell [-].
        Tg: Gas temperature per cell [K].
    """
    er = prob.td.ergun
    if er is None:
        raise ValueError("transport data lack the Ergun constants (data/transport.yaml)")
    case = prob.case.shaft
    gr = prob.grid
    eps = case.bed_porosity
    d = case.pellet.diameter
    mm = prob.td.molar_mass
    inert = case.inert
    kv = er.viscous * (1.0 - eps) ** 2 / (eps**3 * d**2)
    ki = er.inertial * (1.0 - eps) / (eps**3 * d)
    out = np.empty(gr.F.size)
    for i in range(gr.F.size):
        x = {"H2": float(f[i]), "H2O": float(gr.g[i] - f[i]), inert: float(1.0 - gr.g[i])}
        x = {k: max(v, 0.0) for k, v in x.items()}
        mu = tr.mixture_viscosity(prob.td, x, float(Tg[i]))
        m_mix = sum(v * mm[k] for k, v in x.items())
        flux = gr.F[i] / gr.area[i]  # molar flux [mol/(m2 s)]
        out[i] = pr.R_GAS * Tg[i] * flux * (kv * mu + ki * flux * m_mix)
    return out


def pressure_profile(
    prob: en.ThermalProblem, f: FloatArray, Tg: FloatArray, p_top: float
) -> PressureProfile:
    """Pressure field for given gas profiles (exact integration of p^2 per cell).

    Args:
        prob: Thermal problem.
        f: H2 mole fraction per cell [-].
        Tg: Gas temperature per cell [K].
        p_top: Pressure at the top face (gas outlet) [Pa].
    """
    gr = prob.grid
    C = ergun_coefficient(prob, f, Tg)
    p2 = p_top**2 + np.concatenate([[0.0], np.cumsum(2.0 * C * gr.dz)])
    p_faces = np.sqrt(p2)
    p_cells = np.sqrt(p2[:-1] + C * gr.dz)
    if gr.side_cells:
        bustle = float(p_faces[gr.side_cells[0] + 1])
    else:
        bustle = float("nan")
    return PressureProfile(gr.z_faces.copy(), p_faces, p_cells, bustle, float(p_faces[-1]))


def solve_with_pressure(
    prob: en.ThermalProblem,
    p_top: float,
    solve: Callable[[en.ThermalProblem, en.ThermalSolution | None], en.ThermalSolution],
    tol: float = 1.0,
    max_iter: int = 10,
) -> tuple[en.ThermalSolution, PressureProfile, int]:
    """Solution with the Ergun pressure field in the kinetics (fixed point on the pressure).

    Start: the uniform pressure of ``prob``. Each iteration solves with the current cell
    pressures (``solve(problem, previous_solution)``, e.g. a warm start), then recomputes the
    pressure field from the gas profiles. Stops when the cell pressures change by less than
    ``tol``.

    Args:
        prob: Thermal problem (its own ``pressure_cells`` is ignored).
        p_top: Top pressure [Pa] (boundary value).
        solve: Solver callback.
        tol: Convergence threshold on the cell pressures [Pa].
        max_iter: Maximum iterations.

    Returns:
        (solution, pressure profile consistent with it, iterations). The solution's own
        ``converged`` flag must be checked, and the iteration count compared with max_iter.
    """
    p = prob.with_parameters(pressure_cells=None)
    sol = solve(p, None)
    prof = pressure_profile(p, sol.species.x_H2, sol.Tg, p_top)
    for it in range(1, max_iter + 1):
        p = prob.with_parameters(pressure_cells=tuple(float(v) for v in prof.p_cells))
        sol = solve(p, sol)
        new = pressure_profile(p, sol.species.x_H2, sol.Tg, p_top)
        change = float(np.max(np.abs(new.p_cells - prof.p_cells)))
        prof = new
        if change < tol:
            return sol, prof, it
    return sol, prof, max_iter
