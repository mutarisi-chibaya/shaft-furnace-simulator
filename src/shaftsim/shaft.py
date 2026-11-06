"""Isothermal 1D counter-current moving-bed model of the reduction zone of a shaft furnace.

Pellets enter at the top (z = 0) and move down at constant velocity; gas (H2, H2O, inert) enters
at the bottom (z = H) and moves up. Plug flow for both phases, no axial dispersion, uniform
temperature and pressure, constant bed porosity. Each pellet is reduced according to
shaftsim.pellet (grainy pellet model, law of additive reaction times) in the local gas.

Steady-state balances, with z pointing down:
    solid (per pellet):  u_s dX_j/dz = r_j(X, x_H2, x_H2O)                       j = 1, 2, 3
    gas (upward flow):   dF_H2/dz = N_p * 3 n0 * dX/dz,   F_H2O = F_tot - F_inert - F_H2
where N_p is the pellet number flow [1/s], n0 the moles of Fe2O3 per pellet and X the overall
conversion (fraction of removable oxygen removed). One mole of H2 removes one mole of O.

Discretisation (finite volumes, first-order upwind, n cells of height dz; cell 0 at the top):
    (X_i - X_{i-1}) - (dz/u_s) r(X_i, gas_i) = 0          (solid from the cell above)
    (f_{i+1} - f_i) - c (Xbar_i - Xbar_{i-1}) = 0          (gas from the cell below)
with f = F_H2/F_tot, c = 3 n0 N_p / F_tot, X_{-1} = 0 and f_n = inlet. Each cell is evaluated
at its outlet state (upwind for both phases). Because the gas source is written with the same
discrete solid conversion change, the O, H and Fe balances close to round-off for any n.
The nonlinear system is solved by shaftsim.solver.

Whole furnace (Milestone 5): the case may consist of several segments (cylinders or cones, so
the cross-section and the solid velocity vary) and have side gas inlets (the bustle) at
segment boundaries. Everything per cell is collected in a :class:`Grid`: residence time
kappa_i = dz_i (1 - eps) A_i / V_s, total gas flow F_i leaving the cell upward, side
injection S_i entering the cell (with H2 fraction f_s and H2 + H2O fraction g_s). The gas
equation becomes
    F_i f_i = F_{i+1} f_{i+1} + S_i f_s,i - N_O (Xbar_i - Xbar_{i-1})
(N_O the removable oxygen flow), i.e. with the mixed inflow fraction
f_mix,i = (F_{i+1} f_{i+1} + S_i f_s,i) / F_i the same form as above with c_i = N_O / F_i.
Without segments and side inlets the grid is the uniform one above, unchanged.

Units: SI. Mole fractions and conversions are dimensionless.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from functools import cached_property
from pathlib import Path

import numpy as np
import yaml
from numpy.typing import NDArray

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Segment:
    """A vertical section of the furnace: a cylinder or a truncated cone.

    Attributes:
        length: Height of the segment [m].
        r_top: Radius at its top [m].
        r_bottom: Radius at its bottom [m] (equal to r_top for a cylinder).
    """

    length: float
    r_top: float
    r_bottom: float


@dataclass(frozen=True)
class SideInlet:
    """A gas inlet at a segment boundary (e.g. the bustle), counted from the top.

    Attributes:
        depth: Depth of the inlet below the top [m]; must coincide with a segment boundary
            or the bottom. The gas joins the gas rising from below at that face.
        flow: Molar flow [mol/s].
        x: Mole fractions (keys H2, H2O and the inert of the case).
        T: Temperature [K] (used by the energy balance only).
    """

    depth: float
    flow: float
    x: Mapping[str, float]
    T: float


@dataclass(frozen=True)
class ShaftCase:
    """Operating conditions and geometry of an isothermal reduction zone.

    Attributes:
        height: Height of the reduction zone [m].
        area: Cross-section [m2].
        bed_porosity: Void fraction of the bed [-].
        solid_mass_flow: Mass flow of the fresh burden [kg/s].
        pellet: Pellet properties (diameter, density, composition, structure, kinetics).
        gas_flow: Total molar gas flow at the bottom [mol/s].
        x_in: Inlet gas mole fractions (keys H2, H2O and the inert).
        inert: Name of the inert gas species.
        T: Temperature of both phases [K].
        p: Pressure [Pa].
        segments: Optional vertical segments (whole furnace with a cone); if None, one
            cylinder of ``height`` and ``area``. Their lengths must add up to ``height``;
            ``area`` is then the cross-section at the top.
        side_inlets: Optional side gas inlets (bustle). ``gas_flow`` and ``x_in`` are then
            those of the bottom inlet only.
    """

    height: float
    area: float
    bed_porosity: float
    solid_mass_flow: float
    pellet: pl.PelletParameters
    gas_flow: float
    x_in: Mapping[str, float]
    T: float
    p: float
    inert: str = "N2"
    segments: tuple[Segment, ...] | None = None
    side_inlets: tuple[SideInlet, ...] = ()

    @property
    def pellet_mass(self) -> float:
        """Mass of one pellet [kg]."""
        return self.pellet.density * math.pi * self.pellet.diameter**3 / 6.0

    @property
    def pellet_flow(self) -> float:
        """Number of pellets entering per second [1/s]."""
        return self.solid_mass_flow / self.pellet_mass

    @property
    def solid_volume_flow(self) -> float:
        """Volume flow of pellets [m3/s] (pellet volume, not bed volume)."""
        return self.solid_mass_flow / self.pellet.density

    @property
    def total_gas_flow(self) -> float:
        """Gas entering through all inlets [mol/s] (= gas leaving at the top)."""
        return self.gas_flow + sum(s.flow for s in self.side_inlets)

    @property
    def solid_velocity(self) -> float:
        """Downward velocity of the burden in the top cross-section ``area`` [m/s]."""
        solid_volume_flow = self.solid_mass_flow / self.pellet.density
        return solid_volume_flow / ((1.0 - self.bed_porosity) * self.area)

    @property
    def residence_time(self) -> float:
        """Solid residence time in the zone [s]."""
        return self.height / self.solid_velocity

    @property
    def removable_oxygen_flow(self) -> float:
        """Oxygen flow that complete reduction would remove [mol/s]."""
        return pl.O_PER_FE2O3 * pl.initial_moles_fe2o3(self.pellet) * self.pellet_flow


def load_case(
    furnace_yaml: str | Path,
    kinetics_yaml: str | Path,
    T: float | None = None,
) -> ShaftCase:
    """Isothermal reduction-zone case built from case A of Marsigny et al. (2025, Metals 15:862).

    Choices (documented in docs/validation.md):
    - Zone = cylindrical upper part above the bustle (15 - 4.9 = 10.1 m), radius 2.6 m.
    - All gas (bustle + cooling, 2800 mol/s, 93/2/5 % H2/H2O/N2) enters at the bottom of
      the zone; in the real furnace the cooling gas joins from below the bustle.
    - Isothermal at the bustle gas temperature (890 degC) unless T is given; constant pressure
      equal to the mean of the bustle and top pressures.
    - Pellet properties other than diameter and hematite content are those of the Vale pellets
      of data/kinetics.yaml.
    """
    with open(furnace_yaml, encoding="utf-8") as fh:
        f = yaml.safe_load(fh)

    def v(section: str, key: str) -> float:
        return float(f[section][key]["value"])

    pellet, _ = pl.load_kinetics(kinetics_yaml)
    pellet = replace(pellet, diameter=v("solid", "pellet_diameter"), w_Fe2O3=v("solid", "w_Fe2O3"))
    radius = v("geometry", "radius")
    return ShaftCase(
        height=v("geometry", "upper_height") - v("geometry", "bustle_height"),
        area=math.pi * radius**2,
        bed_porosity=v("solid", "bed_porosity"),
        solid_mass_flow=v("solid", "mass_flow"),
        pellet=pellet,
        gas_flow=v("gas", "total_flow"),
        x_in={"H2": v("gas", "x_H2"), "H2O": v("gas", "x_H2O"), "N2": v("gas", "x_N2")},
        T=v("gas", "bustle_temperature") if T is None else T,
        p=0.5 * (v("gas", "bustle_pressure") + v("gas", "top_pressure")),
        inert="N2",
    )


def bed_mass_transfer_coefficient(case: ShaftCase, td: tr.TransportData) -> float:
    """k_g of H2 [m/s] in the bed (Wakao-Kaguei), at the inlet gas composition."""
    if td.wakao_kaguei is None:
        raise ValueError("transport data lack the Wakao-Kaguei constants")
    x = dict(case.x_in)
    u = case.gas_flow * pr.R_GAS * case.T / (case.p * case.area)  # superficial velocity
    m_mix = sum(xi * td.molar_mass[s] for s, xi in x.items())
    rho = case.p * m_mix / (pr.R_GAS * case.T)
    mu = tr.mixture_viscosity(td, x, case.T)
    d_h2 = tr.mixture_diffusivity(td, "H2", x, case.T, case.p, trace_partner="H2O")
    return tr.ranz_marshall_kg(td.wakao_kaguei, case.pellet.diameter, u, rho, mu, d_h2)


def shaft_coefficients(
    case: ShaftCase,
    species: Mapping[str, pr.Species],
    td: tr.TransportData,
    triplets: pl.Triplets | None = None,
) -> tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients]:
    """Pellet characteristic-time coefficients for the (isothermal) shaft conditions.

    The three steps are used at all temperatures with their own equilibrium constants. Below
    the wustite stability limit of the thermodynamic data (~893 K with JANAF Fe0.947O) the
    wustite is metastable: it forms only where the gas is reducing enough for step 2. The
    direct Fe3O4 -> Fe path that real ores follow below ~570 degC is not in the kinetic model
    (Marsigny et al. 2025 state validity above 570 degC); see docs/validation.md.
    """
    k_g = bed_mass_transfer_coefficient(case, td)
    d_mix = tr.mixture_diffusivity(td, "H2", case.x_in, case.T, case.p, trace_partner="H2O")
    return pl.step_coefficients(case.pellet, case.T, case.p, k_g, d_mix, species, td, triplets)


@dataclass(frozen=True)
class Grid:
    """Per-cell geometry and flows of the finite-volume grid (cell 0 at the top).

    Attributes:
        z_faces: Depth of the cell faces below the top, shape (n + 1,) [m].
        z: Cell-centre depths [m].
        dz: Cell heights [m].
        area: Mean cross-section of each cell [m2] (exact mean of pi r^2 over the cell).
        kappa: Solid residence time in each cell, dz (1 - eps) A / V_s [s].
        F: Total gas flow leaving each cell upward [mol/s].
        S: Side gas flow entering each cell [mol/s] (0 for most cells).
        f_side, g_side: H2 and H2 + H2O fractions of that side gas [-].
        g: H2 + H2O fraction of the gas in each cell [-] (fixed by the inlets).
        side_cells: For each side inlet of the case, the index of the cell it enters.
    """

    z_faces: FloatArray
    z: FloatArray
    dz: FloatArray
    area: FloatArray
    kappa: FloatArray
    F: FloatArray
    S: FloatArray
    f_side: FloatArray
    g_side: FloatArray
    g: FloatArray
    side_cells: tuple[int, ...]


#: Limiter regularisation (gradient scale below which the reconstruction fades to first
#: order) for mole/conversion fractions [1/m] and temperatures [K/m]. Numerical choices:
#: without it the limiter has a kink where both gradients vanish (flat zones, e.g. where the
#: gas is buffered at the magnetite/wustite equilibrium), on which Newton's method stalls;
#: where gradients are this small the first-order face error (dz^2 * gradient) is negligible.
LIMITER_EPS_FRACTION = 1e-6
LIMITER_EPS_T = 1e-3


def outflow_faces(
    phi: FloatArray, z: FloatArray, dz: FloatArray, downward: bool, eps: float = 0.0
) -> FloatArray:
    """Value of a transported quantity on the outflow face of each cell (second order).

    MUSCL reconstruction with the van Albada limiter applied to gradients:
    a, b = upstream and downstream gradients (difference / centre spacing, so the scheme
    stays second order where segments of different cell heights meet),
    g = a b (a + b) / (a^2 + b^2 + eps^2), face = phi + g dz / 2 towards the
    outflow. For a linear profile g = a = b (exact); for a monotone profile the face value
    lies between the cell and its downstream neighbour (the limiter is TVD for r > 0), so
    conversions stay in [0, 1]. With eps > 0, g is analytic in a and b (also where both
    vanish) and the scheme fades to first order where |a|, |b| << eps.
    The first and last cells keep first-order faces (face = cell value): each then acts as
    one backward-Euler step (local error O(h^2)), which keeps the global order 2, and the
    outlet values are the cell values exactly as in the first-order scheme.

    Args:
        phi: Cell values, shape (..., n), cell 0 at the top.
        z: Cell-centre depths [m], shape (n,).
        dz: Cell heights [m], shape (n,).
        downward: Flow direction (True: solid, towards larger z; False: gas, upward).
        eps: Limiter regularisation [units of phi per m]; 0 = plain van Albada.

    Returns:
        Outflow face values, same shape as ``phi``.
    """
    phi = np.asarray(phi, dtype=float)
    out = phi.copy()
    if phi.shape[-1] < 3:
        return out
    c = phi[..., 1:-1]
    above = (c - phi[..., :-2]) / (z[1:-1] - z[:-2])  # gradient towards the cell above
    below = (phi[..., 2:] - c) / (z[2:] - z[1:-1])  # gradient towards the cell below
    a, b = (above, below) if downward else (below, above)
    den = a * a + b * b + eps * eps
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(den > 0.0, a * b * (a + b) / den, 0.0)
    half = 0.5 * dz[1:-1]
    out[..., 1:-1] = c + (g * half if downward else -g * half)
    return out


def make_grid(case: ShaftCase, n_cells: int) -> Grid:
    """Build the finite-volume grid of a case.

    Cells are distributed over the segments in proportion to their length (at least one per
    segment) and are uniform within a segment, so segment boundaries and side inlets fall on
    cell faces. A side inlet at depth z_k enters the cell just above the face z_k: it mixes
    with the gas rising from below before that cell.
    """
    segs = case.segments
    if segs is None:
        r = math.sqrt(case.area / math.pi)
        segs = (Segment(case.height, r, r),)
    total = sum(sg.length for sg in segs)
    if not math.isclose(total, case.height, rel_tol=1e-12):
        raise ValueError(f"segment lengths add up to {total} m, not height {case.height} m")
    if n_cells < len(segs):
        raise ValueError("need at least one cell per segment")
    counts = [max(1, round(n_cells * sg.length / total)) for sg in segs]
    counts[-1] = n_cells - sum(counts[:-1])
    if counts[-1] < 1:
        raise ValueError("too few cells for the segments")
    faces = [0.0]
    r_faces = [segs[0].r_top]
    for sg, m in zip(segs, counts, strict=True):
        z0 = faces[-1]
        for k in range(1, m + 1):
            faces.append(z0 + sg.length * k / m)
            r_faces.append(sg.r_top + (sg.r_bottom - sg.r_top) * k / m)
    zf = np.array(faces)
    zf[-1] = case.height
    rf = np.array(r_faces)
    dz = np.diff(zf)
    r_a, r_b = rf[:-1], rf[1:]
    area = math.pi * (r_a**2 + r_a * r_b + r_b**2) / 3.0
    kappa = dz * (1.0 - case.bed_porosity) * area / case.solid_volume_flow

    S = np.zeros(n_cells)
    f_side = np.zeros(n_cells)
    g_side = np.zeros(n_cells)
    side_cells = []
    for si in case.side_inlets:
        k = int(np.argmin(np.abs(zf - si.depth)))
        if k == 0 or not math.isclose(float(zf[k]), si.depth, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"side inlet at {si.depth} m is not on a segment boundary")
        i = k - 1  # the cell just above the face
        if S[i] > 0.0:
            raise ValueError("two side inlets in one cell")
        h2, h2o = si.x["H2"], si.x.get("H2O", 0.0)
        S[i], f_side[i], g_side[i] = si.flow, h2, h2 + h2o
        side_cells.append(i)
    # total flow and H2 + H2O fraction leaving each cell upward, marching up from the bottom
    F = np.empty(n_cells)
    g = np.empty(n_cells)
    F_below = case.gas_flow
    gF_below = case.gas_flow * (case.x_in["H2"] + case.x_in.get("H2O", 0.0))
    for i in range(n_cells - 1, -1, -1):
        F[i] = F_below + S[i]
        gF = gF_below + S[i] * g_side[i]
        g[i] = gF / F[i]
        F_below, gF_below = F[i], gF
    return Grid(
        z_faces=zf,
        z=0.5 * (zf[:-1] + zf[1:]),
        dz=dz,
        area=area,
        kappa=kappa,
        F=F,
        S=S,
        f_side=f_side,
        g_side=g_side,
        g=g,
        side_cells=tuple(side_cells),
    )


@dataclass(frozen=True)
class ShaftProblem:
    """Discrete steady-state problem: everything the residual needs.

    Attributes:
        case: Operating conditions.
        coeffs: Pellet step coefficients (isothermal: the same in every cell).
        n_cells: Number of finite volumes.
        rate_multiplier: Scales all reaction rates (1 = physical; for limiting-case tests).
        cell_coeffs: Optional per-cell coefficients (length n_cells); if given they replace
            ``coeffs`` (non-isothermal shaft). A step that cannot proceed in a cell has an
            infinite a_chem (zero rate).
        availability_width, min_width: Smoothing of the availability constraint (see
            :func:`shaftsim.pellet.step_rates_array`); 0 = the original discontinuous rule.
        scheme: 1 = first-order upwind (donor cell: the outflow of a cell carries the cell
            value); 2 = second-order MUSCL faces (:func:`outflow_faces`). Both are
            conservative. Scheme 2 couples each cell to two neighbours on each side, so it
            is solved by shaftsim.coupled only (the marching species solvers assume 1).
    """

    case: ShaftCase
    coeffs: tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients]
    n_cells: int
    rate_multiplier: float = 1.0
    cell_coeffs: (
        tuple[tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients], ...] | None
    ) = None
    availability_width: float = 0.0
    min_width: float = 0.0
    scheme: int = 1

    def coeffs_of(
        self, i: int
    ) -> tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients]:
        """Step coefficients of cell i (per-cell if given, e.g. non-isothermal; else uniform)."""
        return self.coeffs if self.cell_coeffs is None else self.cell_coeffs[i]

    @cached_property
    def _coeff_arrays(
        self,
    ) -> tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients]:
        """Per-step coefficients with array fields of length n (for vectorised rates)."""
        if self.cell_coeffs is None:
            return self.coeffs
        out = []
        for j in range(3):
            cs = [c[j] for c in self.cell_coeffs]
            out.append(
                pl.StepCoefficients(
                    a_chem=np.array([c.a_chem for c in cs]),  # type: ignore[arg-type]
                    a_diff=np.array([c.a_diff for c in cs]),  # type: ignore[arg-type]
                    a_ext=np.array([c.a_ext for c in cs]),  # type: ignore[arg-type]
                    ratio=np.array([c.ratio for c in cs]),  # type: ignore[arg-type]
                )
            )
        return out[0], out[1], out[2]

    @cached_property
    def grid(self) -> Grid:
        """Per-cell geometry and flows (see :func:`make_grid`)."""
        return make_grid(self.case, self.n_cells)

    @cached_property
    def c_cell(self) -> FloatArray:
        """Removable oxygen flow divided by the gas flow leaving each cell, N_O / F_i [-]."""
        return self.case.removable_oxygen_flow / self.grid.F

    @property
    def dz(self) -> float:
        """Cell height [m] (single-segment grids; ``grid.dz`` in general)."""
        return self.case.height / self.n_cells

    @property
    def c_gas(self) -> float:
        """Removable oxygen flow divided by the bottom gas flow [-] (single-inlet cases)."""
        return self.case.removable_oxygen_flow / self.case.gas_flow

    @property
    def f_in(self) -> float:
        """Inlet H2 flow fraction F_H2/F_tot [-]."""
        return float(self.case.x_in["H2"])

    @property
    def g_in(self) -> float:
        """Bottom-inlet H2 + H2O flow fraction [-] (per cell: ``grid.g``)."""
        return float(self.case.x_in["H2"] + self.case.x_in.get("H2O", 0.0))

    @property
    def f_max(self) -> float:
        """Upper bound of the H2 fraction anywhere: the largest inlet H2 fraction [-]."""
        return float(max([self.f_in] + [s.x["H2"] for s in self.case.side_inlets]))

    def gas_faces(self, f: FloatArray) -> FloatArray:
        """H2 fraction on the upper (outflow) face of each cell [-]: the cell value for
        scheme 1, :func:`outflow_faces` for scheme 2."""
        if self.scheme == 1:
            return f
        gr = self.grid
        return outflow_faces(f, gr.z, gr.dz, downward=False, eps=LIMITER_EPS_FRACTION)

    def solid_faces(self, X: FloatArray) -> FloatArray:
        """Step conversions on the lower (outflow) face of each cell, shape (3, n)."""
        if self.scheme == 1:
            return X
        gr = self.grid
        return outflow_faces(X, gr.z, gr.dz, downward=True, eps=LIMITER_EPS_FRACTION)

    def mixed_inflow(self, f: FloatArray) -> FloatArray:
        """H2 fraction of the gas entering each cell from below, mixed with any side gas [-].

        f_mix,i = (F_{i+1} f_{i+1} + S_i f_s,i) / F_i, with the bottom inlet below the last
        cell; f_{i+1} is the value on the upper face of cell i+1 (its cell value in scheme 1).
        """
        gr = self.grid
        F_below = np.concatenate([gr.F[1:], [self.case.gas_flow]])
        f_below = np.concatenate([self.gas_faces(f)[1:], [self.f_in]])
        return np.asarray((F_below * f_below + gr.S * gr.f_side) / gr.F)

    def oxygen_removed(self, f: FloatArray) -> FloatArray:
        """Oxygen removed from the solid in each cell [mol/s] (= H2 consumed), from the gas."""
        return np.asarray(self.grid.F * (self.mixed_inflow(f) - self.gas_faces(f)))

    def unpack(self, y: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Unknown vector -> (X shape (3, n), f shape (n,)). Ordering per cell: X1, X2, X3, f."""
        Y = y.reshape(self.n_cells, 4)
        return Y[:, :3].T, Y[:, 3]

    def pack(self, X: FloatArray, f: FloatArray) -> FloatArray:
        """Inverse of :meth:`unpack`."""
        return np.column_stack([X.T, f]).ravel()

    def residual(self, y: FloatArray, one_minus: FloatArray | None = None) -> FloatArray:
        """Scaled steady-state residual (dimensionless), same ordering as the unknowns.

        ``one_minus`` optionally gives 1 - X, shape (3, n), at full precision (see
        :func:`shaftsim.pellet.step_rates_array`).
        """
        X, f = self.unpack(y)
        gr = self.grid
        x_h2 = f
        x_h2o = gr.g - f
        r = pl.step_rates_array(
            self._coeff_arrays,
            X,
            x_h2,
            x_h2o,
            self.rate_multiplier,
            one_minus=one_minus,
            availability_width=self.availability_width,
            min_width=self.min_width,
        )
        Xo = self.solid_faces(X)  # outflow (lower) faces; = X in scheme 1
        X_up = np.concatenate([np.zeros((3, 1)), Xo[:, :-1]], axis=1)
        res_s = (Xo - X_up) - gr.kappa * r
        xbar = np.asarray(pl.overall_conversion(Xo))
        xbar_up = np.concatenate([[0.0], xbar[:-1]])
        res_g = (self.mixed_inflow(f) - self.gas_faces(f)) - self.c_cell * (xbar - xbar_up)
        return np.column_stack([res_s.T, res_g]).ravel()


@dataclass(frozen=True)
class ShaftSolution:
    """Steady-state profiles and outlet indicators.

    Attributes:
        z: Cell-centre depth below the top [m].
        X_steps: Step conversions, shape (3, n) [-].
        X: Overall conversion (fraction of removable O removed) [-].
        x_H2, x_H2O: Gas mole fractions in each cell [-].
        metallization: Metallic Fe / total Fe at the solid outlet [-].
        gud: Gas utilisation degree at the top, x_H2O / (x_H2O + x_H2) [-].
        h2_consumption: H2 consumed per tonne of product [kg/t].
        product_flow: Mass flow of reduced burden leaving the zone [kg/s].
        residual_norm: Max-norm of the scaled residual at convergence [-].
        iterations: Nonlinear iterations used.
    """

    z: FloatArray
    X_steps: FloatArray
    X: FloatArray
    x_H2: FloatArray
    x_H2O: FloatArray
    metallization: float
    gud: float
    h2_consumption: float
    product_flow: float
    residual_norm: float
    iterations: int


def make_solution(
    problem: ShaftProblem, y: FloatArray, td: tr.TransportData, residual_norm: float, iters: int
) -> ShaftSolution:
    """Assemble profiles and indicators from a converged unknown vector."""
    case = problem.case
    X, f = problem.unpack(y)
    xbar = np.asarray(pl.overall_conversion(X))
    gr = problem.grid
    x_h2 = f.copy()
    x_h2o = gr.g - f
    m_o = td.molar_mass["H2O"] - td.molar_mass["H2"]
    o_removed = case.removable_oxygen_flow * float(xbar[-1])
    product = case.solid_mass_flow - o_removed * m_o
    h2_in = case.gas_flow * problem.f_in + float(np.sum(gr.S * gr.f_side))
    h2_used = h2_in - float(gr.F[0] * f[0])
    return ShaftSolution(
        z=gr.z.copy(),
        X_steps=X,
        X=xbar,
        x_H2=x_h2,
        x_H2O=x_h2o,
        metallization=float(X[2, -1]),
        gud=float(x_h2o[0] / (x_h2o[0] + x_h2[0])),
        h2_consumption=h2_used * td.molar_mass["H2"] / (product / 1000.0),
        product_flow=product,
        residual_norm=residual_norm,
        iterations=iters,
    )


def element_balances(problem: ShaftProblem, sol: ShaftSolution) -> dict[str, float]:
    """Relative closure (in - out)/in of the Fe, O and H balances over the zone [-].

    Solids: Fe and O from the pellet composition at the inlet (all hematite) and the outlet.
    Gas: H2 and H2O flows at the bottom and side inlets and at the top (outlet).
    """
    case = problem.case
    n0 = pl.initial_moles_fe2o3(case.pellet)
    np_ = case.pellet_flow
    out = pl.solid_moles(case.pellet, sol.X_steps[:, -1])
    fe_in = 2 * n0 * np_
    fe_out = np_ * (
        2 * out["Fe2O3"] + 3 * out["Fe3O4"] + pl.X_FE_WUSTITE_KINETIC * out["Fe0.95O"] + out["Fe"]
    )
    o_solid_in = 3 * n0 * np_
    o_solid_out = np_ * (3 * out["Fe2O3"] + 4 * out["Fe3O4"] + out["Fe0.95O"])
    F = case.gas_flow
    h2o_in, h2_in = F * case.x_in.get("H2O", 0.0), F * case.x_in["H2"]
    for si in case.side_inlets:
        h2o_in += si.flow * si.x.get("H2O", 0.0)
        h2_in += si.flow * si.x["H2"]
    F_top = float(problem.grid.F[0])
    h2o_top, h2_top = F_top * sol.x_H2O[0], F_top * sol.x_H2[0]
    o_in, o_out = o_solid_in + h2o_in, o_solid_out + h2o_top
    h_in, h_out = 2 * (h2_in + h2o_in), 2 * (h2_top + h2o_top)
    return {
        "Fe": (fe_in - fe_out) / fe_in,
        "O": (o_in - o_out) / o_in,
        "H": (h_in - h_out) / h_in,
    }
