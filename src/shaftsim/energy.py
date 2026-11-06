"""Non-isothermal 1D shaft: solid and gas energy balances coupled to the species model.

Energy balances (steady, z downward, per cell i of height dz; upwind as in shaftsim.shaft)
-------------------------------------------------------------------------------------
Enthalpies are absolute (formation-referenced, shaftsim.properties), so reaction heats appear
automatically. Let Hs_i be the enthalpy flow of the burden leaving cell i (all solid species
plus gangue at T_s,i), Hg_i that of the gas leaving cell i upward (at T_g,i), n_i the oxygen
removed in cell i [mol/s] (= H2 consumed = H2O produced), and UA_i = h a A dz:

    solid:  Hs_i - Hs_{i-1} = UA_i (T_g,i - T_s,i) + n_i [H_H2(T_g,i) - H_H2O(T_s,i)]
    gas:    Hg_i - Hg_{i+1} = -UA_i (T_g,i - T_s,i) - n_i [H_H2(T_g,i) - H_H2O(T_s,i)]

The H2 consumed leaves the gas at the gas temperature and the H2O formed leaves the pellet at
the solid temperature, so the heat of reaction is released/absorbed in the solid, as in
REDUCTOR (Marsigny et al. 2025b, eq. 6). Summing over cells, the total enthalpy flow in equals
the total out exactly (discrete conservation).

Heat transfer: Wakao-Kaguei, Nu = 2 + 1.1 Pr^(1/3) Re^0.6 (data/transport.yaml), with the gas
conductivity of shaftsim.transport (modified Eucken, Wilke mixing); specific area
a = 6 (1 - eps_bed) / d_p. Adiabatic wall, no radiation, no axial conduction or dispersion,
constant pressure.

Solution: outer fixed-point iteration between
  (1) the species problem at frozen temperatures (shaftsim.solver; per-cell kinetic
      coefficients evaluated at T_s, external transfer at T_g), and
  (2) the energy problem at frozen species profiles: 2n continuous, piecewise-smooth equations
      in (T_s, T_g) (enthalpies bridged across NIST range boundaries, see enthalpy_ext),
      solved by Newton's method with the analytic block-tridiagonal Jacobian (dH/dT = Cp),
with under-relaxation of the temperatures. Converged when the temperature change is below
1e-6 K and both residuals are below tolerance.

Whole furnace (Milestone 5): with segments and side inlets (shaftsim.shaft.Grid) the gas
balance of a cell also receives the enthalpy of the side gas, Hside_i = S_i h(x_s, T_s):
    gas:    Hg_i - Hg_{i+1} - Hside_i = -UA_i (T_g,i - T_s,i) - n_i [...]
with Hg_i = F_i h(f_i, g_i, T_g,i) for the per-cell flow F_i. The transfer coefficients scale
with the local gas mass flux: Wakao-Kaguei is a + b Re^m X^n, so k_g and h are tabulated as
the two T-dependent parts (a-part, b-part at a reference flux) and the b-part is multiplied by
(G_i / G_ref)^m, exact for the inlet composition. UA_i uses the cell's own area A_i dz_i.

The three reduction steps are used at all temperatures with their own equilibrium constants
(below the wustite stability limit of the data, 893 K, wustite is metastable). An earlier
version switched steps 2 and 3 off below that limit: the switch is discontinuous in T and,
because the wustite reduction is endothermic, cells oscillated across it and the coupled
iteration did not converge (docs/validation.md).

Units: SI (K, Pa, mol/s, W, J/mol).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from functools import cached_property
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import yaml
from numpy.typing import NDArray

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import solver as sv
from shaftsim import transport as tr

FloatArray = NDArray[np.float64]
GANGUE = ("CaO", "SiO2", "MgO", "Al2O3")
GAS = ("H2", "H2O")
#: Fe in the kinetic wustite (Fe0.95O) minus Fe in the thermodynamic wustite (Fe0.947O) per
#: mol O: H(Fe0.95O) := H(Fe0.947O) + 0.003 H(Fe), element-consistent.
_DX_FE = pl.X_FE_WUSTITE_KINETIC - pr.X_WUSTITE


# --------------------------------------------------------------------------- properties


#: Half-width [K] of the window over which H is interpolated linearly across a boundary between
#: two NIST Shomate ranges (see :func:`enthalpy_ext`).
TRANSITION_HALF_WIDTH = 10.0


def _bridged(sp_: pr.Species, T: FloatArray) -> tuple[FloatArray, FloatArray]:
    """NIST H [J/mol] and Cp [J/(mol K)] with every range boundary bridged; T inside the data."""
    h = np.asarray(pr.enthalpy(sp_, T), dtype=float)
    cp = np.asarray(pr.heat_capacity(sp_, T), dtype=float)
    w = TRANSITION_HALF_WIDTH
    for r in sp_.ranges[1:]:
        a, b = r.T_min - w, r.T_min + w
        inside = (T > a) & (T < b)
        if not np.any(inside):
            continue
        ha = float(pr.enthalpy(sp_, a))
        slope = (float(pr.enthalpy(sp_, b)) - ha) / (b - a)
        h = np.where(inside, ha + slope * (T - a), h)
        cp = np.where(inside, slope, cp)
    return h, cp


def enthalpy_ext(sp_: pr.Species, T: FloatArray | float) -> FloatArray:
    """Molar enthalpy [J/mol] for the energy balances: NIST H made continuous, extended linearly.

    Range boundaries: NIST's piecewise Shomate H jumps where two ranges meet (Fe2O3 +0.66 kJ/mol
    at 950 K, SiO2 +0.73 at 847 K (alpha-beta quartz), Fe +0.14/-0.21/-0.36 at 700/1042/1100 K).
    For Fe2O3 the jump carries the enthalpy of the lambda-type magnetic transition that the Cp
    polynomials omit: NIST H agrees with JANAF within 0.04 kJ/mol on both sides, whereas
    integrating Cp alone is off by 0.64 kJ/mol above 950 K. For Fe the jumps are a defect of
    the fit (JANAF has no latent heat there; docs/validation.md, Milestone 1), and integrating
    Cp would not be better (up to 0.47 kJ/mol off). A jump, and a negative one makes H
    non-monotone, can leave a cell energy balance without a root, so Newton's method cycles
    across it. Within +/-TRANSITION_HALF_WIDTH (10 K) of each boundary H is therefore
    interpolated linearly between the NIST values at the window edges: H is continuous and
    strictly increasing (smallest effective Cp 28 J/(mol K), Fe at 1100 K), unchanged outside
    the windows, and a transition enthalpy is spread over 20 K, as real lambda transitions are
    spread.

    Outside [T_min, T_max] (H2O gas below 500 K, solids below 298 K; the burden enters at
    293 K): H(T_lim) + Cp(T_lim) (T - T_lim). The error is small over the few kelvin (solids) to
    ~200 K (H2O) involved; documented in validation.md.
    """
    T = np.asarray(T, dtype=float)
    Tc = np.clip(T, sp_.T_min, sp_.T_max)
    h, cp = _bridged(sp_, Tc)
    return h + cp * (T - Tc)


def cp_ext(sp_: pr.Species, T: FloatArray | float) -> FloatArray:
    """Molar heat capacity [J/(mol K)] consistent with :func:`enthalpy_ext` (its derivative)."""
    T = np.asarray(T, dtype=float)
    return _bridged(sp_, np.clip(T, sp_.T_min, sp_.T_max))[1]


# --------------------------------------------------------------------------- case


@dataclass(frozen=True)
class ThermalCase:
    """Non-isothermal reduction zone.

    Attributes:
        shaft: Geometry, flows, pellet and gas composition (its ``T`` is not used).
        T_solid_in: Burden inlet temperature at the top [K].
        T_gas_in: Gas inlet temperature at the bottom [K] (side inlets carry their own T).
        gangue: Mass fractions of the gangue oxides in the fresh burden [-].
    """

    shaft: sh.ShaftCase
    T_solid_in: float
    T_gas_in: float
    gangue: Mapping[str, float]


def load_thermal_case(furnace_yaml: str | Path, kinetics_yaml: str | Path) -> ThermalCase:
    """Case A of Marsigny et al. (2025b) with the energy balance (see shaft.load_case)."""
    with open(furnace_yaml, encoding="utf-8") as fh:
        f = yaml.safe_load(fh)
    base = sh.load_case(furnace_yaml, kinetics_yaml)
    return ThermalCase(
        shaft=base,
        T_solid_in=float(f["solid_inlet"]["temperature"]["value"]),
        T_gas_in=float(f["gas"]["bustle_temperature"]["value"]),
        gangue={k: float(f["gangue"][k]["value"]) for k in GANGUE},
    )


def load_furnace_case(furnace_yaml: str | Path, kinetics_yaml: str | Path) -> ThermalCase:
    """Whole furnace of case A (Marsigny et al. 2025b): upper cylinder and cooling cone.

    Geometry, top to bottom: reduction zone above the bustle (15 - 4.9 = 10.1 m, radius 2.6 m),
    cylinder below the bustle (4.9 m), cone (10 m, radius 2.6 -> 1.31 m, linear). Gas: bustle
    gas (2300 mol/s, 890 degC) enters at the bustle face as a side inlet, cooling gas
    (500 mol/s, 70 degC) at the bottom; both 93/2/5 % H2/H2O/N2. Burden at 20 degC.

    Pressure: one value for the whole furnace, the length-weighted mean of the piecewise-linear
    profile through the top (150 kPa), bustle (180.4 kPa) and cooling (181.6 kPa) pressures of
    Table 1 (the last two are REDUCTOR outputs). Pressure enters only the transport properties
    (equilibria are equimolar). The mean is 174.6 kPa; the local pressure differs from it by
    at most 14 % (at the top, 150 kPa).
    """
    with open(furnace_yaml, encoding="utf-8") as fh:
        f = yaml.safe_load(fh)

    def v(section: str, key: str) -> float:
        return float(f[section][key]["value"])

    base = sh.load_case(furnace_yaml, kinetics_yaml)
    upper, bustle_h = v("geometry", "upper_height"), v("geometry", "bustle_height")
    cone, r, r_out = (
        v("geometry", "cone_height"),
        v("geometry", "radius"),
        v("geometry", "outlet_radius"),
    )
    zone = upper - bustle_h
    p_top, p_b, p_c = (
        v("gas", "top_pressure"),
        v("gas", "bustle_pressure"),
        v("gas", "cooling_pressure"),
    )
    length = upper + cone
    p_mean = (zone * 0.5 * (p_top + p_b) + (length - zone) * 0.5 * (p_b + p_c)) / length
    x = dict(base.x_in)
    shaft = replace(
        base,
        height=length,
        gas_flow=v("gas", "cooling_flow"),
        p=p_mean,
        segments=(
            sh.Segment(zone, r, r),
            sh.Segment(bustle_h, r, r),
            sh.Segment(cone, r, r_out),
        ),
        side_inlets=(
            sh.SideInlet(zone, v("gas", "bustle_flow"), x, v("gas", "bustle_temperature")),
        ),
    )
    return ThermalCase(
        shaft=shaft,
        T_solid_in=v("solid_inlet", "temperature"),
        T_gas_in=v("gas", "cooling_temperature"),
        gangue={k: v("gangue", k) for k in GANGUE},
    )


# --------------------------------------------------------------------------- problem


@dataclass(frozen=True)
class ThermalProblem:
    """Everything needed to evaluate the discrete energy and species equations.

    ``pressure_cells`` (optional, one value per cell [Pa]) gives a local pressure for the
    pellet kinetics (Milestone 8, shaftsim.pressure); None = the uniform ``case.shaft.p``.
    Pressure enters only through the total gas concentration p/RT and the molecular part of
    the pore diffusivity (both exact rescalings of the tables, see :meth:`cell_coefficients`).
    """

    case: ThermalCase
    species: Mapping[str, pr.Species]
    td: tr.TransportData
    triplets: pl.Triplets | None
    n_cells: int
    rate_multiplier: float = 1.0
    h_multiplier: float = 1.0
    kinetics: pl.Triplets | None = None
    tortuosity_factor: float = 1.0
    transfer_mix: float = 0.0
    availability_width: float = 0.0
    min_width: float = 0.0
    scheme: int = 1
    pressure_cells: tuple[float, ...] | None = None

    #: Fields that change the solution but not the property tables (see with_parameters).
    TABLE_FREE_FIELDS = frozenset(
        {
            "n_cells",
            "rate_multiplier",
            "h_multiplier",
            "kinetics",
            "tortuosity_factor",
            "transfer_mix",
            "availability_width",
            "min_width",
            "scheme",
            "pressure_cells",
        }
    )

    def with_parameters(self, **changes: object) -> ThermalProblem:
        """A copy with other uncertain parameters that reuses this problem's property tables.

        Only the fields in TABLE_FREE_FIELDS may change; they enter as exact factors on the
        tabulated values (Milestone 6): ``kinetics`` replaces the rate constants (a_chem is
        proportional to 1/k_j(T)), ``tortuosity_factor`` multiplies a_diff (proportional to
        the tortuosity, D_eff = eps/tau ...), ``transfer_mix`` blends the gas-solid transfer
        coefficients h and k_g linearly between Wakao-Kaguei (0) and Ranz-Marshall (1).
        """
        bad = set(changes) - self.TABLE_FREE_FIELDS
        if bad:
            raise ValueError(f"these fields change the property tables: {sorted(bad)}")
        new = replace(self, **changes)  # type: ignore[arg-type]
        new.__dict__["_tables"] = self._tables
        return new

    @property
    def dz(self) -> float:
        """Cell height [m] (single-segment grids; ``grid.dz`` in general)."""
        return self.case.shaft.height / self.n_cells

    @cached_property
    def grid(self) -> sh.Grid:
        """Per-cell geometry and flows (shaftsim.shaft.make_grid)."""
        return sh.make_grid(self.case.shaft, self.n_cells)

    @cached_property
    def _flux_factor(self) -> FloatArray:
        """(G_i / G_ref)^m per cell: gas mass flux relative to the table reference (bottom
        inlet flow over the top cross-section), to the Wakao-Kaguei Re exponent [-]."""
        case = self.case.shaft
        wk = self.td.wakao_kaguei
        if wk is None:
            raise ValueError("transport data lack the Wakao-Kaguei constants")
        g_ref = case.gas_flow / case.area
        return np.asarray(((self.grid.F / self.grid.area) / g_ref) ** wk.m)

    @cached_property
    def _flux_factor_rm(self) -> FloatArray:
        """As _flux_factor, with the Ranz-Marshall Re exponent."""
        case = self.case.shaft
        g_ref = case.gas_flow / case.area
        return np.asarray(((self.grid.F / self.grid.area) / g_ref) ** self.td.ranz_marshall.m)

    @cached_property
    def side_enthalpy(self) -> FloatArray:
        """Enthalpy flow of the side gas entering each cell [W] (0 where there is none)."""
        out = np.zeros(self.n_cells)
        for si, i in zip(self.case.shaft.side_inlets, self.grid.side_cells, strict=True):
            g = si.x["H2"] + si.x.get("H2O", 0.0)
            out[i] = self.gas_enthalpy(
                np.array([si.x["H2"]]), np.array([si.T]), si.flow, np.array([g])
            )[0][0]
        return out

    def oxygen_removed(self, f: FloatArray) -> FloatArray:
        """Oxygen removed in each cell [mol/s], from the gas balance (mixed inflow - outflow),
        with the face values of the scheme (:meth:`faces`)."""
        gr = self.grid
        case = self.case.shaft
        fo = self.faces(f, downward=False)
        F_below = np.concatenate([gr.F[1:], [case.gas_flow]])
        f_below = np.concatenate([fo[1:], [case.x_in["H2"]]])
        return np.asarray(F_below * f_below + gr.S * gr.f_side - gr.F * fo)

    def faces(self, phi: FloatArray, downward: bool, temperature: bool = False) -> FloatArray:
        """Outflow-face values of a cell quantity: the cell values for scheme 1, the MUSCL
        reconstruction of :func:`shaftsim.shaft.outflow_faces` for scheme 2 (limiter
        regularisation for fractions or, with ``temperature``, for temperatures)."""
        if self.scheme == 1:
            return phi
        gr = self.grid
        eps = sh.LIMITER_EPS_T if temperature else sh.LIMITER_EPS_FRACTION
        return sh.outflow_faces(phi, gr.z, gr.dz, downward, eps)

    def gangue_flows(self) -> dict[str, float]:
        """Molar flows of the gangue oxides [mol/s]."""
        m = self.case.shaft.solid_mass_flow
        return {k: m * w / self.species[k].molar_mass for k, w in self.case.gangue.items()}  # type: ignore[operator]

    def solid_enthalpy(self, X: FloatArray, Ts: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Enthalpy flow [W] and heat-capacity flow [W/K] of the burden leaving each cell.

        Args:
            X: Step conversions, shape (3, m).
            Ts: Solid temperatures, shape (m,).
        """
        case = self.case.shaft
        n0 = pl.initial_moles_fe2o3(case.pellet) * case.pellet_flow
        x1, x2, x3 = X
        flows = {
            "Fe2O3": n0 * (1 - x1),
            "Fe3O4": 2.0 / 3.0 * n0 * (x1 - x2),
            "W": 2.0 / pl.X_FE_WUSTITE_KINETIC * n0 * (x2 - x3),
            "Fe": 2.0 * n0 * x3,
        }
        H = np.zeros_like(Ts)
        C = np.zeros_like(Ts)
        spc = self.species
        for k in ("Fe2O3", "Fe3O4", "Fe"):
            H += flows[k] * enthalpy_ext(spc[k], Ts)
            C += flows[k] * cp_ext(spc[k], Ts)
        hw = enthalpy_ext(spc[pr.WUSTITE], Ts) + _DX_FE * enthalpy_ext(spc["Fe"], Ts)
        cw = cp_ext(spc[pr.WUSTITE], Ts) + _DX_FE * cp_ext(spc["Fe"], Ts)
        H += flows["W"] * hw
        C += flows["W"] * cw
        for k, n in self.gangue_flows().items():
            H += n * enthalpy_ext(spc[k], Ts)
            C += n * cp_ext(spc[k], Ts)
        return H, C

    def gas_enthalpy(
        self,
        f: FloatArray,
        Tg: FloatArray,
        F: FloatArray | float | None = None,
        g: FloatArray | float | None = None,
    ) -> tuple[FloatArray, FloatArray]:
        """Enthalpy flow [W] and heat-capacity flow [W/K] of a gas stream.

        Args:
            f: H2 fraction [-].
            Tg: Temperature [K].
            F: Total molar flow [mol/s]; default the bottom inlet flow.
            g: H2 + H2O fraction [-] (the rest is the inert); default that of the bottom inlet.
        """
        case = self.case.shaft
        if F is None:
            F = case.gas_flow
        if g is None:
            g = case.x_in["H2"] + case.x_in["H2O"]
        x_inert = 1.0 - g
        spc = self.species
        H = F * (
            f * enthalpy_ext(spc["H2"], Tg)
            + (g - f) * enthalpy_ext(spc["H2O"], Tg)
            + x_inert * enthalpy_ext(spc[case.inert], Tg)
        )
        C = F * (
            f * cp_ext(spc["H2"], Tg)
            + (g - f) * cp_ext(spc["H2O"], Tg)
            + x_inert * cp_ext(spc[case.inert], Tg)
        )
        return H, C

    # ------------------------------------------------------------------ tabulated properties
    #: Temperature grid of the property tables [K] (1 K spacing). Interpolation is linear in
    #: log(value); the error is below 1e-5 relative (tests/test_energy.py).
    T_TABLE: tuple[float, float] = (250.0, 1500.0)

    @cached_property
    def _tables(self) -> dict[str, FloatArray]:
        """Kinetic coefficients vs T_s and transfer coefficients vs T_g, computed once."""
        case = self.case.shaft
        td = self.td
        T = np.arange(self.T_TABLE[0], self.T_TABLE[1] + 0.5, 1.0)
        n = T.size
        out: dict[str, FloatArray] = {
            k: np.empty((3, n)) for k in ("a_chem", "a_diff", "a_diff_kn", "a_ext_kg", "ratio")
        }
        out["T"] = T
        for k in ("kg0", "kg1", "ha0", "ha1", "kgR0", "kgR1", "haR0", "haR1"):
            out[k] = np.empty(n)
        wk = td.wakao_kaguei
        if wk is None:
            raise ValueError("transport data lack the Wakao-Kaguei constants")
        td0 = replace(td, wakao_kaguei=replace(wk, b=0.0))  # the 'a' part alone
        # Ranz-Marshall (single sphere) in place of Wakao-Kaguei, for transfer_mix
        tdR = replace(td, wakao_kaguei=td.ranz_marshall)
        tdR0 = replace(td, wakao_kaguei=replace(td.ranz_marshall, b=0.0))
        # T below the H2O data (500 K): equilibrium constants at 500 K (rates there are
        # negligible anyway, the rate constants still follow Arrhenius at the true T)
        t_eq_min = self.species["H2O"].T_min
        for i, Ti in enumerate(T):
            Ti = float(Ti)
            d_mix = tr.mixture_diffusivity(td, "H2", case.x_in, Ti, case.p, "H2O")
            c = pl.step_coefficients(
                case.pellet,
                Ti,
                case.p,
                1.0,
                d_mix,
                self.species,
                td,
                self.triplets,
                T_eq=max(Ti, t_eq_min),
            )
            # Knudsen part of the pore resistance alone (molecular diffusivity -> infinity):
            # a_diff = molecular part (independent of p) + Knudsen part (~ 1/p), Bosanquet
            c_kn = pl.step_coefficients(
                case.pellet,
                Ti,
                case.p,
                1.0,
                math.inf,
                self.species,
                td,
                self.triplets,
                T_eq=max(Ti, t_eq_min),
            )
            for j in range(3):
                out["a_chem"][j, i] = c[j].a_chem
                out["a_diff"][j, i] = c[j].a_diff
                out["a_diff_kn"][j, i] = c_kn[j].a_diff
                out["a_ext_kg"][j, i] = c[j].a_ext  # a_ext for k_g = 1 m/s
                out["ratio"][j, i] = c[j].ratio
            kg = sh.bed_mass_transfer_coefficient(replace(case, T=Ti), td)
            kg0 = sh.bed_mass_transfer_coefficient(replace(case, T=Ti), td0)
            ha = self._ua_per_dz_direct(Ti) / case.area
            ha0 = self._ua_per_dz_direct(Ti, td0) / case.area
            out["kg0"][i], out["kg1"][i] = kg0, kg - kg0
            out["ha0"][i], out["ha1"][i] = ha0, ha - ha0
            kgR = sh.bed_mass_transfer_coefficient(replace(case, T=Ti), tdR)
            kgR0 = sh.bed_mass_transfer_coefficient(replace(case, T=Ti), tdR0)
            haR = self._ua_per_dz_direct(Ti, tdR) / case.area
            haR0 = self._ua_per_dz_direct(Ti, tdR0) / case.area
            out["kgR0"][i], out["kgR1"][i] = kgR0, kgR - kgR0
            out["haR0"][i], out["haR1"][i] = haR0, haR - haR0
        return out

    def _interp(self, key: str, T: FloatArray, j: int | None = None) -> FloatArray:
        tab = self._tables
        y = tab[key] if j is None else tab[key][j]
        Tc = np.clip(T, tab["T"][0], tab["T"][-1])
        return np.exp(np.interp(Tc, tab["T"], np.log(y)))

    def _ua_per_dz_direct(self, T: float, td: tr.TransportData | None = None) -> float:
        """h a A [W/(K m)] at gas temperature T (Wakao-Kaguei), bottom-inlet composition, flow
        and top cross-section (the table reference)."""
        case = self.case.shaft
        td = self.td if td is None else td
        if td.wakao_kaguei is None:
            raise ValueError("transport data lack the Wakao-Kaguei constants")
        wk = td.wakao_kaguei
        x = dict(case.x_in)
        m_mix = sum(xi * td.molar_mass[s] for s, xi in x.items())
        a = 6.0 * (1.0 - case.bed_porosity) / case.pellet.diameter
        cp = {s: float(cp_ext(self.species[s], T)) for s in x}
        mu = tr.mixture_viscosity(td, x, T)
        lam = tr.mixture_conductivity(td, x, T, cp)
        cp_mass = sum(x[s] * cp[s] for s in x) / m_mix
        u = case.gas_flow * pr.R_GAS * T / (case.p * case.area)
        rho = case.p * m_mix / (pr.R_GAS * T)
        re = rho * u * case.pellet.diameter / mu
        pr_ = cp_mass * mu / lam
        nu = wk.a + wk.b * re**wk.m * pr_**wk.n
        return float(nu * lam / case.pellet.diameter * a * case.area)

    def ua_per_length_ref(self, Tg: FloatArray) -> FloatArray:
        """h a A [W/(K m)] at the table reference flux and area (single-inlet cylinder)."""
        return np.asarray(
            (self._interp("ha0", Tg) + self._interp("ha1", Tg)) * self.case.shaft.area
        )

    def mass_transfer_kg(self, Tg: FloatArray) -> FloatArray:
        """k_g of H2 in each cell [m/s] at the local gas flux (Wakao-Kaguei, blended with
        Ranz-Marshall by ``transfer_mix``)."""
        kg = self._interp("kg0", Tg) + self._interp("kg1", Tg) * self._flux_factor
        lam = self.transfer_mix
        if lam != 0.0:
            kg_r = self._interp("kgR0", Tg) + self._interp("kgR1", Tg) * self._flux_factor_rm
            kg = (1.0 - lam) * kg + lam * kg_r
        return np.asarray(kg)

    def heat_transfer_UA(self, Tg: FloatArray) -> FloatArray:
        """UA of each cell [W/K]: h a A_i dz_i (Wakao-Kaguei at the local gas flux, inlet
        gas composition)."""
        gr = self.grid
        ha = self._interp("ha0", Tg) + self._interp("ha1", Tg) * self._flux_factor
        lam = self.transfer_mix
        if lam != 0.0:
            ha_r = self._interp("haR0", Tg) + self._interp("haR1", Tg) * self._flux_factor_rm
            ha = (1.0 - lam) * ha + lam * ha_r
        return np.asarray(ha * gr.area * gr.dz * self.h_multiplier)

    def cell_coefficients(
        self, Ts: FloatArray, Tg: FloatArray
    ) -> tuple[tuple[pl.StepCoefficients, pl.StepCoefficients, pl.StepCoefficients], ...]:
        """Pellet step coefficients per cell: kinetics and pores at T_s, transfer at T_g."""
        kg = self.mass_transfer_kg(Tg)
        # local pressure: a_chem ~ 1/c_t ~ 1/p; a_diff = molecular part (p-free) + Knudsen part
        # (~ 1/p); a_ext ~ 1/(k_g c_t) is independent of p (k_g ~ 1/p); tables at case p
        scale = (
            None
            if self.pressure_cells is None
            else self.case.shaft.p / np.asarray(self.pressure_cells, dtype=float)
        )
        steps = []
        for j in range(3):
            a_chem = self._interp("a_chem", Ts, j)
            a_diff = self._interp("a_diff", Ts, j)
            if scale is not None:
                a_chem = a_chem * scale
                a_diff = a_diff + self._interp("a_diff_kn", Ts, j) * (scale - 1.0)
            if self.kinetics is not None:
                # a_chem ~ 1/k_j(T): exact rescaling from the table triplets to the new ones
                t_tab = (self.triplets or self.case.shaft.pellet.triplets)[j]
                t_new = self.kinetics[j]
                a_chem = (
                    a_chem * (t_tab.k0 / t_new.k0) * np.exp((t_new.Ea - t_tab.Ea) / (pr.R_GAS * Ts))
                )
            steps.append(
                (
                    a_chem,
                    a_diff * self.tortuosity_factor,
                    self._interp("a_ext_kg", Ts, j) / kg,
                    self._interp("ratio", Ts, j),
                )
            )
        out = []
        for i in range(Ts.size):
            out.append(
                tuple(
                    pl.StepCoefficients(
                        float(steps[j][0][i]),
                        float(steps[j][1][i]),
                        float(steps[j][2][i]),
                        float(steps[j][3][i]),
                    )
                    for j in range(3)
                )
            )
        return tuple(out)  # type: ignore[arg-type]

    def species_problem(self, Ts: FloatArray, Tg: FloatArray) -> sh.ShaftProblem:
        """Species problem at frozen temperature profiles."""
        cc = self.cell_coefficients(Ts, Tg)
        return sh.ShaftProblem(
            self.case.shaft,
            cc[0],
            self.n_cells,
            self.rate_multiplier,
            cc,
            availability_width=self.availability_width,
            min_width=self.min_width,
            scheme=self.scheme,
        )

    # ------------------------------------------------------------------ energy equations

    def energy_residual(
        self, X: FloatArray, f: FloatArray, Ts: FloatArray, Tg: FloatArray, UA: FloatArray
    ) -> FloatArray:
        """Energy residuals [W], ordering per cell (solid, gas)."""
        case = self.case.shaft
        gr = self.grid
        # enthalpy flows on the outflow faces (cell states in scheme 1); the exchange terms
        # use the cell states and cancel between the phases, so the balance stays exact
        Hs, _ = self.solid_enthalpy(self.faces(X, True), self.faces(Ts, True, True))
        Hs_in, _ = self.solid_enthalpy(np.zeros((3, 1)), np.array([self.case.T_solid_in]))
        Hg, _ = self.gas_enthalpy(self.faces(f, False), self.faces(Tg, False, True), gr.F, gr.g)
        f_in = case.x_in["H2"]
        Hg_in, _ = self.gas_enthalpy(np.array([f_in]), np.array([self.case.T_gas_in]))
        Hs_up = np.concatenate([Hs_in, Hs[:-1]])
        Hg_down = np.concatenate([Hg[1:], Hg_in])
        n_o = self.oxygen_removed(f)  # O removed in each cell [mol/s]
        q = UA * (Tg - Ts)
        xfer = n_o * (enthalpy_ext(self.species["H2"], Tg) - enthalpy_ext(self.species["H2O"], Ts))
        r_s = Hs - Hs_up - q - xfer
        r_g = Hg - Hg_down - self.side_enthalpy + q + xfer
        return np.column_stack([r_s, r_g]).ravel()

    def solve_energy(
        self,
        X: FloatArray,
        f: FloatArray,
        Ts: FloatArray,
        Tg: FloatArray,
        tol: float = 1e-6,
        max_iter: int = 50,
    ) -> tuple[FloatArray, FloatArray, float]:
        """Newton's method for (T_s, T_g) at frozen species; returns (Ts, Tg, max |dT| last).

        Scheme 1 only (the analytic Jacobian assumes donor-cell fluxes).
        """
        if self.scheme != 1:
            raise ValueError("solve_energy handles scheme 1 only; use shaftsim.coupled")
        UA = self.heat_transfer_UA(Tg)
        n_o = self.oxygen_removed(f)
        gr = self.grid
        spc = self.species
        dT_max = np.inf
        for _ in range(max_iter):
            r = self.energy_residual(X, f, Ts, Tg, UA)
            _, Cs = self.solid_enthalpy(X, Ts)
            _, Cg = self.gas_enthalpy(f, Tg, gr.F, gr.g)
            cp_h2 = cp_ext(spc["H2"], Tg)
            cp_h2o = cp_ext(spc["H2O"], Ts)
            J = energy_jacobian(Cs, Cg, UA, n_o, cp_h2, cp_h2o)
            dT = spla.spsolve(J, -r)
            Ts = Ts + dT[0::2]
            Tg = Tg + dT[1::2]
            dT_max = float(np.max(np.abs(dT)))
            if dT_max < tol:
                break
        return Ts, Tg, dT_max


def energy_jacobian(
    Cs: FloatArray,
    Cg: FloatArray,
    UA: FloatArray,
    n_o: FloatArray,
    cp_h2: FloatArray,
    cp_h2o: FloatArray,
) -> sp.csc_matrix:
    """Analytic Jacobian of :meth:`ThermalProblem.energy_residual` w.r.t. (T_s, T_g) per cell.

    Block-tridiagonal: the solid row of cell i depends on T_s,i, T_g,i and T_s,i-1 (upwind
    from above); the gas row on T_g,i, T_s,i and T_g,i+1 (upwind from below). UA is held fixed
    (lagged in the Newton iteration).
    """
    n = Cs.size
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n):
        s_, g_ = 2 * i, 2 * i + 1
        rows += [s_, s_]
        cols += [s_, g_]
        vals += [Cs[i] + UA[i] + n_o[i] * cp_h2o[i], -UA[i] - n_o[i] * cp_h2[i]]
        if i > 0:
            rows.append(s_)
            cols.append(2 * (i - 1))
            vals.append(-Cs[i - 1])
        rows += [g_, g_]
        cols += [g_, s_]
        vals += [Cg[i] + UA[i] + n_o[i] * cp_h2[i], -UA[i] - n_o[i] * cp_h2o[i]]
        if i < n - 1:
            rows.append(g_)
            cols.append(2 * (i + 1) + 1)
            vals.append(-Cg[i + 1])
    return sp.csc_matrix((vals, (rows, cols)), shape=(2 * n, 2 * n))


# --------------------------------------------------------------------------- solution


@dataclass(frozen=True)
class ThermalSolution:
    """Converged non-isothermal profiles and indicators.

    Attributes:
        species: Species solution (shaftsim.shaft.ShaftSolution) at the final temperatures.
        Ts, Tg: Solid and gas temperatures per cell [K].
        T_top_gas: Gas temperature leaving the top [K].
        T_solid_out: Burden temperature leaving the bottom [K].
        energy_closure: (enthalpy in - out) / enthalpy in, absolute enthalpies [-].
        energy_closure_sensible: (in - out) / (sensible enthalpy of the inlet gas) [-].
        outer_iterations: Fixed-point iterations.
        converged: Whether all convergence criteria were met.
    """

    species: sh.ShaftSolution
    Ts: FloatArray
    Tg: FloatArray
    T_top_gas: float
    T_solid_out: float
    energy_closure: float
    energy_closure_sensible: float
    outer_iterations: int
    converged: bool


def energy_closure(
    prob: ThermalProblem, X: FloatArray, f: FloatArray, Ts: FloatArray, Tg: FloatArray
) -> tuple[float, float]:
    """Global enthalpy balance of the zone: relative to total inlet and to sensible inlet.

    In: burden at the top, gas at the bottom and side inlets. Out: burden at the bottom, gas
    at the top. The sensible reference is the inlet gas enthalpy relative to 298.15 K.
    """
    case = prob.case.shaft
    gr = prob.grid
    f_in = case.x_in["H2"]
    Hs_in = prob.solid_enthalpy(np.zeros((3, 1)), np.array([prob.case.T_solid_in]))[0][0]
    Hg_in = prob.gas_enthalpy(np.array([f_in]), np.array([prob.case.T_gas_in]))[0][0]
    Hg_in += float(np.sum(prob.side_enthalpy))
    Hs_out = prob.solid_enthalpy(X[:, -1:], Ts[-1:])[0][0]
    Hg_out = prob.gas_enthalpy(f[:1], Tg[:1], gr.F[:1], gr.g[:1])[0][0]
    h_in, h_out = Hs_in + Hg_in, Hs_out + Hg_out
    Hg_298 = prob.gas_enthalpy(np.array([f_in]), np.array([298.15]))[0][0]
    for si in case.side_inlets:
        g = si.x["H2"] + si.x.get("H2O", 0.0)
        Hg_298 += prob.gas_enthalpy(
            np.array([si.x["H2"]]), np.array([298.15]), si.flow, np.array([g])
        )[0][0]
    return (h_in - h_out) / abs(h_in), (h_in - h_out) / abs(Hg_in - Hg_298)


def _interp_profile(values: FloatArray, z_old: FloatArray, z_new: FloatArray) -> FloatArray:
    """Linear interpolation of a cell-centred profile onto other cell centres [same units]."""
    return np.asarray(np.interp(z_new, z_old, values))


def solve_thermal(
    prob: ThermalProblem,
    tol_T: float = 1e-6,
    max_outer: int = 400,
    relax: float = 0.5,
    T_start: tuple[FloatArray, FloatArray] | None = None,
    coarse_cells: int = 25,
    anderson_m: int = 5,
    species_tol: float = 1e-9,
) -> ThermalSolution:
    """Solve the coupled species and energy problem (outer fixed point, see module docstring).

    Numerical choices: the outer map G(T) (species at frozen T, then energy at frozen species)
    is iterated with Anderson acceleration (memory ``anderson_m``, damping ``relax``; history
    reset if the change grows tenfold). Plain damped iteration converges but slowly: the cold
    front where the endothermic wustite reduction slows down moves about one cell per few
    iterations. Grids finer than ``coarse_cells`` therefore start from the interpolated
    solution on half the cells (grid continuation, recursively).

    Args:
        prob: Thermal problem.
        tol_T: Convergence threshold on the temperature change between outer iterations [K].
        max_outer: Maximum outer iterations.
        relax: Damping factor of the (Anderson-accelerated) temperature update.
        T_start: Optional starting (T_s, T_g) profiles.
        coarse_cells: Grid size below which no coarser start is computed.
        anderson_m: Anderson memory (0 = damped fixed-point iteration).
        species_tol: Residual tolerance of the species solver (conversion units). 1e-9 by
            default; Monte Carlo runs use 1e-8, far below any physical effect (the gas balance
            is closed exactly afterwards in any case).

    Returns:
        The solution; ``converged`` must be checked by the caller.
    """
    n = prob.n_cells
    case = prob.case
    X = np.zeros((3, n))
    f = np.full(n, case.shaft.x_in["H2"])
    if T_start is not None:
        Ts, Tg = T_start[0].copy(), T_start[1].copy()
    elif n > coarse_cells:
        coarse_prob = replace(prob, n_cells=n // 2)
        coarse_prob.__dict__["_tables"] = prob._tables  # tables do not depend on the grid
        coarse = solve_thermal(
            coarse_prob, tol_T, max_outer, relax, None, coarse_cells, anderson_m, species_tol
        )
        zc, zf = coarse_prob.grid.z, prob.grid.z
        Ts, Tg = _interp_profile(coarse.Ts, zc, zf), _interp_profile(coarse.Tg, zc, zf)
    else:
        # start: pure counter-current heat exchange (no reaction)
        Ts = np.full(n, case.T_gas_in)
        Tg = np.full(n, case.T_gas_in)
        Ts, Tg, _ = prob.solve_energy(X, f, Ts, Tg)
    y = None

    def picard_map(T: FloatArray) -> tuple[FloatArray, bool]:
        """G(T): species at frozen T, then energy at frozen species. Returns (G(T), ok)."""
        nonlocal y
        sprob = prob.species_problem(T[:n], T[n:])
        y, st = sv.solve(sprob, y0=None, tol=species_tol)
        if not st.converged:
            y, st = sv.solve(sprob, y0=y, tol=species_tol)
        X_, f_ = sprob.unpack(y)
        Ts_new, Tg_new, _ = prob.solve_energy(X_, f_, T[:n].copy(), T[n:].copy())
        return np.concatenate([Ts_new, Tg_new]), st.converged

    x = np.concatenate([Ts, Tg])
    xs: list[FloatArray] = []
    gs: list[FloatArray] = []
    best = np.inf
    converged = False
    outer = 0
    for _k in range(max_outer):
        outer += 1
        gx, ok = picard_map(x)
        g = gx - x
        change = float(np.max(np.abs(g)))
        if change < tol_T and ok:
            x = gx
            converged = True
            break
        if change > 10.0 * best:  # Anderson went astray: restart its history
            xs, gs = [], []
        best = min(best, change)
        xs.append(x.copy())
        gs.append(g.copy())
        if len(xs) > anderson_m + 1:
            xs.pop(0)
            gs.pop(0)
        if len(xs) == 1 or anderson_m == 0:
            x = x + relax * g
        else:
            dG = np.column_stack([gs[i + 1] - gs[i] for i in range(len(gs) - 1)])
            dX = np.column_stack([xs[i + 1] - xs[i] for i in range(len(xs) - 1)])
            gamma, *_ = np.linalg.lstsq(dG, g, rcond=None)
            x = x + relax * g - (dX + relax * dG) @ gamma
    Ts, Tg = x[:n], x[n:]
    # final consistent state: species at the final temperatures, then energy for them
    sprob = prob.species_problem(Ts, Tg)
    y, st = sv.solve(sprob, y0=y, tol=species_tol)
    X, f = sprob.unpack(y)
    Ts, Tg, dT = prob.solve_energy(X, f, Ts, Tg)
    converged = converged and st.converged and dT < tol_T
    sol = sh.make_solution(sprob, y, prob.td, st.residuals[-1], st.iterations)
    closure, closure_s = energy_closure(prob, X, f, Ts, Tg)
    return ThermalSolution(
        species=sol,
        Ts=Ts,
        Tg=Tg,
        T_top_gas=float(Tg[0]),
        T_solid_out=float(Ts[-1]),
        energy_closure=closure,
        energy_closure_sensible=closure_s,
        outer_iterations=outer,
        converged=converged,
    )
