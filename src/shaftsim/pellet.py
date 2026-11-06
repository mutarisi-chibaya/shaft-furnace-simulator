"""Single-pellet reduction kinetics: grainy pellet model with the law of additive reaction times.

Model (Marsigny, Letz, Mirgaux & Patisson 2025, Metals 15:1034, Appendix A, Table A1; based on
Ranzani da Costa et al. 2013 and Wagner 2008). A pellet of diameter d_p is an agglomerate of
dense hematite grains (diameter d_grain). Hematite -> magnetite -> wustite proceed on the
grains; wustite grains break up into dense crystallites (d_cryst) that are reduced to iron.
For each step j (1: Fe2O3->Fe3O4, 2: Fe3O4->Fe0.95O, 3: Fe0.95O->Fe) the characteristic times
of external mass transfer, pore diffusion and chemical reaction are added (law of additive
reaction times, Sohn 1978), giving the instantaneous rate

    dX_j/dt = 1 / (tau_chem f'(X_j) + tau_diff g'(X_j) + tau_ext h'(X_j)),
    f = 1 - (1-X)^(1/3),  g = 1 - 3(1-X)^(2/3) + 2(1-X),  h = X.

X_j is the conversion of oxide j relative to the amount that would be present if all of the
previous oxide had been formed. As in Wagner (2008, Sec. II.5.2), an oxide can only be reduced
when it is present: X_3 <= X_2 <= X_1 (when the constraint is active the step proceeds at the
rate at which its reactant is formed). The overall conversion (fraction of removable oxygen
removed) is X = X_1/9 + 32 X_2/171 + 40 X_3/57.

Characteristic times use b = moles of solid per mole of H2 (3, 19/16, 1); as printed in
Table A1 the chemical time of step 2 uses the factor 2 (2 b = 2.375 would be consistent with
the other entries); the published form is kept because the kinetic constants were fitted with it.
Molar densities: pellet-apparent for external transfer and inter-grain diffusion, true (grain
or crystallite) for the chemical and intra-grain terms, as in Ranzani da Costa et al. (2013),
Table 1.

Departures from the published model, all documented in docs/validation.md:
- Solid-state diffusion and iron sintering are omitted (their constants are not published).
- The 1-D gas discretisation of the TGA furnace is replaced by a single gas state around the
  pellet ("mean" of inlet and outlet composition by default, or static gas).
- Equilibrium constants come from shaftsim.properties (JANAF Fe0.947O), while the
  stoichiometric factors keep the paper's Fe0.95O.
- Pore-scale diffusivities use the bulk gas composition.

Units: SI throughout (s, m, mol, K, Pa). Conversions are dimensionless.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from shaftsim import properties as pr
from shaftsim import transport as tr

FloatArray = NDArray[np.float64]

#: Wustite composition used in the kinetic stoichiometry (paper: Fe0.95O).
X_FE_WUSTITE_KINETIC = 0.95
#: Moles of solid reactant per mole of H2 for steps 1-3 [mol/mol].
B_SOLID_PER_H2 = (3.0, 19.0 / 16.0, 1.0)
#: Fraction of removable oxygen removed by each step (Wagner 2008): 1/9, 32/171, 40/57.
OXYGEN_WEIGHTS = (1.0 / 9.0, 32.0 / 171.0, 40.0 / 57.0)
#: Moles of removable oxygen per mole of Fe2O3 [mol/mol].
O_PER_FE2O3 = 3.0
REACTIONS = (pr.HEMATITE_TO_MAGNETITE, pr.MAGNETITE_TO_WUSTITE, pr.WUSTITE_TO_IRON)


@dataclass(frozen=True)
class KineticTriplet:
    """Arrhenius constants k = k0 exp(-Ea/(R T)); k0 [m/s], Ea [J/mol]."""

    k0: float
    Ea: float

    def k(self, T: float) -> float:
        """Rate constant [m/s] at T [K]."""
        return self.k0 * math.exp(-self.Ea / (pr.R_GAS * T))


@dataclass(frozen=True)
class PelletParameters:
    """Pellet properties and kinetic constants (see data/kinetics.yaml for sources).

    Attributes:
        diameter: Pellet diameter [m].
        density: Initial apparent mass density [kg/m3].
        w_Fe2O3: Hematite mass fraction [-].
        porosity: Initial porosity [-].
        tortuosity: Tortuosity factor [-].
        specific_area: Initial specific surface area [m2/kg].
        d_pores_intergrain: Inter-grain pore diameter [m].
        d_crystallite: Wustite crystallite diameter [m].
        d_pores_crystallite: Pore diameter between crystallites [m].
        porosity_wustite_grain: Porosity inside wustite grains [-].
        true_density: Crystal densities [kg/m3] of Fe2O3, Fe3O4, Fe0.95O, Fe.
        molar_mass: Molar masses [kg/mol] of the same.
        triplets: Kinetic constants of steps 1-3.
    """

    diameter: float
    density: float
    w_Fe2O3: float
    porosity: float
    tortuosity: float
    specific_area: float
    d_pores_intergrain: float
    d_crystallite: float
    d_pores_crystallite: float
    porosity_wustite_grain: float
    true_density: Mapping[str, float]
    molar_mass: Mapping[str, float]
    triplets: tuple[KineticTriplet, KineticTriplet, KineticTriplet]


@dataclass(frozen=True)
class TGAConditions:
    """Thermobalance operating conditions.

    Attributes:
        gas_flow_normal: Gas flow rate at normal conditions [m3(STP)/s].
        x_H2: Inlet H2 mole fraction [-]; the rest is inert unless x_H2O > 0.
        x_H2O: Inlet H2O mole fraction [-].
        inert: Name of the inert gas species.
        pressure: Pressure [Pa].
        furnace_diameter: Inner diameter of the furnace tube [m].
    """

    gas_flow_normal: float
    x_H2: float
    pressure: float
    furnace_diameter: float
    x_H2O: float = 0.0
    inert: str = "He"


def load_kinetics(path: str | Path) -> tuple[PelletParameters, TGAConditions]:
    """Load pellet parameters and the TGA base case from data/kinetics.yaml."""
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)

    def v(section: str, key: str) -> float:
        return float(raw[section][key]["value"])

    kt = raw["kinetic_triplets"]
    triplets = tuple(
        KineticTriplet(float(kt[f"k0_{j}"]["value"]), float(kt[f"Ea_{j}"]["value"]))
        for j in (1, 2, 3)
    )
    pp = PelletParameters(
        diameter=v("pellet", "diameter"),
        density=v("pellet", "density"),
        w_Fe2O3=v("pellet", "w_Fe2O3"),
        porosity=v("pellet", "porosity"),
        tortuosity=v("pellet", "tortuosity"),
        specific_area=v("pellet", "specific_area"),
        d_pores_intergrain=v("structure", "d_pores_intergrain"),
        d_crystallite=v("structure", "d_crystallite"),
        d_pores_crystallite=v("structure", "d_pores_crystallite"),
        porosity_wustite_grain=v("structure", "porosity_wustite_grain"),
        true_density={k: float(e["value"]) for k, e in raw["true_density"].items()},
        molar_mass={k: float(e["value"]) for k, e in raw["molar_mass"].items()},
        triplets=(triplets[0], triplets[1], triplets[2]),
    )
    tga = TGAConditions(
        gas_flow_normal=v("tga_base_case", "gas_flow_normal"),
        x_H2=v("tga_base_case", "x_H2"),
        pressure=v("tga_base_case", "pressure"),
        furnace_diameter=v("tga_base_case", "furnace_diameter"),
    )
    return pp, tga


# --------------------------------------------------------------------------- structure


def initial_moles_fe2o3(pp: PelletParameters) -> float:
    """Moles of Fe2O3 in one pellet [mol]."""
    volume = math.pi * pp.diameter**3 / 6.0
    return pp.density * volume * pp.w_Fe2O3 / pp.molar_mass["Fe2O3"]


def apparent_molar_densities(pp: PelletParameters) -> tuple[float, float, float]:
    """Pellet-apparent molar densities of the reactant of each step [mol/m3 of pellet].

    Fe3O4 and Fe0.95O are counted as if all of the previous oxide had been converted:
    2/3 and 2/0.95 mol per mol of Fe2O3.
    """
    c1 = pp.density * pp.w_Fe2O3 / pp.molar_mass["Fe2O3"]
    return c1, 2.0 / 3.0 * c1, 2.0 / X_FE_WUSTITE_KINETIC * c1


def true_molar_densities(pp: PelletParameters) -> tuple[float, float, float]:
    """Molar densities inside dense grains/crystallites of Fe2O3, Fe3O4, Fe0.95O [mol/m3]."""
    return tuple(  # type: ignore[return-value]
        pp.true_density[s] / pp.molar_mass[s] for s in ("Fe2O3", "Fe3O4", "Fe0.95O")
    )


def grain_diameter(pp: PelletParameters) -> float:
    """Hematite grain diameter from the specific area, dense spheres: d = 6/(rho a) [m]."""
    return 6.0 / (pp.true_density["Fe2O3"] * pp.specific_area)


def layer_porosities(pp: PelletParameters) -> tuple[float, float, float]:
    """Porosity of the magnetite, wustite and iron layers [-].

    Constant pellet volume, single product phase (Marsigny 2025 thesis, eq. 51):
    eps_p = 1 - (1 - eps_r) V_p/V_r, with V the solid volume per mol Fe. Gangue is ignored.
    """
    v = {
        "Fe2O3": pp.molar_mass["Fe2O3"] / pp.true_density["Fe2O3"] / 2.0,
        "Fe3O4": pp.molar_mass["Fe3O4"] / pp.true_density["Fe3O4"] / 3.0,
        "Fe0.95O": pp.molar_mass["Fe0.95O"] / pp.true_density["Fe0.95O"] / X_FE_WUSTITE_KINETIC,
        "Fe": pp.molar_mass["Fe"] / pp.true_density["Fe"],
    }
    eps_m = 1 - (1 - pp.porosity) * v["Fe3O4"] / v["Fe2O3"]
    eps_w = 1 - (1 - eps_m) * v["Fe0.95O"] / v["Fe3O4"]
    eps_fe = 1 - (1 - eps_w) * v["Fe"] / v["Fe0.95O"]
    return eps_m, eps_w, eps_fe


# --------------------------------------------------------------------------- gas


@dataclass(frozen=True)
class GasState:
    """Bulk gas around the pellet: T [K], p [Pa], mole fractions (H2, H2O, inert...)."""

    T: float
    p: float
    x: Mapping[str, float]


@dataclass(frozen=True)
class Environment:
    """Gas state plus the external mass-transfer coefficient k_g [m/s] for H2."""

    gas: GasState
    k_g: float


@dataclass(frozen=True)
class StepTimes:
    """Characteristic times of one reduction step [s]; inf if the step cannot proceed."""

    chem: float
    diff: float
    ext: float


@dataclass(frozen=True)
class StepCoefficients:
    """Temperature-dependent coefficients of the characteristic times of one step.

    tau_chem = a_chem / (x_H2 - x_H2O/ratio), tau_diff = a_diff / (x_H2 - x_eq),
    tau_ext = a_ext / (x_H2 - x_eq), with x_eq = (x_H2 + x_H2O)/(1 + ratio).

    Attributes:
        a_chem, a_diff, a_ext: Coefficients [s].
        ratio: Equilibrium (x_H2O/x_H2) of the step [-].
    """

    a_chem: float
    a_diff: float
    a_ext: float
    ratio: float

    def times(self, x_h2: float, x_h2o: float) -> StepTimes:
        """Characteristic times [s] for the given bulk mole fractions."""
        x_eq = (x_h2 + x_h2o) / (1.0 + self.ratio)
        d_eq = x_h2 - x_eq
        d_chem = x_h2 - x_h2o / self.ratio
        if d_eq <= 0.0 or d_chem <= 0.0:
            return StepTimes(math.inf, math.inf, math.inf)
        return StepTimes(self.a_chem / d_chem, self.a_diff / d_eq, self.a_ext / d_eq)


Triplets = tuple[KineticTriplet, KineticTriplet, KineticTriplet]


def step_coefficients(
    pp: PelletParameters,
    T: float,
    p: float,
    k_g: float,
    d_mix: float,
    species: Mapping[str, pr.Species],
    td: tr.TransportData,
    triplets: Triplets | None = None,
    T_eq: float | None = None,
) -> tuple[StepCoefficients, StepCoefficients, StepCoefficients]:
    """Coefficients of the characteristic times of Marsigny et al. (2025), Table A1.

    Args:
        pp: Pellet parameters.
        T: Temperature [K].
        p: Pressure [Pa].
        k_g: External mass-transfer coefficient of H2 [m/s].
        d_mix: Molecular diffusivity of H2 in the gas [m2/s] (used in the pores).
        species: Thermodynamic data.
        td: Transport data (molar mass of H2 for Knudsen diffusion).
        triplets: Kinetic constants; defaults to ``pp.triplets``.
        T_eq: Temperature [K] for the equilibrium constants (default T); used to stay inside
            the thermodynamic data range (H2O gas from 500 K) at low temperatures.
    """
    trip = pp.triplets if triplets is None else triplets
    ct = p / (pr.R_GAS * T)
    c_app = apparent_molar_densities(pp)
    c_true = true_molar_densities(pp)
    eps_m, eps_w, eps_fe = layer_porosities(pp)
    d_p, d_g, d_cr = pp.diameter, grain_diameter(pp), pp.d_crystallite
    m_h2 = td.molar_mass["H2"]

    def d_eff(eps: float, d_pore: float) -> float:
        return tr.effective_diffusivity(
            eps, pp.tortuosity, d_mix, tr.knudsen_diffusivity(d_pore, T, m_h2)
        )

    d_interg1 = d_eff(eps_m, pp.d_pores_intergrain)
    d_interg2 = d_eff(eps_w, pp.d_pores_intergrain)
    d_intrag2 = d_eff(pp.porosity_wustite_grain, pp.d_pores_crystallite)
    d_interc3 = d_eff(eps_fe, pp.d_pores_intergrain)
    b1, b2, b3 = B_SOLID_PER_H2
    k = [trip[j].k(T) for j in range(3)]
    t_eq = T if T_eq is None else T_eq
    ratio = [float(pr.equilibrium_h2o_h2_ratio(r, species, t_eq)) for r in REACTIONS]
    s1 = StepCoefficients(
        a_chem=c_true[0] * d_g / (2 * b1 * k[0] * ct),
        a_diff=c_app[0] * d_p**2 / (24 * b1 * d_interg1 * ct),
        a_ext=c_app[0] * d_p / (6 * b1 * k_g * ct),
        ratio=ratio[0],
    )
    s2 = StepCoefficients(
        a_chem=c_true[1] * d_g / (2 * k[1] * ct),  # factor 2 as printed (see module docstring)
        a_diff=c_app[1] * d_p**2 / (24 * b2 * d_interg2 * ct)
        + c_true[1] * d_g**2 / (24 * b2 * d_intrag2 * ct),
        a_ext=c_app[1] * d_p / (6 * b2 * k_g * ct),
        ratio=ratio[1],
    )
    s3 = StepCoefficients(
        a_chem=c_true[2] * d_cr / (2 * b3 * k[2] * ct),
        a_diff=c_app[2] * d_p**2 / (24 * b3 * d_interc3 * ct),
        a_ext=c_app[2] * d_p / (6 * b3 * k_g * ct),
        ratio=ratio[2],
    )
    return s1, s2, s3


def characteristic_times(
    pp: PelletParameters,
    env: Environment,
    species: Mapping[str, pr.Species],
    td: tr.TransportData,
) -> tuple[StepTimes, StepTimes, StepTimes]:
    """Characteristic times of the three steps (Marsigny 2025, Table A1) [s]."""
    gas = env.gas
    d_mix = tr.mixture_diffusivity(td, "H2", gas.x, gas.T, gas.p, trace_partner="H2O")
    coeffs = step_coefficients(pp, gas.T, gas.p, env.k_g, d_mix, species, td)
    xa, xc = gas.x.get("H2", 0.0), gas.x.get("H2O", 0.0)
    return coeffs[0].times(xa, xc), coeffs[1].times(xa, xc), coeffs[2].times(xa, xc)


def lart_rate(times: StepTimes, X: float) -> float:
    """dX/dt [1/s] of one step from the law of additive reaction times at conversion X."""
    if X >= 1.0 or math.isinf(times.chem):
        return 0.0
    one = max(1.0 - X, 1e-15)
    f_prime = (1.0 / 3.0) * one ** (-2.0 / 3.0)
    g_prime = 2.0 * (one ** (-1.0 / 3.0) - 1.0)
    return float(1.0 / (times.chem * f_prime + times.diff * g_prime + times.ext))


def step_rates(times: tuple[StepTimes, StepTimes, StepTimes], X: FloatArray) -> FloatArray:
    """dX_j/dt [1/s] for the three steps with the availability constraint X3 <= X2 <= X1."""
    r = np.array([lart_rate(times[j], float(X[j])) for j in range(3)])
    for j in (1, 2):
        if X[j] >= X[j - 1]:
            r[j] = min(r[j], r[j - 1])
    return r


def step_rates_array(
    coeffs: tuple[StepCoefficients, StepCoefficients, StepCoefficients],
    X: FloatArray,
    x_h2: FloatArray,
    x_h2o: FloatArray,
    multiplier: float = 1.0,
    one_minus: FloatArray | None = None,
    availability_width: float = 0.0,
    min_width: float = 0.0,
) -> FloatArray:
    """Vectorised :func:`step_rates` for many cells at once [1/s].

    Args:
        coeffs: Step coefficients (same temperature and pressure in all cells).
        X: Step conversions, shape (3, n).
        x_h2, x_h2o: Bulk mole fractions, shape (n,).
        multiplier: Scales all rates (1 = physical; 0 = no reaction; large = towards local
            equilibrium). For tests and limiting cases only.
        one_minus: Optional 1 - X, shape (3, n), given directly instead of computed from X
            (full precision near complete conversion; shaftsim.coupled passes u^3). The
            results are the same as with X up to rounding.
        availability_width: Width delta [-] over which the availability constraint fades in
            as the gap X_(j-1) - X_j closes (smoothstep); 0 = the original switch, which is
            discontinuous where X_j reaches X_(j-1) with r_j > r_(j-1).
        min_width: Relative width eta of the smooth minimum used in the constraint,
            min(a, b) ~ (a^-p + b^-p)^(-1/p) with p = 1/eta (relative error at most
            eta ln 2, exactly 0 if a or b is 0); 0 = exact min.
            Both widths are used by the coupled Newton solver (shaftsim.coupled), whose
            convergence needs a smooth residual; their effect on the outputs is tested.

    Returns:
        dX_j/dt, shape (3, n), with the availability constraint X3 <= X2 <= X1.
    """
    n = X.shape[1]
    om = 1.0 - X if one_minus is None else np.asarray(one_minus)
    r = np.zeros((3, n))
    for j, c in enumerate(coeffs):
        x_eq = (x_h2 + x_h2o) / (1.0 + c.ratio)
        d_eq = x_h2 - x_eq
        one = np.maximum(om[j], 1e-15 if one_minus is None else 1e-300)
        f_prime = (1.0 / 3.0) * one ** (-2.0 / 3.0)
        g_prime = 2.0 * (one ** (-1.0 / 3.0) - 1.0)
        d_chem = x_h2 - x_h2o / c.ratio
        ok = (d_eq > 0.0) & (d_chem > 0.0) & (om[j] > 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            denom = (c.a_chem / d_chem) * f_prime + (c.a_diff / d_eq) * g_prime + c.a_ext / d_eq
            r[j] = np.where(ok, 1.0 / denom, 0.0)
    r *= multiplier
    for j in (1, 2):
        if availability_width == 0.0 and min_width == 0.0:
            lim = om[j] <= om[j - 1]  # X_j >= X_(j-1): step j limited by step j-1
            r[j] = np.where(lim, np.minimum(r[j], r[j - 1]), r[j])
            continue
        a, b = r[j], r[j - 1]
        if min_width > 0.0:
            # p-norm minimum (a^-p + b^-p)^(-1/p), p = 1/eta, in the overflow-free form
            # lo (1 + (lo/hi)^p)^(-1/p): analytic for a, b > 0, exactly 0 when either rate is 0
            # (a step at equilibrium stays at rest), and lo 2^(-eta) <= m <= lo
            lo, hi = np.minimum(a, b), np.maximum(a, b)
            with np.errstate(divide="ignore", invalid="ignore"):
                q = np.where(hi > 0.0, lo / hi, 0.0)
            m = lo * (1.0 + q ** (1.0 / min_width)) ** (-min_width)
        else:
            m = np.minimum(a, b)
        gap = om[j] - om[j - 1]  # X_(j-1) - X_j: product of step j-1 still available
        if availability_width > 0.0:
            t = np.clip(1.0 - gap / availability_width, 0.0, 1.0)
            w = t * t * (3.0 - 2.0 * t)  # smoothstep: 1 at gap <= 0, 0 at gap >= width
        else:
            w = (gap <= 0.0).astype(float)
        r[j] = a + (m - a) * w
    return r


def overall_conversion(X: FloatArray) -> FloatArray | float:
    """Fraction of removable oxygen removed [-] from the step conversions (last axis or 0)."""
    w = np.asarray(OXYGEN_WEIGHTS)
    Xa = np.asarray(X, dtype=float)
    out = np.tensordot(w, Xa, axes=(0, 0))
    return float(out) if np.ndim(out) == 0 else np.asarray(out)


def solid_moles(pp: PelletParameters, X: FloatArray) -> dict[str, float]:
    """Moles of Fe2O3, Fe3O4, Fe0.95O, Fe in the pellet [mol] for step conversions X."""
    n0 = initial_moles_fe2o3(pp)
    x1, x2, x3 = (float(v) for v in X)
    return {
        "Fe2O3": n0 * (1 - x1),
        "Fe3O4": 2.0 / 3.0 * n0 * (x1 - x2),
        "Fe0.95O": 2.0 / X_FE_WUSTITE_KINETIC * n0 * (x2 - x3),
        "Fe": 2.0 * n0 * x3,
    }


def solid_mass_fractions(pp: PelletParameters, X: FloatArray) -> dict[str, FloatArray]:
    """Mass fractions of the solid [-]: Fe2O3, Fe3O4, Fe0.95O, Fe and gangue.

    Args:
        pp: Pellet parameters (the gangue is the non-hematite part, 1 - w_Fe2O3).
        X: Step conversions, shape (3,) or (3, n).

    Returns:
        Arrays of shape (n,) (or 0-d) that add up to 1.
    """
    n0 = initial_moles_fe2o3(pp)
    x1, x2, x3 = (np.asarray(v, dtype=float) for v in X)
    mm = pp.molar_mass
    mass = {
        "Fe2O3": n0 * (1 - x1) * mm["Fe2O3"],
        "Fe3O4": 2.0 / 3.0 * n0 * (x1 - x2) * mm["Fe3O4"],
        "Fe0.95O": 2.0 / X_FE_WUSTITE_KINETIC * n0 * (x2 - x3) * mm["Fe0.95O"],
        "Fe": 2.0 * n0 * x3 * mm["Fe"],
    }
    pellet_mass = pp.density * math.pi * pp.diameter**3 / 6.0
    mass["gangue"] = np.full_like(x1, pellet_mass * (1.0 - pp.w_Fe2O3))
    total = sum(mass.values())
    return {k: np.asarray(v / total) for k, v in mass.items()}


# --------------------------------------------------------------------------- simulation


@dataclass(frozen=True)
class PelletResult:
    """Time series of an isothermal pellet simulation.

    Attributes:
        t: Time [s].
        X_steps: Step conversions, shape (3, len(t)) [-].
        X: Overall conversion [-].
        x_H2_bulk: H2 mole fraction seen by the pellet [-].
        o_rate: Oxygen removal rate = H2 consumption = H2O production [mol/s].
    """

    t: FloatArray
    X_steps: FloatArray
    X: FloatArray
    x_H2_bulk: FloatArray
    o_rate: FloatArray


def _tga_environment(
    pp: PelletParameters,
    cond: TGAConditions,
    T: float,
    x: Mapping[str, float],
    td: tr.TransportData,
) -> Environment:
    """Gas state and Ranz-Marshall k_g in the thermobalance tube."""
    q = cond.gas_flow_normal * (T / td.normal_temperature) * (td.normal_pressure / cond.pressure)
    u = q / (math.pi * cond.furnace_diameter**2 / 4.0)
    m_mix = sum(xi * td.molar_mass[s] for s, xi in x.items())
    rho = cond.pressure * m_mix / (pr.R_GAS * T)
    mu = tr.mixture_viscosity(td, x, T)
    d_h2 = tr.mixture_diffusivity(td, "H2", x, T, cond.pressure, trace_partner="H2O")
    k_g = tr.ranz_marshall_kg(td.ranz_marshall, pp.diameter, u, rho, mu, d_h2)
    return Environment(GasState(T, cond.pressure, dict(x)), k_g)


def simulate_tga(
    pp: PelletParameters,
    cond: TGAConditions,
    T: float,
    t_end: float,
    species: Mapping[str, pr.Species],
    td: tr.TransportData,
    gas_mode: str = "mean",
    n_out: int = 400,
    triplets: Triplets | None = None,
) -> PelletResult:
    """Isothermal reduction of one pellet in a thermobalance gas flow.

    Args:
        pp: Pellet parameters.
        cond: Gas flow and composition.
        T: Temperature [K] (pellet and gas).
        t_end: Simulated time [s].
        species: Thermodynamic data (shaftsim.properties).
        td: Transport data (shaftsim.transport).
        gas_mode: Gas composition seen by the pellet: "static" (inlet composition, no
            depletion), "mean" (mean of inlet and outlet of the furnace balance) or "outlet"
            (well-mixed furnace). The H2 consumed equals the oxygen removed.
        n_out: Number of output times.
        triplets: Kinetic constants; defaults to ``pp.triplets``.

    Numerical choices: k_g and the pore diffusivities are evaluated once at the inlet gas
    composition (they depend weakly on the H2/H2O split); only the driving forces follow the
    depleted gas. The gas composition is solved at every right-hand-side evaluation from the
    scalar balance r = r_pellet(x_in - theta r/F) by Brent's method (the pellet rate decreases
    with depletion, so the root is unique). Integration: RK45, rtol 1e-6, max_step t_end/500
    so that the switches of the availability constraint are resolved.
    """
    theta = {"static": 0.0, "mean": 0.5, "outlet": 1.0}[gas_mode]
    n0 = initial_moles_fe2o3(pp)
    f_tot = td.normal_pressure * cond.gas_flow_normal / (pr.R_GAS * td.normal_temperature)
    x_in = {"H2": cond.x_H2, "H2O": cond.x_H2O, cond.inert: 1.0 - cond.x_H2 - cond.x_H2O}
    env_in = _tga_environment(pp, cond, T, x_in, td)
    d_mix = tr.mixture_diffusivity(td, "H2", x_in, T, cond.pressure, trace_partner="H2O")
    coeffs = step_coefficients(pp, T, cond.pressure, env_in.k_g, d_mix, species, td, triplets)

    def rates_at(dx: float, X: FloatArray) -> FloatArray:
        xa, xc = x_in["H2"] - dx, x_in["H2O"] + dx
        times = (coeffs[0].times(xa, xc), coeffs[1].times(xa, xc), coeffs[2].times(xa, xc))
        return step_rates(times, X)

    def o_rate(r: FloatArray) -> float:
        return O_PER_FE2O3 * n0 * float(np.dot(OXYGEN_WEIGHTS, r))

    def solve_gas(X: FloatArray) -> tuple[float, float, FloatArray]:
        """Return (oxygen removal rate [mol/s], H2 depletion dx [-], step rates [1/s])."""
        if theta == 0.0:
            r = rates_at(0.0, X)
            return o_rate(r), 0.0, r

        def g(rr: float) -> float:
            return rr - o_rate(rates_at(theta * rr / f_tot, X))

        rr = 0.0
        if g(0.0) < 0.0:
            rr = float(brentq(g, 0.0, f_tot * x_in["H2"] / theta, xtol=1e-16, rtol=1e-10))
        dx = theta * rr / f_tot
        r = rates_at(dx, X)
        return o_rate(r), dx, r

    def rhs(_t: float, y: FloatArray) -> FloatArray:
        return solve_gas(np.clip(y, 0.0, 1.0))[2]

    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(
        rhs,
        (0.0, t_end),
        np.zeros(3),
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9,
        max_step=t_end / 500.0,
    )
    if not sol.success:
        raise RuntimeError(f"pellet integration failed: {sol.message}")
    Xs = np.clip(sol.y, 0.0, 1.0)
    x_h2 = np.empty(len(sol.t))
    o_r = np.empty(len(sol.t))
    for i in range(len(sol.t)):
        rate, dx, _ = solve_gas(Xs[:, i])
        x_h2[i] = x_in["H2"] - dx
        o_r[i] = rate
    return PelletResult(sol.t, Xs, np.asarray(overall_conversion(Xs)), x_h2, o_r)


def load_refit_triplets(path: str | Path) -> Triplets:
    """Kinetic triplets re-fitted by scripts/fit_kinetics.py (data/kinetics_refit.yaml)."""
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)["kinetic_triplets"]
    t = [
        KineticTriplet(float(raw[f"step_{j}"]["k0"]), float(raw[f"step_{j}"]["Ea"]))
        for j in (1, 2, 3)
    ]
    return (t[0], t[1], t[2])


def step_rates_scalar(
    coeffs: tuple[StepCoefficients, StepCoefficients, StepCoefficients],
    x1: float,
    x2: float,
    x3: float,
    x_h2: float,
    x_h2o: float,
    multiplier: float = 1.0,
) -> tuple[float, float, float]:
    """Pure-Python version of :func:`step_rates_array` for one cell [1/s].

    Same mathematics, without numpy overhead: it is called millions of times inside the
    cell-by-cell root finding of shaftsim.solver.
    """
    X = (x1, x2, x3)
    r = [0.0, 0.0, 0.0]
    for j in range(3):
        c = coeffs[j]
        xj = X[j]
        if xj >= 1.0:
            continue
        x_eq = (x_h2 + x_h2o) / (1.0 + c.ratio)
        d_eq = x_h2 - x_eq
        d_chem = x_h2 - x_h2o / c.ratio
        if d_eq <= 0.0 or d_chem <= 0.0:
            continue
        one = max(1.0 - xj, 1e-15)
        f_prime = (1.0 / 3.0) * one ** (-2.0 / 3.0)
        g_prime = 2.0 * (one ** (-1.0 / 3.0) - 1.0)
        r[j] = multiplier / (
            (c.a_chem / d_chem) * f_prime + (c.a_diff / d_eq) * g_prime + c.a_ext / d_eq
        )
    if x2 >= x1:
        r[1] = min(r[1], r[0])
    if x3 >= x2:
        r[2] = min(r[2], r[1])
    return r[0], r[1], r[2]
