"""Gas transport properties and gas-solid mass transfer.

All inputs and outputs are SI: T in K, p in Pa, diameters in m, diffusivities in m2/s,
viscosity in Pa s, mass-transfer coefficient in m/s. Molar masses in kg/mol.

Physics and numerical choices
-----------------------------
- Binary diffusivity: first Chapman-Enskog approximation for a Lennard-Jones 12-6 gas
  (Hirschfelder, Curtiss & Bird 1954, eq. 8.2-44):
      D_ij = (3/16) sqrt(2 pi (k T)^3 / mu_ij) / (p pi sigma_ij^2 Omega11*(T*_ij)),
  with mu_ij the reduced molecular mass, sigma_ij = (sigma_i + sigma_j)/2 and
  eps_ij = sqrt(eps_i eps_j) (Lorentz-Berthelot combining rules), T* = k T / eps_ij.
  This is identical to the familiar 0.0018583 T^1.5 ... form but uses only CODATA constants.
- Pure-gas viscosity: eta = (5/16) sqrt(pi m k T) / (pi sigma^2 Omega22*(T*)).
- Collision integrals Omega11*, Omega22* are interpolated (cubic, in log T* vs log Omega) from
  the Monchick-Mason table in data/collision_integrals.csv. Outside 0.1 <= T* <= 100 an error
  is raised.
- Mixture viscosity: Wilke's rule. Diffusivity of A in a mixture: Blanc's law.
- Pore diffusion: Knudsen D_K = (d/3) sqrt(8 R T / (pi M)), combined with molecular diffusion
  in series (Bosanquet): 1/D = 1/D_mix + 1/D_K, then D_eff = (eps / tau) D. Note: Marsigny et
  al. (2025, Appendix A) print D_eff = (eps/tau)(1/D_mix + 1/D_K), which is dimensionally
  inconsistent; the standard Bosanquet form used here is what that expression must mean.
- External mass transfer: Ranz-Marshall, Sh = a + b Re^m Sc^n (constants in
  data/transport.yaml).
"""

from __future__ import annotations

import csv
import math
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import yaml
from scipy.constants import Avogadro as N_A  # noqa: N814 (physics symbol)
from scipy.constants import Boltzmann as K_B  # noqa: N814 (physics symbol)
from scipy.constants import R as _R
from scipy.interpolate import CubicSpline

R_GAS: float = float(_R)
#: Upper limit of the power-law extrapolation of the collision integrals [-].
TSTAR_EXTRAPOLATION_MAX = 500.0


@dataclass(frozen=True)
class LennardJones:
    """Lennard-Jones 12-6 force constants of one species.

    Attributes:
        sigma: Collision diameter [m].
        eps_over_k: Well depth divided by the Boltzmann constant [K].
    """

    sigma: float
    eps_over_k: float


@dataclass(frozen=True)
class RanzMarshall:
    """Constants of Sh = a + b Re^m Sc^n [-] (also used for the Wakao-Kaguei bed form)."""

    a: float
    b: float
    m: float
    n: float


@dataclass(frozen=True)
class CollisionIntegrals:
    """Reduced collision integrals vs reduced temperature (non-polar 12-6 potential).

    Attributes:
        tstar: Reduced temperatures T* = kT/eps [-], increasing.
        omega11: Omega(1,1)* [-] (diffusion).
        omega22: Omega(2,2)* [-] (viscosity).
    """

    tstar: tuple[float, ...]
    omega11: tuple[float, ...]
    omega22: tuple[float, ...]

    def _interp(self, values: tuple[float, ...], tstar: float) -> float:
        """Cubic interpolation in (log T*, log Omega). Above the table (T* > 100, e.g. He at
        furnace temperatures) the last interval is extended as a power law up to
        T* = TSTAR_EXTRAPOLATION_MAX; Omega varies slowly there (Omega ~ T*^-0.13)."""
        lo, hi = self.tstar[0], self.tstar[-1]
        if lo <= tstar <= hi:
            spline = self._splines[0] if values is self.omega11 else self._splines[1]
            return float(np.exp(spline(math.log(tstar))))
        if hi < tstar <= TSTAR_EXTRAPOLATION_MAX:
            slope = math.log(values[-1] / values[-2]) / math.log(self.tstar[-1] / self.tstar[-2])
            return float(values[-1] * (tstar / hi) ** slope)
        raise ValueError(f"T* = {tstar:.4g} outside [{lo}, {TSTAR_EXTRAPOLATION_MAX}]")

    @cached_property
    def _splines(self) -> tuple[CubicSpline, CubicSpline]:
        """Cubic splines of log Omega vs log T* (built once)."""
        lt = np.log(self.tstar)
        return CubicSpline(lt, np.log(self.omega11)), CubicSpline(lt, np.log(self.omega22))

    def omega11_at(self, tstar: float) -> float:
        """Omega(1,1)*(T*) [-]."""
        return self._interp(self.omega11, tstar)

    def omega22_at(self, tstar: float) -> float:
        """Omega(2,2)*(T*) [-]."""
        return self._interp(self.omega22, tstar)


@dataclass(frozen=True)
class MarreroMason:
    """Constants of the Marrero-Mason (1972) correlation eq. (4.3-1) for one gas pair.

    Attributes:
        A: Pre-factor [atm cm2 s-1 K^-s].
        s: Temperature exponent [-].
        phi0_k: phi0/k [K].
        S, S2: Sutherland-type constants [K], [K^2].
        T_min, T_max: Validity range [K].
    """

    A: float
    s: float
    phi0_k: float
    S: float
    S2: float
    T_min: float
    T_max: float

    def diffusivity(self, T: float, p: float) -> float:
        """D12 [m2/s] at T [K], p [Pa]."""
        if not self.T_min <= T <= self.T_max:
            raise ValueError(f"T = {T} K outside correlation range [{self.T_min}, {self.T_max}]")
        pd_atm_cm2 = math.exp(
            math.log(self.A)
            + self.s * math.log(T)
            - math.log(math.log(self.phi0_k / T) ** 2)
            - self.S / T
            - self.S2 / T**2
        )
        return pd_atm_cm2 * 1e-4 * 101325.0 / p


@dataclass(frozen=True)
class TransportData:
    """Everything needed for transport properties: force constants, molar masses, tables."""

    lj: Mapping[str, LennardJones]
    molar_mass: Mapping[str, float]
    collision: CollisionIntegrals
    ranz_marshall: RanzMarshall
    normal_temperature: float
    normal_pressure: float
    correlations: Mapping[frozenset[str], MarreroMason]
    wakao_kaguei: RanzMarshall | None = None


def load_collision_integrals(path: str | Path) -> CollisionIntegrals:
    """Read data/collision_integrals.csv (comment lines start with '#')."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    return CollisionIntegrals(
        tstar=tuple(float(r["tstar"]) for r in rows),
        omega11=tuple(float(r["omega11"]) for r in rows),
        omega22=tuple(float(r["omega22"]) for r in rows),
    )


def load_transport(
    transport_yaml: str | Path,
    collision_csv: str | Path,
    molar_mass: Mapping[str, float],
) -> TransportData:
    """Load transport data.

    Args:
        transport_yaml: Path to data/transport.yaml.
        collision_csv: Path to data/collision_integrals.csv.
        molar_mass: Molar masses [kg/mol] by species (e.g. from properties.Species).
    """
    with open(transport_yaml, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    lj = {
        k: LennardJones(float(v["sigma"]), float(v["eps_over_k"]))
        for k, v in raw["lennard_jones"].items()
    }
    rm = raw["ranz_marshall"]
    nc = raw["normal_conditions"]
    corr = {
        frozenset(pair.split("-")): MarreroMason(
            A=float(c["A_1e3"]) * 1e-3,
            s=float(c["s"]),
            phi0_k=float(c["phi0_k_1e8"]) * 1e8,
            S=float(c["S"]),
            S2=float(c["S2"]),
            T_min=float(c["T_min"]),
            T_max=float(c["T_max"]),
        )
        for pair, c in (raw.get("marrero_mason") or {}).items()
    }
    return TransportData(
        lj=lj,
        molar_mass=dict(molar_mass),
        collision=load_collision_integrals(collision_csv),
        ranz_marshall=RanzMarshall(float(rm["a"]), float(rm["b"]), float(rm["m"]), float(rm["n"])),
        normal_temperature=float(nc["temperature"]),
        normal_pressure=float(nc["pressure"]),
        correlations=corr,
        wakao_kaguei=(
            RanzMarshall(*(float(raw["wakao_kaguei"][k]) for k in ("a", "b", "m", "n")))
            if "wakao_kaguei" in raw
            else None
        ),
    )


# --------------------------------------------------------------------------- molecular


def binary_diffusivity(td: TransportData, a: str, b: str, T: float, p: float) -> float:
    """Binary diffusion coefficient D_ab [m2/s] at T [K], p [Pa].

    Uses the experimental Marrero-Mason correlation when one exists for the pair, otherwise
    the Chapman-Enskog value (:func:`chapman_enskog_diffusivity`).
    """
    corr = td.correlations.get(frozenset((a, b)))
    if corr is not None:
        return corr.diffusivity(T, p)
    return chapman_enskog_diffusivity(td, a, b, T, p)


def chapman_enskog_diffusivity(td: TransportData, a: str, b: str, T: float, p: float) -> float:
    """First-order Chapman-Enskog binary diffusion coefficient D_ab [m2/s] at T [K], p [Pa]."""
    la, lb = td.lj[a], td.lj[b]
    sigma = 0.5 * (la.sigma + lb.sigma)
    eps_k = math.sqrt(la.eps_over_k * lb.eps_over_k)
    ma, mb = td.molar_mass[a] / N_A, td.molar_mass[b] / N_A
    mu = ma * mb / (ma + mb)
    omega = td.collision.omega11_at(T / eps_k)
    return (
        3.0
        / 16.0
        * math.sqrt(2.0 * math.pi * (K_B * T) ** 3 / mu)
        / (p * math.pi * sigma**2 * omega)
    )


def viscosity(td: TransportData, species: str, T: float) -> float:
    """Chapman-Enskog dynamic viscosity of a pure gas [Pa s] at T [K]."""
    lj = td.lj[species]
    m = td.molar_mass[species] / N_A
    omega = td.collision.omega22_at(T / lj.eps_over_k)
    return 5.0 / 16.0 * math.sqrt(math.pi * m * K_B * T) / (math.pi * lj.sigma**2 * omega)


def mixture_viscosity(td: TransportData, x: Mapping[str, float], T: float) -> float:
    """Wilke's mixing rule for the viscosity of a gas mixture [Pa s].

    eta = sum_i x_i eta_i / sum_j x_j phi_ij,
    phi_ij = [1 + (eta_i/eta_j)^(1/2) (M_j/M_i)^(1/4)]^2 / sqrt(8 (1 + M_i/M_j)).
    """
    names = [s for s, xi in x.items() if xi > 0]
    eta = {s: viscosity(td, s, T) for s in names}
    m = td.molar_mass
    total = 0.0
    for i in names:
        denom = 0.0
        for j in names:
            phi = (1 + math.sqrt(eta[i] / eta[j]) * (m[j] / m[i]) ** 0.25) ** 2 / math.sqrt(
                8 * (1 + m[i] / m[j])
            )
            denom += x[j] * phi
        total += x[i] * eta[i] / denom
    return total


def mixture_diffusivity(
    td: TransportData, a: str, x: Mapping[str, float], T: float, p: float, trace_partner: str
) -> float:
    """Diffusivity of species a in a mixture by Blanc's law [m2/s].

    D_a,mix = (1 - x_a) / sum_{i != a} x_i / D_ai. When a is (nearly) pure, the limit is
    the binary coefficient with `trace_partner` (the species it counter-diffuses with).
    """
    others = {s: xi for s, xi in x.items() if s != a and xi > 0}
    if sum(others.values()) < 1e-12:
        return binary_diffusivity(td, a, trace_partner, T, p)
    return (1 - x.get(a, 0.0)) / sum(
        xi / binary_diffusivity(td, a, s, T, p) for s, xi in others.items()
    )


#: Ratio rho D / eta used for the internal-energy contribution to the conductivity
#: (modified Eucken correction): Svehla (1962), NASA TR R-132, Appendix B, eq. (B3) and p. 9
#: ("close to 1.32 [...] so this value was used for f_int").
F_INT_EUCKEN = 1.32


def thermal_conductivity(td: TransportData, species: str, T: float, cp_molar: float) -> float:
    """Thermal conductivity of a pure gas [W/(m K)], Svehla (1962) eq. (B3).

    lambda = (R/M) [15/4 + f_int (Cp/R - 5/2)] eta, with eta the Chapman-Enskog viscosity,
    Cp the ideal-gas molar heat capacity [J/(mol K)] (from shaftsim.properties) and
    f_int = 1.32. For a monatomic gas Cp/R = 5/2 and only the translational term remains.
    """
    eta = viscosity(td, species, T)
    m = td.molar_mass[species]
    return R_GAS / m * (15.0 / 4.0 + F_INT_EUCKEN * (cp_molar / R_GAS - 2.5)) * eta


def mixture_conductivity(
    td: TransportData, x: Mapping[str, float], T: float, cp_molar: Mapping[str, float]
) -> float:
    """Mixture thermal conductivity [W/(m K)] by Wilke's rule (as Wagner 2008, Sec. III.5.7).

    lambda = sum_i x_i lambda_i / sum_j x_j phi_ij with Wilke's phi_ij from the viscosities.
    """
    names = [s for s, xi in x.items() if xi > 0]
    eta = {s: viscosity(td, s, T) for s in names}
    lam = {s: thermal_conductivity(td, s, T, cp_molar[s]) for s in names}
    m = td.molar_mass
    total = 0.0
    for i in names:
        denom = 0.0
        for j in names:
            phi = (1 + math.sqrt(eta[i] / eta[j]) * (m[j] / m[i]) ** 0.25) ** 2 / math.sqrt(
                8 * (1 + m[i] / m[j])
            )
            denom += x[j] * phi
        total += x[i] * lam[i] / denom
    return total


# --------------------------------------------------------------------------- pores


def knudsen_diffusivity(d_pore: float, T: float, molar_mass: float) -> float:
    """Knudsen diffusion coefficient D_K = (d/3) sqrt(8 R T / (pi M)) [m2/s]."""
    return d_pore / 3.0 * math.sqrt(8.0 * R_GAS * T / (math.pi * molar_mass))


def effective_diffusivity(porosity: float, tortuosity: float, d_mix: float, d_kn: float) -> float:
    """Effective pore diffusivity (eps/tau) / (1/D_mix + 1/D_K) [m2/s] (Bosanquet)."""
    return porosity / tortuosity / (1.0 / d_mix + 1.0 / d_kn)


# --------------------------------------------------------------------------- external


def ranz_marshall_kg(
    rm: RanzMarshall, d_p: float, velocity: float, density: float, viscosity_: float, d_ab: float
) -> float:
    """External mass-transfer coefficient k_g = Sh D / d_p [m/s] for a sphere.

    Args:
        rm: Correlation constants.
        d_p: Pellet diameter [m].
        velocity: Superficial gas velocity past the pellet [m/s].
        density: Gas mass density [kg/m3].
        viscosity_: Gas dynamic viscosity [Pa s].
        d_ab: Diffusivity of the transferred species in the gas [m2/s].
    """
    re = density * velocity * d_p / viscosity_
    sc = viscosity_ / (density * d_ab)
    sh = rm.a + rm.b * float(re) ** rm.m * float(sc) ** rm.n
    return float(sh * d_ab / d_p)
