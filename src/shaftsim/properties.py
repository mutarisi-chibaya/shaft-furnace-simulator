"""Thermodynamic properties of species and equilibrium of the iron-oxide reduction steps.

All values come from NIST Shomate coefficients stored in ``data/species.yaml``; this module
contains no physical data. Units are SI throughout: T in K, Cp and S in J/(mol K),
H and G in J/mol. The only constant used is the molar gas constant R (CODATA, via
``scipy.constants``).

Reference state
---------------
``enthalpy`` returns the *formation-referenced* molar enthalpy, H(T) = dfH(298.15 K) +
[H(T) - H(298.15 K)], so that elements in their standard state have H(298.15 K) = 0.
``entropy`` returns the absolute (third-law) entropy. With these two, the Gibbs energy
G = H - T S is consistent between species, and reaction Gibbs energies follow directly from
sum(nu_i * G_i). Element contributions cancel because every reaction is element-balanced.

Numerical choices
-----------------
- Shomate polynomials are evaluated piecewise. A temperature exactly on a range boundary
  uses the lower range; NIST ranges share endpoints, so both agree to the precision of the
  published coefficients.
- Temperatures outside the tabulated ranges raise ``ValueError``: no silent extrapolation.
- Every function accepts a scalar or a numpy array of temperatures and returns the same shape.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from numpy.typing import ArrayLike, NDArray
from scipy.constants import R as _R_CODATA

FloatArray = NDArray[np.float64]

#: Molar gas constant [J/(mol K)], CODATA exact value via scipy.constants.
R_GAS: float = float(_R_CODATA)

_COEFF_NAMES = ("A", "B", "C", "D", "E", "F", "G", "H")
_RANGE_KEYS = ("T_min", "T_max", *_COEFF_NAMES)


class MissingDataError(ValueError):
    """Raised when a data file has fields that the human still has to fill in."""


@dataclass(frozen=True)
class ShomateRange:
    """Shomate coefficients for one temperature range, in NIST units.

    Attributes:
        T_min, T_max: Range limits [K].
        A..H: NIST Shomate coefficients (t = T/1000; Cp in J/(mol K); H, F in kJ/mol;
            G in J/(mol K)). ``H`` equals dfH at 298.15 K [kJ/mol].
    """

    T_min: float
    T_max: float
    A: float
    B: float
    C: float
    D: float
    E: float
    F: float
    G: float
    H: float


@dataclass(frozen=True)
class Species:
    """A chemical species with its piecewise Shomate data.

    Attributes:
        name: Identifier used everywhere (e.g. ``"Fe3O4"``).
        phase: ``"gas"`` or ``"solid"``.
        elements: Atoms per formula unit, e.g. ``{"Fe": 3, "O": 4}`` (may be
            non-integer for non-stoichiometric phases such as Fe0.947O).
        ranges: Shomate ranges, sorted and contiguous in temperature. All ranges must share
            the same H coefficient (dfH298), as NIST data do.
        molar_mass: Molar mass [kg/mol], or ``None`` if not yet provided.
        source: Key into ``data/SOURCES.md``.
    """

    name: str
    phase: str
    elements: Mapping[str, float]
    ranges: tuple[ShomateRange, ...]
    molar_mass: float | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        if not self.ranges:
            raise ValueError(f"{self.name}: no Shomate ranges")
        for r in self.ranges:
            if not r.T_max > r.T_min:
                raise ValueError(f"{self.name}: range with T_max <= T_min ({r.T_min}, {r.T_max})")
        if len({r.H for r in self.ranges}) > 1:
            raise ValueError(f"{self.name}: Shomate ranges disagree on H (dfH298)")
        for lo, hi in zip(self.ranges[:-1], self.ranges[1:], strict=True):
            if not np.isclose(lo.T_max, hi.T_min, rtol=0.0, atol=1e-6):
                raise ValueError(
                    f"{self.name}: Shomate ranges not contiguous ({lo.T_max} K vs {hi.T_min} K)"
                )

    @property
    def T_min(self) -> float:
        """Lowest temperature covered by the data [K]."""
        return self.ranges[0].T_min

    @property
    def T_max(self) -> float:
        """Highest temperature covered by the data [K]."""
        return self.ranges[-1].T_max


# --------------------------------------------------------------------------- loading


def species_from_dict(name: str, entry: Mapping[str, Any]) -> Species:
    """Build a :class:`Species` from one entry of ``species.yaml``.

    Args:
        name: Species key.
        entry: Mapping with ``phase``, ``elements``, ``shomate`` (list of ranges) and
            optionally ``molar_mass`` [kg/mol] and ``source``.

    Returns:
        The validated species.

    Raises:
        MissingDataError: If any required field is ``null``; the message lists all of them.
    """
    missing: list[str] = []
    for key in ("phase", "elements", "shomate"):
        if entry.get(key) is None:
            missing.append(f"{name}.{key}")
    ranges_raw = entry.get("shomate") or []
    for i, r in enumerate(ranges_raw):
        missing.extend(f"{name}.shomate[{i}].{k}" for k in _RANGE_KEYS if r.get(k) is None)
    if missing:
        raise MissingDataError("Missing values in species data: " + ", ".join(missing))

    ranges = tuple(ShomateRange(**{k: float(r[k]) for k in _RANGE_KEYS}) for r in ranges_raw)
    mm = entry.get("molar_mass")
    return Species(
        name=name,
        phase=str(entry["phase"]),
        elements={str(k): float(v) for k, v in entry["elements"].items()},
        ranges=ranges,
        molar_mass=None if mm is None else float(mm),
        source=entry.get("source"),
    )


def load_species(path: str | Path) -> dict[str, Species]:
    """Load all species from a YAML file in the ``species.yaml`` format.

    Args:
        path: Path to the YAML file.

    Returns:
        Mapping from species name to :class:`Species`.

    Raises:
        MissingDataError: If any value is still ``null``. All missing fields of all species
            are reported at once.
    """
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    out: dict[str, Species] = {}
    missing: list[str] = []
    for name, entry in raw["species"].items():
        try:
            out[name] = species_from_dict(name, entry)
        except MissingDataError as exc:
            missing.append(str(exc).split(": ", 1)[1])
    if missing:
        raise MissingDataError(f"Missing values in {path}: " + ", ".join(missing))
    return out


# --------------------------------------------------------------------------- properties


def _coefficients(sp: Species, T: FloatArray) -> FloatArray:
    """Return the coefficient matrix (len(T), 8) selected per temperature."""
    if np.any(~np.isfinite(T)) or np.any(T < sp.T_min - 1e-9) or np.any(T > sp.T_max + 1e-9):
        raise ValueError(
            f"{sp.name}: temperature outside data range [{sp.T_min}, {sp.T_max}] K: "
            f"min {np.min(T)}, max {np.max(T)}"
        )
    upper = np.array([r.T_max for r in sp.ranges])
    idx = np.minimum(np.searchsorted(upper, T, side="left"), len(sp.ranges) - 1)
    table = np.array([[getattr(r, c) for c in _COEFF_NAMES] for r in sp.ranges])
    return table[idx]


def _prepare(T: ArrayLike) -> tuple[FloatArray, tuple[int, ...]]:
    arr = np.asarray(T, dtype=np.float64)
    return arr.ravel(), arr.shape


def _shape(values: FloatArray, shape: tuple[int, ...]) -> FloatArray | float:
    return float(values[0]) if shape == () else values.reshape(shape)


def heat_capacity(sp: Species, T: ArrayLike) -> FloatArray | float:
    """Molar isobaric heat capacity Cp [J/(mol K)] at temperature T [K]."""
    Tf, shape = _prepare(T)
    c = _coefficients(sp, Tf)
    t = Tf / 1000.0
    cp = c[:, 0] + c[:, 1] * t + c[:, 2] * t**2 + c[:, 3] * t**3 + c[:, 4] / t**2
    return _shape(cp, shape)


def enthalpy(sp: Species, T: ArrayLike) -> FloatArray | float:
    """Formation-referenced molar enthalpy H(T) = dfH298 + (H(T) - H298) [J/mol].

    In NIST form this is A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F (kJ/mol), because the
    ``- H`` term of H - H298 cancels with dfH298 = H.
    """
    Tf, shape = _prepare(T)
    c = _coefficients(sp, Tf)
    t = Tf / 1000.0
    h_kj = (
        c[:, 0] * t
        + c[:, 1] * t**2 / 2
        + c[:, 2] * t**3 / 3
        + c[:, 3] * t**4 / 4
        - c[:, 4] / t
        + c[:, 5]
    )
    return _shape(h_kj * 1000.0, shape)


def sensible_enthalpy(sp: Species, T: ArrayLike) -> FloatArray | float:
    """Molar enthalpy relative to 298.15 K, H(T) - H(298.15 K) [J/mol].

    Uses the NIST form H(T) - H298 = [...] + F - H, i.e. ``enthalpy`` minus the H coefficient
    (dfH298). This works even when 298.15 K is outside the tabulated range (H2O gas data
    start at 500 K), where evaluating H at 298.15 K would be an extrapolation.
    """
    Tf, shape = _prepare(T)
    dfh298 = _coefficients(sp, Tf)[:, 7] * 1000.0
    h = np.asarray(enthalpy(sp, Tf)) - dfh298
    return _shape(h, shape)


def entropy(sp: Species, T: ArrayLike) -> FloatArray | float:
    """Absolute molar entropy S(T) at standard pressure [J/(mol K)]."""
    Tf, shape = _prepare(T)
    c = _coefficients(sp, Tf)
    t = Tf / 1000.0
    s = (
        c[:, 0] * np.log(t)
        + c[:, 1] * t
        + c[:, 2] * t**2 / 2
        + c[:, 3] * t**3 / 3
        - c[:, 4] / (2 * t**2)
        + c[:, 6]
    )
    return _shape(s, shape)


def gibbs(sp: Species, T: ArrayLike) -> FloatArray | float:
    """Molar Gibbs energy G = H - T S at standard pressure [J/mol] (same reference as H)."""
    Tf, shape = _prepare(T)
    g = np.asarray(enthalpy(sp, Tf)) - Tf * np.asarray(entropy(sp, Tf))
    return _shape(g, shape)


# --------------------------------------------------------------------------- reactions


@dataclass(frozen=True)
class Reaction:
    """A reaction defined by stoichiometric coefficients (negative = reactant).

    Attributes:
        name: Identifier, e.g. ``"hematite->magnetite"``.
        stoich: Mapping species name -> stoichiometric coefficient [-].
    """

    name: str
    stoich: Mapping[str, float]


#: Wüstite composition used throughout: Fe(x)O with x = 0.947, the composition of the JANAF
#: wüstite table (Fe-001) from which its data are fitted. Real wüstite varies (roughly
#: x = 0.83-0.95 with T and oxygen potential); a fixed x is a modelling assumption.
#: Must match the ``Fe0.947O`` entry in data/species.yaml (checked by a test).
X_WUSTITE = 0.947
WUSTITE = "Fe0.947O"

#: Reduction steps by hydrogen. Stoichiometry only (pure chemistry, no fitted data).
#: Above the wüstite stability limit: Fe2O3 -> Fe3O4 -> Fe(x)O -> Fe. Below it wüstite is
#: unstable and magnetite reduces directly to iron (MAGNETITE_TO_IRON). Which regime applies
#: follows from the data, see :func:`wustite_is_stable`.
HEMATITE_TO_MAGNETITE = Reaction(
    "hematite->magnetite", {"Fe2O3": -3, "H2": -1, "Fe3O4": 2, "H2O": 1}
)
#: Fe3O4 + (4 - 3/x) H2 -> (3/x) Fe(x)O + (4 - 3/x) H2O
MAGNETITE_TO_WUSTITE = Reaction(
    "magnetite->wustite",
    {
        "Fe3O4": -1,
        "H2": -(4 - 3 / X_WUSTITE),
        WUSTITE: 3 / X_WUSTITE,
        "H2O": 4 - 3 / X_WUSTITE,
    },
)
#: Fe(x)O + H2 -> x Fe + H2O
WUSTITE_TO_IRON = Reaction("wustite->iron", {WUSTITE: -1, "H2": -1, "Fe": X_WUSTITE, "H2O": 1})
MAGNETITE_TO_IRON = Reaction("magnetite->iron", {"Fe3O4": -1, "H2": -4, "Fe": 3, "H2O": 4})

REDUCTION_STEPS: tuple[Reaction, ...] = (
    HEMATITE_TO_MAGNETITE,
    MAGNETITE_TO_WUSTITE,
    WUSTITE_TO_IRON,
)

#: The same steps with idealised stoichiometric FeO, for comparison only. With these data
#: FeO appears stable down to ~460 K, which contradicts the observed ~570 degC limit; this
#: is why the default steps use Fe0.947O (see docs/validation.md).
MAGNETITE_TO_FEO = Reaction("magnetite->FeO", {"Fe3O4": -1, "H2": -1, "FeO": 3, "H2O": 1})
FEO_TO_IRON = Reaction("FeO->iron", {"FeO": -1, "H2": -1, "Fe": 1, "H2O": 1})


def element_imbalance(rxn: Reaction, species: Mapping[str, Species]) -> dict[str, float]:
    """Net atoms of each element produced by one unit of reaction [mol/mol]; zero if balanced."""
    net: dict[str, float] = {}
    for name, nu in rxn.stoich.items():
        for el, n in species[name].elements.items():
            net[el] = net.get(el, 0.0) + nu * n
    return net


def _reaction_sum(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike, prop: Any
) -> FloatArray | float:
    Tf, shape = _prepare(T)
    total = np.zeros_like(Tf)
    for name, nu in rxn.stoich.items():
        total += nu * np.asarray(prop(species[name], Tf))
    return _shape(total, shape)


def reaction_enthalpy(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike
) -> FloatArray | float:
    """Standard reaction enthalpy dH_r(T) = sum(nu_i H_i) [J per mol of reaction].

    Positive means endothermic.
    """
    return _reaction_sum(rxn, species, T, enthalpy)


def reaction_entropy(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike
) -> FloatArray | float:
    """Standard reaction entropy dS_r(T) = sum(nu_i S_i) [J/(K mol of reaction)]."""
    return _reaction_sum(rxn, species, T, entropy)


def reaction_gibbs(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike
) -> FloatArray | float:
    """Standard reaction Gibbs energy dG_r(T) = sum(nu_i G_i) [J per mol of reaction]."""
    return _reaction_sum(rxn, species, T, gibbs)


def equilibrium_constant(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike
) -> FloatArray | float:
    """Equilibrium constant K(T) = exp(-dG_r / (R T)) [-].

    Solids are taken as pure phases (activity 1) and gases as ideal at standard pressure,
    so for the reduction steps K = (p_H2O / p_H2)^n at equilibrium, with n the H2 coefficient.
    """
    Tf, shape = _prepare(T)
    dg = np.asarray(reaction_gibbs(rxn, species, Tf))
    return _shape(np.exp(-dg / (R_GAS * Tf)), shape)


def equilibrium_h2o_h2_ratio(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike
) -> FloatArray | float:
    """Equilibrium gas ratio p_H2O / p_H2 [-] for a reduction step.

    For n H2 consumed and n H2O produced, K = (p_H2O/p_H2)^n, so the ratio is K^(1/n).
    Above this ratio the step cannot proceed; this is why furnaces need excess hydrogen.

    Raises:
        ValueError: If the reaction does not consume H2 and produce H2O in equal amounts.
    """
    n_h2 = -rxn.stoich.get("H2", 0.0)
    n_h2o = rxn.stoich.get("H2O", 0.0)
    if n_h2 <= 0 or not np.isclose(n_h2, n_h2o):
        raise ValueError(f"{rxn.name}: not an H2 -> H2O reduction with equal coefficients")
    Tf, shape = _prepare(T)
    k = np.asarray(equilibrium_constant(rxn, species, Tf))
    return _shape(k ** (1.0 / n_h2), shape)


def equilibrium_h2_fraction(
    rxn: Reaction, species: Mapping[str, Species], T: ArrayLike
) -> FloatArray | float:
    """Equilibrium H2 mole fraction y_H2 = 1 / (1 + p_H2O/p_H2) in a dry H2-H2O gas [-]."""
    Tf, shape = _prepare(T)
    ratio = np.asarray(equilibrium_h2o_h2_ratio(rxn, species, Tf))
    return _shape(1.0 / (1.0 + ratio), shape)


def wustite_disproportionation(species: Mapping[str, Species], wustite: str = WUSTITE) -> Reaction:
    """4 Fe(x)O -> (4x - 3) Fe + Fe3O4, with x read from the species' elements."""
    x = species[wustite].elements["Fe"]
    return Reaction("wustite-disproportionation", {wustite: -4, "Fe": 4 * x - 3, "Fe3O4": 1})


def wustite_is_stable(
    species: Mapping[str, Species], T: ArrayLike, wustite: str = WUSTITE
) -> NDArray[np.bool_] | bool:
    """Whether wüstite is thermodynamically stable relative to Fe + Fe3O4 at T [K].

    Wüstite is stable when its disproportionation (see :func:`wustite_disproportionation`)
    has dG_r > 0. The literature limit is about 570 degC; the value here follows from the
    data (see docs/validation.md for the number obtained).
    """
    dg = np.asarray(reaction_gibbs(wustite_disproportionation(species, wustite), species, T))
    return bool(dg > 0) if dg.shape == () else dg > 0
