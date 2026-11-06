"""Uncertainty quantification of the whole-furnace model (Milestone 6).

Uncertain inputs (data/uncertainty.yaml, sources in data/SOURCES.md):
- kinetics: the six parameters theta = (ln k_j(T_ref), Ea_j), j = 1..3, multivariate normal
  around the Milestone 2 refit with the jackknife covariance (data/kinetics_refit.yaml);
- pellet tortuosity, uniform on [1, 3];
- gas-solid transfer correlation, lambda uniform on [0, 1] between Wakao-Kaguei (0) and
  Ranz-Marshall (1), applied to both h and k_g.

Sampling space: the unit hypercube [0, 1]^8. The first six coordinates are mapped to
independent standard normals z_k (inverse normal CDF, truncated at +/- ``truncation`` sigma)
along the principal axes of the kinetic covariance, theta = mean + V sqrt(L) z; the last two
map linearly to tortuosity and lambda. The inputs are therefore independent, as Sobol indices
require; the kinetics are reported as one group (their principal components have no
physical meaning one by one).

Each sample reuses the property tables of the base problem
(shaftsim.energy.ThermalProblem.with_parameters): the parameters enter as exact factors.

Units: SI internally; outputs in the units named by their keys.
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import yaml
from numpy.typing import NDArray
from scipy.stats import norm, qmc

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr

FloatArray = NDArray[np.float64]
T_REF = 1023.15  # K, reference temperature of the kinetic parameterisation (fit_kinetics.py)
N_KINETIC = 6
NAMES = ("kin_pc1", "kin_pc2", "kin_pc3", "kin_pc4", "kin_pc5", "kin_pc6", "tortuosity", "transfer")
GROUPS = ("kinetics",) * N_KINETIC + ("tortuosity", "transfer correlation")
#: Scalar outputs of a run (see :func:`run_one`).
OUTPUTS = (
    "metallization",
    "h2_consumption_kg_per_t",
    "top_gas_T_C",
    "dri_T_C",
    "peak_solid_T_C",
    "gud",
)


@dataclass(frozen=True)
class UncertaintySpec:
    """Distributions of the uncertain inputs.

    Attributes:
        theta_mean: Kinetic parameters (ln k1, Ea1, ln k2, Ea2, ln k3, Ea3) of the refit;
            ln k in ln(m/s) at T_REF, Ea in J/mol.
        theta_cov: Their covariance (same order and units).
        truncation: Truncation of the kinetic normals [sigma].
        tortuosity_nominal: Tortuosity of the base problem [-].
        tortuosity_range: (low, high) [-].
        transfer_range: (low, high) of lambda [-].
    """

    theta_mean: FloatArray
    theta_cov: FloatArray
    truncation: float
    tortuosity_nominal: float
    tortuosity_range: tuple[float, float]
    transfer_range: tuple[float, float]

    @property
    def principal_axes(self) -> tuple[FloatArray, FloatArray]:
        """(sqrt of eigenvalues, eigenvectors) of the kinetic covariance."""
        w, v = np.linalg.eigh(self.theta_cov)
        return np.sqrt(np.clip(w, 0.0, None)), v

    def zero(self) -> UncertaintySpec:
        """The same spec with no uncertainty (every sample is the nominal case)."""
        return replace(
            self,
            theta_cov=np.zeros_like(self.theta_cov),
            tortuosity_range=(self.tortuosity_nominal, self.tortuosity_nominal),
            transfer_range=(0.0, 0.0),
        )


def load_uncertainty(
    uncertainty_yaml: str | Path, refit_yaml: str | Path, tortuosity_nominal: float
) -> UncertaintySpec:
    """Read data/uncertainty.yaml and the kinetic mean and covariance of the refit."""
    u = yaml.safe_load(Path(uncertainty_yaml).read_text(encoding="utf-8"))
    r = yaml.safe_load(Path(refit_yaml).read_text(encoding="utf-8"))
    cov_key = str(u["kinetics"]["covariance"]["value"])
    if cov_key not in r:
        raise ValueError(f"{refit_yaml} has no '{cov_key}'; rerun scripts/fit_kinetics.py")
    mean = []
    for j in (1, 2, 3):
        # from k0 and Ea (the values every deterministic run uses), not the rounded ln_k_Tref,
        # so that zero uncertainty reproduces the deterministic result exactly
        e = r["kinetic_triplets"][f"step_{j}"]
        k0, ea = float(e["k0"]), float(e["Ea"])
        mean += [math.log(k0) - ea / (pr.R_GAS * T_REF), ea]
    return UncertaintySpec(
        theta_mean=np.array(mean),
        theta_cov=np.array(r[cov_key], dtype=float),
        truncation=float(u["kinetics"]["truncation_sigma"]["value"]),
        tortuosity_nominal=tortuosity_nominal,
        tortuosity_range=(
            float(u["tortuosity"]["low"]["value"]),
            float(u["tortuosity"]["high"]["value"]),
        ),
        transfer_range=(
            float(u["transfer_correlation"]["low"]["value"]),
            float(u["transfer_correlation"]["high"]["value"]),
        ),
    )


def triplets_from_theta(theta: FloatArray) -> pl.Triplets:
    """(ln k_j(T_REF), Ea_j) -> kinetic triplets (k0_j [m/s], Ea_j [J/mol])."""
    out = []
    for j in range(3):
        ln_k, ea = float(theta[2 * j]), float(theta[2 * j + 1])
        out.append(pl.KineticTriplet(math.exp(ln_k + ea / (pr.R_GAS * T_REF)), ea))
    return out[0], out[1], out[2]


@dataclass(frozen=True)
class SampleParameters:
    """Physical parameters of one sample."""

    theta: FloatArray
    tortuosity: float
    transfer_mix: float

    @property
    def kinetics(self) -> pl.Triplets:
        """Kinetic triplets of this sample."""
        return triplets_from_theta(self.theta)


def to_parameters(spec: UncertaintySpec, u: FloatArray) -> SampleParameters:
    """Map a point of the unit hypercube [0, 1]^8 to physical parameters (module docstring)."""
    u = np.asarray(u, dtype=float)
    if u.shape != (len(NAMES),) or np.any(u < 0.0) or np.any(u > 1.0):
        raise ValueError("u must be a point of [0, 1]^8")
    p_lo = norm.cdf(-spec.truncation)
    z = norm.ppf(p_lo + (1.0 - 2.0 * p_lo) * u[:N_KINETIC])
    s, v = spec.principal_axes
    theta = spec.theta_mean + v @ (s * z)
    t_lo, t_hi = spec.tortuosity_range
    m_lo, m_hi = spec.transfer_range
    return SampleParameters(
        theta=theta,
        tortuosity=t_lo + (t_hi - t_lo) * float(u[N_KINETIC]),
        transfer_mix=m_lo + (m_hi - m_lo) * float(u[N_KINETIC + 1]),
    )


def latin_hypercube(n: int, seed: int) -> FloatArray:
    """n points of a Latin hypercube in [0, 1]^8 (scipy.stats.qmc, reproducible by seed)."""
    return np.asarray(qmc.LatinHypercube(d=len(NAMES), seed=seed).random(n))


def sample_problem(
    base: en.ThermalProblem, spec: UncertaintySpec, p: SampleParameters
) -> en.ThermalProblem:
    """The base problem with the parameters of one sample (tables shared)."""
    return base.with_parameters(
        kinetics=p.kinetics,
        tortuosity_factor=p.tortuosity / spec.tortuosity_nominal,
        transfer_mix=p.transfer_mix,
    )


@dataclass(frozen=True)
class RunResult:
    """Outputs of one model run.

    Attributes:
        outputs: Scalar outputs keyed by OUTPUTS (metallisation [-], H2 consumption
            [kg/t DRI], temperatures [degC], GUD [-]).
        Ts, Tg: Solid and gas temperature profiles [degC].
        x_H2: H2 mole fraction profile [-].
        X3: Metallisation profile (wustite -> Fe conversion) [-].
        converged: Whether the solver converged (unconverged runs are kept and counted).
        seconds: Wall time of the run [s].
        retried: Whether the run needed the tighter fallback tolerance (TOL_T_RETRY).
        energy_closure: Relative energy balance closure of the solution [-].
    """

    outputs: dict[str, float]
    Ts: FloatArray
    Tg: FloatArray
    x_H2: FloatArray
    X3: FloatArray
    converged: bool
    seconds: float
    retried: bool = False
    energy_closure: float = float("nan")


#: Solver settings for sampling: looser than the defaults, far below physical effects
#: (docs/validation.md, Milestone 6).
SOLVER = {"tol_T": 1e-4, "species_tol": 1e-8}


#: Fallback temperature tolerance [K]. With ``SOLVER`` about 3 % of the samples stop the
#: outer loop below ``tol_T`` but the final consistency re-solve changes T by slightly more
#: than ``tol_T``, so they are flagged unconverged. These are re-solved from their own end
#: state with the tighter tolerance; the outputs change by ~1e-8 in metallisation.
TOL_T_RETRY = 1e-5


def run_one(prob: en.ThermalProblem) -> RunResult:
    """Solve one sample and collect its outputs (temperatures in degC, see RunResult).

    If the run is flagged unconverged it is re-solved once from its own end state with
    ``tol_T = TOL_T_RETRY``; if still unconverged, the coupled solver finishes from that end
    state (:func:`shaftsim.coupled.finish_coupled`). The result keeps the flag of the last
    solve.
    """
    import time

    t0 = time.perf_counter()
    sol = en.solve_thermal(prob, **SOLVER)  # type: ignore[arg-type]
    retried = not sol.converged
    if retried:
        sol = en.solve_thermal(
            prob,
            tol_T=TOL_T_RETRY,
            species_tol=float(SOLVER["species_tol"]),
            T_start=(sol.Ts, sol.Tg),
        )
    if not sol.converged:  # alternating scheme cycles: coupled finish (Milestone 8)
        from . import coupled as cp

        sol, _ = cp.finish_coupled(prob, sol)
    s = sol.species
    return RunResult(
        outputs={
            "metallization": s.metallization,
            "h2_consumption_kg_per_t": s.h2_consumption,
            "top_gas_T_C": sol.T_top_gas - 273.15,
            "dri_T_C": sol.T_solid_out - 273.15,
            "peak_solid_T_C": float(np.max(sol.Ts)) - 273.15,
            "gud": s.gud,
        },
        Ts=sol.Ts - 273.15,
        Tg=sol.Tg - 273.15,
        x_H2=s.x_H2.copy(),
        X3=s.X_steps[2].copy(),
        converged=sol.converged,
        seconds=time.perf_counter() - t0,
        retried=retried,
        energy_closure=float(sol.energy_closure),
    )


# ----------------------------------------------------------------------------- parallel runs

_WORKER: dict[str, object] = {}


def _init_worker(base: en.ThermalProblem, spec: UncertaintySpec) -> None:
    import warnings

    warnings.simplefilter("ignore")
    _WORKER["base"], _WORKER["spec"] = base, spec


def _run_point(u: FloatArray) -> RunResult:
    base, spec = _WORKER["base"], _WORKER["spec"]
    assert isinstance(base, en.ThermalProblem) and isinstance(spec, UncertaintySpec)
    return run_one(sample_problem(base, spec, to_parameters(spec, u)))


def run_points(
    base: en.ThermalProblem,
    spec: UncertaintySpec,
    points: FloatArray,
    workers: int | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> list[RunResult]:
    """Run the model at each point of the unit hypercube, in parallel processes.

    The base problem (with its property tables, built once here) is sent to every worker.
    Results are returned in the order of ``points``.
    """
    _ = base._tables  # build once, before pickling
    workers = workers or max(1, (os.cpu_count() or 2) // 2)
    out: list[RunResult] = []
    if workers == 1:
        _init_worker(base, spec)
        for k, u in enumerate(points):
            out.append(_run_point(u))
            if progress:
                progress(k + 1, len(points))
        return out
    with ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(base, spec)) as ex:
        for k, r in enumerate(ex.map(_run_point, list(points), chunksize=1)):
            out.append(r)
            if progress:
                progress(k + 1, len(points))
    return out


def output_matrix(results: Sequence[RunResult]) -> FloatArray:
    """(n_runs, len(OUTPUTS)) array of the scalar outputs."""
    return np.array([[r.outputs[k] for k in OUTPUTS] for r in results])


# ----------------------------------------------------------------------------- statistics


def summary(values: FloatArray) -> dict[str, float]:
    """Mean, standard deviation and the 5, 50 and 95 % quantiles of a sample [units of input]."""
    v = np.asarray(values, dtype=float)
    q05, q50, q95 = np.quantile(v, [0.05, 0.5, 0.95])
    return {
        "mean": float(v.mean()),
        "std": float(v.std(ddof=1)) if v.size > 1 else 0.0,
        "p05": float(q05),
        "p50": float(q50),
        "p95": float(q95),
    }


def quantile_bands(profiles: FloatArray, q: tuple[float, ...] = (0.05, 0.5, 0.95)) -> FloatArray:
    """Pointwise quantiles of sampled profiles: shape (len(q), n_cells) from (n_runs, n_cells)."""
    return np.asarray(np.quantile(np.asarray(profiles, dtype=float), q, axis=0))


def running_statistics(values: FloatArray, sizes: Sequence[int]) -> dict[str, FloatArray]:
    """Mean, 5 % and 95 % quantiles of the first n values, for each n in ``sizes``
    (convergence of the Monte Carlo estimates with the number of samples)."""
    v = np.asarray(values, dtype=float)
    out: dict[str, list[float]] = {"n": [], "mean": [], "p05": [], "p95": []}
    for n in sizes:
        s = summary(v[:n])
        out["n"].append(float(n))
        for k in ("mean", "p05", "p95"):
            out[k].append(s[k])
    return {k: np.array(x) for k, x in out.items()}
