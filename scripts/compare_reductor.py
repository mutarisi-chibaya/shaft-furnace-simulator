"""Whole-furnace model vs REDUCTOR case A (Milestone 5): tables and overlay figure.

Usage:
    python scripts/compare_reductor.py [n_cells]

Prints the Table 1 comparison and the profile agreement and writes
figures/reductor_comparison.png. Also used by notebooks/05_reductor_comparison.ipynb.
"""

from __future__ import annotations

import sys
import warnings
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr
from shaftsim import validation as va

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
REF_OUT = ROOT / "tests" / "data" / "reference_model_marsigny2025b_caseA.csv"
REF_PROF = ROOT / "tests" / "data" / "reference_model_marsigny2025b_caseA_profiles.csv"
FIG = ROOT / "figures" / "reductor_comparison.png"

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"

#: Profiles compared (quantity, axis label, scale to display units)
PANELS = (
    ("T_solid", "Solid temperature (°C)", 1.0),
    ("x_H2", "H₂ mole fraction", 1.0),
    ("w_Fe3O4", "Fe₃O₄ mass fraction", 1.0),
    ("w_FeO", "Wüstite mass fraction", 1.0),
    ("w_Fe", "Fe mass fraction", 1.0),
)
#: Regions for the profile scores [m, REDUCTOR z]
REGIONS = {
    "whole furnace": (-9.8, 14.9),
    "above bustle": (4.9, 14.9),
    "below bustle": (-9.8, 4.9),
}


def run(n_cells: int) -> tuple[en.ThermalProblem, en.ThermalSolution]:
    """Solve the whole-furnace case A on n_cells; raises if not converged."""
    warnings.simplefilter("ignore")
    sp = pr.load_species(DATA / "species.yaml")
    td = tr.load_transport(
        DATA / "transport.yaml",
        DATA / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    trip = pl.load_refit_triplets(DATA / "kinetics_refit.yaml")
    case = en.load_furnace_case(DATA / "furnace.yaml", DATA / "kinetics.yaml")
    prob = en.ThermalProblem(case, sp, td, trip, n_cells)
    sol = en.solve_thermal(prob)
    if not sol.converged:
        raise RuntimeError(f"whole-furnace model did not converge on {n_cells} cells")
    return prob, sol


def output_table(prob: en.ThermalProblem, sol: en.ThermalSolution) -> pd.DataFrame:
    """Table 1 quantities: reference, model, difference, relative difference."""
    ref = va.load_reference_outputs(REF_OUT)
    rows = va.compare_outputs(ref, va.model_outputs(prob, sol))
    return pd.DataFrame(
        {
            "key": [c.key for c in rows],
            "REDUCTOR": [c.reference for c in rows],
            "model": [c.model for c in rows],
            "difference": [c.difference for c in rows],
            "relative": [c.relative for c in rows],
        }
    )


def profile_table(prob: en.ThermalProblem, sol: en.ThermalSolution) -> pd.DataFrame:
    """Band agreement of each profile, per region."""
    ref = va.load_reference_profiles(REF_PROF)
    mod = va.model_profiles(prob, sol)
    rows = []
    for q in ("T_solid", "x_H2", "w_Fe2O3", "w_Fe3O4", "w_FeO", "w_Fe"):
        for region, zr in REGIONS.items():
            c = va.compare_profile(ref, mod, q, zr)
            rows.append({"region": region, **asdict(c)})
    return pd.DataFrame(rows)


def plot(prob: en.ThermalProblem, sol: en.ThermalSolution, out: Path = FIG) -> Path:
    """Overlay of the model profiles on the REDUCTOR cross-section averages."""
    ref = va.load_reference_profiles(REF_PROF)
    mod = va.model_profiles(prob, sol)
    t1 = va.load_reference_outputs(REF_OUT)
    fig, axes = plt.subplots(1, len(PANELS), figsize=(15, 6.5), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (q, label, scale) in zip(axes, PANELS, strict=True):
        ax.set_facecolor(SURFACE)
        ax.grid(color=GRID, lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK_2)
        r = ref[ref.quantity == q]
        ax.fill_betweenx(
            r.z,
            r.mean_lo * scale,
            r.mean_hi * scale,
            color=ORANGE,
            alpha=0.18,
            lw=0,
            label="REDUCTOR (colour-band range)",
        )
        ax.plot(r["mean"] * scale, r.z, color=ORANGE, lw=1.6, label="REDUCTOR (band midpoints)")
        ax.plot(mod[q] * scale, mod.z, color=BLUE, lw=2, label="1D model")
        for zl in (4.9, 0.0):
            ax.axhline(zl, color=MUTED, lw=0.7, ls=(0, (4, 3)))
        ax.set_xlabel(label, color=INK_2)
    # Table 1 outlet values at the solid outlet
    axes[0].plot(t1["dri_T"], -10.0, "D", color=ORANGE, ms=6, label="REDUCTOR Table 1 (DRI)")
    axes[2].plot(t1["dri_w_Fe3O4"] / 100, -10.0, "D", color=ORANGE, ms=6)
    axes[3].plot(t1["dri_w_FeO"] / 100, -10.0, "D", color=ORANGE, ms=6)
    axes[4].plot(t1["dri_w_Fe"] / 100, -10.0, "D", color=ORANGE, ms=6)
    axes[0].text(30, 5.1, "bustle", color=MUTED, fontsize=7.5)
    axes[0].text(30, 0.2, "cone", color=MUTED, fontsize=7.5)
    axes[0].set_ylabel("Height z (m)", color=INK_2)
    axes[0].set_ylim(-10.4, 15.2)
    axes[0].legend(frameon=False, fontsize=7.5, loc="upper right", labelcolor=INK)
    fig.text(
        0.01,
        0.005,
        f"Case A of Marsigny et al. (2025, Metals 15:862). REDUCTOR: 2D axisymmetric, "
        f"cross-section averages digitised from Figs. 3, 4, 7-9 (band resolution). "
        f"Model: 1D, {prob.n_cells} cells, no parameter tuned.",
        color=INK_2,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return out


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    prob, sol = run(n)
    pd.set_option("display.width", 160)
    print(f"n = {n}, outer iterations {sol.outer_iterations}, closure {sol.energy_closure:.1e}")
    print(output_table(prob, sol).to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print(profile_table(prob, sol).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    out = plot(prob, sol)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
