"""Axial profiles of the isothermal shaft model (Milestone 3) for case A at two temperatures.

Usage:
    python scripts/plot_shaft_profiles.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import shaft as sh
from shaftsim import solver as sv
from shaftsim import transport as tr

ROOT = Path(__file__).resolve().parents[1]
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def main() -> None:
    warnings.simplefilter("ignore")
    sp = pr.load_species(ROOT / "data" / "species.yaml")
    td = tr.load_transport(
        ROOT / "data" / "transport.yaml",
        ROOT / "data" / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    trip = pl.load_refit_triplets(ROOT / "data" / "kinetics_refit.yaml")
    temps = (1163.15, 973.15)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.5), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for row, T in enumerate(temps):
        case = sh.load_case(ROOT / "data" / "furnace.yaml", ROOT / "data" / "kinetics.yaml", T=T)
        prob = sh.ShaftProblem(case, sh.shaft_coefficients(case, sp, td, trip), 400)
        y, st = sv.solve(prob)
        if not st.converged:
            raise SystemExit(f"not converged at T = {T}")
        s = sh.make_solution(prob, y, td, st.residuals[-1], st.iterations)
        height = case.height - s.z  # height above the gas inlet
        ax1, ax2 = axes[row]
        for ax in (ax1, ax2):
            ax.set_facecolor(SURFACE)
            ax.grid(color=GRID, lw=0.6)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color(GRID)
            ax.tick_params(colors=INK_2)
        for j, (lab, col) in enumerate(
            zip(
                ("Fe₂O₃ → Fe₃O₄", "Fe₃O₄ → FeₓO", "FeₓO → Fe (= metallisation)"),
                (BLUE, ORANGE, AQUA),
                strict=True,
            )
        ):
            ax1.plot(s.X_steps[j], height, color=col, lw=2, label=lab)
        ax1.plot(s.X, height, color=INK, lw=1.2, ls=(0, (4, 3)), label="overall reduction X")
        ax2.plot(s.x_H2, height, color=BLUE, lw=2, label="H₂")
        ax2.plot(s.x_H2O, height, color=ORANGE, lw=2, label="H₂O")
        ax1.set_xlim(-0.02, 1.02)
        ax2.set_xlim(0, 1)
        ax1.set_ylabel("Height above gas inlet (m)", color=INK_2)
        ax1.set_title(
            f"{T - 273.15:.0f} °C: solid conversion (metallisation {s.metallization:.3f})",
            color=INK,
            loc="left",
            fontsize=10,
        )
        ax2.set_title(
            f"{T - 273.15:.0f} °C: gas mole fraction (GUD {s.gud:.3f})",
            color=INK,
            loc="left",
            fontsize=10,
        )
        if row == 0:
            ax1.legend(frameon=False, fontsize=8, loc="lower left", labelcolor=INK)
            ax2.legend(frameon=False, fontsize=8, loc="center right", labelcolor=INK)
    axes[1][0].set_xlabel("Conversion", color=INK_2)
    axes[1][1].set_xlabel("Mole fraction (5 % N₂ not shown)", color=INK_2)
    fig.text(
        0.01,
        0.005,
        "Isothermal 1D model of the reduction zone of case A (Marsigny et al. "
        "2025, Metals 15:862): 10.1 m, 2800 mol/s, 45 kg/s, 400 cells. Not a validation; "
        "the real furnace is not isothermal.",
        color=INK_2,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = ROOT / "figures" / "shaft_isothermal_profiles.png"
    fig.savefig(out, facecolor=SURFACE)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
