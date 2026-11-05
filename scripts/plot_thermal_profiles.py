"""Axial profiles of the non-isothermal shaft model (Milestone 4), case A reduction zone.

Usage:
    python scripts/plot_thermal_profiles.py [n_cells]
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from shaftsim import energy as en
from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr

ROOT = Path(__file__).resolve().parents[1]
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"


def main() -> None:
    warnings.simplefilter("ignore")
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    sp = pr.load_species(ROOT / "data" / "species.yaml")
    td = tr.load_transport(
        ROOT / "data" / "transport.yaml",
        ROOT / "data" / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    trip = pl.load_refit_triplets(ROOT / "data" / "kinetics_refit.yaml")
    case = en.load_thermal_case(ROOT / "data" / "furnace.yaml", ROOT / "data" / "kinetics.yaml")
    sol = en.solve_thermal(en.ThermalProblem(case, sp, td, trip, n))
    ref = en.solve_thermal(en.ThermalProblem(case, sp, td, trip, n, rate_multiplier=0.0))
    if not (sol.converged and ref.converged):
        raise SystemExit("not converged")
    s = sol.species
    H = case.shaft.height
    height = H - s.z
    # cell-outlet values plus the inlet states: burden at the top, gas at the bottom
    hs = np.concatenate([[H], height])
    hg = np.concatenate([height, [0.0]])

    def solid(T: np.ndarray) -> np.ndarray:
        return np.concatenate([[case.T_solid_in], T]) - 273.15

    def gas(T: np.ndarray) -> np.ndarray:
        return np.concatenate([T, [case.T_gas_in]]) - 273.15

    fig, axes = plt.subplots(1, 3, figsize=(13, 6), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(color=GRID, lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK_2)
    ax = axes[0]
    ax.plot(gas(sol.Tg), hg, color=ORANGE, lw=2, label="gas")
    ax.plot(solid(sol.Ts), hs, color=BLUE, lw=2, label="solid")
    ax.plot(
        gas(ref.Tg),
        hg,
        color=ORANGE,
        lw=1,
        ls=(0, (4, 3)),
        alpha=0.7,
        label="gas, no reaction",
    )
    ax.plot(
        solid(ref.Ts),
        hs,
        color=BLUE,
        lw=1,
        ls=(0, (4, 3)),
        alpha=0.7,
        label="solid, no reaction",
    )
    ax.axvline(893.0 - 273.15, color=MUTED, lw=0.8)
    ax.text(893.0 - 273.15 + 5, 5.0, "wüstite limit\n(JANAF data)", color=INK_2, fontsize=7.5)
    ax.set_xlabel("Temperature (°C)", color=INK_2)
    ax.set_ylabel("Height above gas inlet (m)", color=INK_2)
    ax.set_title(
        f"Temperatures (top gas {sol.T_top_gas - 273.15:.0f} °C)",
        color=INK,
        loc="left",
        fontsize=10,
    )
    ax.legend(frameon=False, fontsize=7.5, loc="center left", labelcolor=INK)
    ax = axes[1]
    for j, (lab, col) in enumerate(
        zip(("Fe₂O₃ → Fe₃O₄", "Fe₃O₄ → FeₓO", "FeₓO → Fe"), (BLUE, ORANGE, AQUA), strict=True)
    ):
        ax.plot(np.concatenate([[0.0], s.X_steps[j]]), hs, color=col, lw=2, label=lab)
    ax.plot(np.concatenate([[0.0], s.X]), hs, color=INK, lw=1.2, ls=(0, (4, 3)), label="overall X")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Conversion", color=INK_2)
    ax.set_title(f"Solid (metallisation {s.metallization:.3f})", color=INK, loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="center", labelcolor=INK)
    ax = axes[2]
    ax.plot(s.x_H2, height, color=BLUE, lw=2, label="H₂")
    ax.plot(s.x_H2O, height, color=ORANGE, lw=2, label="H₂O")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Mole fraction (5 % N₂ not shown)", color=INK_2)
    ax.set_title(f"Gas (GUD {s.gud:.3f})", color=INK, loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="center right", labelcolor=INK)
    fig.text(
        0.01,
        0.005,
        f"Non-isothermal 1D model, case A reduction zone (Marsigny et al. 2025, "
        f"Metals 15:862): burden in at 20 °C, 2800 mol/s gas in at 890 °C, {n} cells. "
        "Not a validation (Milestone 5).",
        color=INK_2,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = ROOT / "figures" / "shaft_thermal_profiles.png"
    fig.savefig(out, facecolor=SURFACE)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
