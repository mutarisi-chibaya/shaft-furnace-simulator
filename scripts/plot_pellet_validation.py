"""Plot single-pellet model vs Marsigny et al. (2025) TGA experiments.

Left: temperature series (Fig. 6 data, used to fit the kinetic triplets). Right: H2-He
dilution series at 900 degC (Fig. 9 data, NOT used in the fit). Markers: digitised
experiments; lines: shaftsim.pellet with data/kinetics_refit.yaml.

Usage:
    python scripts/plot_pellet_validation.py
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from shaftsim import pellet as pl
from shaftsim import properties as pr
from shaftsim import transport as tr

ROOT = Path(__file__).resolve().parents[1]
# Sequential blue ramp (reference palette), light -> dark for low -> high values.
RAMP7 = ["#6da7ec", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#0d366b"]
RAMP4 = ["#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def main() -> None:
    sp = pr.load_species(ROOT / "data" / "species.yaml")
    td = tr.load_transport(
        ROOT / "data" / "transport.yaml",
        ROOT / "data" / "collision_integrals.csv",
        {k: v.molar_mass for k, v in sp.items() if v.molar_mass is not None},
    )
    pp, cond = pl.load_kinetics(ROOT / "data" / "kinetics.yaml")
    trip = pl.load_refit_triplets(ROOT / "data" / "kinetics_refit.yaml")
    f6 = pd.read_csv(ROOT / "tests" / "data" / "pellet_reduction_marsigny2025_fig6.csv")
    f9 = pd.read_csv(ROOT / "tests" / "data" / "pellet_reduction_marsigny2025_fig9.csv")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), dpi=150, sharey=True)
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
    for colour, (T_c, g) in zip(RAMP7, f6.groupby("temperature_C"), strict=True):
        g = g.iloc[:: max(len(g) // 40, 1)]
        ax.plot(g.time_min, g.conversion, "o", ms=3, color=colour, alpha=0.6, mec="none")
        r = pl.simulate_tga(pp, cond, float(str(T_c)) + 273.15, 150 * 60.0, sp, td, triplets=trip)
        ax.plot(r.t / 60, r.X, "-", lw=1.8, color=colour, label=f"{T_c} °C")
    ax.set_xlim(0, 150)
    ax.set_title("Temperature series (used for the fit)", color=INK, loc="left", fontsize=10.5)
    ax.set_xlabel("Time (min)", color=INK_2)
    ax.set_ylabel("Conversion X", color=INK_2)
    ax.legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK)

    ax = axes[1]
    for colour, (x_h2, g) in zip(RAMP4, f9.groupby("x_H2_in"), strict=True):
        g = g.iloc[:: max(len(g) // 40, 1)]
        ax.plot(g.time_min, g.conversion, "o", ms=3, color=colour, alpha=0.6, mec="none")
        c = dataclasses.replace(cond, x_H2=float(str(x_h2)))
        r = pl.simulate_tga(pp, c, 1173.15, 40 * 60.0, sp, td, triplets=trip)
        ax.plot(
            r.t / 60, r.X, "-", lw=1.8, color=colour, label=f"{100 * float(str(x_h2)):.0f} % H₂"
        )
    ax.set_xlim(0, 40)
    ax.set_title("H₂–He dilution at 900 °C (held out)", color=INK, loc="left", fontsize=10.5)
    ax.set_xlabel("Time (min)", color=INK_2)
    ax.legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK)
    axes[0].set_ylim(0, 1.02)
    fig.text(
        0.01,
        0.005,
        "Dots: experiments (Marsigny et al. 2025, Figs. 6 and 9, digitised). "
        "Lines: shaftsim pellet model, re-fitted kinetic constants.",
        color=INK_2,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = ROOT / "figures" / "pellet_validation.png"
    fig.savefig(out, facecolor=SURFACE)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
