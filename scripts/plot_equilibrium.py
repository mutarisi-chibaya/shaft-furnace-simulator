"""Plot the Fe-O-H equilibrium (Baur-Glaessner) diagram from data/species.yaml.

Writes figures/equilibrium_h2.png: equilibrium H2 fraction in H2-H2O gas vs temperature for
each reduction step. Above a line the gas can reduce the phase below it.

Usage:
    python scripts/plot_equilibrium.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from shaftsim import properties as pr

ROOT = Path(__file__).resolve().parents[1]

# Reference categorical palette (light mode), fixed order; text in ink tokens.
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"


def main() -> None:
    sp = pr.load_species(ROOT / "data" / "species.yaml")
    dis = pr.wustite_disproportionation(sp)
    t_c = brentq(lambda T: float(pr.reaction_gibbs(dis, sp, T)), 600.0, 1100.0)

    T_lo, T_hi = 500.0, 1500.0  # H2O data start at 500 K; wüstite fit ends at 1600 K
    below = np.linspace(T_lo, t_c, 200)
    above = np.linspace(t_c, T_hi, 300)
    full = np.linspace(T_lo, T_hi, 400)

    def y(rxn: pr.Reaction, T: np.ndarray) -> np.ndarray:
        return np.asarray(pr.equilibrium_h2_fraction(rxn, sp, T))

    def c(T: np.ndarray) -> np.ndarray:
        return T - 273.15

    fig, ax = plt.subplots(figsize=(8, 6.0), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.plot(c(above), y(pr.WUSTITE_TO_IRON, above), color=ORANGE, lw=2, label="Fe$_{0.947}$O / Fe")
    ax.plot(
        c(above),
        y(pr.MAGNETITE_TO_WUSTITE, above),
        color=BLUE,
        lw=2,
        label="Fe$_3$O$_4$ / Fe$_{0.947}$O",
    )
    ax.plot(
        c(below),
        y(pr.MAGNETITE_TO_IRON, below),
        color=AQUA,
        lw=2,
        label="Fe$_3$O$_4$ / Fe (below wüstite limit)",
    )
    ax.plot(
        c(full),
        y(pr.HEMATITE_TO_MAGNETITE, full),
        color=YELLOW,
        lw=2,
        label="Fe$_2$O$_3$ / Fe$_3$O$_4$ (y$_{H2}$ < 0.001)",
    )
    ax.plot(
        c(full),
        y(pr.FEO_TO_IRON, full),
        color=MUTED,
        lw=1.5,
        ls=(0, (4, 3)),
        label="stoichiometric FeO / Fe (not used)",
    )
    ax.axvline((t_c - 273.15), color=GRID, lw=1, zorder=0)
    ax.annotate(
        f"wüstite stable above {t_c - 273.15:.0f} °C",
        ((t_c - 273.15) + 8, 0.03),
        color=INK_2,
        fontsize=8.5,
    )

    # Phase-field labels
    for txt, xy in [
        ("Fe", (1100, 0.85)),
        ("Fe$_{0.947}$O", (1000, 0.42)),
        ("Fe$_3$O$_4$", (520, 0.45)),
        ("Fe$_3$O$_4$", (1150, 0.08)),
    ]:
        ax.text(*xy, txt, color=INK, fontsize=11, ha="center")

    ax.set_xlim((T_lo - 273.15), (T_hi - 273.15))
    ax.set_ylim(0, 1)
    ax.set_xlabel("Temperature (°C)", color=INK_2)
    ax.set_ylabel("Equilibrium H$_2$ fraction, H$_2$/(H$_2$+H$_2$O)", color=INK_2)
    ax.set_title(
        "Fe–O–H equilibrium at 1 bar (NIST/JANAF data)", color=INK, loc="left", fontsize=12
    )
    ax.grid(color=GRID, lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2)
    # Kinks from the NIST Fe Shomate range boundaries (700 K, 1042 K); see docs/validation.md
    for tk in (700.0, 1042.0):
        yk = float(np.asarray(pr.equilibrium_h2_fraction(pr.FEO_TO_IRON, sp, tk)))
        ax.annotate("*", (tk - 273.15, yk + 0.015), color=INK_2, fontsize=10, ha="center")
    ax.text(
        0.0,
        -0.2,
        "* small steps at 427 and 769 °C come from discontinuities in the NIST "
        "Fe Shomate fit (range boundaries 700 and 1042 K)",
        transform=ax.transAxes,
        color=INK_2,
        fontsize=7.5,
    )
    leg = ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2, frameon=False, fontsize=8.5
    )
    for t in leg.get_texts():
        t.set_color(INK)

    out = ROOT / "figures" / "equilibrium_h2.png"
    fig.tight_layout()
    fig.savefig(out, facecolor=SURFACE)
    print(f"wrote {out.relative_to(ROOT)} (wüstite limit {t_c:.1f} K)")


if __name__ == "__main__":
    main()
