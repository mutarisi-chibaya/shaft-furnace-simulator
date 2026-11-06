"""Figures and tables of Milestone 7 from results/cases and results/sensitivity.

Usage:
    python scripts/plot_m7.py

Writes figures/m7_cooling_sweep.png, figures/m7_sensitivity.png and figures/m7_co2.png and
prints the energy/CO2 numbers used in docs/validation.md. Case A energy demand is computed
here on 200 cells (~3 min).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from shaftsim import emissions as em
from shaftsim import operating as op
from shaftsim import scenario as sc

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
sys.path.insert(0, str(ROOT / "scripts"))
import compare_cases as cc  # noqa: E402
import run_sensitivity as rs  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"


def style(ax: Axes) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2)


def cooling_sweep(cells: int) -> list[dict[str, Any]] | None:
    p = cc.OUT / f"cooling_continuation_{cells}cells.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def plot_cooling_sweep() -> None:
    ref = cc.reference_fig10()
    fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ax.plot(list(ref), list(ref.values()), "s", color=INK, ms=7, label="REDUCTOR (2D), Fig. 10")
    r50 = cc.load_results(50)
    if r50 is not None:
        x = list(cc.SWEEP)
        y = [r50[f"B_cool{f:g}"]["metallization"] for f in cc.SWEEP]
        ax.plot(x, y, "o--", color=MUTED, ms=5, lw=1, label="1D model, 50 cells")
    for cells, col in ((200, BLUE),):
        res = cooling_sweep(cells)
        if res is None:
            continue
        # both directions give the same values where they converge (checked here); one line,
        # broken where no converged solution exists
        conv: dict[float, float] = {}
        for r in res:
            if r["converged"]:
                f, m = r["cooling_flow"], 100 * r["metallization"]
                if f in conv and abs(conv[f] - m) > 1e-3:
                    raise SystemExit(f"continuation branches differ at {f} mol/s")
                conv[f] = m
        flows = sorted({r["cooling_flow"] for r in res})
        y = [conv.get(f, np.nan) for f in flows]
        ax.plot(
            flows,
            y,
            "o-",
            color=col,
            ms=5,
            lw=1.6,
            label=f"1D model, {cells} cells (continuation up and down: identical)",
        )
        bad = sorted(
            {r["cooling_flow"] for r in res if not r["converged"]}
            - {r["cooling_flow"] for r in res if r["converged"]}
        )
        if bad:
            ax.axvspan(min(bad) - 12.5, max(bad) + 12.5, color=GRID, alpha=0.6, lw=0)
            ax.text(
                0.5 * (min(bad) + max(bad)),
                ax.get_ylim()[0] + 1,
                "no converged\nsolution",
                ha="center",
                va="bottom",
                color=INK_2,
                fontsize=8,
            )
    ax.set_xlabel("Cooling gas flow (mol/s); bustle flow = 2800 − cooling flow", color=INK_2)
    ax.set_ylabel("Metallisation (%)", color=INK_2)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "m7_cooling_sweep.png", facecolor=SURFACE)
    plt.close(fig)


def plot_sensitivity() -> None:
    # (axis label, offset, factor): display value = (SI value + offset) * factor
    labels: dict[str, tuple[str, float, float]] = {
        "bustle_temperature": ("Bustle gas temperature (°C)", -273.15, 1.0),
        "bustle_flow": ("Bustle gas flow (mol/s)", 0.0, 1.0),
        "pellet_diameter": ("Pellet diameter (mm)", 0.0, 1000.0),
    }
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.6), dpi=150, sharey="row")
    fig.patch.set_facecolor(SURFACE)
    for cells, col, mk in ((50, MUTED, "o--"), (200, BLUE, "o-")):
        res = rs.load(cells)
        if res is None:
            continue
        for j, (key, (lab, offset, factor)) in enumerate(labels.items()):
            rows = sorted((r for r in res if r["input"] == key), key=lambda r: r["value"])
            x = [(r["value"] + offset) * factor for r in rows]
            for i, (out, scale) in enumerate(
                (("metallization", 100.0), ("co2_south_africa_t_per_t", 1.0))
            ):
                y = [scale * r[out] if r["converged"] else np.nan for r in rows]
                axes[i, j].plot(x, y, mk, color=col, ms=4, lw=1.4, label=f"{cells} cells")
            axes[1, j].set_xlabel(lab, color=INK_2)
    for ax in axes.ravel():
        style(ax)
    axes[0, 0].set_ylabel("Metallisation (%)", color=INK_2)
    axes[1, 0].set_ylabel("t CO₂ / t DRI, South African grid", color=INK_2)
    axes[0, 0].legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.text(
        0.01,
        0.005,
        "One input varied at a time over its cited range (data/sensitivity.yaml); "
        "the others at case A. Gaps: no converged solution.",
        color=INK_2,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(FIG / "m7_sensitivity.png", facecolor=SURFACE)
    plt.close(fig)


def plot_co2(base: sc.BaseModel, demand: em.EnergyDemand) -> None:
    d = base.emissions
    fig, ax = plt.subplots(figsize=(7.5, 4.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    g = np.linspace(0, 1000, 101)
    ax.plot(
        g, demand.total * g / 3.6e9, color=BLUE, lw=2, label="H₂ shaft furnace, case A (this model)"
    )
    for gi in d.grids.values():
        y = em.co2_per_dri(demand, gi.value)
        ax.plot(gi.g_per_kWh, y, "o", color=BLUE, ms=7, mec="white", mew=1.5)
        label = f"{gi.name} ({gi.year}): {y:.2f}"
        if gi.g_per_kWh < 200:  # near the origin: label in the free area below the line
            ax.annotate(
                label,
                (gi.g_per_kWh, y),
                xytext=(200, 0.15),
                textcoords="data",
                ha="left",
                va="center",
                color=INK_2,
                fontsize=8.5,
                arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8},
            )
        else:  # the line rises to the right: upper left of the point is clear
            ax.annotate(
                label,
                (gi.g_per_kWh, y),
                xytext=(-8, 8),
                textcoords="offset points",
                ha="right",
                color=INK_2,
                fontsize=8.5,
            )
    ng = d.benchmarks["ng_dri_direct_55_45"]
    ax.axhline(ng, color=MUTED, ls=(0, (5, 3)), lw=1.2)
    ax.text(
        1000,
        ng + 0.04,
        "Natural-gas DRI (H₂:CO 55:45), direct emissions only",
        ha="right",
        color=MUTED,
        fontsize=8,
    )
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, None)
    ax.set_xlabel("Grid CO₂ intensity (g CO₂e/kWh, life-cycle; Ember)", color=INK_2)
    ax.set_ylabel("t CO₂ per t DRI", color=INK_2)
    fig.text(
        0.01,
        0.005,
        f"Electricity {demand.total / 1e6:.1f} GJ/t DRI: electrolysis "
        f"{demand.electrolysis / 1e6:.1f} (18 MJ/Nm³ H₂) + gas heating {demand.heating / 1e6:.1f} "
        f"(η 0.7, no heat recovery) + compression {demand.compression / 1e6:.1f}; "
        "pelletising and EAF excluded.",
        color=INK_2,
        fontsize=7,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIG / "m7_co2.png", facecolor=SURFACE)
    plt.close(fig)


def main() -> None:
    base = sc.load_base(ROOT / "data")
    res = sc.run_scenario(base, op.conditions_of(base.case), 200)
    if not res.converged:
        raise SystemExit("case A did not converge on 200 cells")
    e = res.energy
    ng = base.emissions.benchmarks["ng_dri_direct_55_45"]
    print(
        f"case A, 200 cells: H2 {res.h2_consumption:.2f} kg/t, "
        f"electrolysis {e.electrolysis / 1e6:.2f}, heating {e.heating / 1e6:.2f}, "
        f"compression {e.compression / 1e6:.2f}, total {e.total / 1e6:.2f} GJ/t DRI"
    )
    for g in base.emissions.grids.values():
        print(
            f"  {g.name} {g.g_per_kWh:.0f} g/kWh ({g.year}): "
            f"{em.co2_per_dri(e, g.value):.2f} t CO2/t DRI"
        )
    print(f"  break-even with natural-gas DRI direct ({ng} t/t): {ng / e.total * 3.6e9:.0f} g/kWh")
    plot_cooling_sweep()
    plot_sensitivity()
    plot_co2(base, e)
    print("wrote figures/m7_*.png")


if __name__ == "__main__":
    main()
