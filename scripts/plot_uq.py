"""Figures and tables of Milestone 6 from results/uq/*.npz (made by scripts/run_uq.py).

Usage:
    python scripts/plot_uq.py

Writes figures/uq_profiles.png, figures/uq_outputs.png, figures/uq_convergence.png and
figures/uq_sobol.png, and prints the tables used in docs/validation.md. The table functions
are also used by notebooks/06_uncertainty.ipynb.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from shaftsim import uq

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "uq"
FIG = ROOT / "figures"
sys.path.insert(0, str(ROOT / "scripts"))
import run_uq  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
LABELS = {
    "metallization": "Metallisation (-)",
    "h2_consumption_kg_per_t": "H₂ consumption (kg/t DRI)",
    "top_gas_T_C": "Top gas temperature (°C)",
    "dri_T_C": "DRI temperature (°C)",
    "peak_solid_T_C": "Peak solid temperature (°C)",
    "gud": "Gas utilisation degree (-)",
}
#: REDUCTOR case A (tests/data/reference_model_marsigny2025b_caseA.csv), for orientation
REDUCTOR = {"metallization": 0.933, "top_gas_T_C": 344.0, "dri_T_C": 261.0, "gud": 0.33}
SHOW = ("metallization", "h2_consumption_kg_per_t", "dri_T_C", "top_gas_T_C")


def style(ax: plt.Axes) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2)


def load(name: str) -> dict[str, np.ndarray] | None:
    """Arrays of results/uq/<name>, or None if the file does not exist."""
    p = RES / name
    if not p.exists():
        return None
    with np.load(p) as d:
        return {k: d[k] for k in d.files}


def nominal(cells: int) -> uq.RunResult:
    """Deterministic run with the nominal parameters (same solver settings as the samples)."""
    base, _ = run_uq.setup(cells)
    return uq.run_one(base)


def output_names(d: dict[str, np.ndarray]) -> list[str]:
    return [str(s) for s in d["output_names"]]


def mc_table(mc: dict[str, np.ndarray], nom: uq.RunResult) -> pd.DataFrame:
    """Nominal value and sample statistics of every output (converged runs only)."""
    ok = mc["converged"]
    rows = []
    for k, nm in enumerate(output_names(mc)):
        s = uq.summary(mc["outputs"][ok, k])
        rows.append({"output": nm, "nominal": nom.outputs[nm], **s})
    return pd.DataFrame(rows).set_index("output")


def sobol_table(so: dict[str, np.ndarray]) -> pd.DataFrame:
    """Grouped Sobol indices: one row per (output, group) with S1, ST and 95 % CIs."""
    meta = json.loads(str(so["meta"]))
    rows = []
    for nm, r in meta["indices"].items():
        for i, g in enumerate(meta["groups"]):
            rows.append(
                {
                    "output": nm,
                    "group": g,
                    "S1": r["S1"][i],
                    "S1_conf": r["S1_conf"][i],
                    "ST": r["ST"][i],
                    "ST_conf": r["ST_conf"][i],
                }
            )
    return pd.DataFrame(rows).set_index(["output", "group"])


def grid_table(mc: dict[str, np.ndarray], fine: dict[str, np.ndarray]) -> pd.DataFrame:
    """Same samples on 50 cells and on a finer grid: mean offset, spreads, correlation."""
    n = fine["outputs"].shape[0]
    both = mc["converged"][:n] & fine["converged"]
    rows = []
    for k, nm in enumerate(output_names(mc)):
        a, b = mc["outputs"][:n][both, k], fine["outputs"][both, k]
        rows.append(
            {
                "output": nm,
                "mean offset": float(np.mean(b - a)),
                "std coarse": float(np.std(a, ddof=1)),
                "std fine": float(np.std(b, ddof=1)),
                "std ratio": float(np.std(b, ddof=1) / np.std(a, ddof=1)),
                "corr": float(np.corrcoef(a, b)[0, 1]),
            }
        )
    return pd.DataFrame(rows).set_index("output")


def plot_profiles(mc: dict[str, np.ndarray], nom: uq.RunResult) -> plt.Figure:
    ok = mc["converged"]
    z = 15.0 - mc["z_depth"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 10.5), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    panels = (
        ("Ts", "Solid temperature (°C)", "Solid temperature − nominal (K)", ORANGE),
        ("X3", "Metallisation (wüstite → Fe conversion)", "Metallisation − nominal (-)", BLUE),
        ("x_H2", "H₂ mole fraction", "H₂ mole fraction − nominal (-)", AQUA),
    )
    for j, (key, lab, dlab, col) in enumerate(panels):
        b = uq.quantile_bands(mc[key][ok])
        ref = getattr(nom, key)
        for ax, shift, xl in ((axes[0, j], 0.0, lab), (axes[1, j], ref, dlab)):
            style(ax)
            ax.fill_betweenx(
                z,
                b[0] - shift,
                b[2] - shift,
                color=col,
                alpha=0.25,
                lw=0,
                label="5–95 % of samples",
            )
            ax.plot(b[1] - shift, z, color=col, lw=1.6, label="median")
            ax.plot(ref - shift, z, color=INK, lw=1, ls=(0, (4, 3)), label="nominal")
            for zl in (4.9, 0.0):
                ax.axhline(zl, color=MUTED, lw=0.7, ls=(0, (2, 3)))
            ax.set_xlabel(xl, color=INK_2)
    for row in range(2):
        ax = axes[row, 0]
        ax.set_ylabel("Height z (m)", color=INK_2)
        ax.text(0.02, 5.1, "bustle", color=MUTED, fontsize=7.5, transform=ax.get_yaxis_transform())
        ax.text(0.02, 0.2, "cone", color=MUTED, fontsize=7.5, transform=ax.get_yaxis_transform())
    axes[0, 0].legend(frameon=False, fontsize=7.5, loc="upper right", labelcolor=INK)
    fig.text(
        0.01,
        0.005,
        f"Whole furnace, case A, {ok.sum()} Latin-hypercube samples of kinetics, tortuosity and "
        "transfer correlation (data/uncertainty.yaml); 50 cells. "
        "Bottom row: the same bands minus the nominal run.",
        color=INK_2,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    return fig


def plot_outputs(mc: dict[str, np.ndarray], nom: uq.RunResult) -> plt.Figure:
    ok = mc["converged"]
    names = output_names(mc)
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    for ax, nm in zip(axes, SHOW, strict=True):
        style(ax)
        ax.hist(mc["outputs"][ok, names.index(nm)], bins=30, color=BLUE, alpha=0.75)
        ax.axvline(nom.outputs[nm], color=INK, lw=1.2, ls=(0, (4, 3)), label="nominal")
        if nm in REDUCTOR:
            ax.axvline(REDUCTOR[nm], color=ORANGE, lw=1.6, label="REDUCTOR")
        ax.set_xlabel(LABELS[nm], color=INK_2)
    axes[0].set_ylabel("Samples", color=INK_2)
    axes[0].legend(frameon=False, fontsize=7.5, labelcolor=INK)
    fig.tight_layout()
    return fig


def plot_convergence(mc: dict[str, np.ndarray]) -> plt.Figure:
    ok = mc["converged"]
    names = output_names(mc)
    sizes = [n for n in (16, 32, 64, 128, 256, 512, 1024) if n <= ok.sum()]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    for ax, nm in zip(axes, ("metallization", "dri_T_C"), strict=True):
        style(ax)
        r = uq.running_statistics(mc["outputs"][ok, names.index(nm)], sizes)
        ax.plot(r["n"], r["mean"], "o-", color=BLUE, label="mean")
        ax.plot(r["n"], r["p05"], "s--", color=MUTED, label="5 % / 95 %")
        ax.plot(r["n"], r["p95"], "s--", color=MUTED)
        ax.set_xscale("log", base=2)
        ax.set_xlabel("Number of samples (first n of the Latin hypercube)", color=INK_2)
        ax.set_ylabel(LABELS[nm], color=INK_2)
    axes[0].legend(frameon=False, fontsize=7.5, labelcolor=INK)
    fig.tight_layout()
    return fig


def plot_sobol(so: dict[str, np.ndarray]) -> plt.Figure:
    meta = json.loads(str(so["meta"]))
    groups, ind = meta["groups"], meta["indices"]
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.8), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    x = np.arange(len(groups))
    for ax, nm in zip(axes, SHOW, strict=True):
        style(ax)
        r = ind[nm]
        ax.bar(x - 0.18, r["S1"], 0.36, yerr=r["S1_conf"], color=BLUE, label="first order S1")
        ax.bar(x + 0.18, r["ST"], 0.36, yerr=r["ST_conf"], color=ORANGE, label="total ST")
        ax.set_xticks(x, [g.replace(" correlation", "\ncorrelation") for g in groups], fontsize=8)
        ax.set_title(LABELS[nm], color=INK, fontsize=9, loc="left")
        ax.set_ylim(-0.05, 1.1)
    axes[0].set_ylabel("Sobol index (error bars: 95 % bootstrap)", color=INK_2)
    axes[0].legend(frameon=False, fontsize=7.5, labelcolor=INK)
    fig.tight_layout()
    return fig


def main() -> None:
    mc = load("mc_50cells.npz")
    if mc is None:
        raise SystemExit("results/uq/mc_50cells.npz missing; run scripts/run_uq.py mc")
    nom = nominal(50)
    pd.set_option("display.width", 160)
    ok = mc["converged"]
    print(
        f"MC: {ok.sum()}/{ok.size} converged ({int(mc.get('retried', np.zeros(1)).sum())} with "
        f"the fallback tolerance), {float(mc['runs_per_minute']):.1f} runs/min "
        f"({int(mc['workers'])} workers), median {np.median(mc['seconds']):.1f} s/run"
    )
    print(mc_table(mc, nom).to_string(float_format=lambda v: f"{v:.4g}"))
    for name, fig in (
        ("uq_profiles", plot_profiles(mc, nom)),
        ("uq_outputs", plot_outputs(mc, nom)),
        ("uq_convergence", plot_convergence(mc)),
    ):
        fig.savefig(FIG / f"{name}.png", facecolor=SURFACE)
        plt.close(fig)
    so = load("sobol_50cells.npz")
    if so is not None:
        print(f"\nSobol: {so['converged'].sum()}/{so['converged'].size} converged")
        print(sobol_table(so).to_string(float_format=lambda v: f"{v:.2f}"))
        fig = plot_sobol(so)
        fig.savefig(FIG / "uq_sobol.png", facecolor=SURFACE)
        plt.close(fig)
    for cells in (200, 400):
        fine = load(f"grid_{cells}cells.npz")
        if fine is not None:
            print(f"\nGrid check: first {fine['outputs'].shape[0]} MC samples, {cells} vs 50 cells")
            print(grid_table(mc, fine).to_string(float_format=lambda v: f"{v:.4g}"))
    print("wrote figures/uq_*.png")


if __name__ == "__main__":
    main()
