"""Streamlit app: the 1D hydrogen shaft furnace (case A geometry) at user-chosen conditions.

Run locally:
    streamlit run app/streamlit_app.py

All physics is in shaftsim (scenario.run_scenario); this file only collects inputs, caches
results and draws them. Slider ranges are those of data/sensitivity.yaml and of the REDUCTOR
cases (data/furnace_cases.yaml); the model is not offered outside them.
"""

from __future__ import annotations

import base64
import dataclasses
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import yaml

from shaftsim import emissions as em
from shaftsim import operating as op
from shaftsim import scenario as sc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
T0 = 273.15

# Categorical slots 1-3 of the reference palette (validated all-pairs, light and dark).
PALETTE = {
    "light": {
        "s1": "#2a78d6",
        "s2": "#eb6834",
        "s3": "#1baf7a",
        "ink2": "#52514e",
        "muted": "#8a8984",
    },
    "dark": {
        "s1": "#3987e5",
        "s2": "#d95926",
        "s3": "#199e70",
        "ink2": "#c3c2b7",
        "muted": "#8a8984",
    },
}
#: Presets offered (case D, a 1200 degC burden, lies outside the kinetic data and is left out).
PRESETS = ("A", "B", "C", "E")
PRESET_TEXT = {
    "A": "A: reference (Table 1)",
    "B": "B: 200 mol/s moved from bustle to cooling inlet",
    "C": "C: 32 % N2, bustle gas at 950 °C",
    "E": "E: burden flow 43 kg/s (full metallisation in REDUCTOR)",
}


#: Guided scenarios: (button text, preset, slider overrides, grid key, what to look for).
#: All values lie inside the slider ranges; the notes describe directions, not numbers.
TRY_THIS = {
    "cooler": (
        "Cooler hydrogen",
        "A",
        {"bustle_T": 840.0},
        "south_africa",
        "The hydrogen now enters 50 °C cooler. Watch **metallisation** fall: the last "
        "reduction step (wüstite to iron) needs heat as well as hydrogen.",
    ),
    "cold_gas": (
        "More cold gas at the bottom",
        "B",
        {},
        "south_africa",
        "200 mol/s of gas moves from the hot inlet to the cold one (case B of the published "
        "model). In this 1D model the extra cold gas chills the cone and **metallisation** "
        "falls, while the published 2D model barely changes: a known limit of 1D (see "
        "*Assumptions and limits*).",
    ),
    "clean": (
        "Clean electricity",
        "A",
        {},
        "france",
        "Same furnace, French electricity (mostly nuclear and hydro). **CO₂ per tonne** falls "
        "about seventeen-fold: the footprint of hydrogen iron is set by the electricity.",
    ),
    "reset": ("Back to the reference case", "A", {}, "south_africa", ""),
}

GLOSSARY = """
- **Bustle**: the ring of gas inlets around the middle of the shaft, where the hot hydrogen
  enters.
- **Burden**: the iron ore pellets moving down through the furnace.
- **Pellet**: a small ball of iron ore (here 10–15 mm across).
- **Reduction**: removing the oxygen from the iron ore. Here hydrogen does it and forms water.
- **Hematite, magnetite, wüstite** (Fe₂O₃, Fe₃O₄, FeO): the iron oxides the ore passes
  through, one step at a time, on its way to iron (Fe).
- **DRI** (direct reduced iron): the furnace's product, solid iron that is later melted into
  steel in an electric arc furnace.
- **Metallisation**: the share of the iron in the product that is metal rather than oxide.
- **Top gas**: the gas leaving the top: unused hydrogen and water vapour. The water is
  condensed out and the hydrogen is recycled.
- **Cooling gas**: cold gas fed at the bottom of the cone to cool the iron before it leaves.
- **Gas utilisation**: the share of the top gas (hydrogen plus water vapour) that is water
  vapour, i.e. how much of the hydrogen did its job on one pass.
- **Grid CO₂ intensity**: CO₂ emitted per kilowatt-hour of electricity, counted over the life
  cycle of the power plants.
- **REDUCTOR**: a published two-dimensional model of the same furnace (Marsigny, Mirgaux and
  Patisson 2025), used here as the reference.
- **1D, grid cells**: the model treats each height of the furnace as uniform across its width
  and splits the height into cells. More cells are more accurate but slower.
"""


def apply_try(name: str) -> None:
    """Button callback: set the controls for a guided scenario (runs before the rerun)."""
    _, preset, changes, grid, _ = TRY_THIS[name]
    ss = st.session_state
    ss["preset"] = preset
    ss["_preset_loaded"] = preset
    for k, v in {**preset_values(preset), **changes}.items():
        ss[k] = v
    ss["grid"] = grid
    ss["_tip"] = name
    ss["_tip_state"] = None  # filled in on the next run, cleared when a control changes


def try_this() -> None:
    """Row of guided-scenario buttons."""
    with st.container(key="try"):
        st.markdown("**Try this** (one click sets the controls on the left):")
        cols = st.columns(len(TRY_THIS))
        for col, (name, (text, *_)) in zip(cols, TRY_THIS.items(), strict=True):
            col.button(text, key=f"try_{name}", on_click=apply_try, args=(name,))


# ----------------------------------------------------------------------------- model


@st.cache_resource
def base_model() -> sc.BaseModel:
    return sc.load_base(DATA)


@st.cache_resource
def ranges() -> dict[str, tuple[float, float]]:
    raw = yaml.safe_load((DATA / "sensitivity.yaml").read_text(encoding="utf-8"))
    return {k: (float(v["low"]["value"]), float(v["high"]["value"])) for k, v in raw.items()}


@st.cache_data(show_spinner=False, max_entries=64)
def solve(inputs: tuple[float, ...], p: float, n_cells: int) -> sc.ScenarioResult:
    """Cached run; inputs = (bustle T degC, bustle flow, cooling flow, N2 %, solid flow,
    pellet mm)."""
    t_b, f_b, f_c, n2, solid, dp = inputs
    base = base_model()
    oc = dataclasses.replace(
        op.conditions_of(base.case),
        bustle_T=t_b + T0,
        bustle_flow=f_b,
        cooling_flow=f_c,
        x_gas={"H2": round(0.98 - n2 / 100, 10), "H2O": 0.02, "N2": n2 / 100},
        solid_flow=solid,
        pellet_diameter=dp / 1000,
        p=p,
    )
    return sc.run_scenario(base, oc, n_cells)


def preset_values(name: str) -> dict[str, float]:
    oc = base_model().reference_cases[name]
    return {
        "bustle_T": round(oc.bustle_T - T0, 1),
        "bustle_flow": oc.bustle_flow,
        "cooling_flow": oc.cooling_flow,
        "n2": round(100 * oc.x_gas["N2"], 1),
        "solid_flow": oc.solid_flow,
        "pellet": round(oc.pellet_diameter * 1000, 1),
    }


# ----------------------------------------------------------------------------- charts


def colors() -> dict[str, str]:
    theme = getattr(getattr(st.context, "theme", None), "type", None)
    return PALETTE["dark" if theme == "dark" else "light"]


def label_positions(series: list[np.ndarray]) -> list[int]:
    """For each series, the index (middle 70 % of the height) where it is farthest from the
    others: where a direct label collides least."""
    n = series[0].size
    lo, hi = int(0.15 * n), int(0.85 * n)
    span = max(float(np.ptp(np.concatenate(series))), 1e-12)
    min_sep = max(1, int(0.12 * n))  # labels at least 12 % of the height apart
    idx = np.arange(n)
    out: list[int] = []
    for i, s in enumerate(series):
        others = [o for j, o in enumerate(series) if j != i]
        if not others:
            out.append((lo + hi) // 2)
            continue
        gap = np.min([np.abs(s - o) for o in others], axis=0) / span
        if gap[lo:hi].max() < 0.05:  # curves nearly coincide: spread the labels in height
            out.append(lo + (hi - lo) * (2 * i + 1) // (2 * len(series)))
            continue
        ok = (idx >= lo) & (idx < hi)
        for k in out:
            ok &= np.abs(idx - k) >= min_sep
        score = np.where(ok, gap, -np.inf)
        out.append(int(np.argmax(score)) if np.isfinite(score.max()) else (lo + hi) // 2)
    return out


def profile_chart(
    z: np.ndarray,
    series: list[tuple[str, np.ndarray, str]],
    xlabel: str,
    fmt: str,
    xrange: tuple[float, float] | None = None,
) -> go.Figure:
    c = colors()
    allx = np.concatenate([x for _, x, _ in series])
    mid = 0.5 * (float(allx.min()) + float(allx.max()))
    fig = go.Figure()
    for name, x, slot in series:
        fig.add_trace(
            go.Scatter(
                x=x,
                y=z,
                name=name,
                mode="lines",
                line={"color": c[slot], "width": 2},
                hovertemplate=f"{name}: %{{x:{fmt}}}<extra></extra>",
            )
        )
    for (name, x, _), k in zip(series, label_positions([s for _, s, _ in series]), strict=True):
        right = x[k] > mid  # keep the label inside the plot
        fig.add_annotation(
            x=x[k],
            y=z[k],
            text=name,
            showarrow=False,
            xanchor="right" if right else "left",
            xshift=-6 if right else 6,
            yshift=9,
            font={"color": c["ink2"], "size": 12},
        )
    for zl, txt in ((4.9, "bustle"), (0.0, "cone")):
        fig.add_hline(
            y=zl,
            line={"color": c["muted"], "width": 1, "dash": "dot"},
            annotation_text=txt,
            annotation_position="bottom left",
            annotation_font={"color": c["muted"], "size": 11},
        )
    fig.update_layout(
        height=520,
        margin={"l": 10, "r": 10, "t": 40, "b": 10},
        hovermode="y unified",
        legend={"orientation": "h", "y": 1.08, "x": 0},
        xaxis_title=xlabel,
        yaxis_title="Height z (m)",
    )
    if xrange is not None:
        fig.update_xaxes(range=list(xrange))
    return fig


def co2_chart(demand: em.EnergyDemand, data: em.EmissionsData, custom: float | None) -> go.Figure:
    c = colors()
    gi = np.linspace(0.0, 1000.0, 101)  # g CO2e/kWh
    y = demand.total * gi / 3.6e9
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=gi,
            y=y,
            mode="lines",
            line={"color": c["s1"], "width": 2},
            name="This scenario",
            hovertemplate="%{x:.0f} g/kWh: %{y:.2f} t CO₂/t DRI<extra></extra>",
        )
    )
    points = [(g.name, g.g_per_kWh, g.year) for g in data.grids.values()]
    if custom is not None:
        points.append(("Your value", custom, None))
    for name, g, year in points:
        yy = demand.total * g / 3.6e9
        fig.add_trace(
            go.Scatter(
                x=[g],
                y=[yy],
                mode="markers",
                showlegend=False,
                marker={"size": 10, "color": c["s1"], "line": {"width": 2, "color": "white"}},
                hovertemplate=f"{name}: %{{y:.2f}} t CO₂/t DRI<extra></extra>",
            )
        )
        # the line rises to the right: upper-left and lower-right of a point are clear
        left_edge = g < 200.0
        fig.add_annotation(
            x=g,
            y=yy,
            text=f"{name}{f' ({year})' if year else ''}: {yy:.2f}",
            showarrow=False,
            xanchor="left" if left_edge else "right",
            xshift=10 if left_edge else -10,
            yshift=-10 if left_edge else 10,
            font={"color": c["ink2"], "size": 12},
        )
    ng = data.benchmarks["ng_dri_direct_55_45"]
    fig.add_hline(
        y=ng,
        line={"color": c["muted"], "width": 1.5, "dash": "dash"},
        annotation_text="Natural-gas DRI, direct emissions only",
        annotation_position="top right",
        annotation_font={"color": c["muted"], "size": 11},
    )
    fig.update_layout(
        height=380,
        margin={"l": 10, "r": 10, "t": 20, "b": 10},
        showlegend=False,
        xaxis_title="Grid CO₂ intensity (g CO₂e/kWh, life-cycle)",
        yaxis_title="t CO₂ per t DRI",
        yaxis_rangemode="tozero",
    )
    return fig


# ----------------------------------------------------------------------------- page


ASSUMPTIONS = """
**What the model is.** A steady-state, one-dimensional (plug-flow) model of a hydrogen
direct-reduction shaft furnace: three-step reduction (hematite → magnetite → wüstite → iron)
in grainy pellets with kinetics refitted to published thermogravimetry, gas–solid heat and
mass transfer (Wakao–Kaguei), and solid and gas energy balances with NIST/JANAF
thermodynamics. Mass and energy balances close to round-off (shown below each run).

**Geometry and fixed inputs.** Case A of Marsigny, Mirgaux & Patisson (2025): 15 m cylinder
(radius 2.6 m) above a 10 m cone, bustle 4.9 m above the cone, gas at 2 % H₂O. Pellet
structure and kinetics are those of the Vale pellets of the kinetic study.

**Limits you should know before reading the numbers.**
- **Grid.** The app solves on 50 or 100 cells for speed. That is not grid-converged: the
  converged case A (grid study to 3200 cells) has a metallisation about 0.04 lower (0.03 at
  100 cells), a DRI about 40 K hotter and a top gas about 3 K cooler. Differences between
  scenarios are usually less affected than the absolute values, except near the cooling-gas
  flow where the cooling front moves into the cone (see Solver below).
- **1D.** No radial gradients, so no H₂-starved core and no wall effects. Against the 2D model
  REDUCTOR, the outlet values agree within a few percent for case A, but the reduction front
  sits 1.5–2 m lower and the cooling zone behaves differently, so the agreement is partly
  compensation.
- **Not modelled:** the local pressure (one mean pressure from the REDUCTOR pressures of the
  chosen preset is used; an Ergun pressure field exists in the code but changes case A by
  less than 0.001), heat losses, sticking and fines, burden descent other than plug flow, carbon.
- **Kinetics from one source** (seven TGA curves of one pellet type, 600–900 °C). Pellet size
  is varied with the same kinetics; bustle gas above 900 °C is a mild extrapolation.
- **Solver.** Where the cooling front sits inside the cone (e.g. cooling gas around
  625–750 mol/s) the usual iteration cycles and a fully coupled solver takes over, which
  takes longer. If neither converges, the app shows an error instead of numbers.
- **Energy and CO₂.** Electricity for electrolysis (18 MJ/Nm³ H₂), bustle gas heating
  (heater efficiency 0.7, no heat recovery, so an upper bound), as in the authors' gas-loop
  model, and gas compression (isothermal, efficiency 0.7, Shao et al. 2021). Pellet making,
  the EAF and plant construction are not included.
  Grid intensities are Ember's life-cycle values. The natural-gas line is a direct-emission
  figure per tonne of DRI and excludes upstream methane and electricity.

Sources: `data/SOURCES.md`. Validation and limitations in detail: `docs/validation.md`.
"""


def overview() -> None:
    """Plain-language summary for visitors without a metallurgy background (read in 30 s).

    The findings restate results documented in docs/validation.md (Milestones 5, 7 and 8).
    """
    with st.container(border=True, key="overview"):
        pic, text = st.columns([5, 4])
        svg = (ROOT / "figures" / "furnace_schematic.svg").read_bytes()
        pic.markdown(
            '<img alt="Schematic of the hydrogen shaft furnace" style="width:100%" '
            f'src="data:image/svg+xml;base64,{base64.b64encode(svg).decode()}">',
            unsafe_allow_html=True,
        )
        a = b = c = text
        a.markdown(
            "**What this is**\n\n"
            "A computer model of a furnace that turns iron ore into iron using **hydrogen "
            "instead of coal**, the key step towards low-carbon steel. Inside the furnace "
            "the hydrogen takes the oxygen out of the ore and leaves as water vapour."
        )
        b.markdown(
            "**What it answers**\n\n"
            "- How much of the ore becomes usable iron\n"
            "- How much hydrogen and electricity that takes\n"
            "- How much CO₂ that electricity causes, in South Africa and elsewhere"
        )
        c.markdown(
            "**Main findings**\n\n"
            "1. Without any tuning, its prediction of how much of the iron becomes metal is "
            "within 3 percentage points of a published research model.\n"
            "2. On South Africa's coal-heavy grid, iron made this way causes about **5× the "
            "CO₂** of iron made with natural gas. It only pays off with clean electricity.\n"
            "3. How hot the hydrogen is, and how much is fed in, matter most; the size of "
            "the ore pellets hardly matters."
        )
        st.caption(
            "Change the settings on the left, or use the buttons below; the results update. "
            "A step-by-step guide is in docs/tutorial.md; the words used are explained in the "
            "glossary at the bottom of the page."
        )


def plain_summary(res: sc.ScenarioResult, co2: float, grid_label: str, ng: float) -> str:
    """One sentence on the headline results in everyday words."""
    ratio = co2 / ng if ng > 0 else float("nan")
    if co2 >= ng:
        compare = f"{ratio:.1f}× the direct emissions of natural-gas iron ({ng:.2f} t)"
    else:
        compare = f"less than the direct emissions of natural-gas iron ({ng:.2f} t)"
    return (
        f"**In plain words:** {100 * res.metallization:.1f} % of the iron leaves the furnace as "
        f"metal (the published reference model: {REDUCTOR_A_METALLISATION:.1f} % for its "
        f"reference case). With {grid_label}, each tonne of iron causes **{co2:.1f} t of CO₂**, "
        f"{compare}."
    )


#: REDUCTOR metallisation of case A [%] (tests/data, Marsigny et al. 2025b, Table 1).
REDUCTOR_A_METALLISATION = float(
    pd.read_csv(ROOT / "tests" / "data" / "reference_model_marsigny2025b_cases.csv", comment="#")
    .set_index(["case", "key"])
    .loc[("A", "metallization"), "value"]
)


def main() -> None:
    st.set_page_config(page_title="H₂ shaft furnace", layout="wide")
    st.title("Making iron with hydrogen: a furnace simulator")
    st.caption(
        "A physics-based model of a hydrogen direct-reduction shaft furnace. Every parameter is "
        "cited in the repository."
    )
    overview()
    try_this()
    base = base_model()
    rng = ranges()

    with st.sidebar:
        st.header("Operating conditions")
        st.caption("What goes into the furnace. Hover over ⓘ to see what each control means.")
        preset = st.selectbox(
            "Start from REDUCTOR case",
            PRESETS,
            format_func=lambda k: PRESET_TEXT[k],
            key="preset",
            help="Four operating cases of a published 2D model of this furnace (REDUCTOR, "
            "Marsigny et al. 2025b). Picking one sets all the controls below.",
        )
        if st.session_state.get("_preset_loaded") != preset:
            for k, v in preset_values(preset).items():
                st.session_state[k] = v
            st.session_state["_preset_loaded"] = preset
        lo, hi = rng["bustle_temperature"]
        st.slider(
            "Bustle gas temperature (°C)",
            lo - T0,
            hi - T0,
            step=5.0,
            format="%.0f",
            key="bustle_T",
            help="How hot the hydrogen is when it enters through the ring of inlets (the "
            "bustle) around the middle of the furnace. Range of the published studies.",
        )
        lo, hi = rng["bustle_flow"]
        st.slider(
            "Bustle gas flow (mol/s)",
            lo,
            hi,
            step=10.0,
            format="%.0f",
            key="bustle_flow",
            help="How much hot gas (mostly hydrogen) enters at the bustle, in moles per second. "
            "Range of the published studies.",
        )
        st.slider(
            "Cooling gas flow (mol/s)",
            500.0,
            800.0,
            step=10.0,
            format="%.0f",
            key="cooling_flow",
            help="Cold gas fed at the bottom of the cone; it cools the iron before it leaves. "
            "Range of Fig. 10 of Marsigny et al. (2025b).",
        )
        st.slider(
            "N₂ in the feed gas (%)",
            5.0,
            32.0,
            step=1.0,
            format="%.0f",
            key="n2",
            help="Nitrogen dilutes the hydrogen and does not react. The rest of the gas is "
            "hydrogen and 2 % water vapour. Range of cases A–C.",
        )
        st.slider(
            "Burden flow (kg/s)",
            43.0,
            45.0,
            step=0.5,
            format="%.1f",
            key="solid_flow",
            help="How fast iron ore pellets (the 'burden') are fed at the top, in kilograms "
            "per second. Range of cases A–E.",
        )
        lo, hi = rng["pellet_diameter"]
        st.slider(
            "Pellet diameter (mm)",
            lo * 1000,
            hi * 1000,
            step=0.5,
            format="%.1f",
            key="pellet",
            help="Size of the iron ore pellets. Range of the sensitivity study.",
        )
        st.header("Electricity")
        grid_keys = list(base.emissions.grids)
        grid_choice = st.radio(
            "Grid CO₂ intensity",
            [*grid_keys, "custom"],
            key="grid",
            help="CO₂ emitted per kWh of electricity, over the life cycle of the power plants "
            "(Ember, 2025). The electricity makes the hydrogen and heats the gas.",
            format_func=lambda k: (
                f"{base.emissions.grids[k].name} ({base.emissions.grids[k].g_per_kWh:.0f} g/kWh, "
                f"{base.emissions.grids[k].year})"
                if k != "custom"
                else "Custom"
            ),
        )
        custom = (
            st.slider("Custom intensity (g CO₂e/kWh)", 0.0, 1000.0, 300.0, step=10.0, format="%.0f")
            if grid_choice == "custom"
            else None
        )
        st.header("Numerics")
        n_cells = st.radio(
            "Grid cells",
            (50, 100),
            horizontal=True,
            help="The furnace height is split into cells: more cells are more accurate but "
            "slower. 50 cells: ~15 s per run; 100 cells: ~40 s. Neither is grid-converged.",
        )

    ss = st.session_state
    inputs = (
        float(ss.bustle_T),
        float(ss.bustle_flow),
        float(ss.cooling_flow),
        float(ss.n2),
        float(ss.solid_flow),
        float(ss.pellet),
    )
    p = base.reference_cases[preset].p
    with st.spinner("Solving the furnace (10–40 s; results are cached)…"):
        res = solve(inputs, p, int(n_cells))
    if not res.converged:
        st.error(
            "The solver did not converge for these conditions; the values below are not a solution."
        )

    grid = base.emissions.grids.get(grid_choice) if grid_choice != "custom" else None
    intensity = grid.value if grid is not None else float(custom or 0.0) / 3.6e9
    co2 = em.co2_per_dri(res.energy, intensity)
    ng = base.emissions.benchmarks["ng_dri_direct_55_45"]
    st.subheader("Results for the settings on the left")
    tip = ss.get("_tip")
    state = (inputs, grid_choice)
    if tip and ss.get("_tip_state") is None:
        ss["_tip_state"] = state  # the controls as the button set them
    if tip and ss.get("_tip_state") == state and TRY_THIS[tip][4]:
        st.info("**What to look for:** " + TRY_THIS[tip][4])
    elif tip and ss.get("_tip_state") != state:
        ss["_tip"] = None  # a control was changed by hand: the note no longer applies
    kpis = st.container(key="kpis")
    cols = kpis.columns(6)
    cols[0].metric(
        "Metallisation (%)",
        f"{100 * res.metallization:.1f}",
        help="Share of the iron that leaves the furnace as metal; the rest is still oxide. "
        "Higher is better.",
    )
    cols[1].metric(
        "H₂ (kg/t DRI)",
        f"{res.h2_consumption:.1f}",
        help="Hydrogen used per tonne of product. DRI = direct reduced iron, the furnace's "
        "product, which is melted into steel afterwards.",
    )
    cols[2].metric(
        "Top gas (°C)",
        f"{res.top_gas_T - T0:.0f}",
        help="Temperature of the gas leaving the top of the furnace (it is recycled).",
    )
    cols[3].metric(
        "DRI outlet (°C)",
        f"{res.dri_T - T0:.0f}",
        help="Temperature of the iron product leaving the bottom of the furnace.",
    )
    cols[4].metric(
        "Electricity (GJ/t DRI)",
        f"{res.energy.total / 1e6:.2f}",
        help="Electricity per tonne of iron: making the hydrogen by electrolysis of water, "
        "heating the gas and compressing it. 1 GJ = 278 kWh.",
    )
    cols[5].metric(
        "CO₂ (t/t DRI)",
        f"{co2:.2f}",
        help="CO₂ caused by that electricity on the grid chosen on the left (life-cycle "
        f"values). For comparison, iron made with natural gas emits about {ng:.2f} t directly.",
    )
    if grid is None:
        grid_label = "the custom electricity grid"
    elif grid.name.lower() == "world":
        grid_label = "the world-average electricity grid"
    else:
        grid_label = f"{grid.name}'s electricity grid"
    st.markdown(plain_summary(res, co2, grid_label, ng))
    st.caption(
        f"{n_cells} cells, {res.seconds:.0f} s. Gas utilisation {res.gud:.3f}; DRI "
        f"{res.dri_flow * 3.6:.1f} t/h. Balance closure: energy {abs(res.energy_closure):.0e}, "
        f"elements {res.element_closure:.0e}. Not grid-converged: see Assumptions and limits."
    )

    if inputs == tuple(float(v) for v in preset_values(preset).values()):
        ref = pd.read_csv(
            ROOT / "tests" / "data" / "reference_model_marsigny2025b_cases.csv", comment="#"
        )
        r = ref[ref["case"] == preset].set_index("key")["value"]
        st.container(key="reductor").info(
            f"REDUCTOR (2D) for case {preset}: metallisation {r['metallization']:.1f} %, "
            f"top gas {r['top_gas_T']:.0f} °C, DRI {r['dri_T']:.0f} °C, GUD {r['gud']:.2f} "
            "(Marsigny et al. 2025b, Table 1)."
        )

    st.subheader("CO₂ per tonne of DRI and the electricity grid")
    left, right = st.container(key="co2").columns([3, 2])
    left.plotly_chart(co2_chart(res.energy, base.emissions, custom), use_container_width=True)
    breakeven = ng / res.energy.total * 3.6e9
    right.markdown(
        f"""
Electricity per tonne of DRI: **{res.energy.electrolysis / 1e6:.2f} GJ** for electrolysis,
**{res.energy.heating / 1e6:.2f} GJ** for gas heating and
**{res.energy.compression / 1e6:.2f} GJ** for gas compression.

- On South Africa's grid ({base.emissions.grids["south_africa"].g_per_kWh:.0f} g/kWh in
  {base.emissions.grids["south_africa"].year}) hydrogen DRI emits
  **{em.co2_per_dri(res.energy, base.emissions.grids["south_africa"].value):.2f} t CO₂/t DRI**,
  about {em.co2_per_dri(res.energy, base.emissions.grids["south_africa"].value) / ng:.0f}× the
  direct emissions of natural-gas DRI ({ng:.2f} t/t).
- Hydrogen DRI beats natural-gas DRI only below about **{breakeven:.0f} g CO₂e/kWh**, i.e. on
  a largely low-carbon supply (dedicated renewables or a clean grid), not on today's average grid.
- For context: blast furnace steel emits about
  {base.emissions.benchmarks["bf_bof_crude_steel"]:.1f} t CO₂ per tonne of *crude steel*
  (IEA 2020), a different basis that includes the steel plant.
"""
    )

    st.subheader("Inside the furnace (for specialists)")
    st.caption(
        "How temperature, gas composition and the three reduction steps change from the top "
        "of the furnace (height 15 m) down to the bottom of the cone (−10 m). The hot gas enters "
        "at the bustle; cold gas enters at the bottom."
    )
    z = res.height
    c1, c2, c3 = st.container(key="profiles").columns(3)
    c1.plotly_chart(
        profile_chart(
            z, [("Solid", res.Ts - T0, "s2"), ("Gas", res.Tg - T0, "s1")], "Temperature (°C)", ".0f"
        ),
        use_container_width=True,
    )
    c2.plotly_chart(
        profile_chart(
            z, [("H₂", res.x_H2, "s1"), ("H₂O", res.x_H2O, "s2")], "Mole fraction (-)", ".3f"
        ),
        use_container_width=True,
    )
    c3.plotly_chart(
        profile_chart(
            z,
            [
                ("Fe₂O₃ → Fe₃O₄", res.X_steps[0], "s1"),
                ("Fe₃O₄ → FeO", res.X_steps[1], "s2"),
                ("FeO → Fe", res.X_steps[2], "s3"),
            ],
            "Conversion of each step (-)",
            ".3f",
            xrange=(-0.02, 1.02),
        ),
        use_container_width=True,
    )
    st.caption(
        "Solid and gas temperatures nearly coincide: gas–solid heat transfer is fast compared "
        "with the flows. Conversion: fraction of each reduction step completed."
    )

    with st.expander("Assumptions and limits", expanded=False):
        st.markdown(ASSUMPTIONS)
    with st.expander("Glossary: the words used here"):
        st.markdown(GLOSSARY)
    with st.expander("Profile data (table)"):
        st.dataframe(
            pd.DataFrame(
                {
                    "height_m": z,
                    "T_solid_C": res.Ts - T0,
                    "T_gas_C": res.Tg - T0,
                    "x_H2": res.x_H2,
                    "x_H2O": res.x_H2O,
                    "X_hematite": res.X_steps[0],
                    "X_magnetite": res.X_steps[1],
                    "X_wustite": res.X_steps[2],
                }
            ).round(4),
            use_container_width=True,
        )


main()
