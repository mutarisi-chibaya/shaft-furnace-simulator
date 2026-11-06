# Hydrogen Shaft Furnace Simulator – build plan

**Total time:** about 6–7 weeks of evenings
**What makes it impressive:** it is a real reactor model, the kind of furnace modelling Ex Mente does for clients. It couples kinetics, heat transfer and counter-current flow, solves a hard numerical problem, targets the most important decarbonisation technology in steelmaking, and is validated against published work.

## CV entry (once built)

> **Hydrogen Shaft Furnace Simulator** (Python, SciPy, Streamlit)
> - Built a 1D counter-current moving-bed model of hydrogen direct reduction of iron ore, coupling shrinking-core kinetics with gas–solid heat and mass transfer
> - Reproduced a published shaft furnace model within X% on metallisation and gas composition
> - Quantified prediction uncertainty from kinetic parameters with Monte Carlo and Sobol sensitivity analysis, identifying which measurements matter most
> - Built an interactive app for sensitivity studies and CO2-per-tonne estimates under different electricity grids, including South Africa's

Fill in X with your real validation result.

## Folder layout

```
shaft-furnace-simulator/
  pyproject.toml
  README.md  CLAUDE.md  PLAN.md
  data/
    species.yaml          # Shomate coefficients from NIST (you fill in)
    kinetics.yaml         # rate constants, diffusivities (you fill in from literature)
    furnace.yaml          # geometry and operating conditions of the reference case
    SOURCES.md
  src/shaftsim/
    properties.py         # Cp, H, equilibrium constants for the reduction steps
    pellet.py             # single-pellet shrinking-core model
    transport.py          # heat/mass transfer correlations, optional Ergun pressure drop
    shaft.py              # 1D counter-current reactor model
    solver.py             # numerical solution to steady state
    analysis.py           # sensitivity studies, CO2 accounting
  app/streamlit_app.py
  tests/  tests/data/  notebooks/  figures/  docs/validation.md
```

---

## Milestone 1 – Thermodynamics of the reduction steps (week 1)

**Science.** Reduction happens in steps above about 570 °C:

```
3 Fe2O3 + H2 = 2 Fe3O4 + H2O
Fe3O4 + H2 = 3 FeO + H2O
FeO + H2 = Fe + H2O
```

Each step stops when the gas reaches its equilibrium H2O/H2 ratio, K(T). That limit is why real furnaces need excess hydrogen.

**Your manual task:** fill `data/species.yaml` from the NIST WebBook (Fe, FeO, Fe3O4, Fe2O3, H2, H2O; add CO, CO2, CH4 if you want a natural-gas comparison later). Put a few JANAF rows into `tests/data/janaf_reference.csv`.

**Build:** `properties.py` with Cp, H, G per species and K(T) per reduction step.
**Tests:** properties within 0.1% of JANAF; K values match hand calculations.

---

## Milestone 2 – Single-pellet kinetic model (weeks 1–2)

**Science.** In the shrinking-core model, a pellet reduces from the outside in. The rate of each step is limited by three resistances in series: gas film mass transfer, diffusion through the reduced outer layer, and chemical reaction at the interface. The driving force is the distance from equilibrium, for example pH2 − pH2O / K. Three-interface versions model all three steps at once.

**Your manual task:** find a paper that reports isothermal reduction curves (reduction degree vs. time) for hematite pellets in hydrogen, with the rate parameters it used. Put the parameters in `data/kinetics.yaml` and digitise a few curves into `tests/data/pellet_reduction.csv`. (A free tool like WebPlotDigitizer helps.)

**Build:** `pellet.py` returning the reaction rate of each step for given gas composition, temperature and current conversion; an ODE integration of a pellet in constant gas.

**Decision (2026-11-05):** instead of a three-interface shrinking core model, the pellet model is the grainy pellet model with the law of additive reaction times (GM-LART) of Marsigny et al. (2025, *Metals* 15:1034). Reason: it is the only open-access source found that gives isothermal pure-H2 pellet curves over 600–900 °C together with the model and fitted constants, so the implementation can be checked against it; it is also cheap enough for the shaft model, and the same group publishes the shaft-furnace model planned for Milestone 5. Its three resistances (external transfer, pore diffusion, interface reaction) are the same physics as the SCM. See `docs/validation.md`.
**Tests:** reproduces the digitised curves within a tolerance you justify; rate is zero at equilibrium; conversion never exceeds 1.

---

## Milestone 3 – Isothermal 1D shaft (weeks 2–3)

**Science.** Divide the shaft height into cells. Solids enter at the top and move down; gas enters at the bottom and moves up. Each cell has solid conservation for each oxide and gas conservation for H2 and H2O, linked by the pellet rates from Milestone 2. Because the two streams flow in opposite directions, boundary conditions sit at opposite ends. That is the hard numerical part.

**Recommended method:** finite volumes with upwind differencing, marched in pseudo-time until steady state. It is more robust than shooting methods and easier to explain. Ask Claude Code to compare it with `scipy.integrate.solve_bvp` and justify the choice.

**Build:** `shaft.py`, `solver.py`; outputs axial profiles of gas composition and solid conversion; metallisation at the outlet; H2 utilisation.
**Tests:** element balances (Fe, O, H) close to 1e-6; infinitely fast kinetics reproduces the equilibrium limit from Milestone 1; zero kinetics gives no reduction; results converge as the grid is refined.

---

## Milestone 4 – Energy balance and heat transfer (weeks 3–4)

**Science.** Hot gas heats the cold pellets as they descend; the reactions absorb or release heat (hydrogen reduction of iron oxides is overall endothermic, so the gas cools). Add an energy equation for each phase, linked by a gas–solid heat transfer coefficient from a packed-bed correlation, plus reaction heats from Milestone 1. Optionally add Ergun pressure drop.

**Your manual task:** take the correlation(s) and pellet/bed properties from your chosen reference paper into `data/` with citations.

**Build:** non-isothermal shaft model; temperature profiles of gas and solid.
**Tests:** energy balance closes to 1e-6 of the inlet enthalpy flow; with reactions switched off, gas and solid temperatures approach each other as expected for counter-current heat exchange.

---

## Milestone 5 – Validation against a published model (week 5)

This is the step that makes the project credible to an engineering firm.

**Your manual task:** choose one published shaft-furnace model of direct reduction that reports its input conditions and its results (outlet metallisation, gas composition or axial profiles). Searching for counter-current moving-bed direct reduction models (for example the work by Parisi and Laborde, or by Valipour) is a good starting point; check the paper yourself. Put its inputs in `data/furnace.yaml` and digitise its results into `tests/data/reference_model.csv`.

**Build:** a notebook that runs the reference case and overlays your profiles on the published ones; a short discussion of where and why they differ.
**Write:** `docs/validation.md` with the comparison table and your error numbers.

---

## Milestone 6 – Uncertainty quantification (week 6)

This milestone answers a gap that a December 2025 review of hydrogen shaft furnace modelling names directly: reduction kinetics are uncertain, and models rarely say how confident their predictions are. ([Processes 14(1):138](https://www.mdpi.com/2227-9717/14/1/138))

**Your manual task:** from the literature, record a plausible range (or spread) for each uncertain parameter in `data/uncertainty.yaml`: rate constants, activation energies, effective diffusivities, heat-transfer coefficient. Cite where each range comes from.

**Build:**
1. Monte Carlo over the uncertain parameters (Latin hypercube sampling with `scipy.stats.qmc`), reporting the spread of outlet metallisation, H2 consumption and peak temperatures, with confidence bands on the axial profiles.
2. Global sensitivity analysis (Sobol indices, e.g. with the SALib package) showing which parameters matter most. This tells a plant which lab measurement would reduce uncertainty the most.
3. Speed: if one run is slow, profile it and speed it up (vectorisation, sparse Jacobians). Report runs per minute.

**Tests:** with zero uncertainty the Monte Carlo collapses to the single deterministic result; sampled parameters stay inside their ranges; results are reproducible with a fixed random seed.

---

## Milestone 7 – Analysis, app and polish (week 7)

**Analysis:** sensitivity of metallisation and H2 consumption to gas inlet temperature, gas flow and pellet size. CO2 per tonne of DRI as a function of electricity grid intensity (kWh per kg H2 and grid factors as named, cited parameters). Show the South African case honestly.

**App:** Streamlit with sliders for operating conditions, plots of axial temperature, gas composition and conversion profiles, outlet KPIs, and an "Assumptions and limits" panel. Deploy free on Streamlit Community Cloud.

**Polish:** limitations section (1D, steady state, no radial effects, no sticking or fines, kinetics from one source), GitHub Actions (ruff, mypy, pytest), README with the best two figures and the app link.

## Milestone 8 – Numerics and pressure (approved 2026-11-06)

**Solver:** fully coupled Newton (species and energy in one sparse system) with coupled pseudo-transient globalisation; non-smooth points of the rate law smoothed only over a width whose effect on outputs is tested below 1e-6. Must converge where the alternating iteration cycles (case B sweep 625–750 mol/s on 200 cells).

**Scheme:** second-order upwind with a smooth limiter (van Albada), conservative so the balances stay exact; observed order measured against the reaction-free shooting reference.

**Pressure:** Ergun pressure drop (constants cited), local pressure in the transport properties; bustle and cooling pressures compared with REDUCTOR Table 1 (cases A–E), not tuned.

**Energy (conditional):** gas-loop heat recovery and compression from Shao et al. (2021) if openly accessible.

**Then:** rerun the Milestone 5 and 7 comparisons on converged grids; update tests (regression values with reasons), figures, notebooks, docs, app.

## Prompt for each milestone

```
Read CLAUDE.md and PLAN.md. We are on Milestone <N>. Summarise what this milestone requires, list the data files I must provide and check whether they exist, then propose an implementation plan. Wait for my approval before writing code.
```
