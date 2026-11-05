# Hydrogen Shaft Furnace Simulator

A physics-based simulator of a hydrogen direct-reduction shaft furnace, the technology at the centre of green steelmaking. Iron ore pellets move down the shaft while hot hydrogen flows up and removes the oxygen, producing direct reduced iron (DRI) with water as the main by-product instead of CO2.

The simulator predicts what happens inside the furnace, along its full height: temperatures, gas composition and how far the ore has been reduced at each point.

## Why this matters

- **Industry needs fast, validated models for scale-up.** The EU-funded MaxH2DR project (2022–2026), with Tata Steel Nederland as a partner, aims to develop validated digital twins of industrial hydrogen direct-reduction plants, because pilot trials alone cannot support scale-up. ([Łukasiewicz-GIT](https://git.lukasiewicz.gov.pl/en/?p=53997))
- **The modelling gaps are documented.** A December 2025 review identifies limited predictive accuracy, computational cost, uncertain reduction kinetics and scarce plant validation data as open problems, and names fast surrogate models and uncertainty quantification as future directions. ([Processes 14(1):138](https://www.mdpi.com/2227-9717/14/1/138))
- **It is happening in Southern Africa.** ArcelorMittal South Africa is studying green-hydrogen DRI at its Saldanha works using the Midrex process ([Mining Weekly](https://www.miningweekly.com/article/study-under-way-into-producing-green-direct-reduced-iron-at-mothballed-saldanha-works-2023-10-11)), and HyIron in Namibia is building a hydrogen DRI plant ([SteelOrbis](https://www.steelorbis.com/steel-news/latest-news/hyiron-produces-first-green-hydrogen-in-southern-africa-1385155.htm)).

This simulator is a fast screening model: it runs thousands of operating scenarios in minutes, with quantified uncertainty, before committing to expensive CFD or pilot trials.

## What it does

- Models the stepwise reduction of hematite → magnetite → wüstite → iron in each pellet using a shrinking-core kinetic model
- Solves a 1D steady-state counter-current moving-bed model with coupled mass and energy balances
- Predicts metallisation, hydrogen utilisation and specific hydrogen consumption (kg H2 per tonne DRI)
- Runs sensitivity studies (gas temperature, gas flow, pellet size) and estimates CO2 per tonne for different electricity grids, including South Africa's
- Quantifies how uncertainty in the kinetic parameters spreads into the predictions (Monte Carlo, global sensitivity analysis)
- Interactive Streamlit app showing axial profiles as you change operating conditions
- Validated against a published shaft furnace model and single-pellet experimental data

## Status

Under development. See `PLAN.md`.

- **Milestone 1 – thermodynamics:** done. `src/shaftsim/properties.py`, data generated from NIST/JANAF by `scripts/fetch_thermo_data.py`.
- **Milestone 2 – single-pellet kinetics:** done. `src/shaftsim/pellet.py` (grainy pellet model, law of additive reaction times) and `src/shaftsim/transport.py`. Kinetic constants re-fitted to the TGA experiments of Marsigny et al. (2025) and checked on held-out data.
- **Milestone 3 – isothermal 1D shaft:** done. `src/shaftsim/shaft.py` (counter-current finite volumes, exact element balances) and `src/shaftsim/solver.py` (discrete shooting, coupled Gauss–Seidel, pseudo-transient Newton). Reference geometry and flows from case A of Marsigny et al. (2025).
- **Milestone 4 – energy balance and heat transfer:** done. `src/shaftsim/energy.py` (solid and gas energy equations with absolute enthalpies, Wakao–Kaguei heat transfer, energy closure to round-off; Anderson-accelerated coupling with the species solver). Case A reduction zone: top gas ~345 °C, metallisation ~0.70 (not yet compared with the reference model; see Milestone 5).
- **Milestone 5 – comparison with a published model:** done. Whole furnace (cylinder, cooling cone, bustle side inlet) compared with REDUCTOR case A (Marsigny et al. 2025): outlet metallisation ~90–92 % vs 93.3 %, top gas ~360 vs 344 °C, DRI ~300–310 vs 261 °C, with no parameter tuned. The axial path differs (radial effects a 1D model cannot capture), so the outlet agreement is not a validation of the profiles. `notebooks/05_reductor_comparison.ipynb`, `figures/reductor_comparison.png`.
- **Milestone 6 – uncertainty quantification:** done. `src/shaftsim/uq.py` with `scripts/run_uq.py` (parallel, ~25 runs/min on 6 cores at 50 cells). 512 Latin-hypercube samples and grouped Sobol indices (1280 runs) for the kinetic constants (jackknife covariance of the refit), pellet tortuosity (1–3) and the gas–solid transfer correlation. Spread: metallisation SD ≈ 0.002–0.004, DRI temperature SD ≈ 2–4 K. Tortuosity carries ~80 % of the variance (an upper bound, as it is coupled to the kinetic fit), the kinetics ~20 % and the transfer correlation < 6 %. The parameter spread is smaller than both the grid error at 50 cells and the 1D-vs-2D gap to REDUCTOR. `notebooks/06_uncertainty.ipynb`, `figures/uq_*.png`.

Results and limitations are in `docs/validation.md`.

## Development

```
python -m venv .venv
.venv\Scripts\activate            # Windows; source .venv/bin/activate elsewhere
pip install -e ".[dev]"
pytest && ruff check . && mypy
```

To rebuild the data from the primary sources: `python scripts/fetch_thermo_data.py`, `python scripts/fetch_references.py` (downloads the open-access papers into `references/`), `python scripts/extract_collision_integrals.py`, `python scripts/digitize_marsigny_fig6.py`, `python scripts/fit_kinetics.py`.

## Data sources

All kinetic, thermodynamic and transport parameters come from cited literature, listed in `data/SOURCES.md`.
