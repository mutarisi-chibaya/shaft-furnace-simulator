# Project rules for Claude Code

Purpose: a physics-based 1D simulator of a hydrogen direct-reduction shaft furnace, built as a serious portfolio project. Work milestone by milestone as described in `PLAN.md`.

## Structure
- Installable package `shaftsim` (src layout, `pyproject.toml`).
- Parameters in `data/` (YAML), every source listed in `data/SOURCES.md`.
- Tests in `tests/`, reference/validation values in `tests/data/` (provided by the human).
- Streamlit app in `app/`, notebooks in `notebooks/`, figures in `figures/`, write-up in `docs/validation.md`.

## Conventions
- Python 3.11+. Tools: pytest, ruff, mypy. Dependencies: numpy, scipy, pandas, matplotlib, pyyaml, streamlit, plotly, SALib.
- Units: SI internally (J, mol, K, Pa, m, s, kg). degC and tonnes only for display.
- Every public function has type hints, a docstring stating units, and a test.
- Small, pure modules; no global state. Numerical choices are explained in docstrings.

## Hard rules
- NEVER invent physical, kinetic or thermodynamic parameters. All come from `data/*.yaml`, which the human fills in from cited sources. If something is missing, stop and say exactly what is needed.
- Validation targets come only from CSV files the human provides in `tests/data/`.
- Every model must close its mass and energy balances; test this.
- Always propose a plan for a milestone and wait for approval before implementing it.
- Never overstate what the simulator can do. Document assumptions and limitations honestly.
