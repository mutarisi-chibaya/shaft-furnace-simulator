# Milestone 8 handoff (numerics and pressure) — work in progress

Status on 2026-11-06: **step 1 of 5 in progress.** Read `CLAUDE.md`, then `PLAN.md`
(Milestone 8, approved by the user), then this file. The user approved the plan, so no
re-approval is needed to continue it; a change of plan does need approval.

## Session 2 update (2026-11-06): band solved — read this first

The coupled solver now converges in the band (case B, 200 cells, 700 mol/s: 115 iterations
from the prolonged 100-cell start; sweep results in `build/m8_sweep.log`). Four causes, found
and fixed in this order (scripts `m8_stall.py`, `m8_cell.py`, `m8_step.py`, `m8_trace.py`,
`m8_zone.py` in `build/m8_scratch/`):

1. **Zero slope at the lower bound.** A u-row reads (1 - u^3) - X_up - kappa r(u), r ~ u^2,
   so dR/du ~ 0 at u = U_MIN: after a Newton overshoot to complete conversion, u could not
   come back although R > 0 (stall at 7e5 x tol, cell 100). Fix: `_release_bound` resets such
   u to cbrt(u^3 + R) before each step.
2. **Smooth min was negative when one rate is 0.** (a + b - sqrt(...))/2 gives -a eta^2/4 for
   b = 0, so a step behind one at equilibrium got a negative rate and wanted X < 0 (clipped,
   stall). Fix: p-norm min lo (1 + (lo/hi)^p)^(-1/p), p = 1/eta (0 if either is 0).
3. **FD step too large for the transition.** FD_STEP in u was 1e-7 (dX ~ 3e-7 near X = 0),
   ~10 % Jacobian error inside a 1e-5 smoothstep: Newton converged linearly, one cell per
   iteration (a front walking along the column). FD_STEP u = 1e-9 -> quadratic, 11 iterations.
4. Added (useful, not decisive): backtracking line search (`line_search`), more retries,
   return of the best state, per-iteration `callback`.

Tried and **not needed** (removed again): width continuation, secant predictor in the width,
and a softplus on the equilibrium kink of the rate law (at 1e-6 it shifted metallisation by
3.6e-6, above the 1e-6 target).

**Sweep, case B, 200 cells** (`build/m8_sweep.log`; start = alternating on 100 cells,
prolonged): coupled converges at 600, 625, 650, 675, 700, 725, 750 mol/s (solve 28-54 s,
96-175 iterations); alternating fails at 650/675/700 (400 outer iterations, nonsense end
states) but converges at 600/625/725/750, where coupled - alternating metallisation =
1.9e-7, 1.7e-7, 4.7e-7, **9.1e-7** (smoothing effect; 750 is close to the 1e-6 target) and
DRI/top temperatures agree to 1e-4 K. **775 mol/s: coupled fails from the 100-cell start**
(stuck at 2.5e-4 / 1.7 K) while alternating converges in 18 s: the coupled solver is a
finisher, not globally robust. Production therefore keeps the alternating solver first.

**Integration (done):** `scenario.run_scenario` and `uq.run_one`: alternating -> tolerance
retry -> `coupled.finish_coupled` from the end state if still unconverged. Converged results
are unchanged bit for bit. Not yet: `scripts/run_sensitivity.py` and
`scripts/run_cooling_continuation.py` (own solver calls), rerun of the Milestone 7 sweep and
the n.c. entries in docs/validation.md, timing of case A on 200 cells.

Smoothing effect (default `DEFAULT_SMOOTHING = (1e-5, 1e-5)`): case A 50 cells: met -4.2e-8,
DRI +9.1e-5 K, top -4.2e-6 K (1e-6: 4e-9 / 7e-7 K). Band 700 mol/s: met 1e-3 -> 1e-4 -> 1e-5:
0.8801277, 0.8800729, 0.8800691 (about 14x smaller per decade, so ~3e-7 from the limit).
Without smoothing the band point still stalls (no exact root of the discontinuous system).
Case A 50 cells, zero smoothing: coupled = alternating to 7e-11 (met), 7e-9 K.

New: `coupled.thermal_solution`, `coupled.finish_coupled`, `coupled.solve_thermal_coupled`
(alternating start on half the cells, coupled finish), `tests/test_coupled.py` (9 tests,
13 s).

## Goal (approved plan, short)

1. **Fully coupled Newton solver** (species + energy in one system) that converges where the
   alternating solver `energy.solve_thermal` cycles: case B sweep, cooling flow 625–750 mol/s,
   200 cells (`docs/validation.md`, Milestone 7).
2. **Second-order scheme** with a smooth limiter (van Albada), conservative (balances exact).
3. **Ergun pressure drop** (constants 150 / 1.75 quoted from Wagner 2008), local pressure in
   the transport properties (exact 1/p scaling, like the Milestone 6 factors); predicted bustle
   and cooling pressures compared with REDUCTOR Table 1, cases A–E (not tuned).
4. **Gas-loop energy terms** (heat recovery, compression) only if Shao et al. (2021),
   *Metall. Mater. Trans. B* 52:451, is openly accessible; otherwise skip and document.
5. **Rerun** Milestone 5 and 7 comparisons on converged grids; update tests, figures,
   notebooks 05 and 07, docs, README, app.

Acceptance (from the plan): all tests pass (regression values pinned to the old scheme
updated with reasons); coupled solver equals the old one where both converge (case A,
1e-8); balance closure 1e-12; observed order about 2 on the reaction-free reference;
smoothing effect on outputs below 1e-6; 625–750 mol/s converges on 200 cells; 100 cells
within 0.005 of the extrapolated metallisation; case A on 200 cells under 2 min. **If a target
is missed, report it honestly; do not relax it silently.**

## Done so far (uncommitted; 208 tests pass; ruff and mypy clean)

| File | Change |
|---|---|
| `src/shaftsim/coupled.py` (new) | Coupled solver. Unknowns per cell (X1, X2, X3, f, Ts, Tg); solved for core radii u = (1 - X)^(1/3); residual = exactly the existing species residual + energy residual (gas row negated, energy rows / heat-capacity flow -> kelvin); coloured finite-difference Jacobian (18 evaluations, block-tridiagonal); pseudo-transient continuation with sign-aware diagonal and SER (at least doubling dtau after a decrease); projection onto bounds (u >= U_MIN = 1e-8); active set freezing core radii at the bound. Public: `pack`, `unpack`, `to_internal`, `from_internal`, `residual`, `residual_internal`, `jacobian`, `project`, `solve_coupled(prob, y0, ..., smoothing=(delta, eta))`, `CoupledStats`. **No tests yet.** |
| `src/shaftsim/pellet.py` | `step_rates_array(..., one_minus=None, availability_width=0.0, min_width=0.0)`: optional exact 1 - X (precision near complete conversion) and a continuous availability constraint (smoothstep over width delta, smooth min with relative width eta). Defaults reproduce the old behaviour exactly. |
| `src/shaftsim/shaft.py` | `ShaftProblem.residual(y, one_minus=None)`; new fields `availability_width`, `min_width` **at the end** of the dataclass (positional construction elsewhere). |
| `src/shaftsim/energy.py` | `ThermalProblem` fields `availability_width`, `min_width` (in `TABLE_FREE_FIELDS`), passed to `species_problem`. |
| `PLAN.md` | Milestone 8 section added. |

## Findings so far (with numbers)

1. **Coupled solver works from a good start.** Case A, 50 cells, started from the alternating
   solution: converged in 7 iterations, 1.1 s, species residual 8e-11, energy 1.9e-7 K;
   metallisation equal to 1e-8 and temperatures to 2e-7 K. From a crude start (linear T, no
   reaction) it does not converge (weak global convergence) -> use grid continuation starts.
2. **First cause of failure (fixed): infinite slope at complete conversion.** The rate falls
   like (1 - X)^(2/3); Newton in X overshoots past 1 and stalls (first stall: X1 = 1 in cell 6,
   X2 = 1 in cell 22). Fixed with the u-variables + exact `one_minus` + `U_MIN` + active set.
3. **Second cause (open): the availability switch is discontinuous.** At the band state
   (case B sweep, 700 mol/s, 200 cells, start = 100-cell solution prolonged) the coupled
   solver without smoothing stalls at species 1.9e-6, energy 8e-3 K (metallisation 0.8807,
   DRI 70 °C). Diagnosis (`build/m8_scratch/m8_diag5.py`): in cell 19 near the top, where
   X2 = X3 = 0, the X3 residual jumps by 2.6e-4 for any step size (the same at FD step 1e-2 and
   1e-6): step 3 is capped by step 2 while X3 >= X2 and uncapped as soon as X2 > X3.
4. **Smoothing implemented, effect small, convergence not yet achieved.** Case A, 50 cells,
   coupled runs started from the exact solution, 300 iterations, **not converged to tolerance**
   (so these numbers are indicative only): delta = eta = 1e-6 -> metallisation +2.9e-9, DRI
   -2.9e-6 K; 1e-5 -> -5.9e-8, +5.6e-5 K; 1e-4 -> -4.2e-7, +4.0e-4 K. At the exact solution,
   delta = 1e-6 alone changes the species residual by 2.8e-6 (some cells near complete
   reduction have 0 < X2 - X3 < 1e-6, so the transition zone is occupied).
5. **Band with smoothing fails.** 700 mol/s, 200 cells, start = prolonged 100-cell solution:
   with (1e-4, 1e-4) 11 accepted / 12 rejected steps, residual species 1.8e-4, energy 3.4 K;
   with 1e-5 and 1e-6 worse (species ~8e-2, energy ~65 K) — the loop exits after 12
   consecutive rejections.
6. Earlier (Milestone 7, still true): damped, Anderson and a pseudo-transient energy-only
   prototype of the alternating scheme all end in a period-2 cycle in the band; outside the
   band, continuation up and down gives identical solutions (no hysteresis).

## Suggested next steps for the solver (time-box about 2 hours, then decide)

1. **Continuation in the smoothing width**: solve with delta = eta = 1e-3, then 1e-4, 1e-5,
   1e-6, each from the previous solution. A steep transition (width 1e-6 in X) is hard to
   reach directly.
2. **Jacobian accuracy**: the FD step in u is 1e-7, comparable to a 1e-6 transition in X; check
   FD against directional derivatives at smoothed states (as in `m8_diag5.py`); consider
   analytic derivatives of the rate law for the species rows, or FD steps scaled to delta.
3. **Globalisation**: add a backtracking line search on the Newton direction (now only dtau
   is reduced); allow more than 12 retries or restart from the best state.
4. **Better starts in the band**: continuation in the cooling flow with the coupled solver from
   600 mol/s (converged with the old solver) upward, instead of prolonging the 100-cell
   solution.
5. When it converges: compare with the old solver at converged points (case A, 50/200
   cells) to 1e-8, measure the smoothing effect on converged solutions (must be < 1e-6), add
   `tests/test_coupled.py`, then wire it into `solve_thermal` (alternating start, coupled
   finish) or as a separate entry point used by `scenario.run_scenario`, `uq.run_one` and the
   scripts.
6. **If it still fails after the time box**: keep the coupled solver as an optional polishing
   step (it is exact and fast from good starts), document the band as unsolved, and continue
   with steps 2–5 of the plan (Ergun does not depend on the solver).

## How to reproduce the experiments

Scripts copied to `build/m8_scratch/` (not committed; they hard-code `ROOT` and use
`scripts/run_uq.py:setup`):

| Script | What it does |
|---|---|
| `m8_try.py check 50` | Case A: alternating vs coupled (from the alternating solution and from a crude start) |
| `m8_try.py band 200 700` | Band point from the prolonged 100-cell solution, no smoothing |
| `m8_widths.py A` / `m8_widths.py 700` | Smoothing widths: effect on case A, convergence in the band |
| `m8_diag2.py` | Saves the stalled band state to `build/m8_band_state.npy`; prints the worst cells |
| `m8_diag5.py` | One Newton iteration at the stalled state; FD vs directional-derivative check (finds the jump) |
| `m8_dbg.py` | Residual with and without smoothing at the exact case A solution |
| `proto_ptc.py`, `diag_b2.py` | Milestone 7 prototypes (energy-only pseudo-transient; outer-iteration trajectory) |

Run with `.venv/Scripts/python -W ignore build/m8_scratch/<script> <args>`. A 200-cell band
run takes 1–3 min; the 100-cell start about 1 min.

## Gotchas

- **Shell:** bash heredocs in this environment mangle quotes and backslashes. Write patch
  scripts with the file-writing tool and run them; do not inline multi-line Python with quotes.
- `ShaftProblem` is constructed positionally (`energy.species_problem` passes 5 arguments):
  keep new fields at the end.
- Use `ThermalProblem.with_parameters(...)` to change table-free fields (shares the property
  tables); `dataclasses.replace` drops the cached tables (about 5 s to rebuild) unless you copy
  `__dict__["_tables"]` as `solve_thermal` does for coarse grids.
- Smoothing widths default to 0 everywhere = the original model; the existing 208 tests depend
  on that.
- Windows file names are case-insensitive (`app_x.png` overwrote `app_X.png` once).
- Full test suite: about 3.5 min (`pytest -q`); `tests/test_app.py` needs Streamlit (~30 s).
- Never invent parameters (`CLAUDE.md`); new constants (Ergun) need a quote in
  `data/*.yaml` checked by `tests/test_provenance.py` and a row in `data/SOURCES.md`.

## Not part of Milestone 8 (documented future work)

Two-zone or 2D model, broader kinetics (Wagner 2008 and Ranzani 2011 data are in
`references/`), plant data, wall heat losses.
