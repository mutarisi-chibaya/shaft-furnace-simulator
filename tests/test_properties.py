"""Tests for shaftsim.properties.

Two kinds of tests:

1. Mathematical consistency, using SYNTHETIC coefficients. These numbers are not physical
   data and are never used outside this file; they only check that the code implements
   the Shomate equations and the thermodynamic identities correctly (dH/dT = Cp,
   dS/dT = Cp/T, Hess's law, K = exp(-dG/RT)).
2. Validation against JANAF (tests/data/janaf_reference.csv) using the real NIST
   coefficients in data/species.yaml. These skip until the human provides both files.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from shaftsim import properties as pr
from shaftsim.properties import (
    MAGNETITE_TO_IRON,
    MAGNETITE_TO_WUSTITE,
    REDUCTION_STEPS,
    WUSTITE_TO_IRON,
    MissingDataError,
    Reaction,
    Species,
    species_from_dict,
)

from .conftest import DATA_DIR, TEST_DATA_DIR

# ----------------------------------------------------------------------------- fixtures


def _synthetic_range(t_min: float, t_max: float, seed: int) -> dict[str, float]:
    """Arbitrary, smooth, NON-PHYSICAL Shomate coefficients for identity tests."""
    rng = np.random.default_rng(seed)
    a, b, c, d, e, f, g, h = rng.uniform(-1.0, 1.0, 8) * [20, 10, 5, 1, 0.5, 100, 50, 100]
    return {
        "T_min": t_min,
        "T_max": t_max,
        "A": 30 + a,
        "B": b,
        "C": c,
        "D": d,
        "E": e,
        "F": f,
        "G": 150 + g,
        "H": h,
    }


def _synthetic_species(name: str, elements: dict[str, float], seed: int) -> Species:
    return species_from_dict(
        name,
        {
            "phase": "solid",
            "elements": elements,
            "shomate": [
                _synthetic_range(298.0, 900.0, seed),
                # same H (dfH298) in every range, as in NIST data
                _synthetic_range(900.0, 2000.0, seed + 100)
                | {"H": _synthetic_range(298.0, 900.0, seed)["H"]},
            ],
        },
    )


@pytest.fixture
def synth() -> dict[str, Species]:
    formulas: dict[str, dict[str, float]] = {
        "H2": {"H": 2},
        "H2O": {"H": 2, "O": 1},
        "Fe": {"Fe": 1},
        "FeO": {"Fe": 1, "O": 1},
        "Fe3O4": {"Fe": 3, "O": 4},
        "Fe2O3": {"Fe": 2, "O": 3},
        "Fe0.947O": {"Fe": 0.947, "O": 1},
    }
    return {n: _synthetic_species(n, el, i) for i, (n, el) in enumerate(formulas.items())}


T_GRID = np.array([300.0, 500.0, 850.0, 950.0, 1200.0, 1500.0, 1900.0])

# ----------------------------------------------------------------------------- identities


def test_enthalpy_derivative_is_heat_capacity(synth: dict[str, Species]) -> None:
    sp = synth["Fe3O4"]
    dT = 1e-3
    dh = np.asarray(pr.enthalpy(sp, T_GRID + dT)) - np.asarray(pr.enthalpy(sp, T_GRID - dT))
    np.testing.assert_allclose(dh / (2 * dT), pr.heat_capacity(sp, T_GRID), rtol=1e-6)


def test_entropy_derivative_is_cp_over_t(synth: dict[str, Species]) -> None:
    sp = synth["H2O"]
    dT = 1e-3
    ds = np.asarray(pr.entropy(sp, T_GRID + dT)) - np.asarray(pr.entropy(sp, T_GRID - dT))
    np.testing.assert_allclose(
        ds / (2 * dT), np.asarray(pr.heat_capacity(sp, T_GRID)) / T_GRID, rtol=1e-6
    )


def test_gibbs_is_h_minus_ts(synth: dict[str, Species]) -> None:
    sp = synth["Fe"]
    expected = np.asarray(pr.enthalpy(sp, T_GRID)) - T_GRID * np.asarray(pr.entropy(sp, T_GRID))
    np.testing.assert_allclose(pr.gibbs(sp, T_GRID), expected, rtol=1e-12)


def test_shomate_formula_by_hand(synth: dict[str, Species]) -> None:
    """Evaluate one point by hand, directly from the NIST formulas."""
    sp = synth["FeO"]
    r = sp.ranges[1]
    T = 1200.0
    t = T / 1000
    cp = r.A + r.B * t + r.C * t**2 + r.D * t**3 + r.E / t**2
    h = (r.A * t + r.B * t**2 / 2 + r.C * t**3 / 3 + r.D * t**4 / 4 - r.E / t + r.F) * 1000
    s = r.A * math.log(t) + r.B * t + r.C * t**2 / 2 + r.D * t**3 / 3 - r.E / (2 * t**2) + r.G
    assert pr.heat_capacity(sp, T) == pytest.approx(cp, rel=1e-12)
    assert pr.enthalpy(sp, T) == pytest.approx(h, rel=1e-12)
    assert pr.entropy(sp, T) == pytest.approx(s, rel=1e-12)


def test_sensible_enthalpy_is_nist_form(synth: dict[str, Species]) -> None:
    """H - H298 = enthalpy - 1000 * H coefficient, without evaluating at 298.15 K."""
    sp = synth["H2"]
    expected = np.asarray(pr.enthalpy(sp, T_GRID)) - 1000 * sp.ranges[0].H
    np.testing.assert_allclose(pr.sensible_enthalpy(sp, T_GRID), expected, rtol=1e-12)


def test_ranges_must_share_formation_enthalpy() -> None:
    r1, r2 = _synthetic_range(298, 800, 1), _synthetic_range(800, 1500, 1)
    r2["H"] = r1["H"] + 1.0
    with pytest.raises(ValueError, match="disagree on H"):
        species_from_dict("bad", {"phase": "solid", "elements": {"Fe": 1}, "shomate": [r1, r2]})


def test_boundary_uses_lower_range(synth: dict[str, Species]) -> None:
    sp = synth["Fe"]
    lower_only = Species("x", "solid", sp.elements, (sp.ranges[0],))
    assert pr.heat_capacity(sp, 900.0) == pr.heat_capacity(lower_only, 900.0)


def test_scalar_and_array_shapes(synth: dict[str, Species]) -> None:
    sp = synth["Fe"]
    assert isinstance(pr.heat_capacity(sp, 1000.0), float)
    grid = np.full((2, 3), 1000.0)
    assert np.asarray(pr.gibbs(sp, grid)).shape == (2, 3)


@pytest.mark.parametrize("T", [250.0, 2500.0, float("nan")])
def test_out_of_range_raises(synth: dict[str, Species], T: float) -> None:
    with pytest.raises(ValueError, match="outside data range"):
        pr.heat_capacity(synth["Fe"], T)


def test_non_contiguous_ranges_raise() -> None:
    with pytest.raises(ValueError, match="not contiguous"):
        species_from_dict(
            "bad",
            {
                "phase": "solid",
                "elements": {"Fe": 1},
                "shomate": [_synthetic_range(298, 800, 1), _synthetic_range(900, 1500, 1)],
            },
        )


def test_missing_values_are_all_reported(tmp_path: Path) -> None:
    f = tmp_path / "s.yaml"
    f.write_text(
        "species:\n"
        "  A:\n    phase: gas\n    elements: {H: 2}\n    shomate:\n"
        "      - {T_min: 298, T_max: 1000, A: null, B: 1, C: 1, D: 1, E: 1, F: 1, G: 1, H: 1}\n"
        "  B:\n    phase: gas\n    elements: {H: 2}\n    shomate: null\n",
        encoding="utf-8",
    )
    with pytest.raises(MissingDataError) as exc:
        pr.load_species(f)
    assert "A.shomate[0].A" in str(exc.value)
    assert "B.shomate" in str(exc.value)


# ----------------------------------------------------------------------------- reactions


ALL_STEPS = [*REDUCTION_STEPS, MAGNETITE_TO_IRON, pr.MAGNETITE_TO_FEO, pr.FEO_TO_IRON]


@pytest.mark.parametrize("rxn", ALL_STEPS)
def test_reactions_are_element_balanced(rxn: Reaction) -> None:
    """Uses the formulas in data/species.yaml (available before coefficients are filled)."""
    raw = yaml.safe_load((DATA_DIR / "species.yaml").read_text(encoding="utf-8"))["species"]
    species = {
        n: Species(n, e["phase"], e["elements"], (pr.ShomateRange(1, 2, *[0.0] * 8),))
        for n, e in raw.items()
    }
    assert all(abs(v) < 1e-12 for v in pr.element_imbalance(rxn, species).values())


def test_wustite_constant_matches_data() -> None:
    raw = yaml.safe_load((DATA_DIR / "species.yaml").read_text(encoding="utf-8"))["species"]
    assert raw[pr.WUSTITE]["elements"] == {"Fe": pr.X_WUSTITE, "O": 1}


def test_disproportionation_is_balanced(synth: dict[str, Species]) -> None:
    for w in ("FeO", "Fe0.947O"):
        rxn = pr.wustite_disproportionation(synth, w)
        assert all(abs(v) < 1e-12 for v in pr.element_imbalance(rxn, synth).values())


def test_equilibrium_constant_definition(synth: dict[str, Species]) -> None:
    T = 1100.0
    dg = float(pr.reaction_gibbs(WUSTITE_TO_IRON, synth, T))
    assert pr.equilibrium_constant(WUSTITE_TO_IRON, synth, T) == pytest.approx(
        math.exp(-dg / (pr.R_GAS * T)), rel=1e-12
    )


def test_gibbs_helmholtz_consistency(synth: dict[str, Species]) -> None:
    """dG_r = dH_r - T dS_r for a reaction."""
    for rxn in REDUCTION_STEPS:
        dg = np.asarray(pr.reaction_gibbs(rxn, synth, T_GRID))
        dh = np.asarray(pr.reaction_enthalpy(rxn, synth, T_GRID))
        ds = np.asarray(pr.reaction_entropy(rxn, synth, T_GRID))
        np.testing.assert_allclose(dg, dh - T_GRID * ds, rtol=1e-10, atol=1e-6)


def test_hess_law(synth: dict[str, Species]) -> None:
    """Fe3O4 + 4 H2 -> 3 Fe + 4 H2O equals (Fe3O4 -> wüstite) + (3/x) (wüstite -> Fe)."""
    k_direct = np.log(np.asarray(pr.equilibrium_constant(MAGNETITE_TO_IRON, synth, T_GRID)))
    k_steps = np.log(
        np.asarray(pr.equilibrium_constant(MAGNETITE_TO_WUSTITE, synth, T_GRID))
    ) + 3 / pr.X_WUSTITE * np.log(
        np.asarray(pr.equilibrium_constant(WUSTITE_TO_IRON, synth, T_GRID))
    )
    np.testing.assert_allclose(k_direct, k_steps, rtol=1e-9, atol=1e-9)


def test_equilibrium_ratio_and_fraction(synth: dict[str, Species]) -> None:
    T = 1000.0
    k = float(pr.equilibrium_constant(MAGNETITE_TO_IRON, synth, T))
    ratio = float(pr.equilibrium_h2o_h2_ratio(MAGNETITE_TO_IRON, synth, T))
    assert ratio == pytest.approx(k**0.25, rel=1e-12)
    y = float(pr.equilibrium_h2_fraction(MAGNETITE_TO_IRON, synth, T))
    assert y == pytest.approx(1 / (1 + ratio), rel=1e-12)
    assert 0.0 < y < 1.0


def test_equilibrium_ratio_rejects_non_reduction(synth: dict[str, Species]) -> None:
    with pytest.raises(ValueError):
        pr.equilibrium_h2o_h2_ratio(Reaction("x", {"FeO": -4, "Fe": 1, "Fe3O4": 1}), synth, 1000)


def test_wustite_stability_shape(synth: dict[str, Species]) -> None:
    assert isinstance(pr.wustite_is_stable(synth, 1000.0), bool)
    assert np.asarray(pr.wustite_is_stable(synth, T_GRID)).shape == T_GRID.shape


# ----------------------------------------------------------------------------- JANAF validation

JANAF_CSV = TEST_DATA_DIR / "janaf_reference.csv"

# Tolerances (see docs/validation.md, Milestone 1, for the measured errors).
# - S and H - H298: 0.1 % (PLAN.md). H - H298 also gets an absolute floor of 0.05 kJ/mol,
#   because it is ~0 near 298 K where a relative test is meaningless.
# - Cp: 0.2 %. Shomate polynomials are least-squares fits to the JANAF Cp; their residuals
#   reach 0.13 % (FeO). Cp itself does not enter K; H and S (its integrals) do.
# - K: 1 % in K = 0.01 R T in dG_r (~83 J/mol at 1000 K).
REL_TOL_S_H = 1e-3
REL_TOL_CP = 2e-3
ABS_TOL_H_KJ = 0.05
REL_TOL_K = 1e-2

# KNOWN DEVIATION of the NIST Fe (alpha-delta) Shomate fit around the Curie point (1042 K).
# NIST's 700-1042 and 1042-1100 K ranges do not join continuously (H jumps of 141, -212 and
# -365 J/mol at 700, 1042 and 1100 K, where JANAF has no latent heat), and inside
# 700-1100 K the fit deviates from JANAF by up to 0.37 kJ/mol in H and 0.48 % in S. Cp also
# deviates by up to 0.49 % at 1300-1500 K. These are properties of the published
# coefficients, not of this code; they are kept visible here instead of being refitted.
FE_WINDOW_K = (700.0, 1100.0)  # (exclusive, inclusive)
FE_REL_TOL_S = 5e-3
FE_ABS_TOL_H_KJ = 0.4
FE_REL_TOL_CP = 5e-3
# Resulting error in dG_r per mol of Fe in the reaction (G error of Fe = dH - T dS).
FE_ABS_TOL_DG_PER_FE = 400.0  # J/mol


def _in_fe_window(T: float) -> bool:
    return FE_WINDOW_K[0] < T <= FE_WINDOW_K[1]


@pytest.fixture(scope="module")
def janaf() -> pd.DataFrame:
    if not JANAF_CSV.exists():
        pytest.skip("tests/data/janaf_reference.csv not provided yet (see tests/data/README.md)")
    return pd.read_csv(JANAF_CSV)


def test_properties_match_janaf(real_species: dict[str, Species], janaf: pd.DataFrame) -> None:
    for row in janaf.to_dict("records"):
        sp = real_species[str(row["species"])]
        T = float(row["T_K"])
        tol_cp, tol_s, tol_h = REL_TOL_CP, REL_TOL_S_H, ABS_TOL_H_KJ
        if sp.name == "Fe":
            tol_cp = FE_REL_TOL_CP
            if _in_fe_window(T):
                tol_s, tol_h = FE_REL_TOL_S, FE_ABS_TOL_H_KJ
        assert pr.heat_capacity(sp, T) == pytest.approx(row["Cp_J_molK"], rel=tol_cp), (sp.name, T)
        assert pr.entropy(sp, T) == pytest.approx(row["S_J_molK"], rel=tol_s), (sp.name, T)
        assert float(pr.sensible_enthalpy(sp, T)) / 1000 == pytest.approx(
            row["H_minus_H298_kJ_mol"], rel=REL_TOL_S_H, abs=tol_h
        ), (sp.name, T)


def test_fe_deviation_is_still_present(
    real_species: dict[str, Species], janaf: pd.DataFrame
) -> None:
    """Guards the documentation: if the Fe data change and the deviation disappears,
    the wider tolerance above and the note in docs/validation.md must be removed."""
    fe = janaf[(janaf.species == "Fe") & (janaf.T_K > 700) & (janaf.T_K <= 1100)]
    worst = max(
        abs(float(pr.sensible_enthalpy(real_species["Fe"], T)) / 1000 - h)
        for T, h in zip(fe.T_K, fe.H_minus_H298_kJ_mol, strict=True)
    )
    assert worst > ABS_TOL_H_KJ


def test_wustite_stability_matches_janaf(
    real_species: dict[str, Species], janaf: pd.DataFrame
) -> None:
    """Sign of dG for 4 Fe(x)O -> (4x-3) Fe + Fe3O4 agrees with JANAF at every tabulated T,
    for both Fe0.947O and stoichiometric FeO."""
    dfg = janaf.pivot_table(index="T_K", columns="species", values="dfG_kJ_mol")
    for w in ("Fe0.947O", "FeO"):
        rxn = pr.wustite_disproportionation(real_species, w)
        for T_label, v in dfg[list(rxn.stoich)].dropna().iterrows():
            T = float(str(T_label))
            janaf_stable = sum(n * v[c] for c, n in rxn.stoich.items()) > 0
            assert pr.wustite_is_stable(real_species, T, w) == janaf_stable, (w, T)


def test_real_sensible_enthalpy_zero_at_298(real_species: dict[str, Species]) -> None:
    """NIST F coefficients give H(298.15 K) - H298 close to 0. Measured offsets: Fe3O4 -7 and
    Fe2O3 -26 J/mol (the Fe2O3 offset is also seen against JANAF at every T). Must stay
    inside the enthalpy tolerance used against JANAF."""
    for sp in real_species.values():
        if sp.T_min <= 298.15:
            assert abs(float(pr.sensible_enthalpy(sp, 298.15))) < 1000 * ABS_TOL_H_KJ, sp.name


def test_equilibrium_constants_match_janaf(
    real_species: dict[str, Species], janaf: pd.DataFrame
) -> None:
    """K from NIST Shomate vs K from JANAF delta_f G (independent path)."""
    dfg = janaf.pivot_table(index="T_K", columns="species", values="dfG_kJ_mol")
    checked = 0
    for rxn in ALL_STEPS:
        cols = list(rxn.stoich)
        sub = dfg[cols].dropna() if set(cols) <= set(dfg.columns) else dfg.iloc[0:0]
        for T_label, values in sub.iterrows():
            T = float(str(T_label))
            dg_ref = 1000 * sum(rxn.stoich[c] * values[c] for c in cols)
            k_ref = math.exp(-dg_ref / (pr.R_GAS * T))
            n_fe = abs(rxn.stoich.get("Fe", 0.0))
            if n_fe and _in_fe_window(T):
                dg = float(pr.reaction_gibbs(rxn, real_species, T))
                assert abs(dg - dg_ref) <= n_fe * FE_ABS_TOL_DG_PER_FE, (rxn.name, T)
            else:
                k = float(pr.equilibrium_constant(rxn, real_species, T))
                assert k == pytest.approx(k_ref, rel=REL_TOL_K), (rxn.name, T)
            checked += 1
    if checked == 0:
        pytest.skip("janaf_reference.csv has no temperature with dfG for all species of a step")


def test_real_wustite_step_is_endothermic(real_species: dict[str, Species]) -> None:
    """Sanity check on real data, not a fitted target: FeO + H2 -> Fe + H2O is endothermic
    at typical shaft temperatures (well established; only the sign is tested)."""
    assert float(pr.reaction_enthalpy(WUSTITE_TO_IRON, real_species, 1000.0)) > 0
