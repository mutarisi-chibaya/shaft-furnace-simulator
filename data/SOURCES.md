# Data sources

Every number in `data/*.yaml` and `tests/data/*.csv` must have an entry here.
The `source` field in the YAML files refers to the key in the first column.

## Thermodynamic data (Milestone 1)

`data/species.yaml` and `tests/data/janaf_reference.csv` are **generated** by
`scripts/fetch_thermo_data.py`, which downloads the pages below and extracts the numbers
programmatically (nothing is typed by hand). Rerun it to reproduce or audit the files.
Accessed 2026-11-04.

| Key | Used for | Reference | URL |
|-----|----------|-----------|-----|
| NIST-H2 | Shomate coefficients, H2(g), 298–6000 K | NIST Chemistry WebBook, SRD 69 (coefficients cite Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1 |
| NIST-H2O | Shomate coefficients, H2O(g), 500–6000 K | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=1 |
| NIST-Fe | Shomate coefficients, Fe(cr), α–δ phase, 298–1809 K (4 ranges; the overlapping "γ phase" column is not used) | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C7439896&Units=SI&Mask=2 |
| NIST-FeO | Shomate coefficients, stoichiometric FeO(cr), 298–1650 K (comparison only) | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1345251&Units=SI&Mask=2 |
| NIST-Fe3O4 | Shomate coefficients, magnetite, 298–3000 K | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1309382&Units=SI&Mask=2 |
| NIST-Fe2O3 | Shomate coefficients, hematite, 298–2500 K | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317608&Units=SI&Mask=2 |
| JANAF-Fe-001-fit | Shomate coefficients for wüstite Fe0.947O(cr), 298.15–1600 K, **fitted** by `scripts/fetch_thermo_data.py` to JANAF table Fe-001 (A–E least squares on Cp; F, G, H from the 298.15 K values) | Chase, M.W. (1998) NIST-JANAF Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9 | https://janaf.nist.gov/tables/Fe-001.txt |
| NIST-CaO | Shomate coefficients, CaO(cr), 298–3200 K (gangue, burden heat capacity, Milestone 4) | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1305788&Units=SI&Mask=2 |
| NIST-SiO2 | Shomate coefficients, quartz SiO2(cr), 298–1996 K (gangue) | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C14808607&Units=SI&Mask=2 |
| NIST-MgO | Shomate coefficients, MgO(cr), 298–3105 K (gangue) | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1309484&Units=SI&Mask=2 |
| NIST-Al2O3 | Shomate coefficients, α-Al2O3 (corundum), 298–2327 K (gangue; the γ, δ, κ fits on the same page are filtered out) | NIST Chemistry WebBook, SRD 69 (Chase 1998) | https://webbook.nist.gov/cgi/cbook.cgi?ID=C1344281&Units=SI&Mask=2 |
| JANAF-H-050 | Validation rows, H2 (ref) | Chase (1998) | https://janaf.nist.gov/tables/H-050.txt |
| JANAF-H-064 | Validation rows, H2O(g) | Chase (1998) | https://janaf.nist.gov/tables/H-064.txt |
| JANAF-Fe-004 | Validation rows, Fe(cr) α–δ | Chase (1998) | https://janaf.nist.gov/tables/Fe-004.txt |
| JANAF-Fe-018 | Validation rows, FeO(cr) | Chase (1998) | https://janaf.nist.gov/tables/Fe-018.txt |
| JANAF-Fe-032 | Validation rows, Fe3O4(cr) | Chase (1998) | https://janaf.nist.gov/tables/Fe-032.txt |
| JANAF-Fe-030 | Validation rows, Fe2O3(cr) | Chase (1998) | https://janaf.nist.gov/tables/Fe-030.txt |
| JANAF-Fe-001 | Validation rows, Fe0.947O(cr) (same table the fit uses; see docs/validation.md) | Chase (1998) | https://janaf.nist.gov/tables/Fe-001.txt |

Molar masses are the NIST WebBook "Molecular weight" values / 1000. For Fe0.947O:
0.947 M(Fe) + [M(FeO) − M(Fe)].

Note: the NIST Shomate coefficients are themselves fits to the JANAF tables (Chase 1998), so
the JANAF comparison checks the fits and this code, not two independent measurements.

## Constants

Molar gas constant R: CODATA exact value via `scipy.constants.R`.

## Pellet kinetics and transport (Milestone 2)

References are downloaded by `scripts/fetch_references.py` into `references/` (not committed).
Values in `data/kinetics.yaml` carry a verbatim quote that `tests/test_provenance.py` checks
against the PDF text layer. Values from scanned pages (no usable text layer) say so and are
cross-checked numerically where possible. Accessed 2026-11-05.

| Key | Used for | Reference | URL |
|-----|----------|-----------|-----|
| Marsigny2025-T2 / -T3 / -S2.1 / -S3.5.3 | Published kinetic triplets (Table 2), base case and pellet properties (Table 3, Sec. 2.1), tortuosity (Sec. 3.5.3), model form (Appendix A, Table A1, read from the page image), TGA curves (Fig. 6, Fig. 9, digitised by `scripts/digitize_marsigny_fig6.py`), repeatability of t50/t90 (Sec. 2.2) | Marsigny, A.; Letz, J.-B.; Mirgaux, O.; Patisson, F. (2025) Kinetics of the Reduction of Iron Ore Pellets with Hydrogen: A Parametric Experimental and Modeling Study. *Metals* 15, 1034. doi:10.3390/met15091034. CC BY 4.0 | https://hal.science/hal-05287190/document |
| Marsigny2025-thesis-eq12 | Ranz–Marshall constants (eq. 12, p. 34, read from the page image); layer-porosity method (Article 2 appendix, eq. 51) | Marsigny, A. (2025) Experimental and modeling study of the iron ore direct reduction by a H2-rich gas. PhD thesis, Université de Lorraine. CC BY 4.0 | https://theses.hal.science/tel-05795515/document |
| Wagner2008-II.4.2 / -T2 | Grain/crystallite structure (Sec. II.4.2), true densities and molar masses (Tableau 2, from COACH/THERMODATA), combination of the three steps (Sec. II.5.2) | Wagner, D. (2008) Étude expérimentale et modélisation de la réduction du minerai de fer par l'hydrogène. PhD thesis, INPL Nancy | https://theses.hal.science/tel-00280689/document |
| Ranzani2013 | Structure of the GM-LART model, choice of molar densities per term (Table 1) | Ranzani da Costa, A.; Wagner, D.; Patisson, F. (2013) Modelling a new, low CO2 emissions, hydrogen steelmaking process. *J. Clean. Prod.* 46, 27–35 (HAL author copy) | https://hal.science/hal-00943356/document |
| Svehla1962-TI | Lennard-Jones σ and ε/k of H2, H2O, He (p. 22), N2 (p. 23), read from the page image | Svehla, R.A. (1962) Estimated Viscosities and Thermal Conductivities of Gases at High Temperatures. NASA TR R-132 (public domain) | https://ntrs.nasa.gov/citations/19630012982 |
| MarreroMason1972 | Experimental correlations for D(He–H2), D(H2–N2): Table 12, eq. (4.3-1), journal p. 35, read from the page image | Marrero, T.R.; Mason, E.A. (1972) Gaseous Diffusion Coefficients. *J. Phys. Chem. Ref. Data* 1, 3–118 | https://srd.nist.gov/JPCRD/jpcrd1.pdf |
| Cantera-MMCollisionInt | Collision integrals Ω(2,2)*, A* (δ* = 0), extracted by `scripts/extract_collision_integrals.py` into `data/collision_integrals.csv` | Cantera v3.0.0, `src/transport/MMCollisionInt.cpp` (tables of Monchick & Mason 1961), BSD-3 | https://github.com/Cantera/cantera/blob/v3.0.0/src/transport/MMCollisionInt.cpp |
| NIST-He, NIST-N2, JANAF-He-001, JANAF-N-023 | Shomate coefficients, molar masses and JANAF checks of the inert gases | as above (NIST WebBook / NIST-JANAF) | https://webbook.nist.gov |

`data/kinetics_refit.yaml` is generated by `scripts/fit_kinetics.py`: the kinetic triplets
re-fitted to the digitised Fig. 6 curves (see docs/validation.md, Milestone 2).

## Shaft furnace reference case (Milestones 3-5)

| Key | Used for | Reference | URL |
|-----|----------|-----------|-----|
| Marsigny2025b-S2.1.1 / -S3.1 | Geometry, gas flow and composition, bustle temperature, pellet diameter, bed porosity, hematite content (text quotes in `data/furnace.yaml`) | Marsigny, A.; Mirgaux, O.; Patisson, F. (2025) Operating Solutions to Improve the Direct Reduction of Iron Ore by Hydrogen in a Shaft Furnace. *Metals* 15, 862. doi:10.3390/met15080862. CC BY 4.0 | https://hal.science/hal-05287188/document |
| Marsigny2025b-T1-image | Table 1 (case A inputs and outputs: gas flows, pressures, solid flow, top gas, DRI), read from the page image (page 7) because the text layer scrambles the columns. The outputs are in `tests/data/reference_model_marsigny2025b_caseA.csv`, written by `scripts/extract_reductor_table1.py`, which checks every value against the text layer and the internal consistency of the table (composition vs metallisation, GUD, flows) | same | same |
| Marsigny2025b-S3.1 (cone, cooling gas) | Cooling cone, 10 m long, outlet radius 1.31 m; cooling gas at 70 °C (Milestone 5, quotes in `data/furnace.yaml`) | same | same |
| Marsigny2025b-Figs3-9 | Case A maps of solid temperature, H2 fraction and the iron-phase mass fractions (Figs. 3, 4, 6–9), digitised by `scripts/digitize_reductor_maps.py` into cross-section averages: `tests/data/reference_model_marsigny2025b_caseA_profiles.csv` (colour-band resolution; check figure `figures/digitization_check_reductor.png`) | same | same |
| Wagner2008-III.5.3 | Wakao–Kaguei packed-bed correlation Nu = 2 + 1.1 Pr^(1/3) Re^0.6 (p. 89), used for mass transfer by analogy | Wagner (2008), see above | see above |
| Wagner2008-III.5.4 | Ergun equation for the bed pressure drop, constants 150 (viscous) and 1.75 (inertial), spherical pellets (pp. 90–91; Ergun 1952, ref. [60] there) — Milestone 8 | Wagner (2008), see above | see above |

The radius of 2.6 m ("2.6 m wide" in the text) was checked against Fig. 2 of the paper (radial axis 0–2.6 m, axial solid mass flux ~2 kg m⁻² s⁻¹ = 45 kg/s over π·2.6² m²).

## Energy balance and heat transfer (Milestone 4)

| Key | Used for | Reference | URL |
|-----|----------|-----------|-----|
| Svehla1962-B3 | Pure-gas thermal conductivity, modified Eucken correction λ = (R/M)[15/4 + f_int(Cp/R − 5/2)]η with f_int = 1.32 (Appendix B, eq. B3, and p. 9), read from the page image | Svehla (1962), NASA TR R-132, see above | https://ntrs.nasa.gov/citations/19630012982 |
| Svehla1962-TIII-H2 | `tests/data/svehla1962_table3_H2.csv`: viscosity and conductivity of H2 at 800/1000/1200 K (Table III, report p. 74), read from the page image; implementation check only (same theory and constants) | same | same |
| Wagner2008-III.5.7 | Wilke's rule for the mixture thermal conductivity (p. 92) | Wagner (2008), see above | see above |
| Wagner2008-III.5.3 | Wakao–Kaguei Nu = 2 + 1.1 Pr^(1/3) Re^0.6 for gas–solid heat transfer (the same entry as for mass transfer, above) | Wagner (2008), see above | see above |

The gangue mass fractions and the burden inlet temperature (20 °C) are from Marsigny2025b
(`data/furnace.yaml`, with quotes). The burden heat capacity uses the NIST Shomate data
above; the energy model bridges the enthalpy steps between NIST ranges over ±10 K
(`shaftsim.energy.enthalpy_ext`, a numerical choice, no new data).

## Uncertainty quantification (Milestone 6)

`data/uncertainty.yaml` lists the uncertain inputs; quotes are checked by `tests/test_provenance.py`.

| Key | Used for | Reference | URL |
|-----|----------|-----------|-----|
| Marsigny2025-S3.5.3 | Tortuosity range 1–3 ("Physically, tortuosity should be in the range of 1–3 [34]; hence, we chose a moderate value of 2") | Marsigny et al. (2025), *Metals* 15, 1034, see above | https://hal.science/hal-05287190/document |
| Marsigny2025-Fig6-refit | Kinetic parameter covariance: jackknife over the seven TGA temperature curves (leave one curve out) of the Milestone 2 refit, `data/kinetics_refit.yaml` (generated by `scripts/fit_kinetics.py` from the digitised Fig. 6) | derived from Marsigny et al. (2025), Fig. 6 | see above |
| Wagner2008-III.5.3, Marsigny2025-thesis-eq12 | Ends of the transfer-correlation range: Wakao–Kaguei (packed bed) and Ranz–Marshall (single sphere), both already in `data/transport.yaml`. No published accuracy figure for Wakao–Kaguei was found (searched 2026-11-05), so the spread between these two published choices is used | see above | see above |
| numerical-choice | Truncation of the kinetic normals at ±3σ (a sampling choice, not a physical parameter) | — | — |

## Operating cases, energy and CO2 (Milestone 7)

`data/furnace_cases.yaml` and `data/grid_intensity.yaml` are generated by scripts; `data/emissions.yaml` quotes are checked by `tests/test_provenance.py`.

| Key | Used for | Reference | URL |
|-----|----------|-----------|-----|
| Marsigny2025b-T1-image | Operating conditions of cases A–E (`data/furnace_cases.yaml`) and their outputs (`tests/data/reference_model_marsigny2025b_cases.csv`), read from the page image of Table 1 and checked against the PDF text layer by `scripts/extract_reductor_cases.py` | Marsigny, Mirgaux & Patisson (2025), *Metals* 15, 862, see above | https://hal.science/hal-05287188/document |
| Marsigny2025b-Fig10 | Metallisation vs cooling gas flow in case B (`tests/data/reference_model_marsigny2025b_fig10.csv`), digitised by `scripts/digitize_reductor_fig10.py` and checked against Table 1 | as above | as above |
| Marsigny2025b-S2.2 | Electrolysis electricity 18 MJ/Nm3 H2, heater efficiency 0.7, condenser exit temperature 70 °C (gas-loop model, Sec. 2.2) | as above | as above |
| DIN1343-normal-state | Definition of the normal cubic metre: 0 °C, 101.325 kPa (needed to convert MJ/Nm3 to J/mol). A convention, not a measured parameter | DIN 1343:1990, *Reference conditions, normal conditions, normal volume* | — |
| Marsigny2025-thesis-T20 | Energy breakdown of case A (electrolysis 9.73, heating 2.28, compression 0.25, total 12.26 GJ/t DRI; `tests/data/reference_energy_marsigny2025_thesis_t20.csv`) and direct CO2 of natural-gas shaft furnaces (0.29 and 0.45 t CO2/t DRI), Table 20, p. 125 | Marsigny, A. (2025) PhD thesis, Université de Lorraine, see above. CC BY 4.0 | https://theses.hal.science/tel-05795515/document |
| Marsigny2025-thesis-IEA2020 | Route intensities per tonne of crude steel: BF-BOF 2.2, natural-gas DRI-EAF 1.4 t CO2/t (the thesis citing IEA, *Iron and Steel Technology Roadmap*, 2020) | as above | as above |
| Ember2026-yearly | Grid emissions intensity (life-cycle, gCO2e/kWh) of South Africa, the world and France, 2015–2025 (`data/grid_intensity.yaml`, made by `scripts/fetch_grid_intensity.py`; file last modified 3 Nov 2026) | Ember, *Yearly electricity data*. CC BY 4.0 | https://ember-energy.org/data/yearly-electricity-data/ |
| Marsigny2025-thesis-dp | Pellet diameter range of the sensitivity study: "porous pellets in blast furnaces and shaft furnaces (dp = 10-15 mm)" (`data/sensitivity.yaml`) | Marsigny, A. (2025) PhD thesis, see above | https://theses.hal.science/tel-05795515/document |
