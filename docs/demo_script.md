# Demo video script: hydrogen shaft furnace simulator

About 4½ minutes. Each scene gives what to show on screen and the narration. Narration lines
start with `>`; they are written as they should be spoken (numbers and units spelled out where
a voice would otherwise stumble) and are read by `scripts/make_voiceover.py` to produce the
voice-over (South African English male voice). All numbers are what the app shows on its
default 50-cell grid.

Sources of the spoken facts: blast furnace 2.2 t CO₂ per tonne of crude steel (IEA 2020, via
`data/emissions.yaml`); Saldanha and HyIron (links in `README.md`); coal 81.4 % of South
Africa's electricity generation in 2025 (Ember yearly data, the file downloaded by
`scripts/fetch_grid_intensity.py`); everything else from the app itself.

Automatic version: with the app running on port 8599, `python scripts/record_demo.py` drives
the app through these scenes and writes `build/demo/demo.mp4`, with the narration as captions.
If `build/voiceover/scene_XX.mp3` exist (from `scripts/make_voiceover.py`), each scene is timed
to its audio and the voice-over is mixed in.

Recording by hand instead? Before recording: run the app (`streamlit run app/streamlit_app.py`), open it at 1600×1000 or
larger, select case A and run each scenario once so the results are cached and appear
instantly on camera.

## Scene 1: Title (0:00–0:15)

**Screen:** title card: "Hydrogen shaft furnace simulator — a physics-based model of green
iron making". Then fade to the app with case A.

> Making steel from iron ore in a blast furnace emits about two tonnes of C-O-two for every tonne of steel. One way to cut that is to make the iron with hydrogen instead of coal. This is a simulator of the furnace at the heart of that process.

## Scene 2: Why it matters (0:15–0:45)

**Screen:** the furnace schematic (`figures/furnace_schematic.svg`) while the furnace is described, then a card with the projects.

> In a hydrogen shaft furnace, iron ore pellets move down while hot hydrogen flows up and strips the oxygen out. The product is direct reduced iron, and the main by-product is water. Projects are under way close to home: ArcelorMittal South Africa is studying hydrogen-based iron at Saldanha, and HyIron is building a plant in Namibia. Before anyone builds one, they need fast and honest models.

## Scene 3: What the simulator is (0:45–1:15)

**Screen:** the app's "Assumptions and limits" panel briefly, then `data/SOURCES.md` scrolling.

> This is a one-dimensional, steady-state model of the whole furnace. It tracks the three reduction steps inside the pellets, heat and mass transfer between gas and solid, and the energy balance of both phases. Every physical parameter comes from a cited source and is checked by a script, and the mass and energy balances close to machine precision.

## Scene 4: The reference case (1:15–2:00)

**Screen:** case A selected. Point at the six tiles from left to right, then at the blue
REDUCTOR box, then trace the three profile charts from the top of the furnace to the bottom.

> Here is the reference case from the literature. Ninety-four point six percent metallisation, forty-nine kilograms of hydrogen per tonne of iron, top gas at three hundred and sixty-three degrees. The published two-dimensional model, REDUCTOR, reports ninety-three point three percent for the same furnace.
> The profiles show where things happen. Hematite is reduced near the top. Magnetite and wüstite react around the bustle, where the hot gas enters. The last of the wüstite is reduced in the cone at the bottom.

## Scene 5: Turning the main dial (2:00–2:40)

**Screen:** drag the bustle gas temperature from 890 to 840 °C; wait for the update; then
to 920 °C.

> The bustle gas temperature is the main dial. Drop it by fifty degrees and metallisation falls to eighty percent. Raise it to nine hundred and twenty and the iron is fully metallised. The last reduction step needs heat as well as hydrogen, and the hot gas is where that heat comes from.

## Scene 6: Where one dimension is not enough (2:40–3:20)

**Screen:** select case B. Point at the DRI outlet tile (74 °C) and the metallisation tile
(87.6). Then open the conversion chart and point at the cone.

> Now case B, where more of the gas goes in at the bottom as cooling gas. Here the model and the published two-dimensional model disagree. In one dimension, the extra cold gas chills the whole cone, the last wüstite stops reacting, and metallisation drops to eighty-seven point six percent. The two-dimensional model keeps a hot core and stays near ninety-four. The app shows the disagreement, and the documentation explains it.

## Scene 7: The South African question (3:20–4:10)

**Screen:** back to case A. Show the CO₂ chart; select World, then France, then Custom and slide
to about 130 g/kWh.

> Hydrogen iron is only as clean as the electricity behind it. Electrolysis, gas heating and compression need about twelve and a half gigajoules of electricity per tonne. On South Africa's grid, which is still more than eighty percent coal, that is two point four tonnes of C-O-two per tonne of iron. That is about five times the direct emissions of a natural gas furnace. On France's low-carbon grid it is zero point one four tonnes. The break-even with natural gas is about one hundred and thirty grams per kilowatt-hour, so for South Africa, green iron means dedicated renewable power, not the average grid.

## Scene 8: How far to trust it (4:10–4:35)

**Screen:** the "Assumptions and limits" panel, scrolled slowly; then `docs/validation.md`.

> A model is only useful if you know its limits. This one is one-dimensional and steady. Its fast grid overestimates metallisation by about four points, and a few extreme operating points still do not converge. When a run does not converge, the app shows an error instead of a number. All of this is documented, with the tests behind it.

## Scene 9: Close (4:35–4:50)

**Screen:** the repository page with the README figures; app link on screen.

> The code, the data sources, the validation and the app are all in the repository, ready to run. Thanks for watching.
