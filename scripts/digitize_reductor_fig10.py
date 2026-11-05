"""Metallisation vs cooling gas flow in case B of REDUCTOR (Marsigny, Mirgaux & Patisson 2025,
Metals 15:862, Fig. 10), digitised from the embedded image.

Usage:
    python scripts/digitize_reductor_fig10.py

Fig. 10 is a bar chart (cooling gas flow 500, 600, 700, 800 mol/s at a constant total gas flow
of 2800 mol/s, case B otherwise). The script extracts the embedded JPEG of page 14 with
pdfimages (poppler), finds
- the y-axis ticks (90-94 %) as short dark marks left of the axis; a straight-line fit gives
  the pixel -> % calibration (the residual of the fit is reported);
- each bar's top edge as the centre of the dark outline above the purple fill;
then checks the 500 and 700 mol/s bars against Table 1 (cases A and B: 93.3 and 93.7 %),
which were not used for the calibration, and writes
tests/data/reference_model_marsigny2025b_fig10.csv.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "references" / "marsigny2025_metals15_862.pdf"
OUT = ROOT / "tests" / "data" / "reference_model_marsigny2025b_fig10.csv"
PAGE = 14
TICKS = (94.0, 93.0, 92.0, 91.0, 90.0)  # top to bottom
FLOWS = (500.0, 600.0, 700.0, 800.0)
TABLE1 = {500.0: 93.3, 700.0: 93.7}  # cases A and B
CHECK_TOL = 0.05  # % metallisation


def extract_image() -> np.ndarray:
    exe = shutil.which("pdfimages")
    if exe is None:
        raise SystemExit("pdfimages (poppler) not found")
    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "img"
        subprocess.run(
            [exe, "-png", "-f", str(PAGE), "-l", str(PAGE), str(PDF), str(prefix)], check=True
        )
        files = sorted(Path(tmp).glob("img-*.png"))
        if len(files) != 1:
            raise SystemExit(f"expected one image on page {PAGE}, found {len(files)}")
        return np.asarray(Image.open(files[0]).convert("RGB")).astype(int)


def runs(idx: np.ndarray) -> list[np.ndarray]:
    """Split sorted indices into runs of consecutive values."""
    return np.split(idx, np.where(np.diff(idx) > 1)[0] + 1) if idx.size else []


def main() -> None:
    if not PDF.exists():
        raise SystemExit(f"{PDF.name} not found; run scripts/fetch_references.py")
    a = extract_image()
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    dark = (r + g + b) < 200
    purple = (b > 150) & (r > 70) & (r < 140) & (g < 120)

    axis = int(np.argmax(dark[:, : a.shape[1] // 3].sum(0)))
    tick_rows = [float(np.mean(t)) for t in runs(np.where(dark[:, axis - 5])[0])]
    if len(tick_rows) != len(TICKS):
        raise SystemExit(f"found {len(tick_rows)} ticks, expected {len(TICKS)}")
    slope, icpt = np.polyfit(tick_rows, TICKS, 1)
    resid = np.max(np.abs(np.polyval([slope, icpt], tick_rows) - TICKS))
    print(f"axis at column {axis}; ticks at rows {tick_rows}; calibration residual {resid:.4f} %")

    bars = [c for c in runs(np.where(purple.sum(0) > 50)[0]) if c.size > 40]
    if len(bars) != len(FLOWS):
        raise SystemExit(f"found {len(bars)} bars, expected {len(FLOWS)}")
    values = []
    for flow, cols in zip(FLOWS, bars, strict=True):
        inner = cols[5:-5]  # away from the side outlines
        tops = []
        for c in inner:
            fill_top = int(np.where(purple[:, c])[0].min())
            edge = np.where(dark[fill_top - 8 : fill_top, c])[0] + fill_top - 8
            if edge.size == 0:
                raise SystemExit(f"no outline above bar {flow:g}")
            tops.append(float(np.mean(edge)))
        row = float(np.median(tops))
        val = float(np.polyval([slope, icpt], row))
        values.append(val)
        print(f"  {flow:g} mol/s: top row {row:.1f} (spread {np.ptp(tops):.1f} px) -> {val:.3f} %")

    for flow, ref in TABLE1.items():
        got = values[FLOWS.index(flow)]
        status = "ok" if abs(got - ref) <= CHECK_TOL else "FAILED"
        print(f"check vs Table 1 at {flow:g} mol/s: {got:.3f} vs {ref} % {status}")
        if status != "ok":
            raise SystemExit("digitisation check failed")

    px = 1.0 / abs(slope)  # pixels per %
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write(
            "# REDUCTOR, case B with varied cooling gas flow (total gas flow 2800 mol/s): burden\n"
            "# metallisation (Marsigny, Mirgaux & Patisson 2025, Metals 15:862, Fig. 10,\n"
            "# CC BY 4.0), digitised by scripts/digitize_reductor_fig10.py (do not edit).\n"
            f"# Resolution: 1 px = {1 / px:.4f} %; checked against Table 1 at 500 and 700\n"
            f"# mol/s (within {CHECK_TOL} %).\n"
        )
        fh.write("cooling_flow_mol_s,bustle_flow_mol_s,metallization_pct\n")
        for flow, val in zip(FLOWS, values, strict=True):
            fh.write(f"{flow:g},{2800 - flow:g},{val:.2f}\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
