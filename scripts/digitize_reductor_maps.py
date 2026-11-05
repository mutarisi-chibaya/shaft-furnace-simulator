"""Digitise the case A column of the REDUCTOR maps (Marsigny et al. 2025b, Metals 15:862,
Figs. 3, 4 and 6-9) into cross-section-averaged axial profiles.

Usage:
    python scripts/digitize_reductor_maps.py

Method (everything is detected from the images except the colour-bar labels, which are read
from the figures and listed in FIGURES below):
1. The figure images are extracted from the PDF with ``pdfimages`` (poppler). Each map is
   split over two images (z = 15 .. about -2.7 m and about -2.5 .. -10 m).
2. Axes: the z axis of case A is the thick dark column at the left; the major ticks (long dark
   marks to its right) are z = 15, 10, 5, 0 in the upper image and -5 in the lower one, whose
   bottom axis line is z = -10. z(y) is a least-squares line through these. The radial scale
   is the same as the vertical one (equal aspect: the coloured region of the cylinder is
   2.6 m wide at this scale, and the R = 2 tick of the lower image is checked).
3. Colour bar: the bar at the right is cut into uniform colour blocks (top to bottom); there
   must be one block more than labels. Labels sit on block boundaries, so block k means a
   value between label k-1 and label k; the open top and bottom blocks are bounded by the
   physical limits given below (mass fractions: the stoichiometric maximum of the phase for
   this burden, from data/species.yaml and data/furnace.yaml).
4. Every pixel of the case A panel is assigned to the nearest colour block (RGB distance
   < 60, otherwise ignored: JPEG edges, axis), and its value is the block midpoint (band
   edges give the lower and upper bounds).
5. Cross-section average of each row: sum(v r) / sum(r) over the classified pixels (r the
   pixel-centre radius), i.e. area-weighted for an axisymmetric field; plus the centre value
   (r < 0.15 m) and the wall value (outer 0.15 m).

Output: tests/data/reference_model_marsigny2025b_caseA_profiles.csv (one row per quantity and
z, every 0.1 m) and figures/digitization_check_reductor.png. The resolution is the colour
band: +/- 50 K for temperature, +/- 0.025 for the H2 fraction, +/- 0.05 for mass fractions.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from numpy.typing import NDArray
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "references" / "marsigny2025_metals15_862.pdf"
OUT = ROOT / "tests" / "data" / "reference_model_marsigny2025b_caseA_profiles.csv"
FIG = ROOT / "figures" / "digitization_check_reductor.png"
IntArray = NDArray[np.int64]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class MapFigure:
    """One map figure: the PDF page, image indices and colour-bar labels (top to bottom)."""

    number: int
    quantity: str
    unit: str
    page: int
    images: tuple[int, int]  # (upper, lower) index among the images of that page
    labels: tuple[float, ...]
    lower: float  # physical bound below the lowest label
    upper: float  # physical bound above the highest label


def _steps(hi: float, lo: float, step: float) -> tuple[float, ...]:
    n = round((hi - lo) / step)
    return tuple(round(hi - k * step, 10) for k in range(n + 1))


def stoichiometric_maxima() -> dict[str, float]:
    """Largest mass fraction each iron phase can reach in this burden (all iron in it)."""
    sp = yaml.safe_load((ROOT / "data" / "species.yaml").read_text(encoding="utf-8"))["species"]
    fu = yaml.safe_load((ROOT / "data" / "furnace.yaml").read_text(encoding="utf-8"))
    m = {k: float(sp[k]["molar_mass"]) for k in ("Fe2O3", "Fe3O4", "FeO", "Fe")}
    w = float(fu["solid"]["w_Fe2O3"]["value"])
    gangue = sum(float(v["value"]) for v in fu["gangue"].values())
    fe = 2.0 * w / m["Fe2O3"]  # mol Fe per kg of fresh burden
    out = {}
    for k, fe_per in (("Fe2O3", 2.0), ("Fe3O4", 3.0), ("FeO", 1.0), ("Fe", 1.0)):
        mass = fe / fe_per * m[k]
        out[k] = mass / (mass + gangue)
    return out


def figures() -> tuple[MapFigure, ...]:
    mx = stoichiometric_maxima()
    t_in = yaml.safe_load((ROOT / "data" / "furnace.yaml").read_text(encoding="utf-8"))[
        "solid_inlet"
    ]["temperature"]["value"]
    frac = _steps(0.9, 0.1, 0.1)
    return (
        MapFigure(3, "T_solid", "degC", 10, (0, 1), _steps(1200, 100, 100), t_in - 273.15, 1300),
        MapFigure(4, "x_H2", "-", 10, (2, 3), _steps(0.95, 0.30, 0.05), 0.25, 1.0),
        MapFigure(6, "w_Fe2O3", "-", 11, (2, 3), frac, 0.0, mx["Fe2O3"]),
        MapFigure(7, "w_Fe3O4", "-", 12, (0, 1), frac, 0.0, mx["Fe3O4"]),
        MapFigure(8, "w_FeO", "-", 12, (2, 3), frac, 0.0, mx["FeO"]),
        MapFigure(9, "w_Fe", "-", 13, (0, 1), frac, 0.0, mx["Fe"]),
    )


def extract_images(page: int, tmp: Path) -> list[Path]:
    exe = shutil.which("pdfimages")
    if exe is None:
        raise SystemExit("pdfimages (poppler) is required")
    prefix = tmp / f"p{page}"
    subprocess.run([exe, "-j", "-f", str(page), "-l", str(page), str(PDF), str(prefix)], check=True)
    return sorted(tmp.glob(f"p{page}-*.jpg"))


def load(path: Path) -> IntArray:
    return np.asarray(Image.open(path).convert("RGB")).astype(np.int64)


def axis_column(a: IntArray) -> tuple[float, int]:
    """Centre of the thick z axis of case A and the first column right of it."""
    dark = a.sum(axis=2) < 150
    cnt = dark[:, :350].sum(axis=0)
    cols = np.where(cnt > 0.9 * cnt.max())[0]
    run = cols[cols >= cols[-1] - 15]
    run = run[np.concatenate([[True], np.diff(run) == 1])]
    return float(run.mean()), int(run.max()) + 1


def major_ticks(a: IntArray, right: int) -> list[float]:
    """Rows of the long tick marks right of the axis (length >= 20 px)."""
    dark = a.sum(axis=2) < 150
    runs = np.zeros(a.shape[0], dtype=int)
    for y in range(a.shape[0]):
        k = 0
        while right + k < a.shape[1] and dark[y, right + k]:
            k += 1
        runs[y] = k
    ys = np.where(runs >= 20)[0]
    groups = np.split(ys, np.where(np.diff(ys) > 2)[0] + 1)
    return [float(g.mean()) for g in groups if g.size]


def colour_blocks(a: IntArray, n_blocks: int) -> FloatArray:
    """Median RGB of the colour-bar blocks, top to bottom."""
    sat = (a.max(axis=2) - a.min(axis=2)) > 60
    mid = a.shape[0] // 3
    xs = np.where(sat[mid] | sat[a.shape[0] // 2])[0]
    segs = np.split(xs, np.where(np.diff(xs) > 3)[0] + 1)
    bar = segs[-1]
    xc = int(bar.mean())
    col = a[:, xc - 3 : xc + 4].mean(axis=1)
    ok = (col.max(axis=1) - col.min(axis=1)) > 60
    blocks: list[list[FloatArray]] = []
    prev: FloatArray | None = None
    for y in range(a.shape[0]):
        if not ok[y]:
            prev = None
            continue
        c = col[y]
        if prev is None or np.abs(c - prev).sum() > 40:
            blocks.append([])
        blocks[-1].append(c)
        prev = c
    meds = [np.median(np.array(b), axis=0) for b in blocks if len(b) >= 5]
    if len(meds) != n_blocks:
        raise SystemExit(f"colour bar: found {len(meds)} blocks, expected {n_blocks}")
    return np.array(meds)


def band_values(fig: MapFigure) -> tuple[FloatArray, FloatArray, FloatArray]:
    """(mid, low, high) per colour block, top block first."""
    edges = (fig.upper, *fig.labels, fig.lower)
    hi = np.array(edges[:-1], dtype=float)
    lo = np.array(edges[1:], dtype=float)
    return 0.5 * (hi + lo), lo, hi


def classify(panel: IntArray, colours: FloatArray) -> tuple[IntArray, FloatArray]:
    d = np.linalg.norm(panel[..., None, :] - colours[None, None, :, :], axis=-1)
    return np.argmin(d, axis=-1), np.min(d, axis=-1)


@dataclass(frozen=True)
class PanelCheck:
    """What the check figure needs from one image: pixels, axis and colour calibration."""

    image: IntArray
    x_axis: float
    slope: float
    intercept: float
    px_per_m: float
    colours: FloatArray


def digitise(fig: MapFigure, paths: list[Path]) -> tuple[pd.DataFrame, dict[str, object]]:
    up, low = (load(paths[i]) for i in fig.images)
    colours = colour_blocks(up, len(fig.labels) + 1)
    mid, lo, hi = band_values(fig)
    rows = []
    check: dict[str, object] = {}
    scale = None
    for part, a in (("upper", up), ("lower", low)):
        xc, right = axis_column(a)
        ticks = major_ticks(a, right)
        labels: tuple[float, ...] = (15.0, 10.0, 5.0, 0.0)
        if part == "lower":
            # bottom axis line (z = -10): the long horizontal dark run near the bottom
            dark = a.sum(axis=2) < 150
            line = np.where(dark[:, right : right + 120].sum(axis=1) > 100)[0]
            ticks = [*[t for t in ticks if t < line.min() - 5][:1], float(line.mean())]
            labels = (-5.0, -10.0)
        if len(ticks) != len(labels):
            raise SystemExit(f"Fig. {fig.number} {part}: ticks {ticks} vs labels {labels}")
        slope, icpt = np.polyfit(np.array(ticks), np.array(labels), 1)  # z = slope y + icpt
        px_per_m = -1.0 / slope
        scale = px_per_m if scale is None else scale
        sat = (a.max(axis=2) - a.min(axis=2)) > 60
        if part == "upper":
            # the two images abut: the seam is the last row of the upper one (its last two
            # rows are dropped, JPEG edge)
            z_seam = float(slope * (a.shape[0] - 1) + icpt)
        for y in range(a.shape[0]):
            z = slope * y + icpt
            if part == "upper" and (z > 14.95 or z < z_seam + 0.05):
                continue
            if part == "lower" and (z >= z_seam + 0.05 or z < -9.95):
                continue
            xs = right
            while xs < a.shape[1] and not sat[y, xs]:
                xs += 1
                if xs > right + 10:
                    break
            xe = xs
            while xe < a.shape[1] and sat[y, xe]:
                xe += 1
            if xe - xs < 5:
                continue
            seg = a[y, right:xe][None, :, :]
            k, dist = classify(seg, colours)
            k, dist = k[0], dist[0]
            r = (np.arange(right, xe) + 0.5 - xc) / px_per_m
            use = dist < 60.0
            if use.sum() < 5:
                continue
            w = r[use]
            v_mid, v_lo, v_hi = mid[k[use]], lo[k[use]], hi[k[use]]
            mean = float(np.sum(v_mid * w) / np.sum(w))
            rows.append(
                {
                    "z": z,
                    "mean": mean,
                    "mean_lo": float(np.sum(v_lo * w) / np.sum(w)),
                    "mean_hi": float(np.sum(v_hi * w) / np.sum(w)),
                    "centre": float(np.median(v_mid[w < 0.15])) if np.any(w < 0.15) else np.nan,
                    "wall": float(np.median(v_mid[w > w.max() - 0.15])),
                    "r_edge": float((xe - xc) / px_per_m),
                    "classified": float(use.mean()),
                }
            )
        check[part] = PanelCheck(a, xc, float(slope), float(icpt), px_per_m, colours)
    df = pd.DataFrame(rows).sort_values("z", ascending=False)
    # one row per 0.1 m: the best-classified image row within +/- 0.15 m (tick marks hide
    # part of the rows at whole metres), the nearest one on ties
    grid = np.round(np.arange(14.9, -9.85, -0.1), 1)
    zz, cl = df.z.to_numpy(), df.classified.to_numpy()
    pick = []
    for g in grid:
        near = np.where(np.abs(zz - g) <= 0.15)[0]
        if near.size == 0:
            raise SystemExit(f"Fig. {fig.number}: no image row within 0.15 m of z = {g} m")
        best = near[np.lexsort((np.abs(zz[near] - g), -cl[near]))[0]]
        pick.append(int(best))
    out = df.iloc[pick].copy()
    out["z"] = grid
    out.insert(0, "quantity", fig.quantity)
    out.insert(1, "unit", fig.unit)
    out.insert(2, "figure", fig.number)
    check["scale"] = scale
    return out, check


def main() -> None:
    figs = figures()
    frames = []
    checks = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        pages = {f.page: extract_images(f.page, tmp) for f in figs}
        for f in figs:
            df, ck = digitise(f, pages[f.page])
            frames.append(df)
            checks.append((f, df, ck))
            print(
                f"Fig. {f.number} {f.quantity}: {len(df)} rows, px/m {ck['scale']:.2f}, "
                f"edge radius at z = 10 m: {df.r_edge.iloc[np.argmin(np.abs(df.z - 10))]:.2f} m, "
                f"at z = -9 m: {df.r_edge.iloc[np.argmin(np.abs(df.z + 9))]:.2f} m, "
                f"classified {df.classified.mean():.3f}"
            )
    data = pd.concat(frames)
    header = (
        "# Case A of REDUCTOR (Marsigny, Mirgaux & Patisson 2025, Metals 15:862, CC BY 4.0),\n"
        "# Figs. 3, 4, 6-9, digitised by scripts/digitize_reductor_maps.py (do not edit).\n"
        "# z: height above the bottom of the upper cylinder (bustle at 4.9 m, top 15 m, solid\n"
        "# outlet -10 m). mean: cross-section (area-weighted) average of the colour-band\n"
        "# midpoints; mean_lo/mean_hi: same with the band edges; centre: r < 0.15 m; wall:\n"
        "# outer 0.15 m; r_edge: detected radius of the coloured region [m]; classified:\n"
        "# fraction of panel pixels matched to a colour band. T in degC, mass fractions of\n"
        "# the solid, H2 mole fraction of the gas.\n"
    )
    cols = ["quantity", "unit", "figure", "z", "mean", "mean_lo", "mean_hi", "centre", "wall"]
    cols += ["r_edge", "classified"]
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write(header)
        data[cols].to_csv(fh, index=False, float_format="%.5g")
    print(f"wrote {OUT.relative_to(ROOT)}")

    # check figure: reconstructed band map next to the original panel, per quantity
    fig, axes = plt.subplots(2, len(checks), figsize=(2.2 * len(checks), 9), dpi=120)
    for j, (f, _, ck) in enumerate(checks):
        pc = ck["upper"]
        assert isinstance(pc, PanelCheck)
        x1 = int(pc.x_axis + 2.7 * pc.px_per_m)
        y0 = int((15.2 - pc.intercept) / pc.slope)
        y1 = int((-2.3 - pc.intercept) / pc.slope)
        crop = pc.image[y0:y1, int(pc.x_axis) - 5 : x1]
        k, dist = classify(crop, pc.colours)
        rec = np.where(dist < 60, pc.colours[k].transpose(2, 0, 1), 255.0)
        axes[0, j].imshow(crop.astype(np.uint8))
        axes[0, j].set_title(f"Fig. {f.number}\noriginal", fontsize=8)
        axes[1, j].imshow(np.clip(rec.transpose(1, 2, 0), 0, 255).astype(np.uint8))
        axes[1, j].set_title("classified", fontsize=8)
        for ax in axes[:, j]:
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle("REDUCTOR case A maps (z = 15 to -2.3 m): original vs colour-band classification")
    fig.tight_layout()
    fig.savefig(FIG)
    print(f"wrote {FIG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
