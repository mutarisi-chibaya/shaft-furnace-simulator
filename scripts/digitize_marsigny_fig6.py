"""Digitise reduction curves from Marsigny et al. (2025), Metals 15:1034 (Figures 6 and 9).

The figures are raster images inside the PDF, so the curves are recovered from the pixels,
reproducibly, instead of by hand:

1. The image is extracted from the PDF page with ``pdfimages`` (poppler).
2. Axes: the plot frame (long dark lines) is detected; its edges are the axis limits read
   from the tick labels (Fig. 6: 0-150 min, Fig. 9: 0-40 min; conversion 0-1). The time
   calibration is cross-checked against the detected tick marks (must agree within 1 % of
   the axis range).
3. Each curve has its own matplotlib colour (given in the legend). Pixels are classified by
   nearest colour within a distance threshold; the legend box is excluded by masking a region
   no curve enters.
4. Experimental curves are solid, simulated ones dotted, in the same colour. They are
   separated by connected-component size: solid lines form long components, dots small ones.
   Only the EXPERIMENTAL curves are exported; the dotted ones are drawn on the check image but
   their extraction is unreliable where curves are close (dots and anti-aliased edges mix).
5. For each pixel column the median row of the solid pixels gives one point. Points at
   X > 0.985 are dropped (curves merge with the top frame line there).

Outputs:
    tests/data/pellet_reduction_marsigny2025_fig6.csv  (temperature series, base gas)
    tests/data/pellet_reduction_marsigny2025_fig9.csv  (H2 dilution series, 900 degC)
    figures/digitization_check_fig6.png, figures/digitization_check_fig9.png

Usage:
    python scripts/digitize_marsigny_fig6.py
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "references" / "marsigny2025_metals15_1034.pdf"
COLOUR_TOL = 60.0  # max RGB distance to accept a pixel as a curve colour
SOLID_MIN_EXTENT = 40  # px; components with a larger bounding box are solid lines

NAVY, BLUE, CYAN = (0, 0, 128), (0, 0, 255), (0, 255, 255)
LIGHTGREEN, DARKGREEN, ORANGE, RED = (144, 238, 144), (0, 100, 0), (255, 165, 0), (255, 0, 0)


@dataclass(frozen=True)
class FigureSpec:
    name: str
    page: int
    image_index: int  # index among the images on that page
    t_max_min: float
    x_ticks_inner: tuple[float, ...]  # tick labels strictly inside the frame
    curves: dict[float, tuple[int, int, int]]  # curve label -> RGB
    label_column: str
    legend_mask: tuple[float, float]  # exclude t > a [min] and X < b
    output: str


FIGURES = (
    FigureSpec(
        name="fig6",
        page=10,
        image_index=0,
        t_max_min=150.0,
        x_ticks_inner=(20, 40, 60, 80, 100, 120, 140),
        curves={
            600: NAVY,
            650: BLUE,
            700: CYAN,
            750: LIGHTGREEN,
            800: DARKGREEN,
            850: ORANGE,
            900: RED,
        },
        label_column="temperature_C",
        legend_mask=(110.0, 0.45),
        output="pellet_reduction_marsigny2025_fig6.csv",
    ),
    FigureSpec(
        name="fig9",
        page=13,
        image_index=0,
        t_max_min=40.0,
        x_ticks_inner=(5, 10, 15, 20, 25, 30, 35),
        curves={0.20: NAVY, 0.55: BLUE, 0.80: CYAN, 1.00: LIGHTGREEN},
        label_column="x_H2_in",
        legend_mask=(25.0, 0.30),
        output="pellet_reduction_marsigny2025_fig9.csv",
    ),
)


def extract_image(spec: FigureSpec, tmp: Path) -> np.ndarray:
    exe = shutil.which("pdfimages")
    if exe is None:
        raise SystemExit("pdfimages (poppler) is required")
    subprocess.run(
        [exe, "-png", "-f", str(spec.page), "-l", str(spec.page), str(PDF), str(tmp / "img")],
        check=True,
    )
    files = sorted(tmp.glob("img*.png"))
    img = mpimg.imread(files[spec.image_index])
    arr = np.asarray(img[..., :3], dtype=float)
    return arr * 255.0 if arr.max() <= 1.0 else arr


def find_frame(img: np.ndarray) -> tuple[int, int, int, int]:
    """Return (left, right, top, bottom) pixel indices of the plot frame."""
    dark = img.sum(axis=2) < 200
    cols = np.where(dark.sum(axis=0) > 0.6 * img.shape[0])[0]
    rows = np.where(dark.sum(axis=1) > 0.6 * img.shape[1])[0]
    return int(cols.min()), int(cols.max()), int(rows.min()), int(rows.max())


def tick_columns(img: np.ndarray, bottom: int, left: int, right: int) -> np.ndarray:
    """Columns of x tick marks just below the bottom frame line (matplotlib: outward)."""
    dark = img.sum(axis=2) < 200
    band = dark[bottom + 3 : bottom + 9, left + 3 : right - 3]
    cols = np.where(band.sum(axis=0) >= 4)[0] + left + 3
    groups = np.split(cols, np.where(np.diff(cols) > 2)[0] + 1)
    return np.array([g.mean() for g in groups if len(g)])


def digitize(spec: FigureSpec) -> None:
    with tempfile.TemporaryDirectory() as td:
        img = extract_image(spec, Path(td))
    left, right, top, bottom = find_frame(img)

    def to_t(col: np.ndarray) -> np.ndarray:
        return (col - left) / (right - left) * spec.t_max_min

    def to_x(row: np.ndarray) -> np.ndarray:
        return (bottom - row) / (bottom - top)

    ticks = tick_columns(img, bottom, left, right)
    inner = ticks[(ticks > left + 10) & (ticks < right - 10)]
    expected = np.array(spec.x_ticks_inner, dtype=float)
    if len(inner) != len(expected):
        raise SystemExit(f"{spec.name}: expected {len(expected)} inner ticks, found {len(inner)}")
    err = float(np.max(np.abs(to_t(inner) - expected)))
    print(f"{spec.name}: x-calibration max tick error {err:.3f} min ({len(inner)} ticks)")
    if err > 0.01 * spec.t_max_min:
        raise SystemExit("axis calibration failed the tick check")

    rr, cc = np.mgrid[0 : img.shape[0], 0 : img.shape[1]]
    inside = (cc > left + 2) & (cc < right - 2) & (rr > top + 2) & (rr < bottom - 2)
    legend = (to_t(cc) > spec.legend_mask[0]) & (to_x(rr) < spec.legend_mask[1])
    dists = {
        k: np.sqrt(((img - np.array(c, float)) ** 2).sum(axis=2)) for k, c in spec.curves.items()
    }
    nearest = np.min(np.stack(list(dists.values())), axis=0)

    rows_out: list[list[object]] = []
    fig, ax = plt.subplots(figsize=(8, 7), dpi=110)
    extent = (
        float(to_t(np.array(0.0))),
        float(to_t(np.array(float(img.shape[1])))),
        float(to_x(np.array(float(img.shape[0])))),
        float(to_x(np.array(0.0))),
    )
    ax.imshow(img.astype(np.uint8), extent=extent, alpha=0.35, aspect="auto")
    for label, rgb in spec.curves.items():
        mask = (dists[label] < COLOUR_TOL) & (dists[label] <= nearest) & inside & ~legend
        lab, _ = ndimage.label(mask, structure=np.ones((3, 3)))
        solid = np.zeros_like(mask)
        dotted = np.zeros_like(mask)
        for i, s in enumerate(ndimage.find_objects(lab), start=1):
            h, w = s[0].stop - s[0].start, s[1].stop - s[1].start
            target = solid if max(h, w) >= SOLID_MIN_EXTENT else dotted
            target[s] |= lab[s] == i
        for kind, m in (("experimental", solid), ("simulated", dotted)):
            ys, xs = np.nonzero(m)
            pts = [
                (
                    float(to_t(np.array(float(c)))),
                    float(to_x(np.array(float(np.median(ys[xs == c]))))),
                )
                for c in np.unique(xs)
            ]
            pts_a = np.array(pts) if pts else np.zeros((0, 2))
            if len(pts_a):
                pts_a = pts_a[pts_a[:, 1] <= 0.985]
            if kind == "experimental":
                rows_out += [[label, kind, round(t, 3), round(x, 4)] for t, x in pts_a]
            if len(pts_a):
                ax.plot(
                    pts_a[:, 0],
                    pts_a[:, 1],
                    "-" if kind == "experimental" else ":",
                    color=np.array(rgb) / 255,
                    lw=1,
                    label=f"{label} {kind} ({len(pts_a)} pts)",
                )
    ax.set_xlim(0, spec.t_max_min)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Conversion")
    ax.set_title(
        f"Digitisation check: Marsigny et al. 2025, {spec.name} (image faded behind)", fontsize=10
    )
    ax.legend(fontsize=6, loc="lower right", ncol=2)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / f"digitization_check_{spec.name}.png")
    plt.close(fig)

    out = ROOT / "tests" / "data" / spec.output
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([spec.label_column, "curve", "time_min", "conversion"])
        w.writerows(rows_out)
    print(
        f"{spec.name}: frame px {left}-{right} x {top}-{bottom}; wrote {out.name} "
        f"({len(rows_out)} points)"
    )


def main() -> None:
    for spec in FIGURES:
        digitize(spec)


if __name__ == "__main__":
    main()
