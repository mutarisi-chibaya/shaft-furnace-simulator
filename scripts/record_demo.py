"""Record the demo video of docs/demo_script.md by driving the running app with Playwright.

Usage (app running locally):
    streamlit run app/streamlit_app.py --server.port 8599 --server.headless true
    python scripts/record_demo.py [--url http://localhost:8599]

Needs `pip install playwright`, a local Chrome (Playwright "chrome" channel) and ffmpeg.

- Captions: the narration of each scene is shown sentence by sentence at the bottom.
- Timing: if build/voiceover/scene_XX.mp3 exist (scripts/make_voiceover.py), each scene lasts
  as long as its audio and the audio is mixed in; otherwise sentences are timed at 145 words
  per minute and the video is silent with captions.
- A rehearsal pass without recording runs every scenario first, so the results come from the
  app's cache during the recorded pass.
Output: build/demo/demo.mp4 (build/ is not committed).
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path

from playwright.sync_api import Locator, Page, sync_playwright
from playwright.sync_api import TimeoutError as PlaywrightTimeout

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build" / "demo"
VOICE = ROOT / "build" / "voiceover"
sys.path.insert(0, str(ROOT / "scripts"))
import make_voiceover  # noqa: E402

WPM = 145.0
VIEW = {"width": 1600, "height": 1000}
TIMEOUT_MS = 240_000

OVERLAY_JS = """
() => {
  if (document.getElementById('demo-caption')) return;
  const c = document.createElement('div');
  c.id = 'demo-caption';
  Object.assign(c.style, {position: 'fixed', left: '50%', bottom: '28px',
    transform: 'translateX(-50%)', maxWidth: '1150px', padding: '12px 22px',
    background: 'rgba(11,11,11,0.82)', color: '#fcfcfb', font: '500 21px/1.4 system-ui, sans-serif',
    borderRadius: '10px', zIndex: 2147483001, textAlign: 'center', display: 'none'});
  document.body.appendChild(c);
  const k = document.createElement('div');
  k.id = 'demo-card';
  Object.assign(k.style, {position: 'fixed', inset: '0', background: '#fcfcfb', color: '#0b0b0b',
    zIndex: 2147483000, display: 'none', flexDirection: 'column', justifyContent: 'center',
    padding: '0 160px', font: '400 26px/1.5 system-ui, sans-serif'});
  document.body.appendChild(k);
  const s = document.createElement('style');
  s.textContent = '.demo-hl { outline: 3px solid #eb6834 !important; outline-offset: 6px;'
    + ' border-radius: 8px; transition: outline-color .3s; }';
  document.head.appendChild(s);
}
"""


class Demo:
    """Browser actions used by the scenes."""

    def __init__(self, page: Page) -> None:
        self.page = page

    def wait_results(self) -> None:
        """Wait for the app's rerun to start (status widget) and finish, then for the KPIs."""
        status = '[data-testid="stStatusWidget"]'
        try:
            self.page.wait_for_selector(status, state="visible", timeout=2500)
        except PlaywrightTimeout:
            pass  # rerun already finished (cached) or nothing to rerun
        self.page.wait_for_selector(status, state="hidden", timeout=TIMEOUT_MS)
        self.page.wait_for_selector('[data-testid="stMetric"]', timeout=TIMEOUT_MS)
        self.page.wait_for_function(
            "() => !document.querySelector('[data-testid=\"stSpinner\"]')", timeout=TIMEOUT_MS
        )
        self.page.wait_for_timeout(1200)
        self.page.evaluate(OVERLAY_JS)

    def caption(self, text: str | None) -> None:
        self.page.evaluate(
            "t => { const c = document.getElementById('demo-caption');"
            " c.textContent = t || ''; c.style.display = t ? 'block' : 'none'; }",
            text,
        )

    def card(self, html: str | None) -> None:
        self.page.evaluate(
            "h => { const k = document.getElementById('demo-card');"
            " k.innerHTML = h || ''; k.style.display = h ? 'flex' : 'none'; }",
            html,
        )

    def highlight(self, selector: str | None, nth: int = 0) -> None:
        self.page.evaluate(
            "([s, n]) => { document.querySelectorAll('.demo-hl').forEach("
            "e => e.classList.remove('demo-hl'));"
            " if (!s) return; const e = document.querySelectorAll(s)[n];"
            " if (e) { e.classList.add('demo-hl'); e.scrollIntoView({block: 'center',"
            " behavior: 'smooth'}); } }",
            [selector, nth],
        )

    def scroll_top(self) -> None:
        self.page.evaluate(
            '() => document.querySelectorAll(\'section, [data-testid="stMain"],'
            ' [data-testid="stAppViewContainer"]\').forEach(e => e.scrollTo({top: 0,'
            " behavior: 'smooth'}))"
        )

    def preset(self, label_start: str) -> None:
        self.page.locator('[data-testid="stSidebar"] [data-testid="stSelectbox"]').click()
        self.page.get_by_role("option").filter(has_text=label_start).first.click()
        self.wait_results()

    def radio(self, text_start: str) -> None:
        self.page.locator('[data-testid="stSidebar"] label').filter(
            has_text=text_start
        ).first.click()
        self.wait_results()

    def slider(self, label: str) -> Locator:
        return self.page.locator('[data-testid="stSlider"]').filter(has_text=label).first

    def highlight_slider(self, label: str) -> None:
        self.highlight(None)
        el = self.slider(label)
        el.evaluate("e => { e.classList.add('demo-hl'); e.scrollIntoView({block: 'center'}); }")

    def slider_to(self, label: str, value: float, lo: float, hi: float, step: float) -> None:
        """Drag a sidebar slider to `value` (one rerun), then correct with arrow keys.

        Streamlit's slider (React Aria): a horizontal track holding a thumb that wraps a
        native range input, whose value is the slider value.
        """
        box = self.slider(label)
        box.scroll_into_view_if_needed()  # measure where the drag will happen
        self.page.wait_for_timeout(300)
        rng = box.locator('input[type="range"]')
        thumb = rng.locator("xpath=../..")
        track = box.locator('[role="group"] > div[data-orientation="horizontal"]').first
        tr, tb = track.bounding_box(), thumb.bounding_box()
        assert tr is not None and tb is not None
        x = tr["x"] + tr["width"] * (value - lo) / (hi - lo)
        y = tb["y"] + tb["height"] / 2
        self.page.mouse.move(tb["x"] + tb["width"] / 2, y)
        self.page.mouse.down()
        self.page.mouse.move(x, y, steps=25)
        self.page.mouse.up()
        self.wait_results()
        for _ in range(3):  # correct a miss: all key steps at once, then one rerun
            now = float(rng.evaluate("e => e.value"))
            n = round(abs(value - now) / step)
            if n == 0:
                break
            rng.focus()
            for _ in range(n):
                rng.press("ArrowRight" if now < value else "ArrowLeft", delay=20)
            self.wait_results()


def sentences(paragraphs: list[str]) -> list[str]:
    text = " ".join(paragraphs)
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]


def audio_seconds(path: Path) -> float | None:
    probe = shutil.which("ffprobe")
    if not path.exists() or probe is None:
        return None
    out = subprocess.run(
        [probe, "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0", str(path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return float(out.stdout.strip())


def timed(sents: list[str], total: float | None) -> list[tuple[str, float]]:
    """Seconds per sentence: proportional to words, scaled to the scene audio if given."""
    words = [len(s.split()) for s in sents]
    base = [60.0 * w / WPM for w in words]
    scale = (total / sum(base)) if total else 1.0
    return [(s, b * scale) for s, b in zip(sents, base, strict=True)]


def diagram_html() -> str:
    """The furnace schematic (figures/furnace_schematic.svg) as a full-height card."""
    svg = (ROOT / "figures" / "furnace_schematic.svg").read_text(encoding="utf-8")
    svg = svg.replace("<svg ", '<svg style="height:66vh;width:auto" ', 1)
    return (
        "<div style='display:flex;flex-direction:column;align-items:center'>"
        "<div style='font:700 44px/1.15 system-ui,sans-serif;margin-bottom:12px'>"
        f"How the furnace works</div>{svg}</div>"
    )


def card_html(title: str, lines: list[str]) -> str:
    items = "".join(f"<li style='margin:10px 0'>{x}</li>" for x in lines)
    return (
        f"<div style='font:700 52px/1.15 system-ui,sans-serif;margin-bottom:28px'>{title}</div>"
        f"<ul style='color:#52514e;padding-left:28px'>{items}</ul>"
    )


#: Per scene: actions keyed by sentence index (run when that sentence starts).
Actions = dict[int, Callable[[Demo], None]]


def scene_actions() -> list[Actions]:
    def bustle(v: float) -> Callable[[Demo], None]:
        return lambda d: d.slider_to("Bustle gas temperature", v, 815.0, 950.0, 5.0)

    return [
        {  # 1 title
            0: lambda d: d.card(
                card_html(
                    "Hydrogen shaft furnace simulator",
                    [
                        "A physics-based model of green iron making",
                        "1D furnace model · cited data · tested · open source",
                    ],
                )
            ),
        },
        {  # 2 why it matters: the furnace first (sentences 1-2), then the projects
            0: lambda d: d.card(diagram_html()),
            2: lambda d: d.card(
                card_html(
                    "Why it matters",
                    [
                        "Hydrogen replaces coal as the reducing agent: the by-product is water",
                        "ArcelorMittal South Africa: study of hydrogen DRI at Saldanha",
                        "HyIron, Namibia: hydrogen DRI plant under construction",
                        "Scale-up needs fast, honest models",
                    ],
                )
            ),
        },
        {  # 3 what it is
            0: lambda d: (d.card(None), d.scroll_top(), d.highlight(".st-key-overview")),
            1: lambda d: d.highlight(".st-key-profiles"),
            2: lambda d: (
                d.page.get_by_text("Assumptions and limits", exact=True).click(),
                d.highlight('[data-testid="stExpander"]', 0),
            ),
        },
        {  # 4 reference case
            0: lambda d: (
                d.page.get_by_text("Assumptions and limits", exact=True).click(),
                d.scroll_top(),
                d.highlight(".st-key-kpis"),
            ),
            2: lambda d: d.highlight(".st-key-reductor"),
            3: lambda d: d.highlight(".st-key-profiles"),
        },
        {  # 5 bustle temperature
            0: lambda d: (d.scroll_top(), d.highlight_slider("Bustle gas temperature")),
            1: lambda d: (bustle(840.0)(d), d.highlight('[data-testid="stMetric"]', 0)),
            2: lambda d: (bustle(920.0)(d), d.highlight('[data-testid="stMetric"]', 0)),
        },
        {  # 6 case B
            0: lambda d: (
                d.highlight(None),
                d.preset("B:"),
                d.highlight('[data-testid="stMetric"]', 3),
            ),
            2: lambda d: d.highlight('[data-testid="stMetric"]', 0),
            3: lambda d: d.highlight(".st-key-reductor"),
        },
        {  # 7 South Africa and the grid
            0: lambda d: (
                d.highlight(None),
                d.preset("A:"),
                d.highlight(".st-key-co2"),
            ),
            1: lambda d: d.highlight('[data-testid="stMetric"]', 4),
            2: lambda d: d.highlight('[data-testid="stMetric"]', 5),
            4: lambda d: (d.radio("France"), d.highlight(".st-key-co2")),
            5: lambda d: (
                d.radio("Custom"),
                d.slider_to("Custom intensity", 130.0, 0.0, 1000.0, 10.0),
                d.highlight(".st-key-co2"),
            ),
        },
        {  # 8 limits
            0: lambda d: (
                d.highlight(None),
                d.radio("South Africa"),
                d.page.get_by_text("Assumptions and limits", exact=True).click(),
                d.highlight('[data-testid="stExpander"]', 0),
            ),
        },
        {  # 9 close
            0: lambda d: (
                d.highlight(None),
                d.card(
                    card_html(
                        "Hydrogen shaft furnace simulator",
                        [
                            "Code, data sources and validation: the GitHub repository",
                            "Run it: streamlit run app/streamlit_app.py",
                            "Plain-language guide: docs/tutorial.md",
                            "Built with Python · NumPy · SciPy · Streamlit",
                        ],
                    )
                ),
            ),
        },
    ]


def run(
    url: str, record: bool, plan: list[tuple[str, list[tuple[str, float]]]]
) -> tuple[Path | None, float]:
    """Play the scenes; return (raw video, seconds of page loading before scene 1)."""
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome", headless=True)
        opts: dict[str, object] = {"viewport": VIEW, "color_scheme": "light"}
        if record:
            OUT.mkdir(parents=True, exist_ok=True)
            opts["record_video_dir"] = str(OUT / "raw")
            opts["record_video_size"] = VIEW
        ctx = browser.new_context(**opts)  # type: ignore[arg-type]
        page = ctx.new_page()
        t_video = time.perf_counter()  # the recording starts with the page
        page.goto(url)
        d = Demo(page)
        d.wait_results()
        lead = time.perf_counter() - t_video
        for (title, sents), actions in zip(plan, scene_actions(), strict=True):
            print(f"  {title}", flush=True)
            for i, (sentence, secs) in enumerate(sents):
                t0 = time.perf_counter()
                d.caption(sentence)
                if i in actions:
                    actions[i](d)
                if record:
                    rest = secs - (time.perf_counter() - t0)
                    if rest > 0:
                        page.wait_for_timeout(int(1000 * rest))
            if record:
                page.wait_for_timeout(500)
        d.caption(None)
        video = page.video
        ctx.close()
        browser.close()
        return (Path(video.path()) if (record and video) else None), lead


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://localhost:8599")
    ap.add_argument("--no-rehearsal", action="store_true")
    a = ap.parse_args()
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise SystemExit("ffmpeg not found")
    sc = make_voiceover.scenes(make_voiceover.SCRIPT)
    audio = [VOICE / f"scene_{i:02d}.mp3" for i in range(1, len(sc) + 1)]
    durations = [audio_seconds(f) for f in audio]
    with_voice = all(x is not None for x in durations)
    plan = [
        (title, timed(sentences(ps), dur if with_voice else None))
        for (title, ps), dur in zip(sc, durations, strict=True)
    ]
    total = sum(s for _, ss in plan for _, s in ss)
    print(
        f"{len(plan)} scenes, {total / 60:.1f} min, "
        f"{'with voice-over' if with_voice else 'captions only'}"
    )
    if not a.no_rehearsal:
        print("rehearsal (fills the app cache):")
        run(a.url, False, plan)
    print("recording:")
    raw, lead = run(a.url, True, plan)
    assert raw is not None
    final = OUT / "demo.mp4"
    # cut the page loading before scene 1 from the video input only (the voice-over starts at
    # scene 1); input seeking is frame-accurate when re-encoding
    cmd = [ffmpeg, "-y", "-ss", f"{lead:.2f}", "-i", str(raw)]
    if with_voice:
        lst = OUT / "audio_list.txt"
        lst.write_text("".join(f"file '{f.as_posix()}'\n" for f in audio), encoding="utf-8")
        cmd += [
            "-f",
            "concat",
            "-safe",
            "0",
            "-i",
            str(lst),
            "-map",
            "0:v",
            "-map",
            "1:a",
            "-c:a",
            "aac",
            "-b:a",
            "128k",
            "-shortest",
        ]
    cmd += [
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "20",
        "-movflags",
        "+faststart",
        str(final),
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    print(f"wrote {final.relative_to(ROOT)} ({final.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
