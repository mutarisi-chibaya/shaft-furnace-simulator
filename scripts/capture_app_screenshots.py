"""Screenshots of the running Streamlit app for docs/tutorial.md.

Usage (app running locally):
    streamlit run app/streamlit_app.py --server.port 8599 --server.headless true
    python scripts/capture_app_screenshots.py [--url http://localhost:8599]

Needs `pip install playwright` and a local Chrome (used through Playwright's "chrome"
channel; no browser download). Writes figures/tutorial/*.png in the light theme.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "tutorial"
TIMEOUT_MS = 240_000  # a cold 100-cell run takes ~40 s; allow for slow machines


def wait_for_results(page: Page) -> None:
    """Wait until the spinner is gone and the KPI tiles are shown."""
    page.wait_for_selector('[data-testid="stMetric"]', timeout=TIMEOUT_MS)
    page.wait_for_function(
        "() => !document.querySelector('[data-testid=\"stSpinner\"]')", timeout=TIMEOUT_MS
    )
    page.wait_for_timeout(2500)  # let Plotly finish drawing


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://localhost:8599")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome", headless=True)
        ctx = browser.new_context(
            viewport={"width": 1600, "height": 1000}, color_scheme="light", device_scale_factor=1
        )
        page = ctx.new_page()
        page.goto(a.url)
        wait_for_results(page)

        page.screenshot(path=str(OUT / "01_overview.png"))
        page.locator('[data-testid="stSidebar"]').screenshot(path=str(OUT / "02_sidebar.png"))
        page.locator(".st-key-kpis").screenshot(path=str(OUT / "03_kpis.png"))
        page.locator(".st-key-profiles").screenshot(path=str(OUT / "04_profiles.png"))
        co2 = page.get_by_text("CO₂ per tonne of DRI and the electricity grid")
        co2.scroll_into_view_if_needed()
        page.locator(".st-key-co2").screenshot(path=str(OUT / "05_co2.png"))
        page.get_by_text("Assumptions and limits", exact=True).click()
        page.wait_for_timeout(800)
        page.locator('[data-testid="stExpander"]').nth(0).screenshot(
            path=str(OUT / "06_limits.png")
        )
        browser.close()
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
