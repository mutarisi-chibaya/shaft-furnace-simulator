"""Headless test of the Streamlit app (Milestone 7): it runs, shows the KPIs and profiles, and
switching the preset loads that case's operating conditions. About 30 s (two solves)."""

from __future__ import annotations

import pytest

from .conftest import ROOT

st_testing = pytest.importorskip("streamlit.testing.v1")


def test_app_runs_case_a_and_switches_preset() -> None:
    at = st_testing.AppTest.from_file(str(ROOT / "app" / "streamlit_app.py"), default_timeout=300)
    at.run()
    assert not at.exception, at.exception
    labels = [m.label for m in at.metric]
    assert labels == [
        "Metallisation (%)",
        "H₂ (kg/t DRI)",
        "Top gas (°C)",
        "DRI outlet (°C)",
        "Electricity (GJ/t DRI)",
        "CO₂ (t/t DRI)",
    ]
    # case A on 50 cells: the regression value of tests/test_furnace.py (94.56 %)
    assert at.metric[0].value == "94.6"
    assert len(at.get("plotly_chart")) == 4
    assert any("REDUCTOR (2D) for case A" in i.value for i in at.info)
    assert not at.error

    at.selectbox(key="preset").set_value("C").run()
    assert not at.exception, at.exception
    assert at.slider(key="n2").value == 32.0
    assert at.slider(key="bustle_T").value == 950.0
    assert any("REDUCTOR (2D) for case C" in i.value for i in at.info)


def test_try_this_buttons_set_controls_and_explain() -> None:
    """Layer 2: a guided-scenario button sets the controls and shows what to look for; the
    note disappears once a control is changed by hand. The glossary is on the page."""
    at = st_testing.AppTest.from_file(str(ROOT / "app" / "streamlit_app.py"), default_timeout=300)
    at.run()
    at.button(key="try_cooler").click().run()
    assert not at.exception, at.exception
    assert at.slider(key="bustle_T").value == 840.0
    assert any("What to look for" in i.value for i in at.info)
    at.slider(key="bustle_T").set_value(850.0).run()
    assert not any("What to look for" in i.value for i in at.info)
    at.button(key="try_clean").click().run()
    assert at.radio(key="grid").value == "france"
    assert any("Glossary" in e.label for e in at.expander)
