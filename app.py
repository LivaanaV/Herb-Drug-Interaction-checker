"""Streamlit entry point for the HDI Safety Checker.

This serves the same static site as index.html (engine.js + data.js), so all
prediction logic stays in engine.js. The scripts are inlined into the page so
Streamlit can render it without a separate web server. Nothing here changes
a number or a rule.

Run locally:   streamlit run app.py
Deploy:        Streamlit Community Cloud, main file = app.py
"""
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

ROOT = Path(__file__).resolve().parent
PAGE_HEIGHT = 1200


@st.cache_data(show_spinner=False)
def build_page() -> str:
    html = (ROOT / "index.html").read_text(encoding="utf-8")
    scripts = (
        ('<script src="data.js"></script>', (ROOT / "data.js").read_text(encoding="utf-8")),
        ('<script src="engine.js"></script>', (ROOT / "engine.js").read_text(encoding="utf-8")),
    )
    for tag, code in scripts:
        if tag not in html:
            raise RuntimeError(f"index.html no longer contains {tag!r}; update app.py")
        if "</script" in code.lower():
            raise RuntimeError(f"{tag} contains a closing script tag and cannot be inlined")
        html = html.replace(tag, "<script>\n" + code + "\n</script>", 1)
    return html


st.set_page_config(page_title="HDI Safety Checker", layout="wide")
components.html(build_page(), height=PAGE_HEIGHT, scrolling=True)
