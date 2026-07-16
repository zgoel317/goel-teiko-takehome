"""Loblaw Bio immune-cell analysis dashboard."""
import os
import sys

# `streamlit run dashboard/app.py` puts dashboard/ on sys.path, not the repo
# root, so project imports (load_data, analysis.*, dashboard.*) would fail both
# locally (`make dashboard`) and on Streamlit Community Cloud. Put the repo root
# (this file's grandparent) on the path before any project import.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st

st.set_page_config(page_title="Loblaw Bio — Immune Cell Analysis",
                   layout="wide")

from dashboard.components import configure_plotly
from dashboard.sections import overview, response, subset

configure_plotly()

st.sidebar.title("Loblaw Bio")
st.sidebar.caption("Immune cell population analysis — miraclib trial")
with st.sidebar.expander("About & methodology"):
    st.markdown(
        "Relative cell-population frequencies from `cell-count.csv`. "
        "Response analysis compares miraclib responders vs non-responders "
        "(melanoma, PBMC) using a linear mixed-effects model (primary), a "
        "Mann-Whitney U cross-check on per-subject means, Benjamini-Hochberg "
        "FDR correction, and bootstrap effect-size CIs."
    )

st.title("Immune Cell Population Analysis")

tab_overview, tab_response, tab_subset = st.tabs(
    ["Overview", "Response Analysis", "Subset Explorer"])
with tab_overview:
    overview.render()
with tab_response:
    response.render()
with tab_subset:
    subset.render()
