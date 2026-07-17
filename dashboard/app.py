"""Loblaw Bio immune-cell analysis dashboard."""
import os
import sys

# streamlit run dashboard/app.py puts dashboard/ on sys.path, not the repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st

st.set_page_config(page_title="Loblaw Bio — Immune Cell Analysis",
                   layout="wide")

from dashboard.components import configure_plotly
from dashboard.sections import part2, part3, part4

configure_plotly()

st.sidebar.title("Loblaw Bio")
st.sidebar.caption("Immune cell population analysis — miraclib trial "
                   "(source: cell-count.csv)")

st.title("Immune Cell Population Analysis")

tab2, tab3, tab4 = st.tabs([
    "Part 2 · Cell-Type Frequencies",
    "Part 3 · Statistical Analysis",
    "Part 4 · Subset Analysis",
])
with tab2:
    part2.render()
with tab3:
    part3.render()
with tab4:
    part4.render()
