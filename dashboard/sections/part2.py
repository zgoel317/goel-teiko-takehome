"""Part 2 · Cell-Type Frequencies: the per-sample relative-frequency summary."""
import plotly.express as px
import streamlit as st

from dashboard.data import load_dataset_summary, load_frequencies
from dashboard.components import kpi_row


def render():
    st.subheader("Part 2 · Relative cell-population frequencies")
    st.write("For each sample, the relative frequency of each immune cell "
             "population as a percentage of that sample's total cell count.")

    summary = load_dataset_summary()
    kpi_row([("Samples", summary["samples"]),
             ("Subjects", summary["subjects"]),
             ("Projects", summary["projects"])])

    freq = load_frequencies()
    mean_comp = (freq.groupby("population", as_index=False)["percentage"].mean()
                 .sort_values("percentage", ascending=False))
    fig = px.bar(mean_comp, x="population", y="percentage",
                 title="Mean relative frequency by population (all samples)",
                 labels={"percentage": "Mean %", "population": "Population"})
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("**Per-sample frequency table**")
    st.dataframe(freq, use_container_width=True, hide_index=True)
