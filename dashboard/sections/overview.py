"""Overview tab (Part 2): dataset KPIs + relative-frequency table + composition."""
import plotly.express as px
import streamlit as st

from dashboard.data import load_dataset_summary, load_frequencies_annotated
from dashboard.components import kpi_row

_ALL = "All"


def _select(label, values):
    return st.selectbox(label, [_ALL] + sorted(values))


def render():
    st.subheader("Data overview — relative cell-population frequencies")
    summary = load_dataset_summary()
    kpi_row([("Samples", summary["samples"]),
             ("Subjects", summary["subjects"]),
             ("Projects", summary["projects"])])

    data = load_frequencies_annotated()
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        condition = _select("Condition", data["condition"].dropna().unique())
    with c2:
        treatment = _select("Treatment", data["treatment"].dropna().unique())
    with c3:
        sample_type = _select("Sample type", data["sample_type"].dropna().unique())
    with c4:
        populations = st.multiselect("Populations",
                                     sorted(data["population"].unique()))

    view = data
    if condition != _ALL:
        view = view[view["condition"] == condition]
    if treatment != _ALL:
        view = view[view["treatment"] == treatment]
    if sample_type != _ALL:
        view = view[view["sample_type"] == sample_type]
    if populations:
        view = view[view["population"].isin(populations)]

    if view.empty:
        st.info("No samples match these filters.")
        return

    mean_comp = (view.groupby("population", as_index=False)["percentage"].mean()
                 .sort_values("percentage", ascending=False))
    fig = px.bar(mean_comp, x="population", y="percentage",
                 title="Mean relative frequency by population (%)",
                 labels={"percentage": "Mean %", "population": "Population"})
    st.plotly_chart(fig, use_container_width=True)

    st.dataframe(
        view[["sample", "total_count", "population", "count", "percentage"]]
        .reset_index(drop=True),
        use_container_width=True, hide_index=True,
    )
