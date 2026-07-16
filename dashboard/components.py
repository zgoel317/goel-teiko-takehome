"""Shared dashboard presentation helpers: palette, plot config, KPI tiles,
dataframe styling, and the plain-language findings summary."""
import plotly.express as px
import streamlit as st

# Okabe-Ito colour-blind-safe categorical palette (5 populations / 2 groups).
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]


def configure_plotly():
    """Apply the shared palette + a clean template to every Plotly Express chart."""
    px.defaults.color_discrete_sequence = PALETTE
    px.defaults.template = "plotly_white"


def kpi_row(items):
    """Render a row of st.metric tiles. items: list of (label, value)."""
    cols = st.columns(len(items))
    for col, (label, value) in zip(cols, items):
        col.metric(label, value)


def style_comparison(df):
    """Style the Part 3 statistics table: highlight significant_primary rows."""
    def _highlight(row):
        color = "background-color: #FFF3CD" if row.get("significant_primary") else ""
        return [color] * len(row)
    return df.style.apply(_highlight, axis=1).format(precision=3)


def finding_summary(comparison, alpha=0.05):
    """Plain-language markdown summary of the responder comparison result."""
    primary = comparison[comparison["significant_primary"]]
    if primary.empty:
        return ("**No cell population** shows a statistically significant "
                "difference in relative frequency between responders and "
                f"non-responders (primary (mixed-model) FDR-adjusted q ≥ {alpha:.2f}).")
    lines = []
    for _, r in primary.iterrows():
        direction = "higher" if r["mean_diff"] > 0 else "lower"
        line = (f"**{r['population']}** is {direction} in responders "
                f"({r['mean_diff']:+.2f} pp; mixed-model q={r['q_mixed']:.3f})")
        if r["concordant"]:
            line += (" — **confirmed** by both the mixed model and the "
                     "non-parametric test.")
        else:
            line += (f" — but **not concordant** with the conservative "
                     f"per-subject test (q={r['q_mw']:.3f}); treat as a "
                     f"**candidate**, not confirmed.")
        lines.append(line)
    return "  \n".join(lines)
