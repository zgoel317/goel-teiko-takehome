"""Shared dashboard presentation helpers: palette, plot config, KPI tiles,
dataframe styling, and the plain-language findings summary."""
import pandas as pd
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


def format_stats_table(comparison):
    """Paper-quality per-population display table from a compare_responders result."""
    df = comparison
    return pd.DataFrame({
        "Population": df["population"].values,
        "Δ resp−non (pp)": df["mean_diff"].round(2).values,
        "95% CI (pp)": [f"[{lo:.2f}, {hi:.2f}]"
                        for lo, hi in zip(df["boot_ci_low"], df["boot_ci_high"])],
        "Mixed p": df["mixed_p"].round(4).values,
        "FDR q": df["q_mixed"].round(4).values,
        "MWU q": df["q_mw"].round(4).values,
        "Effect size": df["rank_biserial"].round(3).values,
        "Significant": df["significant_primary"].map({True: "yes", False: "no"}).values,
        "Concordant": df["concordant"].map({True: "yes", False: "no"}).values,
    })


def style_stats_table(display_df):
    """Highlight rows that are significant (Significant == 'yes')."""
    def _highlight(row):
        color = "background-color: #FFF3CD" if row.get("Significant") == "yes" else ""
        return [color] * len(row)
    return display_df.style.apply(_highlight, axis=1)


def subset_coverage(metadata):
    """Baseline (t=0) PBMC subset sizes per condition x treatment, flagging the
    degenerate subsets (no responder/non-responder label)."""
    base = metadata[(metadata["sample_type"] == "PBMC") & (metadata["time"] == 0)]
    rows = []
    for (condition, treatment), grp in base.groupby(["condition", "treatment"]):
        labelled = int(grp["response"].isin(["yes", "no"]).sum())
        rows.append({
            "condition": condition,
            "treatment": treatment,
            "samples": len(grp),
            "responder_labelled": labelled,
            "degenerate": labelled == 0,
        })
    return (pd.DataFrame(rows)
            .sort_values(["condition", "treatment"])
            .reset_index(drop=True))


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
