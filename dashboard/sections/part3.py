"""Part 3 · Statistical Analysis: a fixed responder vs non-responder report on
the required cohort (melanoma + miraclib + PBMC), presented as a paper."""
import streamlit as st

from dashboard.data import load_cohort, load_comparison, load_boxplot_figure
from dashboard.components import finding_summary, format_stats_table, style_stats_table

_CONDITION, _TREATMENT, _SAMPLE_TYPE = "melanoma", "miraclib", "PBMC"


def render():
    st.subheader("Part 3 · Statistical analysis: responders vs non-responders")

    st.markdown("#### Question")
    st.write("Do immune cell population relative frequencies differ between "
             "miraclib **responders** and **non-responders** in melanoma PBMC "
             "samples — and could any population predict response to miraclib?")

    cohort = load_cohort(_CONDITION, _TREATMENT, _SAMPLE_TYPE)
    n_subjects = cohort["subject_id"].nunique()
    counts = cohort.drop_duplicates("subject_id")["response"].value_counts()
    st.markdown("#### Cohort")
    st.write(
        f"Melanoma + miraclib + PBMC: **{n_subjects} subjects** "
        f"({int(counts.get('yes', 0))} responders, "
        f"{int(counts.get('no', 0))} non-responders), each sampled at "
        "timepoints 0, 7, and 14.")

    st.markdown("#### Methods")
    st.markdown(
        "- **Primary — linear mixed-effects model** "
        "(`percentage ~ response + C(time)`, random intercept per subject): "
        "uses all timepoints while respecting the repeated measures per subject "
        "(avoids pseudoreplication). Where the between-subject variance is ~0 the "
        "model reduces to — and is computed as — ordinary least squares.\n"
        "- **Cross-check — Mann-Whitney U** on per-subject mean %: "
        "distribution-free, one independent observation per subject.\n"
        "- **Multiple comparisons — Benjamini-Hochberg FDR** across the five "
        "populations (reported as q).\n"
        "- **Effect size — bootstrap 95% CI** of the responder − non-responder "
        "difference in mean %.\n"
        "- A population is a **candidate** if significant under the primary "
        "model, and **confirmed** only if also significant under the "
        "non-parametric cross-check."
    )

    result = load_comparison(_CONDITION, _TREATMENT, _SAMPLE_TYPE)

    st.markdown("#### Results")
    st.info(finding_summary(result))
    st.plotly_chart(
        load_boxplot_figure(_CONDITION, _TREATMENT, _SAMPLE_TYPE),
        use_container_width=True)
    st.markdown("**Per-population statistics**")
    st.caption("Highlighted rows are significant under the primary mixed-effects "
               "model (FDR q < 0.05).")
    st.dataframe(style_stats_table(format_stats_table(result)),
                 use_container_width=True, hide_index=True)

    st.markdown("#### Interpretation")
    confirmed = result.loc[result["concordant"], "population"].tolist()
    candidates = result.loc[result["significant_primary"]
                            & ~result["concordant"], "population"].tolist()
    group_var = float(result["group_var"].abs().max())
    lines = []
    if confirmed:
        lines.append(f"**Confirmed:** {', '.join(confirmed)} differ between "
                     "responders and non-responders under both the mixed model "
                     "and the non-parametric cross-check.")
    if candidates:
        lines.append(f"**Candidate(s):** {', '.join(candidates)} reach "
                     "significance under the primary model but not the "
                     "conservative per-subject test — worth validating, not an "
                     "established biomarker.")
    if not confirmed and not candidates:
        lines.append("No population shows a significant difference between "
                     "responders and non-responders.")
    lines.append(f"The between-subject (random-effect) variance is "
                 f"~{group_var:.2g} across populations, so the repeated "
                 "timepoints add little within-subject correlation here.")
    st.markdown("\n\n".join(lines))
