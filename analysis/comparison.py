"""Part 3: compare cell-population frequencies between miraclib responders and
non-responders (melanoma, PBMC). Builds on Part 2's frequency summary table."""
import warnings

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import plotly.express as px

from load_data import POPULATION_COLUMNS
from analysis.frequencies import sample_frequencies

POPULATIONS = list(POPULATION_COLUMNS)

_METADATA_QUERY = """
SELECT sa.sample_id                         AS sample,
       sa.subject_id                        AS subject_id,
       sa.sample_type                       AS sample_type,
       sa.time_from_treatment_start         AS time,
       su.condition                         AS condition,
       su.treatment                         AS treatment,
       su.response                          AS response,
       su.sex                               AS sex
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
"""


def cohort_frequencies(conn, condition="melanoma", treatment="miraclib",
                       sample_type="PBMC", timepoint=None):
    """Per-sample relative frequencies for the response-comparison cohort.

    Reuses Part 2's sample_frequencies (single source of truth) and joins the
    sample/subject metadata, then filters to the cohort and to responders /
    non-responders only (response in {'yes','no'}).
    """
    freq = sample_frequencies(conn)
    meta = pd.read_sql_query(_METADATA_QUERY, conn)
    df = freq.merge(meta, on="sample", how="inner")
    mask = (
        (df["condition"] == condition)
        & (df["treatment"] == treatment)
        & (df["sample_type"] == sample_type)
        & (df["response"].isin(["yes", "no"]))
    )
    df = df[mask]
    if timepoint is not None:
        df = df[df["time"] == timepoint]
    keep = ["sample", "total_count", "population", "count", "percentage",
            "subject_id", "time", "response", "sex"]
    return df[keep].reset_index(drop=True)


def subject_means(cohort):
    """Collapse to one value per (subject, population): mean % across timepoints.

    Restores independence (one row per subject) for the non-parametric and
    bootstrap steps. sort=False preserves the population order of the input
    (canonical, as produced by cohort_frequencies) instead of pandas' default
    alphabetical sort.
    """
    return (
        cohort.groupby(["subject_id", "response", "population"], sort=False)["percentage"]
        .mean()
        .reset_index()
    )


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR-adjusted q-values (order-preserving, clipped [0,1])."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p downward
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(q_sorted, 0, 1)
    return q


def mannwhitney_by_population(cohort, populations=POPULATIONS):
    """Two-sided Mann-Whitney U on per-subject means, responder vs non-responder."""
    sm = subject_means(cohort)
    rows = []
    for pop in populations:
        d = sm[sm["population"] == pop]
        r = d[d["response"] == "yes"]["percentage"].to_numpy()
        nr = d[d["response"] == "no"]["percentage"].to_numpy()
        u, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        n1, n2 = len(r), len(nr)
        # rank-biserial correlation: +1 when responders stochastically dominate
        rank_biserial = 2.0 * u / (n1 * n2) - 1.0
        rows.append(dict(population=pop, mw_U=float(u), mw_p=float(p),
                         rank_biserial=rank_biserial,
                         n_responders=n1, n_nonresponders=n2,
                         median_resp=float(np.median(r)),
                         median_nonresp=float(np.median(nr))))
    return pd.DataFrame(rows)


_RESPONSE_TERM = "response[T.yes]"


def _fit_population_effect(d, formula):
    """Estimate the responder effect for one population.

    Primary: linear mixed model with a per-subject random intercept. When the
    random-effect variance collapses to ~0 the covariance matrix inverted during
    the solver's gradient step becomes near-singular, and whether numpy's LAPACK
    backend raises ``LinAlgError`` on it is platform-dependent (it does not on
    macOS/Accelerate but does on some Linux/OpenBLAS builds — e.g. Streamlit
    Cloud and CI). In that limit the mixed model reduces to ordinary least
    squares, so we fall back to OLS (reporting ``group_var`` as 0) rather than
    crash. Returns (coef, ci_low, ci_high, p, group_var).

    Suppresses the two benign diagnostics that fire as the variance collapses
    (ConvergenceWarning and the "Random effects covariance is singular"
    UserWarning).
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        warnings.filterwarnings(
            "ignore", message="Random effects covariance is singular")
        try:
            fit = smf.mixedlm(formula, d, groups=d["subject_id"]).fit(method="lbfgs")
            ci = fit.conf_int()
            return (float(fit.params[_RESPONSE_TERM]),
                    float(ci.loc[_RESPONSE_TERM, 0]),
                    float(ci.loc[_RESPONSE_TERM, 1]),
                    float(fit.pvalues[_RESPONSE_TERM]),
                    float(fit.cov_re.iloc[0, 0]))
        except (np.linalg.LinAlgError, ValueError):
            ols = smf.ols(formula, d).fit()
            ci = ols.conf_int()
            return (float(ols.params[_RESPONSE_TERM]),
                    float(ci.loc[_RESPONSE_TERM, 0]),
                    float(ci.loc[_RESPONSE_TERM, 1]),
                    float(ols.pvalues[_RESPONSE_TERM]),
                    0.0)


def mixedlm_by_population(cohort, populations=POPULATIONS):
    """Per-population responder effect via a mixed-effects model (OLS fallback).

    `percentage ~ response (+ time)`, random intercept per subject. Records the
    random-effect variance (group_var). See `_fit_population_effect` for the
    fallback behaviour when the random-effect covariance is singular."""
    rows = []
    for pop in populations:
        d = cohort[cohort["population"] == pop]
        formula = "percentage ~ response"
        if d["time"].nunique() > 1:
            formula += " + C(time)"
        coef, ci_low, ci_high, mixed_p, group_var = _fit_population_effect(d, formula)
        rows.append(dict(population=pop, coef=coef, ci_low=ci_low,
                         ci_high=ci_high, mixed_p=mixed_p, group_var=group_var))
    return pd.DataFrame(rows)


def bootstrap_diff_ci(cohort, populations=POPULATIONS, n_boot=2000, seed=0):
    """Bootstrap 95% CI for the responder - non-responder difference in
    per-subject mean %, per population. Seeded for reproducibility."""
    rng = np.random.default_rng(seed)
    sm = subject_means(cohort)
    rows = []
    for pop in populations:
        d = sm[sm["population"] == pop]
        r = d[d["response"] == "yes"]["percentage"].to_numpy()
        nr = d[d["response"] == "no"]["percentage"].to_numpy()
        mean_diff = r.mean() - nr.mean()
        boots = np.empty(n_boot)
        for b in range(n_boot):
            rb = rng.choice(r, len(r), replace=True).mean()
            nrb = rng.choice(nr, len(nr), replace=True).mean()
            boots[b] = rb - nrb
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append(dict(population=pop, mean_diff=float(mean_diff),
                         boot_ci_low=float(lo), boot_ci_high=float(hi)))
    return pd.DataFrame(rows)


def compare_responders(conn, condition="melanoma", treatment="miraclib",
                       sample_type="PBMC", alpha=0.05, seed=0):
    """Full responder vs non-responder comparison, one row per population.

    Combines the mixed-effects model (primary), Mann-Whitney U on per-subject
    means (non-parametric cross-check), bootstrap effect-size CIs, and
    Benjamini-Hochberg FDR across the five populations.
    """
    cohort = cohort_frequencies(conn, condition, treatment, sample_type)

    mm = mixedlm_by_population(cohort, POPULATIONS)
    mw = mannwhitney_by_population(cohort, POPULATIONS)
    boot = bootstrap_diff_ci(cohort, POPULATIONS, seed=seed)

    sm = subject_means(cohort)
    group_means = (
        sm.pivot_table(index="population", columns="response",
                       values="percentage", aggfunc="mean")
        .rename(columns={"yes": "mean_resp", "no": "mean_nonresp"})
        .reset_index()
    )

    res = (mm.merge(mw, on="population")
              .merge(boot, on="population")
              .merge(group_means, on="population"))

    res["q_mixed"] = benjamini_hochberg(res["mixed_p"].to_numpy())
    res["q_mw"] = benjamini_hochberg(res["mw_p"].to_numpy())
    res["significant_primary"] = res["q_mixed"] < alpha
    res["concordant"] = (res["q_mixed"] < alpha) & (res["q_mw"] < alpha)

    # canonical population order + tidy column order
    res["population"] = pd.Categorical(res["population"], categories=POPULATIONS,
                                       ordered=True)
    res = res.sort_values("population").reset_index(drop=True)
    res["population"] = res["population"].astype(str)
    columns = ["population", "n_responders", "n_nonresponders", "mean_resp",
               "mean_nonresp", "mean_diff", "boot_ci_low", "boot_ci_high",
               "coef", "ci_low", "ci_high", "mixed_p", "q_mixed", "group_var",
               "mw_U", "mw_p", "q_mw", "rank_biserial", "median_resp",
               "median_nonresp", "significant_primary", "concordant"]
    return res[columns]


def responder_boxplots(cohort, populations=POPULATIONS):
    """One boxplot panel per population comparing responders vs non-responders,
    built from per-subject mean % (individual points overlaid)."""
    sm = subject_means(cohort)
    sm = sm[sm["population"].isin(populations)]
    fig = px.box(
        sm, x="response", y="percentage", color="response",
        facet_col="population", facet_col_wrap=len(populations),
        points="all", boxmode="overlay",
        category_orders={"population": list(populations),
                         "response": ["no", "yes"]},
        labels={"percentage": "Relative frequency (%)", "response": "Response"},
    )
    # strip the "population=" prefix plotly adds to each facet title
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(
        showlegend=False,
        title="Cell-population relative frequency: responders vs non-responders",
    )
    return fig
