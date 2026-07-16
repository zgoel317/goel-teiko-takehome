import numpy as np
import pandas as pd
import plotly.graph_objects as go

from analysis.comparison import responder_boxplots


def _cohort():
    rng = np.random.default_rng(0)
    rows = []
    for resp, is_r in [("yes", 1), ("no", 0)]:
        for s in range(10):
            for t in (0, 7, 14):
                for pop in ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]:
                    rows.append(dict(subject_id=f"{resp}{s}", response=resp, time=t,
                                     population=pop,
                                     percentage=20 + (5 * is_r if pop == "cd4_t_cell" else 0)
                                     + rng.normal(0, 1)))
    return pd.DataFrame(rows)


def test_returns_plotly_figure_with_data():
    fig = responder_boxplots(_cohort())
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0


def test_covers_all_five_populations():
    fig = responder_boxplots(_cohort())
    text = str(fig.to_dict())
    for pop in ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]:
        assert pop in text


def test_both_response_groups_present():
    fig = responder_boxplots(_cohort())
    text = str(fig.to_dict())
    assert "yes" in text and "no" in text


def test_plot_uses_per_subject_means_not_raw_samples():
    # subject_means collapses the 3 timepoints: 10 subjects x 2 groups x 5 pops = 100
    # points, NOT the 300 raw sample rows. Guards the statistical-honesty requirement.
    fig = responder_boxplots(_cohort())
    n_points = sum(len(t.y) for t in fig.data if getattr(t, "y", None) is not None)
    assert n_points == 100


def test_populations_param_filters_panels():
    fig = responder_boxplots(_cohort(), populations=["b_cell", "cd8_t_cell"])
    text = str(fig.to_dict())
    assert "cd4_t_cell" not in text and "nk_cell" not in text
    n_points = sum(len(t.y) for t in fig.data if getattr(t, "y", None) is not None)
    assert n_points == 40  # 10 subjects x 2 groups x 2 populations
