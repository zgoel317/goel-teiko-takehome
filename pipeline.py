"""Runs the full data pipeline: build the database (Part 1) and generate every
required output table and plot (Parts 2-4) into outputs/.

Run directly:
    python3 pipeline.py
"""
from pathlib import Path

from load_data import build_database, CSV_PATH, DB_PATH
from db.connection import get_connection
from analysis.frequencies import sample_frequencies
from analysis.comparison import (
    compare_responders,
    cohort_frequencies,
    responder_boxplots,
)
from analysis.subsets import (
    baseline_subset,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)
from analysis.questions import avg_bcells_melanoma_male_responders_at_baseline

ROOT = Path(__file__).resolve().parent
OUTPUTS_DIR = ROOT / "outputs"


def run_pipeline(csv_path=CSV_PATH, db_path=DB_PATH, outdir=OUTPUTS_DIR):
    """Build the DB and write all Part 2-4 artifacts. Returns {name: path}."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    build_database(csv_path=csv_path, db_path=db_path)  # Part 1
    conn = get_connection(db_path)
    written = {}
    try:
        # Part 2
        freq_path = outdir / "frequencies.csv"
        sample_frequencies(conn).to_csv(freq_path, index=False)
        written["frequencies"] = freq_path

        # Part 3
        comp_path = outdir / "responder_comparison.csv"
        compare_responders(conn).to_csv(comp_path, index=False)
        written["responder_comparison"] = comp_path

        box_path = outdir / "responder_boxplots.html"
        responder_boxplots(cohort_frequencies(conn)).write_html(box_path)
        written["responder_boxplots"] = box_path

        # Part 4
        subset = baseline_subset(conn)
        subset_path = outdir / "baseline_subset.csv"
        subset.to_csv(subset_path, index=False)
        written["baseline_subset"] = subset_path

        spp_path = outdir / "baseline_samples_per_project.csv"
        samples_per_project(subset).to_csv(spp_path, index=False)
        written["baseline_samples_per_project"] = spp_path

        sbr_path = outdir / "baseline_subjects_by_response.csv"
        subjects_by_response(subset).to_csv(sbr_path, index=False)
        written["baseline_subjects_by_response"] = sbr_path

        sbs_path = outdir / "baseline_subjects_by_sex.csv"
        subjects_by_sex(subset).to_csv(sbs_path, index=False)
        written["baseline_subjects_by_sex"] = sbs_path

        answer = avg_bcells_melanoma_male_responders_at_baseline(conn)
        answer_path = outdir / "final_answer.txt"
        answer_path.write_text(f"{answer}\n")
        written["final_answer"] = answer_path
    finally:
        conn.close()

    return {name: str(path) for name, path in written.items()}


def main():
    written = run_pipeline()
    print("Pipeline complete. Wrote:")
    for name, path in written.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
