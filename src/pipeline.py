"""Regenerate the complete Kaprekar dataset, tables, summary, and figures."""

from __future__ import annotations

from pathlib import Path

from .analysis import generate_analysis
from .experiment import generate_results, validate_results, write_results
from .visualization import generate_figures


def main() -> None:
    project_root = Path.cwd()
    data_path = project_root / "data" / "kaprekar_results.csv"
    summary_path = project_root / "data" / "analysis_summary.json"
    tables_dir = project_root / "tables"
    figures_dir = project_root / "figures"

    records = generate_results()
    validate_results(records)
    write_results(records, data_path)
    summary = generate_analysis(data_path, tables_dir, summary_path)
    generate_figures(data_path, tables_dir, figures_dir)
    print(
        f"Pipeline complete: {summary['total_states']:,} states, "
        f"{summary['unique_one_step_outputs']} one-step outputs, "
        f"maximum valid distance {summary['maximum_distance']}."
    )


if __name__ == "__main__":
    main()

