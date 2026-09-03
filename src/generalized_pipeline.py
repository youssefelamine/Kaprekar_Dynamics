"""Run the exact generalized Kaprekar census, figures, reports, and proof."""

from __future__ import annotations

import argparse
from pathlib import Path

from .generalized_analysis import run_generalized_census
from .generalized_report import write_generalized_report
from .generalized_visualization import generate_generalized_figures
from .proof_certificate import write_pair_certificate, write_proof_report
from .unified_report import write_unified_report


def parse_integer_selection(text: str) -> tuple[int, ...]:
    """Parse an inclusive ``start:end`` range or a comma-separated integer list."""

    text = text.strip()
    if ":" in text:
        parts = text.split(":")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError("ranges must use start:end")
        start, end = map(int, parts)
        if end < start:
            raise argparse.ArgumentTypeError("range end must not precede range start")
        return tuple(range(start, end + 1))
    try:
        values = tuple(sorted(set(int(part.strip()) for part in text.split(",") if part.strip())))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("values must be integers") from exc
    if not values:
        raise argparse.ArgumentTypeError("selection cannot be empty")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bases", type=parse_integer_selection, default=parse_integer_selection("2:16"))
    parser.add_argument("--digits", type=parse_integer_selection, default=parse_integer_selection("2:6"))
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    return parser.parse_args()


def run_pipeline(project_root: Path, bases: tuple[int, ...], digit_widths: tuple[int, ...]) -> dict:
    """Execute the additive generalized pipeline and return its aggregate summary."""

    data_dir = project_root / "data" / "generalized"
    tables_dir = project_root / "tables" / "generalized"
    figures_dir = project_root / "figures" / "generalized"
    report_dir = project_root / "report"

    summary = run_generalized_census(bases, digit_widths, data_dir, tables_dir)
    write_pair_certificate(tables_dir / "kaprekar_6174_pair_certificate.csv")
    generate_generalized_figures(tables_dir / "system_summary.csv", figures_dir)
    write_proof_report(
        report_dir / "6174_finite_proof.md",
        "../tables/generalized/kaprekar_6174_pair_certificate.csv",
    )
    write_generalized_report(
        data_dir / "generalized_summary.json",
        tables_dir / "system_summary.csv",
        report_dir / "generalized_kaprekar_report.md",
    )
    write_unified_report(project_root)
    for chapter in (
        report_dir / "6174_finite_proof.md",
        report_dir / "generalized_kaprekar_report.md",
    ):
        chapter.unlink()
    return summary


def main() -> None:
    args = parse_args()
    summary = run_pipeline(args.project_root.resolve(), args.bases, args.digits)
    print(
        f"Generalized pipeline complete: {summary['system_count']} systems, "
        f"{summary['ordered_states_represented']:,} ordered states represented by "
        f"{summary['weighted_class_records']:,} exact weighted classes."
    )


if __name__ == "__main__":
    main()
