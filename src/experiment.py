"""Exhaustively enumerate the 10,000 four-digit Kaprekar states."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Iterable

from .kaprekar import STATE_COUNT, analyze_start, normalize_state


CSV_FIELDS = (
    "start",
    "start_state",
    "is_repdigit",
    "distinct_digits",
    "digit_multiset",
    "x",
    "y",
    "descending",
    "ascending",
    "first_result",
    "first_result_state",
    "trajectory",
    "final_attractor",
    "final_attractor_state",
    "iterations_to_attractor",
    "distance_to_6174",
    "distance_to_0000",
    "reaches_6174",
    "reaches_0000",
    "other_cycle_detected",
    "cycle_length",
    "cycle",
)


def _format_path(states: Iterable[int], *, close_cycle: bool = False) -> str:
    values = list(states)
    if close_cycle and values:
        values.append(values[0])
    return " -> ".join(normalize_state(state) for state in values)


def serialize_record(record: dict[str, Any]) -> dict[str, Any]:
    """Convert list/tuple fields from :func:`analyze_start` to CSV text."""

    serialized = dict(record)
    serialized["trajectory"] = _format_path(record["trajectory"])
    serialized["cycle"] = _format_path(record["cycle"], close_cycle=True)
    serialized["final_attractor_state"] = normalize_state(record["final_attractor"])
    for nullable_distance in ("distance_to_6174", "distance_to_0000"):
        if serialized[nullable_distance] is None:
            serialized[nullable_distance] = ""
    return serialized


def generate_results() -> list[dict[str, Any]]:
    """Analyze every state from 0000 through 9999."""

    return [serialize_record(analyze_start(state)) for state in range(STATE_COUNT)]


def write_results(records: list[dict[str, Any]], output_path: Path) -> None:
    """Write exhaustive records to *output_path* with stable column ordering."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows({field: record[field] for field in CSV_FIELDS} for record in records)


def validate_results(records: list[dict[str, Any]]) -> None:
    """Fail fast if enumeration or classification is internally inconsistent."""

    if len(records) != STATE_COUNT:
        raise AssertionError(f"expected {STATE_COUNT} rows, generated {len(records)}")
    starts = [int(record["start"]) for record in records]
    if starts != list(range(STATE_COUNT)):
        raise AssertionError("results do not cover each state exactly once in numeric order")
    for record in records:
        if record["start_state"] != normalize_state(int(record["start"])):
            raise AssertionError(f"bad normalization for state {record['start']}")
        outcomes = sum(bool(record[key]) for key in ("reaches_6174", "reaches_0000", "other_cycle_detected"))
        if outcomes != 1:
            raise AssertionError(f"non-exclusive outcome classification for {record['start_state']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/kaprekar_results.csv"),
        help="destination CSV (default: data/kaprekar_results.csv)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = generate_results()
    validate_results(records)
    write_results(records, args.output)
    print(f"Wrote {len(records):,} exhaustive state records to {args.output}")


if __name__ == "__main__":
    main()
