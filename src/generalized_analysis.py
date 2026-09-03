"""Exact symmetry-weighted analysis for generalized Kaprekar systems."""

from __future__ import annotations

import csv
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any, Iterable

from .kaprekar import (
    DIGIT_SYMBOLS,
    canonical_cycle,
    digits_in_base,
    format_state_general,
    integer_from_digits,
    kaprekar_step_general,
)


CLASS_FIELDS = (
    "base",
    "digits",
    "digit_multiset",
    "class_weight",
    "is_repdigit",
    "successor",
    "successor_state",
    "attractor_id",
    "cycle_length",
    "cycle_member_count",
    "cycle_member_state",
    "noncycle_member_count",
    "noncycle_transient_depth",
    "all_members_same_depth",
)

SYSTEM_FIELDS = (
    "base",
    "digits",
    "ordered_state_count",
    "repdigit_state_count",
    "valid_state_count",
    "permutation_class_count",
    "unique_output_count",
    "attractor_count",
    "fixed_point_count",
    "nontrivial_cycle_count",
    "maximum_cycle_length",
    "maximum_transient_depth_all",
    "maximum_transient_depth_valid",
    "valid_mean_transient_depth",
    "valid_median_transient_depth",
    "valid_population_standard_deviation",
    "valid_attractor_count",
    "unanimous_valid_attractor",
    "largest_valid_basin_size",
    "largest_valid_basin_percentage",
)

ATTRACTOR_FIELDS = (
    "base",
    "digits",
    "attractor_id",
    "cycle",
    "cycle_length",
    "kind",
    "basin_size_all",
    "basin_percentage_all",
    "basin_size_valid",
    "basin_percentage_valid",
    "repdigit_basin_size",
    "maximum_transient_depth_all",
    "maximum_transient_depth_valid",
)

DEPTH_FIELDS = (
    "base",
    "digits",
    "attractor_id",
    "transient_depth",
    "state_count_all",
    "state_percentage_all",
    "state_count_valid",
    "state_percentage_valid",
)


@dataclass
class SystemAnalysis:
    """All exact records produced for one ``(base, digits)`` system."""

    summary: dict[str, Any]
    class_records: list[dict[str, Any]]
    attractor_records: list[dict[str, Any]]
    depth_records: list[dict[str, Any]]


def multiset_weight(multiset: tuple[int, ...]) -> int:
    """Return the number of ordered fixed-width states represented by a multiset."""

    denominator = math.prod(math.factorial(count) for count in Counter(multiset).values())
    return math.factorial(len(multiset)) // denominator


def multiset_successor(multiset: tuple[int, ...], base: int) -> int:
    """Apply the sort-and-subtract operation directly to a sorted digit multiset."""

    ascending = integer_from_digits(multiset, base)
    descending = integer_from_digits(tuple(reversed(multiset)), base)
    return descending - ascending


def _cycle_id(cycle: tuple[int, ...], digits: int, base: int) -> str:
    return " -> ".join(format_state_general(state, digits, base) for state in cycle)


def _decompose_functional_graph(
    successors: dict[int, int],
) -> tuple[dict[int, tuple[tuple[int, ...], int]], tuple[tuple[int, ...], ...]]:
    """Return ``node -> (terminal cycle, depth)`` and all cycles."""

    information: dict[int, tuple[tuple[int, ...], int]] = {}
    cycles: set[tuple[int, ...]] = set()

    for start in sorted(successors):
        if start in information:
            continue
        path: list[int] = []
        local_index: dict[int, int] = {}
        current = start
        while current not in information and current not in local_index:
            if current not in successors:
                raise AssertionError("the image of a Kaprekar map must be closed under the map")
            local_index[current] = len(path)
            path.append(current)
            current = successors[current]

        if current in local_index:
            cycle_start = local_index[current]
            cycle = canonical_cycle(tuple(path[cycle_start:]))
            cycles.add(cycle)
            for node in path[cycle_start:]:
                information[node] = (cycle, 0)
            prefix = path[:cycle_start]
        else:
            prefix = path

        for node in reversed(prefix):
            successor = successors[node]
            terminal_cycle, successor_depth = information[successor]
            information[node] = (terminal_cycle, successor_depth + 1)

    return information, tuple(sorted(cycles))


def _weighted_statistics(distribution: Counter[int]) -> tuple[float, float, float]:
    total = sum(distribution.values())
    if total <= 0:
        raise ValueError("a weighted distribution cannot be empty")
    mean = sum(depth * count for depth, count in distribution.items()) / total
    variance = sum((depth - mean) ** 2 * count for depth, count in distribution.items()) / total

    positions = ((total - 1) // 2, total // 2)
    median_values: list[int] = []
    cumulative = 0
    position_index = 0
    for depth, count in sorted(distribution.items()):
        cumulative += count
        while position_index < len(positions) and positions[position_index] < cumulative:
            median_values.append(depth)
            position_index += 1
    return mean, sum(median_values) / 2, math.sqrt(variance)


def analyze_system(base: int, digits: int) -> SystemAnalysis:
    """Analyze one generalized system exactly through weighted permutation classes."""

    if not 2 <= base <= 36:
        raise ValueError("base must be between 2 and 36")
    if digits < 2:
        raise ValueError("digits must be at least 2")

    multisets = list(combinations_with_replacement(range(base), digits))
    class_successors = {multiset: multiset_successor(multiset, base) for multiset in multisets}
    outputs = set(class_successors.values())
    output_successors = {
        output: kaprekar_step_general(output, digits=digits, base=base) for output in outputs
    }
    if not set(output_successors.values()).issubset(outputs):
        raise AssertionError("unique Kaprekar outputs are not closed under iteration")
    node_information, cycles = _decompose_functional_graph(output_successors)

    cycle_node_by_multiset: dict[tuple[int, ...], int] = {}
    for cycle in cycles:
        for node in cycle:
            key = tuple(sorted(digits_in_base(node, digits, base)))
            if key in cycle_node_by_multiset:
                raise AssertionError("a digit multiset contains more than one cycle node")
            cycle_node_by_multiset[key] = node

    all_depths: Counter[int] = Counter()
    valid_depths: Counter[int] = Counter()
    all_depths_by_attractor: dict[str, Counter[int]] = defaultdict(Counter)
    valid_depths_by_attractor: dict[str, Counter[int]] = defaultdict(Counter)
    repdigit_basin_counts: Counter[str] = Counter()
    class_records: list[dict[str, Any]] = []

    for multiset in multisets:
        weight = multiset_weight(multiset)
        successor = class_successors[multiset]
        cycle, successor_depth = node_information[successor]
        attractor_id = _cycle_id(cycle, digits, base)
        cycle_member = cycle_node_by_multiset.get(multiset)
        cycle_member_count = int(cycle_member is not None)
        noncycle_count = weight - cycle_member_count
        noncycle_depth = successor_depth + 1
        is_repdigit = len(set(multiset)) == 1

        if cycle_member is not None and cycle_member not in cycle:
            raise AssertionError("cycle-member multiset was assigned to a different attractor")
        if cycle_member_count:
            all_depths[0] += 1
            all_depths_by_attractor[attractor_id][0] += 1
            if not is_repdigit:
                valid_depths[0] += 1
                valid_depths_by_attractor[attractor_id][0] += 1
        if noncycle_count:
            all_depths[noncycle_depth] += noncycle_count
            all_depths_by_attractor[attractor_id][noncycle_depth] += noncycle_count
            if not is_repdigit:
                valid_depths[noncycle_depth] += noncycle_count
                valid_depths_by_attractor[attractor_id][noncycle_depth] += noncycle_count
        if is_repdigit:
            repdigit_basin_counts[attractor_id] += weight

        class_records.append(
            {
                "base": base,
                "digits": digits,
                "digit_multiset": "".join(DIGIT_SYMBOLS[value] for value in multiset),
                "class_weight": weight,
                "is_repdigit": is_repdigit,
                "successor": successor,
                "successor_state": format_state_general(successor, digits, base),
                "attractor_id": attractor_id,
                "cycle_length": len(cycle),
                "cycle_member_count": cycle_member_count,
                "cycle_member_state": (
                    format_state_general(cycle_member, digits, base) if cycle_member is not None else ""
                ),
                "noncycle_member_count": noncycle_count,
                "noncycle_transient_depth": noncycle_depth if noncycle_count else "",
                "all_members_same_depth": not (cycle_member_count and noncycle_count),
            }
        )

    ordered_state_count = base**digits
    valid_state_count = ordered_state_count - base
    if sum(all_depths.values()) != ordered_state_count:
        raise AssertionError("weighted class sizes do not recover the ordered state space")
    if sum(valid_depths.values()) != valid_state_count:
        raise AssertionError("weighted valid class sizes are inconsistent")

    valid_mean, valid_median, valid_pstdev = _weighted_statistics(valid_depths)
    valid_basin_sizes = {
        attractor_id: sum(counts.values())
        for attractor_id, counts in valid_depths_by_attractor.items()
        if sum(counts.values()) > 0
    }
    largest_valid_basin = max(valid_basin_sizes.values())

    attractor_records: list[dict[str, Any]] = []
    for cycle in cycles:
        attractor_id = _cycle_id(cycle, digits, base)
        all_distribution = all_depths_by_attractor[attractor_id]
        valid_distribution = valid_depths_by_attractor[attractor_id]
        all_basin = sum(all_distribution.values())
        valid_basin = sum(valid_distribution.values())
        attractor_records.append(
            {
                "base": base,
                "digits": digits,
                "attractor_id": attractor_id,
                "cycle": attractor_id + " -> " + format_state_general(cycle[0], digits, base),
                "cycle_length": len(cycle),
                "kind": "fixed point" if len(cycle) == 1 else "non-trivial cycle",
                "basin_size_all": all_basin,
                "basin_percentage_all": 100 * all_basin / ordered_state_count,
                "basin_size_valid": valid_basin,
                "basin_percentage_valid": 100 * valid_basin / valid_state_count,
                "repdigit_basin_size": repdigit_basin_counts[attractor_id],
                "maximum_transient_depth_all": max(all_distribution),
                "maximum_transient_depth_valid": max(valid_distribution) if valid_distribution else "",
            }
        )

    depth_records: list[dict[str, Any]] = []
    for attractor_id in sorted(all_depths_by_attractor):
        all_distribution = all_depths_by_attractor[attractor_id]
        valid_distribution = valid_depths_by_attractor[attractor_id]
        for depth in range(max(all_distribution) + 1):
            all_count = all_distribution[depth]
            valid_count = valid_distribution[depth]
            depth_records.append(
                {
                    "base": base,
                    "digits": digits,
                    "attractor_id": attractor_id,
                    "transient_depth": depth,
                    "state_count_all": all_count,
                    "state_percentage_all": 100 * all_count / ordered_state_count,
                    "state_count_valid": valid_count,
                    "state_percentage_valid": 100 * valid_count / valid_state_count,
                }
            )

    summary = {
        "base": base,
        "digits": digits,
        "ordered_state_count": ordered_state_count,
        "repdigit_state_count": base,
        "valid_state_count": valid_state_count,
        "permutation_class_count": len(multisets),
        "unique_output_count": len(outputs),
        "attractor_count": len(cycles),
        "fixed_point_count": sum(len(cycle) == 1 for cycle in cycles),
        "nontrivial_cycle_count": sum(len(cycle) > 1 for cycle in cycles),
        "maximum_cycle_length": max(map(len, cycles)),
        "maximum_transient_depth_all": max(all_depths),
        "maximum_transient_depth_valid": max(valid_depths),
        "valid_mean_transient_depth": valid_mean,
        "valid_median_transient_depth": valid_median,
        "valid_population_standard_deviation": valid_pstdev,
        "valid_attractor_count": len(valid_basin_sizes),
        "unanimous_valid_attractor": len(valid_basin_sizes) == 1,
        "largest_valid_basin_size": largest_valid_basin,
        "largest_valid_basin_percentage": 100 * largest_valid_basin / valid_state_count,
    }
    return SystemAnalysis(summary, class_records, attractor_records, depth_records)


def _package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for package in ("pandas", "numpy", "matplotlib", "networkx"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def run_generalized_census(
    bases: Iterable[int],
    digit_widths: Iterable[int],
    data_dir: Path,
    tables_dir: Path,
) -> dict[str, Any]:
    """Run and export a deterministic collection of generalized systems."""

    bases = tuple(sorted(set(int(base) for base in bases)))
    digit_widths = tuple(sorted(set(int(digits) for digits in digit_widths)))
    if not bases or not digit_widths:
        raise ValueError("at least one base and digit width are required")
    data_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    class_path = data_dir / "weighted_classes.csv"
    system_summaries: list[dict[str, Any]] = []
    attractor_records: list[dict[str, Any]] = []
    depth_records: list[dict[str, Any]] = []
    total_class_records = 0

    with class_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CLASS_FIELDS)
        writer.writeheader()
        for base in bases:
            for digits in digit_widths:
                analysis = analyze_system(base, digits)
                writer.writerows(analysis.class_records)
                total_class_records += len(analysis.class_records)
                system_summaries.append(analysis.summary)
                attractor_records.extend(analysis.attractor_records)
                depth_records.extend(analysis.depth_records)

    def write_rows(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    write_rows(tables_dir / "system_summary.csv", SYSTEM_FIELDS, system_summaries)
    write_rows(tables_dir / "attractors_and_basins.csv", ATTRACTOR_FIELDS, attractor_records)
    write_rows(tables_dir / "depth_distributions.csv", DEPTH_FIELDS, depth_records)

    summary = {
        "base_values": list(bases),
        "digit_widths": list(digit_widths),
        "system_count": len(system_summaries),
        "ordered_states_represented": sum(row["ordered_state_count"] for row in system_summaries),
        "weighted_class_records": total_class_records,
        "attractor_records": len(attractor_records),
        "systems": system_summaries,
    }
    (data_dir / "generalized_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    metadata = {
        "algorithm": "exact symmetry-weighted digit-multiset enumeration",
        "base_values": list(bases),
        "digit_widths": list(digit_widths),
        "python": sys.version.split()[0],
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "packages": _package_versions(),
    }
    (data_dir / "reproducibility_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary

