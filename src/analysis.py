"""Statistical, permutation, and graph analysis of exhaustive results."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from statistics import fmean, median, pstdev
from typing import Any, Iterable

import networkx as nx
import pandas as pd


TEXT_COLUMNS = {
    "start_state": "string",
    "digit_multiset": "string",
    "first_result_state": "string",
    "trajectory": "string",
    "final_attractor_state": "string",
    "cycle": "string",
}


def load_results(path: Path) -> pd.DataFrame:
    """Load and validate the exhaustive experiment CSV."""

    frame = pd.read_csv(path, dtype=TEXT_COLUMNS, keep_default_na=True)
    if len(frame) != 10_000 or frame["start"].nunique() != 10_000:
        raise ValueError("the exhaustive dataset must contain 10,000 unique starts")
    if sorted(frame["start"].astype(int)) != list(range(10_000)):
        raise ValueError("the exhaustive dataset does not cover 0000..9999 exactly")
    for column in ("start_state", "digit_multiset", "first_result_state", "final_attractor_state"):
        frame[column] = frame[column].str.zfill(4)
    for column in ("is_repdigit", "reaches_6174", "reaches_0000", "other_cycle_detected"):
        if frame[column].dtype != bool:
            frame[column] = frame[column].map({"True": True, "False": False})
    return frame.sort_values("start").reset_index(drop=True)


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def _parse_path(text: str) -> list[int]:
    return [int(part.strip()) for part in text.split("->")]


def _numeric_summary(label: str, values: Iterable[int]) -> dict[str, Any]:
    data = [int(value) for value in values]
    return {
        "scope": label,
        "count": len(data),
        "minimum": min(data),
        "maximum": max(data),
        "mean": fmean(data),
        "median": median(data),
        "population_standard_deviation": pstdev(data),
    }


def _overall_classification(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    scopes = {
        "four-digit states 0000-9999": frame,
        "ordinary integers 1000-9999": frame[frame["start"] >= 1000],
    }
    for scope, subset in scopes.items():
        size = len(subset)
        categories = {
            "all states": size,
            "valid non-repdigit states": int((~subset["is_repdigit"]).sum()),
            "repdigits": int(subset["is_repdigit"].sum()),
            "basin of 6174": int(subset["reaches_6174"].sum()),
            "basin of 0000": int(subset["reaches_0000"].sum()),
            "other-cycle basin": int(subset["other_cycle_detected"].sum()),
        }
        for category, count in categories.items():
            rows.append(
                {
                    "scope": scope,
                    "category": category,
                    "count": count,
                    "percentage_of_scope": 100 * count / size,
                }
            )
    return pd.DataFrame(rows)


def _iteration_distribution(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    scopes = {
        "all valid states 0000-9999": frame[~frame["is_repdigit"]],
        "valid ordinary integers 1000-9999": frame[(frame["start"] >= 1000) & ~frame["is_repdigit"]],
    }
    largest_depth = int(max(subset["distance_to_6174"].max() for subset in scopes.values()))
    for scope, subset in scopes.items():
        counts = subset["distance_to_6174"].astype(int).value_counts()
        for depth in range(largest_depth + 1):
            count = int(counts.get(depth, 0))
            rows.append(
                {
                    "scope": scope,
                    "transformations_to_6174": depth,
                    "count": count,
                    "percentage": 100 * count / len(subset),
                }
            )
    return pd.DataFrame(rows)


def _cycle_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cycle, group in frame.groupby("cycle", sort=True):
        cycle_length = int(group["cycle_length"].iloc[0])
        cycle_nodes = [part.strip() for part in str(cycle).split("->")][:-1]
        rows.append(
            {
                "cycle": cycle,
                "cycle_nodes": " ".join(cycle_nodes),
                "cycle_length": cycle_length,
                "kind": "fixed point" if cycle_length == 1 else "non-trivial cycle",
                "basin_size": len(group),
                "maximum_basin_depth": int(group["iterations_to_attractor"].max()),
            }
        )
    return pd.DataFrame(rows).sort_values(["cycle_length", "cycle"]).reset_index(drop=True)


def _validate_graph(frame: pd.DataFrame, cycles: pd.DataFrame) -> dict[str, int]:
    graph = nx.DiGraph()
    graph.add_edges_from(zip(frame["start"].astype(int), frame["first_result"].astype(int)))
    cyclic_components: set[frozenset[int]] = set()
    for component in nx.strongly_connected_components(graph):
        if len(component) > 1:
            cyclic_components.add(frozenset(component))
        else:
            node = next(iter(component))
            if graph.has_edge(node, node):
                cyclic_components.add(frozenset(component))
    tabulated = {
        frozenset(int(node) for node in row.cycle_nodes.split())
        for row in cycles.itertuples(index=False)
    }
    if cyclic_components != tabulated:
        raise AssertionError("trajectory cycles disagree with graph strongly connected components")
    return {
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "weakly_connected_components": nx.number_weakly_connected_components(graph),
        "strongly_connected_components": nx.number_strongly_connected_components(graph),
        "cyclic_strong_components": len(cyclic_components),
    }


def _predecessor_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = frame["first_result"].value_counts().reindex(range(10_000), fill_value=0)
    all_counts = pd.DataFrame(
        {
            "state": range(10_000),
            "state_string": [f"{state:04d}" for state in range(10_000)],
            "direct_predecessor_count": counts.to_numpy(dtype=int),
        }
    )
    all_counts["is_one_step_output"] = all_counts["direct_predecessor_count"] > 0
    top = (
        all_counts.sort_values(["direct_predecessor_count", "state"], ascending=[False, True])
        .head(20)
        .reset_index(drop=True)
    )
    top.insert(0, "rank", range(1, len(top) + 1))
    return all_counts, top


def _permutation_classes(frame: pd.DataFrame) -> pd.DataFrame:
    lookup = frame.set_index("start")
    rows: list[dict[str, Any]] = []
    for multiset, group in frame.groupby("digit_multiset", sort=True):
        first_results = group["first_result"].unique()
        if len(first_results) != 1:
            raise AssertionError(f"permutation class {multiset} has multiple first results")
        first_result = int(first_results[0])
        successor = lookup.loc[first_result]
        distances = sorted(group["distance_to_6174"].dropna().astype(int).unique())
        attractors = sorted(group["final_attractor_state"].unique())
        rows.append(
            {
                "digit_multiset": multiset,
                "class_size": len(group),
                "members": " ".join(group["start_state"]),
                "distinct_digits": int(group["distinct_digits"].iloc[0]),
                "is_repdigit_class": bool(group["is_repdigit"].all()),
                "x": int(group["x"].iloc[0]),
                "y": int(group["y"].iloc[0]),
                "first_result": first_result,
                "first_result_state": f"{first_result:04d}",
                "distance_after_first_to_6174": (
                    "" if pd.isna(successor["distance_to_6174"]) else int(successor["distance_to_6174"])
                ),
                "member_distance_values_to_6174": " ".join(map(str, distances)),
                "member_distance_min_to_6174": distances[0] if distances else "",
                "member_distance_max_to_6174": distances[-1] if distances else "",
                "all_members_same_distance": len(distances) <= 1,
                "member_iterations_to_attractor_min": int(group["iterations_to_attractor"].min()),
                "member_iterations_to_attractor_max": int(group["iterations_to_attractor"].max()),
                "final_attractors": " ".join(attractors),
            }
        )
    return pd.DataFrame(rows)


def _xy_table(frame: pd.DataFrame) -> pd.DataFrame:
    lookup = frame.set_index("start")
    rows: list[dict[str, Any]] = []
    for (x, y), group in frame.groupby(["x", "y"], sort=True):
        outputs = group["first_result"].unique()
        if len(outputs) != 1:
            raise AssertionError(f"digit-difference pair {(x, y)} has multiple outputs")
        output = int(outputs[0])
        output_row = lookup.loc[output]
        rows.append(
            {
                "x": int(x),
                "y": int(y),
                "input_state_count": len(group),
                "permutation_class_count": group["digit_multiset"].nunique(),
                "kaprekar_output": output,
                "kaprekar_output_state": f"{output:04d}",
                "output_distance_to_6174": (
                    "" if pd.isna(output_row["distance_to_6174"]) else int(output_row["distance_to_6174"])
                ),
                "output_final_attractor": output_row["final_attractor_state"],
            }
        )
    result = pd.DataFrame(rows)
    if len(result) != result["kaprekar_output"].nunique():
        raise AssertionError("(x,y) pairs do not map one-to-one to Kaprekar outputs")
    return result


def _trajectory_frequency(frame: pd.DataFrame, predecessor_counts: pd.DataFrame) -> pd.DataFrame:
    containing_counter: Counter[int] = Counter()
    post_first_counter: Counter[int] = Counter()
    valid = frame[~frame["is_repdigit"]]
    for text in valid["trajectory"]:
        path = _parse_path(str(text))
        containing_counter.update(path[:-1])
        # A fixed point has [n, n].  Otherwise discard only the closing repeat.
        post_first = path[1:] if len(path) == 2 else path[1:-1]
        post_first_counter.update(post_first)

    result = predecessor_counts.copy()
    result["valid_trajectories_containing_state"] = [containing_counter[state] for state in range(10_000)]
    result["valid_post_first_trajectory_occurrences"] = [post_first_counter[state] for state in range(10_000)]
    return result


def _one_step_outputs(frame: pd.DataFrame, predecessor_counts: pd.DataFrame) -> pd.DataFrame:
    output_rows = predecessor_counts[predecessor_counts["direct_predecessor_count"] > 0].copy()
    lookup = frame.set_index("start")
    output_rows["distance_to_6174"] = [lookup.loc[state, "distance_to_6174"] for state in output_rows["state"]]
    output_rows["final_attractor_state"] = [lookup.loc[state, "final_attractor_state"] for state in output_rows["state"]]
    return output_rows.reset_index(drop=True)


def generate_analysis(input_path: Path, tables_dir: Path, summary_path: Path) -> dict[str, Any]:
    """Create every summary table and a machine-readable analysis summary."""

    tables_dir.mkdir(parents=True, exist_ok=True)
    frame = load_results(input_path)
    valid = frame[~frame["is_repdigit"]].copy()
    ordinary_valid = valid[valid["start"] >= 1000].copy()

    overall = _overall_classification(frame)
    distribution = _iteration_distribution(frame)
    cycles = _cycle_table(frame)
    graph_summary = _validate_graph(frame, cycles)
    predecessor_counts, top_predecessors = _predecessor_tables(frame)
    permutation_classes = _permutation_classes(frame)
    xy_states = _xy_table(frame)
    trajectory_frequency = _trajectory_frequency(frame, predecessor_counts)
    one_step_outputs = _one_step_outputs(frame, predecessor_counts)

    summaries = pd.DataFrame(
        [
            _numeric_summary("all valid states 0000-9999", valid["distance_to_6174"].astype(int)),
            _numeric_summary(
                "valid ordinary integers 1000-9999", ordinary_valid["distance_to_6174"].astype(int)
            ),
        ]
    )

    max_depth = int(valid["distance_to_6174"].max())
    maximum_states = valid[valid["distance_to_6174"] == max_depth][
        [
            "start",
            "start_state",
            "digit_multiset",
            "x",
            "y",
            "first_result",
            "first_result_state",
            "distance_to_6174",
            "trajectory",
        ]
    ].copy()
    maximum_class_counts = (
        maximum_states.groupby("digit_multiset", as_index=False)
        .agg(
            state_count=("start", "size"),
            x=("x", "first"),
            y=("y", "first"),
            first_result_state=("first_result_state", "first"),
            example_state=("start_state", "first"),
            example_trajectory=("trajectory", "first"),
        )
        .sort_values("digit_multiset")
    )

    fixed_points = frame.loc[frame["start"] == frame["first_result"], "start_state"].tolist()
    valid_outputs = set(valid["first_result"].astype(int))
    valid_path_states = set(
        state
        for text in valid["trajectory"]
        for state in _parse_path(str(text))[:-1]
    )
    valid_post_first_states = set(
        trajectory_frequency.loc[
            trajectory_frequency["valid_post_first_trajectory_occurrences"] > 0, "state"
        ].astype(int)
    )

    reductions = pd.DataFrame(
        [
            {"stage": "all four-digit states", "count": len(frame), "notes": "0000 through 9999"},
            {
                "stage": "sorted-digit permutation classes",
                "count": frame["digit_multiset"].nunique(),
                "notes": "combinations with repetition of 10 digits taken 4 at a time",
            },
            {
                "stage": "digit-difference (x,y) pairs",
                "count": frame[["x", "y"]].drop_duplicates().shape[0],
                "notes": "0 <= y <= x <= 9",
            },
            {
                "stage": "unique outputs after one step",
                "count": frame["first_result"].nunique(),
                "notes": "includes 0000 from repdigits",
            },
            {
                "stage": "unique outputs from valid states",
                "count": len(valid_outputs),
                "notes": "excludes the repdigit output 0000",
            },
            {"stage": "attractors", "count": len(cycles), "notes": "cycles discovered exhaustively"},
        ]
    )

    table_map = {
        "overall_state_space.csv": overall,
        "iteration_distribution.csv": distribution,
        "summary_statistics.csv": summaries,
        "maximum_distance_states.csv": maximum_states,
        "maximum_distance_permutation_classes.csv": maximum_class_counts,
        "cycles_and_basins.csv": cycles,
        "predecessor_counts.csv": predecessor_counts,
        "top_predecessor_states.csv": top_predecessors,
        "state_space_reduction.csv": reductions,
        "permutation_classes.csv": permutation_classes,
        "xy_reduced_states.csv": xy_states,
        "trajectory_state_frequency.csv": trajectory_frequency,
        "unique_one_step_outputs.csv": one_step_outputs,
    }
    for filename, table in table_map.items():
        _write_csv(table, tables_dir / filename)

    valid_summary = _numeric_summary("all valid states 0000-9999", valid["distance_to_6174"].astype(int))
    ordinary_summary = _numeric_summary(
        "valid ordinary integers 1000-9999", ordinary_valid["distance_to_6174"].astype(int)
    )
    full_distribution = distribution[distribution["scope"] == "all valid states 0000-9999"]
    ordinary_distribution = distribution[
        distribution["scope"] == "valid ordinary integers 1000-9999"
    ]
    summary: dict[str, Any] = {
        "total_states": len(frame),
        "ordinary_integer_states": int((frame["start"] >= 1000).sum()),
        "valid_states": len(valid),
        "ordinary_valid_states": len(ordinary_valid),
        "repdigits": int(frame["is_repdigit"].sum()),
        "reaches_6174": int(frame["reaches_6174"].sum()),
        "reaches_0000": int(frame["reaches_0000"].sum()),
        "other_cycle_basin_size": int(frame["other_cycle_detected"].sum()),
        "fixed_points": fixed_points,
        "cycles": cycles.to_dict(orient="records"),
        "graph": graph_summary,
        "valid_iteration_statistics": valid_summary,
        "ordinary_valid_iteration_statistics": ordinary_summary,
        "valid_iteration_distribution": {
            str(int(row.transformations_to_6174)): int(row.count)
            for row in full_distribution.itertuples(index=False)
        },
        "ordinary_valid_iteration_distribution": {
            str(int(row.transformations_to_6174)): int(row.count)
            for row in ordinary_distribution.itertuples(index=False)
        },
        "maximum_distance": max_depth,
        "maximum_distance_state_count": len(maximum_states),
        "maximum_distance_permutation_class_count": len(maximum_class_counts),
        "permutation_class_count": len(permutation_classes),
        "valid_permutation_class_count": int((~permutation_classes["is_repdigit_class"]).sum()),
        "permutation_classes_with_nonuniform_start_distance": int(
            (~permutation_classes["all_members_same_distance"]).sum()
        ),
        "xy_pair_count": len(xy_states),
        "unique_one_step_outputs": int(frame["first_result"].nunique()),
        "unique_valid_one_step_outputs": len(valid_outputs),
        "unique_states_in_valid_trajectories_including_starts": len(valid_path_states),
        "unique_post_first_states_in_valid_trajectories": len(valid_post_first_states),
        "top_predecessors": top_predecessors.head(10).to_dict(orient="records"),
    }

    # Convert NumPy scalar values before serializing.
    def json_default(value: Any) -> Any:
        if hasattr(value, "item"):
            return value.item()
        if isinstance(value, float) and math.isnan(value):
            return None
        raise TypeError(f"cannot serialize {type(value)}")

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, default=json_default) + "\n", encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/kaprekar_results.csv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("tables"))
    parser.add_argument("--summary", type=Path, default=Path("data/analysis_summary.json"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = generate_analysis(args.input, args.tables_dir, args.summary)
    print(
        "Generated analysis for "
        f"{summary['total_states']:,} states; maximum valid distance = {summary['maximum_distance']}"
    )


if __name__ == "__main__":
    main()
