"""Publication-quality visualizations for the Kaprekar experiment."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Hashable, Iterable

# Keep Matplotlib's cache inside the writable, disposable project tree.
os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

from .analysis import load_results


BLUE = "#1769aa"
LIGHT_BLUE = "#64b5f6"
ORANGE = "#ef6c00"
DARK = "#263238"
GRID = "#cfd8dc"


def configure_style() -> None:
    """Apply a restrained, report-friendly Matplotlib style."""

    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 14,
            "axes.labelsize": 11,
            "axes.titleweight": "bold",
            "axes.edgecolor": DARK,
            "axes.labelcolor": DARK,
            "xtick.color": DARK,
            "ytick.color": DARK,
            "grid.color": GRID,
            "grid.alpha": 0.55,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, figures_dir: Path, stem: str) -> None:
    """Save one figure as a high-resolution PNG and a vector PDF."""

    figures_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(figures_dir / f"{stem}.{suffix}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _ring_positions(
    identifiers: Iterable[Hashable],
    basins: Iterable[str],
    depths: Iterable[int],
    *,
    second_basin_x: float = 10.0,
) -> dict[Hashable, tuple[float, float]]:
    records = sorted(zip(identifiers, basins, depths), key=lambda row: (row[1], row[2], str(row[0])))
    grouped: dict[tuple[str, int], list[Hashable]] = {}
    for identifier, basin, depth in records:
        grouped.setdefault((str(basin), int(depth)), []).append(identifier)

    positions: dict[Hashable, tuple[float, float]] = {}
    for (basin, depth), group in grouped.items():
        center_x = 0.0 if basin == "6174" else second_basin_x
        if depth == 0:
            for identifier in group:
                positions[identifier] = (center_x, 0.0)
            continue
        radius = float(depth) if basin == "6174" else 1.25 * depth
        offset = 0.43 * depth
        for index, identifier in enumerate(group):
            angle = 2 * math.pi * index / len(group) + offset
            positions[identifier] = (center_x + radius * math.cos(angle), radius * math.sin(angle))
    return positions


def figure_iteration_distribution(frame: pd.DataFrame, figures_dir: Path) -> None:
    valid = frame[~frame["is_repdigit"]]
    counts = valid["distance_to_6174"].astype(int).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    bars = ax.bar(counts.index, counts.values, color=BLUE, edgecolor="white", width=0.78)
    ax.bar_label(bars, labels=[f"{value:,}" for value in counts.values], padding=3, fontsize=9)
    ax.set_xlabel("Kaprekar transformations required to first reach 6174")
    ax.set_ylabel("Number of valid starting states")
    ax.set_title("Iteration distribution for all 9,990 valid four-digit states")
    ax.set_xticks(range(int(counts.index.max()) + 1))
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, max(counts.values) * 1.15)
    fig.text(
        0.5,
        -0.01,
        "Distance counts transformations; the initial state 6174 has distance zero.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_1_iteration_distribution")


def figure_full_graph(frame: pd.DataFrame, figures_dir: Path) -> None:
    positions = _ring_positions(
        frame["start"].astype(int),
        frame["final_attractor_state"].astype(str),
        frame["iterations_to_attractor"].astype(int),
    )
    segments_6174: list[list[tuple[float, float]]] = []
    segments_0000: list[list[tuple[float, float]]] = []
    for row in frame.itertuples(index=False):
        if int(row.start) == int(row.first_result):
            continue
        segment = [positions[int(row.start)], positions[int(row.first_result)]]
        (segments_6174 if row.final_attractor_state == "6174" else segments_0000).append(segment)

    fig, ax = plt.subplots(figsize=(11.5, 8.0))
    ax.add_collection(LineCollection(segments_6174, colors=BLUE, linewidths=0.16, alpha=0.07, zorder=1))
    ax.add_collection(LineCollection(segments_0000, colors=ORANGE, linewidths=0.8, alpha=0.65, zorder=2))
    valid = frame[frame["final_attractor_state"] == "6174"]
    repdigits = frame[frame["final_attractor_state"] == "0000"]
    ax.scatter(
        [positions[int(state)][0] for state in valid["start"]],
        [positions[int(state)][1] for state in valid["start"]],
        s=1.8,
        c=valid["iterations_to_attractor"],
        cmap="Blues",
        alpha=0.65,
        linewidths=0,
        zorder=3,
    )
    ax.scatter(
        [positions[int(state)][0] for state in repdigits["start"]],
        [positions[int(state)][1] for state in repdigits["start"]],
        s=20,
        color=ORANGE,
        edgecolor="white",
        linewidth=0.4,
        zorder=4,
    )
    ax.scatter(*positions[6174], s=130, color="#0d47a1", edgecolor="white", zorder=6)
    ax.scatter(*positions[0], s=130, color="#e65100", edgecolor="white", zorder=6)
    ax.text(*positions[6174], "6174", ha="center", va="center", color="white", fontsize=8, weight="bold", zorder=7)
    ax.text(*positions[0], "0000", ha="center", va="center", color="white", fontsize=8, weight="bold", zorder=7)
    ax.set_title("Full 10,000-node Kaprekar state-transition graph")
    ax.text(0, -7.8, "Basin of 6174: 9,990 states", ha="center", color=BLUE, weight="bold")
    ax.text(10, -1.8, "Basin of 0000: 10 states", ha="center", color=ORANGE, weight="bold")
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.margins(0.08)
    ax.axis("off")
    fig.text(
        0.5,
        0.015,
        "Nodes are placed on rings by distance to their attractor; faint directed links show deterministic successors.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_2a_full_state_graph")


def figure_reduced_graph(frame: pd.DataFrame, figures_dir: Path) -> None:
    output_states = sorted(frame["first_result"].astype(int).unique())
    reduced = frame.set_index("start").loc[output_states].reset_index()
    positions = _ring_positions(
        reduced["start"].astype(int),
        reduced["final_attractor_state"].astype(str),
        reduced["iterations_to_attractor"].astype(int),
        second_basin_x=8.0,
    )
    fig, ax = plt.subplots(figsize=(10.5, 7.5))
    for row in reduced.itertuples(index=False):
        if int(row.start) == int(row.first_result):
            continue
        ax.add_patch(
            FancyArrowPatch(
                positions[int(row.start)],
                positions[int(row.first_result)],
                arrowstyle="-|>",
                mutation_scale=6,
                shrinkA=4.5,
                shrinkB=5.5,
                color="#78909c",
                linewidth=0.75,
                alpha=0.65,
                zorder=1,
            )
        )
    valid = reduced[reduced["final_attractor_state"] == "6174"]
    exceptional = reduced[reduced["final_attractor_state"] == "0000"]
    scatter = ax.scatter(
        [positions[int(state)][0] for state in valid["start"]],
        [positions[int(state)][1] for state in valid["start"]],
        c=valid["iterations_to_attractor"],
        cmap="viridis",
        s=85,
        edgecolor="white",
        linewidth=0.7,
        zorder=3,
    )
    ax.scatter(
        [positions[int(state)][0] for state in exceptional["start"]],
        [positions[int(state)][1] for state in exceptional["start"]],
        color=ORANGE,
        s=100,
        edgecolor="white",
        zorder=3,
    )
    for row in reduced.itertuples(index=False):
        x_pos, y_pos = positions[int(row.start)]
        ax.annotate(f"{int(row.start):04d}", (x_pos, y_pos), xytext=(3, 2), textcoords="offset points", fontsize=6)
    colorbar = fig.colorbar(scatter, ax=ax, shrink=0.72, pad=0.02)
    colorbar.set_label("Distance to attractor")
    ax.set_title("Reduced graph of the 55 possible one-step outputs")
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.margins(0.12)
    ax.axis("off")
    fig.text(
        0.5,
        0.015,
        "Every one of the 10,000 inputs enters this closed 55-state image after one transformation.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_2b_reduced_55_state_graph")


def figure_grouped_distance_graph(frame: pd.DataFrame, figures_dir: Path) -> None:
    grouped = (
        frame.groupby(["final_attractor_state", "iterations_to_attractor"]).size().rename("count").reset_index()
    )
    fig, ax = plt.subplots(figsize=(13.5, 5.6))

    def draw_basin(attractor: str, y: float, color: str, start_x: float) -> None:
        basin = grouped[grouped["final_attractor_state"] == attractor].set_index("iterations_to_attractor")
        depths = sorted(basin.index.astype(int), reverse=True)
        coordinates: list[float] = []
        for index, depth in enumerate(depths):
            count = int(basin.loc[depth, "count"])
            x = start_x + 1.48 * index
            coordinates.append(x)
            box = FancyBboxPatch(
                (x - 0.55, y - 0.40),
                1.10,
                0.80,
                boxstyle="round,pad=0.04,rounding_size=0.18",
                facecolor=color,
                edgecolor="white",
                linewidth=1.5,
                zorder=3,
            )
            ax.add_patch(box)
            label = f"depth {depth}\n{count:,} state{'s' if count != 1 else ''}"
            ax.text(x, y, label, ha="center", va="center", color="white", fontsize=8.5, weight="bold", zorder=4)
        for x1, x2 in zip(coordinates, coordinates[1:]):
            ax.add_patch(
                FancyArrowPatch(
                    (x1 + 0.56, y),
                    (x2 - 0.56, y),
                    arrowstyle="-|>",
                    mutation_scale=15,
                    linewidth=1.5,
                    color=DARK,
                    alpha=0.75,
                    zorder=2,
                )
            )

    draw_basin("6174", 1.05, BLUE, 0.0)
    draw_basin("0000", -1.05, ORANGE, 8.88)
    ax.text(-0.78, 1.05, "6174 basin", ha="right", va="center", color=BLUE, weight="bold", fontsize=11)
    ax.text(8.1, -1.05, "0000 basin", ha="right", va="center", color=ORANGE, weight="bold", fontsize=11)
    ax.set_xlim(-1.8, 11.1)
    ax.set_ylim(-2.15, 2.15)
    ax.set_title("State-transition graph aggregated by distance to each attractor")
    ax.axis("off")
    fig.text(
        0.5,
        0.035,
        "Every non-fixed edge reduces the remaining basin distance by exactly one transformation.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_2c_graph_grouped_by_distance")


def figure_permutation_graph(frame: pd.DataFrame, tables_dir: Path, figures_dir: Path) -> None:
    classes = pd.read_csv(
        tables_dir / "permutation_classes.csv",
        dtype={"digit_multiset": "string", "first_result_state": "string", "final_attractors": "string"},
    )
    classes["digit_multiset"] = classes["digit_multiset"].str.zfill(4)
    classes["first_result_state"] = classes["first_result_state"].str.zfill(4)
    classes["successor_class"] = classes["first_result_state"].map(lambda text: "".join(sorted(str(text))))
    classes["class_depth"] = classes["member_iterations_to_attractor_min"].astype(int)
    positions = _ring_positions(
        classes["digit_multiset"],
        classes["final_attractors"].astype(str),
        classes["class_depth"],
        second_basin_x=9.0,
    )
    segments = [
        [positions[str(row.digit_multiset)], positions[str(row.successor_class)]]
        for row in classes.itertuples(index=False)
        if str(row.digit_multiset) != str(row.successor_class)
    ]
    fig, ax = plt.subplots(figsize=(11.0, 7.8))
    ax.add_collection(LineCollection(segments, colors="#607d8b", linewidths=0.28, alpha=0.14, zorder=1))
    valid = classes[classes["final_attractors"] == "6174"]
    repdigit = classes[classes["final_attractors"] == "0000"]
    ax.scatter(
        [positions[str(key)][0] for key in valid["digit_multiset"]],
        [positions[str(key)][1] for key in valid["digit_multiset"]],
        s=8 + valid["class_size"] * 0.7,
        c=valid["class_depth"],
        cmap="Blues",
        alpha=0.8,
        linewidth=0,
        zorder=3,
    )
    ax.scatter(
        [positions[str(key)][0] for key in repdigit["digit_multiset"]],
        [positions[str(key)][1] for key in repdigit["digit_multiset"]],
        s=22,
        color=ORANGE,
        edgecolor="white",
        linewidth=0.3,
        zorder=3,
    )
    ax.scatter(*positions["1467"], s=110, color="#0d47a1", edgecolor="white", zorder=5)
    ax.scatter(*positions["0000"], s=110, color="#e65100", edgecolor="white", zorder=5)
    ax.annotate("class 1467\n(contains 6174)", positions["1467"], xytext=(8, 8), textcoords="offset points", fontsize=8)
    ax.annotate("class 0000", positions["0000"], xytext=(8, 8), textcoords="offset points", fontsize=8)
    ax.set_title("Permutation-class quotient graph (715 sorted-digit multisets)")
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.margins(0.1)
    ax.axis("off")
    fig.text(
        0.5,
        0.015,
        "Node area reflects class size; class edges apply one Kaprekar step and then re-sort the output digits.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_2d_permutation_class_graph")


def figure_basin_depth(frame: pd.DataFrame, figures_dir: Path) -> None:
    depths = range(int(frame["iterations_to_attractor"].max()) + 1)
    valid_counts = (
        frame[frame["final_attractor_state"] == "6174"]["iterations_to_attractor"]
        .value_counts()
        .reindex(depths, fill_value=0)
    )
    zero_counts = (
        frame[frame["final_attractor_state"] == "0000"]["iterations_to_attractor"]
        .value_counts()
        .reindex(depths, fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(9.3, 5.5))
    width = 0.38
    x = np.arange(len(list(depths)))
    bars_valid = ax.bar(x - width / 2, valid_counts.values, width, color=BLUE, label="6174 basin")
    bars_zero = ax.bar(x + width / 2, zero_counts.values, width, color=ORANGE, label="0000 basin")
    ax.bar_label(bars_valid, labels=[f"{v:,}" if v else "" for v in valid_counts.values], padding=2, fontsize=8)
    ax.bar_label(bars_zero, labels=[f"{v:,}" if v else "" for v in zero_counts.values], padding=2, fontsize=8)
    ax.set_xticks(x, list(depths))
    ax.set_xlabel("Graph distance to the basin attractor")
    ax.set_ylabel("Number of states")
    ax.set_title("Basin depth across the complete 10,000-state graph")
    ax.legend(frameon=False)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, max(valid_counts.values) * 1.15)
    save_figure(fig, figures_dir, "figure_3_basin_depth")


def figure_first_step_compression(tables_dir: Path, figures_dir: Path) -> None:
    reductions = pd.read_csv(tables_dir / "state_space_reduction.csv")
    selected_names = [
        "all four-digit states",
        "sorted-digit permutation classes",
        "digit-difference (x,y) pairs",
        "unique outputs after one step",
        "attractors",
    ]
    selected = reductions.set_index("stage").loc[selected_names].reset_index()
    labels = [
        "All states",
        "Digit-multiset classes",
        "Feasible (x, y) pairs",
        "One-step outputs",
        "Attractors",
    ]
    fig, ax = plt.subplots(figsize=(9.8, 5.8))
    y = np.arange(len(selected))
    colors = [DARK, "#455a64", LIGHT_BLUE, BLUE, ORANGE]
    bars = ax.barh(y, selected["count"], color=colors, edgecolor="white", height=0.65)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Number of states or structural classes (logarithmic scale)")
    ax.set_title("First-step compression of the four-digit state space")
    ax.grid(axis="x", which="both")
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    for bar, count in zip(bars, selected["count"]):
        ax.text(float(count) * 1.12, bar.get_y() + bar.get_height() / 2, f"{int(count):,}", va="center", weight="bold")
    ax.set_xlim(1, 25_000)
    fig.text(
        0.5,
        0.015,
        "The 55 feasible difference pairs map one-to-one to the 55 possible first-step outputs.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_4_first_step_compression")


def figure_xy_space(tables_dir: Path, figures_dir: Path) -> None:
    xy = pd.read_csv(tables_dir / "xy_reduced_states.csv", keep_default_na=True)
    convergent = xy[xy["output_final_attractor"].astype(str).str.zfill(4) == "6174"].copy()
    exceptional = xy.drop(convergent.index)
    fig, ax = plt.subplots(figsize=(8.2, 6.7))
    sizes = 35 + 5.5 * np.sqrt(convergent["input_state_count"])
    scatter = ax.scatter(
        convergent["x"],
        convergent["y"],
        c=convergent["output_distance_to_6174"],
        s=sizes,
        cmap="viridis_r",
        norm=Normalize(vmin=0, vmax=int(convergent["output_distance_to_6174"].max())),
        edgecolor="white",
        linewidth=0.8,
        alpha=0.95,
    )
    if not exceptional.empty:
        ax.scatter(
            exceptional["x"],
            exceptional["y"],
            s=95,
            marker="X",
            color=ORANGE,
            edgecolor="white",
            linewidth=0.8,
            label="(0, 0) maps to 0000",
            zorder=4,
        )
    for row in xy.itertuples(index=False):
        ax.annotate(str(row.kaprekar_output_state).zfill(4), (row.x, row.y), xytext=(3, 3), textcoords="offset points", fontsize=6)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("Distance of K(n) to 6174")
    ax.set_xticks(range(10))
    ax.set_yticks(range(10))
    ax.set_xlim(-0.55, 9.55)
    ax.set_ylim(-0.55, 9.55)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x=a-d$")
    ax.set_ylabel(r"$y=b-c$")
    ax.set_title(r"Reduced digit-difference state space: $0\leq y\leq x\leq 9$")
    ax.grid()
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left")
    fig.text(
        0.5,
        0.01,
        "Labels show K(n); marker area increases with the number of original states represented.",
        ha="center",
        fontsize=9,
        color=DARK,
    )
    save_figure(fig, figures_dir, "figure_5_xy_reduced_state_space")


def figure_predecessors(tables_dir: Path, figures_dir: Path) -> None:
    top = pd.read_csv(tables_dir / "top_predecessor_states.csv", dtype={"state_string": "string"})
    top["state_string"] = top["state_string"].str.zfill(4)
    plot_data = top.iloc[::-1]
    fig, ax = plt.subplots(figsize=(9.2, 7.0))
    colors = [ORANGE if state == "6174" else BLUE for state in plot_data["state_string"]]
    bars = ax.barh(plot_data["state_string"], plot_data["direct_predecessor_count"], color=colors, edgecolor="white")
    ax.bar_label(bars, labels=[f"{value:,}" for value in plot_data["direct_predecessor_count"]], padding=3, fontsize=8)
    ax.set_xlabel("Number of direct predecessors among 10,000 states")
    ax.set_ylabel("Kaprekar output state")
    ax.set_title("States with the largest direct predecessor counts")
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, max(plot_data["direct_predecessor_count"]) * 1.14)
    ax.legend(
        handles=[Patch(facecolor=BLUE, label="Transient state"), Patch(facecolor=ORANGE, label="Fixed point 6174")],
        frameon=False,
        loc="lower right",
    )
    save_figure(fig, figures_dir, "figure_6_predecessor_distribution")


def generate_figures(input_path: Path, tables_dir: Path, figures_dir: Path) -> None:
    """Generate all required state-space figures in PNG and PDF formats."""

    configure_style()
    frame = load_results(input_path)
    figure_iteration_distribution(frame, figures_dir)
    figure_full_graph(frame, figures_dir)
    figure_reduced_graph(frame, figures_dir)
    figure_grouped_distance_graph(frame, figures_dir)
    figure_permutation_graph(frame, tables_dir, figures_dir)
    figure_basin_depth(frame, figures_dir)
    figure_first_step_compression(tables_dir, figures_dir)
    figure_xy_space(tables_dir, figures_dir)
    figure_predecessors(tables_dir, figures_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/kaprekar_results.csv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_figures(args.input, args.tables_dir, args.figures_dir)
    print(f"Generated PNG and PDF figures in {args.figures_dir}")


if __name__ == "__main__":
    main()
