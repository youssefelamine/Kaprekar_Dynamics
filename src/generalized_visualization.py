"""Comparative figures for generalized Kaprekar systems."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"{stem}.{suffix}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _heatmap(
    summary: pd.DataFrame,
    column: str,
    title: str,
    colorbar_label: str,
    output_dir: Path,
    stem: str,
    *,
    value_format: str = ".0f",
    cmap: str = "viridis",
) -> None:
    pivot = summary.pivot(index="digits", columns="base", values=column).sort_index()
    values = pivot.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(13.0, 5.2))
    image = ax.imshow(values, aspect="auto", cmap=cmap)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label(colorbar_label)
    ax.set_xticks(range(len(pivot.columns)), pivot.columns)
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    ax.set_xlabel("Numeric base")
    ax.set_ylabel("Fixed digit width")
    ax.set_title(title, weight="bold", pad=12)

    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values[row_index, column_index]
            red, green, blue, _ = image.cmap(image.norm(value))
            luminance = 0.299 * red + 0.587 * green + 0.114 * blue
            color = "#111111" if luminance > 0.58 else "white"
            ax.text(
                column_index,
                row_index,
                format(value, value_format),
                ha="center",
                va="center",
                fontsize=7.5,
                color=color,
                weight="bold",
            )
    _save(fig, output_dir, stem)


def _base_ten_comparison(summary: pd.DataFrame, output_dir: Path) -> None:
    decimal = summary[summary["base"] == 10].sort_values("digits")
    if decimal.empty:
        return
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), sharex=True)
    digits = decimal["digits"]

    axes[0, 0].bar(digits, decimal["valid_attractor_count"], color="#1769aa")
    axes[0, 0].set_ylabel("Valid attractors")
    axes[0, 0].set_title("Attractor count")

    axes[0, 1].plot(digits, decimal["maximum_transient_depth_valid"], marker="o", color="#6a1b9a")
    axes[0, 1].set_ylabel("Transformations")
    axes[0, 1].set_title("Maximum valid transient depth")

    axes[1, 0].plot(digits, decimal["valid_mean_transient_depth"], marker="o", color="#00897b")
    axes[1, 0].set_ylabel("Transformations")
    axes[1, 0].set_title("Mean valid transient depth")

    axes[1, 1].bar(digits, decimal["largest_valid_basin_percentage"], color="#ef6c00")
    axes[1, 1].set_ylabel("Percent of valid states")
    axes[1, 1].set_ylim(0, 105)
    axes[1, 1].set_title("Largest valid basin share")

    for ax in axes.flat:
        ax.set_xticks(digits)
        ax.set_xlabel("Decimal digit width")
        ax.grid(axis="y", alpha=0.3)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Base-10 Kaprekar dynamics across widths 2–6", fontsize=16, weight="bold")
    fig.tight_layout()
    _save(fig, output_dir, "figure_g6_base10_width_comparison")


def generate_generalized_figures(system_summary_path: Path, output_dir: Path) -> None:
    """Generate comparative PNG and PDF figures from the system summary."""

    summary = pd.read_csv(system_summary_path)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 14})
    _heatmap(
        summary,
        "attractor_count",
        "Kaprekar attractors across bases and digit widths",
        "Number of terminal cycles",
        output_dir,
        "figure_g1_attractor_count_heatmap",
        cmap="magma_r",
    )
    _heatmap(
        summary,
        "nontrivial_cycle_count",
        "Non-trivial cycle counts across generalized systems",
        "Cycles of length greater than one",
        output_dir,
        "figure_g2_nontrivial_cycle_heatmap",
        cmap="plasma_r",
    )
    _heatmap(
        summary,
        "maximum_transient_depth_valid",
        "Maximum valid transient depth",
        "Kaprekar transformations",
        output_dir,
        "figure_g3_maximum_depth_heatmap",
        cmap="viridis_r",
    )
    _heatmap(
        summary,
        "maximum_cycle_length",
        "Maximum terminal cycle length",
        "States in longest cycle",
        output_dir,
        "figure_g4_maximum_cycle_length_heatmap",
        cmap="cividis_r",
    )
    _heatmap(
        summary,
        "largest_valid_basin_percentage",
        "Share of valid states in the largest basin",
        "Percent of valid states",
        output_dir,
        "figure_g5_largest_basin_share_heatmap",
        value_format=".1f",
        cmap="YlGnBu",
    )
    _base_ten_comparison(summary, output_dir)
