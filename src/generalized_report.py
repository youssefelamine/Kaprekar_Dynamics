"""Write the generalized scientific report from generated census tables."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def _systems_list(frame: pd.DataFrame, columns: list[str], limit: int = 10) -> str:
    selected = frame.sort_values(["base", "digits"]).head(limit)
    return ", ".join(
        f"base {int(row.base)}, width {int(row.digits)}"
        for row in selected.itertuples(index=False)
    )


def _base_ten_table(summary: pd.DataFrame) -> str:
    rows = [
        "| Width | Valid states | Valid attractors | Non-trivial cycles | Maximum depth | Mean depth | Largest basin |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary[summary["base"] == 10].sort_values("digits").itertuples(index=False):
        rows.append(
            f"| {int(row.digits)} | {int(row.valid_state_count):,} | "
            f"{int(row.valid_attractor_count)} | {int(row.nontrivial_cycle_count)} | "
            f"{int(row.maximum_transient_depth_valid)} | {row.valid_mean_transient_depth:.4f} | "
            f"{row.largest_valid_basin_percentage:.2f}% |"
        )
    return "\n".join(rows)


def write_generalized_report(
    summary_path: Path,
    system_summary_path: Path,
    report_path: Path,
) -> None:
    """Create a report whose numerical statements are read from generated outputs."""

    aggregate = json.loads(summary_path.read_text(encoding="utf-8"))
    systems = pd.read_csv(system_summary_path)
    unanimous_count = int(systems["unanimous_valid_attractor"].sum())
    nontrivial_systems = systems[systems["nontrivial_cycle_count"] > 0]
    maximum_attractors = int(systems["attractor_count"].max())
    max_attractor_systems = systems[systems["attractor_count"] == maximum_attractors]
    maximum_depth = int(systems["maximum_transient_depth_valid"].max())
    max_depth_systems = systems[systems["maximum_transient_depth_valid"] == maximum_depth]
    maximum_cycle_length = int(systems["maximum_cycle_length"].max())
    max_cycle_systems = systems[systems["maximum_cycle_length"] == maximum_cycle_length]

    report = rf"""# Generalized Kaprekar Dynamics Across Bases 2–16 and Widths 2–6

## Abstract

This extension compares {aggregate['system_count']} fixed-width Kaprekar systems formed by bases {min(aggregate['base_values'])}–{max(aggregate['base_values'])} and digit widths {min(aggregate['digit_widths'])}–{max(aggregate['digit_widths'])}. Exact symmetry weighting represents {aggregate['ordered_states_represented']:,} ordered states with {aggregate['weighted_class_records']:,} sorted-digit multiset classes; no sampling is used. For each system, the experiment discovers every terminal cycle, computes weighted basin sizes and transient-depth distributions, and distinguishes repdigits from valid non-repdigit states. The results show that the single-attractor behavior familiar from decimal width four is not universal: {unanimous_count} of {aggregate['system_count']} systems have one valid attractor, while {len(nontrivial_systems)} contain at least one non-trivial cycle. The largest observed attractor count is {maximum_attractors}, the largest valid transient depth is {maximum_depth}, and the longest terminal cycle contains {maximum_cycle_length} states. These values are generated from the exported census and should not be generalized beyond the surveyed grid.

## 1. Research design

For base $b$ and width $d$, the state space contains $b^d$ zero-padded strings. Exactly $b$ are repdigits, leaving $b^d-b$ valid states. Direct enumeration over the selected grid would visit {aggregate['ordered_states_represented']:,} states. Instead, the experiment groups strings by sorted digit multiset. A multiset with digit multiplicities $m_i$ represents exactly

```text
d! / (m₀! m₁! ⋯ mₖ!).
```

ordered strings, all of which have the same first successor.

The analyzer enumerates every multiset, constructs the closed functional graph of unique outputs, discovers cycles generically, and propagates terminal cycles and depths backward. Basin statistics are then recovered exactly from multinomial weights. If a class contains a cycle node, that one ordering is counted at depth zero and the remaining permutations are counted after their shared first edge. This prevents quotienting from hiding the distinction between a fixed or periodic ordering and its non-periodic permutations.

## 2. Survey-wide results

- Systems analyzed: **{aggregate['system_count']}**.
- Ordered states represented: **{aggregate['ordered_states_represented']:,}**.
- Weighted permutation classes: **{aggregate['weighted_class_records']:,}**.
- Attractor records discovered: **{aggregate['attractor_records']:,}**.
- Systems with one valid attractor: **{unanimous_count}**.
- Systems containing non-trivial cycles: **{len(nontrivial_systems)}**.
- Largest attractor count: **{maximum_attractors}**, attained by {_systems_list(max_attractor_systems, ['base', 'digits'])}.
- Largest valid transient depth: **{maximum_depth}**, attained by {_systems_list(max_depth_systems, ['base', 'digits'])}.
- Longest terminal cycle: **{maximum_cycle_length} states**, attained by {_systems_list(max_cycle_systems, ['base', 'digits'])}.

![Generalized attractor count](../figures/generalized/figure_g1_attractor_count_heatmap.png)

**Figure G1.** Number of terminal cycles, including the zero fixed point, in every surveyed system.

![Non-trivial cycles](../figures/generalized/figure_g2_nontrivial_cycle_heatmap.png)

**Figure G2.** Number of cycles of length greater than one.

![Maximum transient depth](../figures/generalized/figure_g3_maximum_depth_heatmap.png)

**Figure G3.** Maximum transformations required for a valid state to enter its terminal cycle.

![Maximum cycle length](../figures/generalized/figure_g4_maximum_cycle_length_heatmap.png)

**Figure G4.** Length of the longest terminal cycle in each system.

![Largest basin share](../figures/generalized/figure_g5_largest_basin_share_heatmap.png)

**Figure G5.** Percentage of valid ordered states belonging to the largest valid basin.

## 3. Decimal widths 2–6

{_base_ten_table(systems)}

![Base-10 comparison](../figures/generalized/figure_g6_base10_width_comparison.png)

**Figure G6.** Attractor and convergence statistics for decimal widths two through six. Width four is unusual in having unanimous convergence of valid states to the fixed point 6174; neighboring widths must be described in terms of multiple cycles or fixed points when the census reports them.

## 4. Finite proof certificate for 6174

The companion [finite proof chapter](6174_finite_proof.md) proves the reduction from 10,000 four-digit decimal states to 55 digit-difference pairs and checks their complete transition graph. The 54 nonzero pairs all reach `(6,2)`, the pair representing 6174. Their maximum pair distance is six, which yields the sharp seven-transformation bound for original states.

## 5. Interpretation

The survey demonstrates that a Kaprekar “constant” is not the generic outcome of changing base or width. Every finite system eventually reaches a cycle, but the number, length, and basin balance of those cycles vary. The largest-basin heatmap distinguishes unanimous systems from systems in which several attractors compete, while the depth heatmap separates basin structure from convergence speed.

Symmetry weighting is exact because digit order is destroyed by the first transformation. It is also what makes the broader comparison practical: {aggregate['weighted_class_records']:,} classes replace {aggregate['ordered_states_represented']:,} ordered starts without changing any basin total or depth frequency.

## 6. Limitations

- The census is complete only for the selected bases and widths.
- Exact class weighting depends on the standard fixed-width sort-and-subtract convention.
- Heatmaps summarize systems and cannot display every cycle; the attractor CSV is the authoritative detailed record.
- The proof certificate establishes the classical decimal four-digit theorem, not a universal theorem for the generalized grid.

## 7. Reproducibility and data

Run `python3 -m src.generalized_pipeline` from the project root. The pipeline regenerates:

- [`weighted_classes.csv`](../data/generalized/weighted_classes.csv): every weighted digit-multiset class;
- [`generalized_summary.json`](../data/generalized/generalized_summary.json): aggregate and per-system summaries;
- [`system_summary.csv`](../tables/generalized/system_summary.csv): one row per surveyed system;
- [`attractors_and_basins.csv`](../tables/generalized/attractors_and_basins.csv): every discovered terminal cycle;
- [`depth_distributions.csv`](../tables/generalized/depth_distributions.csv): exact weighted depth counts;
- [`kaprekar_6174_pair_certificate.csv`](../tables/generalized/kaprekar_6174_pair_certificate.csv): the checked 55-pair proof table.

## References

The historical and mathematical bibliography is shared with the [classical 6174 report](kaprekar_report.md), particularly Devlin and Zeng on four-digit base-dependent maximum distances and Kay and Downes-Ward on fixed points and cycles in generalized bases.
"""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
