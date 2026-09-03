"""Validation tests for generalized weighted Kaprekar analysis and proof."""

from __future__ import annotations

import csv
import math
import unittest
from collections import Counter
from pathlib import Path

from src.generalized_analysis import analyze_system, multiset_weight
from src.generalized_pipeline import parse_integer_selection
from src.kaprekar import (
    analyze_general_state,
    digits_in_base,
    format_state_general,
    integer_from_digits,
    kaprekar_components_general,
    kaprekar_step_general,
    trajectory_general,
)
from src.proof_certificate import EXPECTED_MAXIMUM_PAIRS, generate_pair_certificate


def independent_reference_step(n: int, digits: int, base: int) -> int:
    values = [0] * digits
    remaining = n
    for index in range(digits - 1, -1, -1):
        remaining, values[index] = divmod(remaining, base)
    low = 0
    high = 0
    for digit in sorted(values):
        low = low * base + digit
    for digit in sorted(values, reverse=True):
        high = high * base + digit
    return high - low


class GeneralizedCoreTests(unittest.TestCase):
    def test_fixed_width_base_conversion_and_formatting(self) -> None:
        self.assertEqual(digits_in_base(0xAF, digits=4, base=16), (0, 0, 10, 15))
        self.assertEqual(format_state_general(0xAF, digits=4, base=16), "00AF")
        self.assertEqual(integer_from_digits((0, 0, 10, 15), base=16), 0xAF)

    def test_components_and_leading_zeros(self) -> None:
        self.assertEqual(kaprekar_components_general(1, digits=4, base=2), (8, 1, 7))
        self.assertEqual(format_state_general(7, digits=4, base=2), "0111")

    def test_reference_agreement_on_representative_systems(self) -> None:
        for base, digits in ((2, 6), (3, 4), (5, 3), (10, 3), (16, 2)):
            for state in range(base**digits):
                self.assertEqual(
                    kaprekar_step_general(state, digits, base),
                    independent_reference_step(state, digits, base),
                    msg=f"base={base}, digits={digits}, state={state}",
                )

    def test_generic_cycle_detection(self) -> None:
        result = analyze_general_state(9, digits=2, base=10)
        self.assertEqual(result["cycle_length"], 5)
        self.assertEqual(result["transient_depth"], 0)
        self.assertEqual(trajectory_general(9, digits=2, base=10), [9, 81, 63, 27, 45, 9])

    def test_selection_parser(self) -> None:
        self.assertEqual(parse_integer_selection("2:5"), (2, 3, 4, 5))
        self.assertEqual(parse_integer_selection("10,2,10,4"), (2, 4, 10))


class WeightedAnalysisTests(unittest.TestCase):
    def test_multinomial_class_weight(self) -> None:
        self.assertEqual(multiset_weight((0, 1, 2, 3)), 24)
        self.assertEqual(multiset_weight((0, 0, 1, 2)), 12)
        self.assertEqual(multiset_weight((7, 7, 7, 7)), 1)

    def test_weighted_totals_for_representative_systems(self) -> None:
        for base, digits in ((2, 2), (3, 4), (7, 3), (10, 4)):
            analysis = analyze_system(base, digits)
            self.assertEqual(sum(row["class_weight"] for row in analysis.class_records), base**digits)
            self.assertEqual(
                sum(row["class_weight"] for row in analysis.class_records if row["is_repdigit"]),
                base,
            )
            self.assertEqual(analysis.summary["valid_state_count"], base**digits - base)
            self.assertEqual(
                analysis.summary["permutation_class_count"],
                math.comb(base + digits - 1, digits),
            )

    def test_cycle_class_depth_splitting(self) -> None:
        analysis = analyze_system(10, 2)
        class_09 = next(row for row in analysis.class_records if row["digit_multiset"] == "09")
        self.assertEqual(class_09["class_weight"], 2)
        self.assertEqual(class_09["cycle_member_count"], 1)
        self.assertEqual(class_09["cycle_member_state"], "09")
        self.assertEqual(class_09["noncycle_member_count"], 1)
        self.assertEqual(class_09["noncycle_transient_depth"], 1)
        self.assertFalse(class_09["all_members_same_depth"])

    def test_weighted_analysis_matches_brute_force(self) -> None:
        for base, digits in ((3, 3), (4, 3), (5, 2)):
            weighted = analyze_system(base, digits)
            brute_depths: Counter[int] = Counter()
            brute_cycles: set[tuple[int, ...]] = set()
            for state in range(base**digits):
                result = analyze_general_state(state, digits, base)
                brute_cycles.add(result["cycle"])
                if not result["is_repdigit"]:
                    brute_depths[result["transient_depth"]] += 1
            exported_depths: Counter[int] = Counter()
            for row in weighted.depth_records:
                exported_depths[int(row["transient_depth"])] += int(row["state_count_valid"])
            self.assertEqual(exported_depths, brute_depths)
            self.assertEqual(weighted.summary["attractor_count"], len(brute_cycles))
            self.assertEqual(weighted.summary["maximum_transient_depth_valid"], max(brute_depths))

    def test_decimal_four_digit_summary_matches_existing_dataset(self) -> None:
        dataset = Path("data/kaprekar_results.csv")
        if not dataset.exists():
            self.skipTest("classical dataset has not been generated")
        with dataset.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        valid = [row for row in rows if row["is_repdigit"] == "False"]
        depths = Counter(int(float(row["distance_to_6174"])) for row in valid)
        weighted = analyze_system(10, 4)
        weighted_depths: Counter[int] = Counter()
        for row in weighted.depth_records:
            weighted_depths[int(row["transient_depth"])] += int(row["state_count_valid"])
        self.assertEqual(weighted_depths, depths)
        self.assertEqual(weighted.summary["permutation_class_count"], 715)
        self.assertEqual(weighted.summary["unique_output_count"], 55)

    def test_selected_grid_has_expected_class_count(self) -> None:
        count = sum(math.comb(base + digits - 1, digits) for base in range(2, 17) for digits in range(2, 7))
        self.assertEqual(count, 244_999)


class ProofCertificateTests(unittest.TestCase):
    def test_pair_certificate(self) -> None:
        rows = generate_pair_certificate()
        self.assertEqual(len(rows), 55)
        self.assertEqual(len({row["kaprekar_output"] for row in rows}), 55)
        self.assertEqual(sum(row["attractor_pair"] == "(6,2)" for row in rows), 54)
        self.assertEqual(max(row["pair_distance_to_attractor"] for row in rows), 6)
        maximum_pairs = {
            (row["x"], row["y"]) for row in rows if row["is_maximum_nonzero_pair"]
        }
        self.assertEqual(maximum_pairs, EXPECTED_MAXIMUM_PAIRS)
        fixed = next(row for row in rows if (row["x"], row["y"]) == (6, 2))
        self.assertEqual(fixed["kaprekar_output"], 6174)
        self.assertEqual((fixed["next_x"], fixed["next_y"]), (6, 2))


if __name__ == "__main__":
    unittest.main()

