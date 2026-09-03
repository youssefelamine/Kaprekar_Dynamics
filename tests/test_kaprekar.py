"""Automated tests for the four-digit Kaprekar implementation."""

from __future__ import annotations

import itertools
import unittest

from src.kaprekar import (
    analyze_start,
    digit_differences,
    kaprekar_components,
    kaprekar_step,
    kaprekar_step_general,
    normalize_state,
    trajectory,
)


def reference_step(n: int) -> int:
    """Independent literal implementation used for exhaustive cross-checking."""

    text = str(n).rjust(4, "0")
    low = int("".join(sorted(text)))
    high = int("".join(sorted(text, reverse=True)))
    return high - low


class FourDigitKaprekarTests(unittest.TestCase):
    def test_known_3524_trajectory(self) -> None:
        self.assertEqual(kaprekar_step(3524), 3087)
        self.assertEqual(kaprekar_step(3087), 8352)
        self.assertEqual(kaprekar_step(8352), 6174)
        self.assertEqual(kaprekar_step(6174), 6174)
        self.assertEqual(trajectory(3524), [3524, 3087, 8352, 6174, 6174])

    def test_leading_zero_behavior(self) -> None:
        self.assertEqual(normalize_state(0), "0000")
        self.assertEqual(normalize_state(999), "0999")
        self.assertEqual(kaprekar_components(1000), (1000, 1, 999))
        self.assertEqual(kaprekar_step(1000), 999)
        self.assertEqual(normalize_state(kaprekar_step(1000)), "0999")

    def test_repdigits(self) -> None:
        for state in (0, 1111, 2222, 3333, 4444, 5555, 6666, 7777, 8888, 9999):
            with self.subTest(state=state):
                self.assertEqual(kaprekar_step(state), 0)

    def test_permutations_have_the_same_first_step(self) -> None:
        outputs = {
            kaprekar_step(int("".join(permutation)))
            for permutation in set(itertools.permutations("3524"))
        }
        self.assertEqual(outputs, {3087})

    def test_algebraic_difference_formula_exhaustively(self) -> None:
        for state in range(10_000):
            x, y = digit_differences(state)
            self.assertEqual(kaprekar_step(state), 999 * x + 90 * y)

    def test_against_independent_reference_for_all_states(self) -> None:
        for state in range(10_000):
            self.assertEqual(kaprekar_step(state), reference_step(state), msg=f"state={state:04d}")

    def test_every_state_reaches_a_detected_cycle(self) -> None:
        for state in range(10_000):
            path = trajectory(state)
            self.assertIn(path[-1], path[:-1], msg=f"state={state:04d}")

    def test_analysis_counts_transformations_not_list_entries(self) -> None:
        result = analyze_start(3524)
        self.assertEqual(result["distance_to_6174"], 3)
        self.assertEqual(result["iterations_to_attractor"], 3)
        self.assertTrue(result["reaches_6174"])
        self.assertFalse(result["other_cycle_detected"])
        fixed = analyze_start(6174)
        self.assertEqual(fixed["distance_to_6174"], 0)
        self.assertEqual(fixed["iterations_to_attractor"], 0)

    def test_input_validation(self) -> None:
        for invalid in (-1, 10_000):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                kaprekar_step(invalid)
        for invalid in (3.5, "3524", True):
            with self.subTest(invalid=invalid), self.assertRaises(TypeError):
                kaprekar_step(invalid)  # type: ignore[arg-type]
        with self.assertRaises(RuntimeError):
            trajectory(3524, max_iterations=1)


class GeneralizedKaprekarTests(unittest.TestCase):
    def test_general_decimal_function_agrees_at_width_four(self) -> None:
        for state in range(10_000):
            self.assertEqual(kaprekar_step_general(state, digits=4), kaprekar_step(state))

    def test_three_digit_example(self) -> None:
        self.assertEqual(kaprekar_step_general(210, digits=3), 198)
        self.assertEqual(kaprekar_step_general(211, digits=3), 99)
        self.assertEqual(kaprekar_step_general(495, digits=3), 495)

    def test_binary_width_is_preserved(self) -> None:
        # 0001_2 -> 1000_2 - 0001_2 = 0111_2.
        self.assertEqual(kaprekar_step_general(1, digits=4, base=2), 7)


if __name__ == "__main__":
    unittest.main()
