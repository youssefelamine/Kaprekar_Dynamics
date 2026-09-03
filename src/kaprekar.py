"""Core transformations for Kaprekar's four-digit routine.

All public four-digit functions accept integer states in the closed interval
0..9999.  Integers are only a storage representation: every transformation
operates on the corresponding zero-padded four-character digit string.
"""

from __future__ import annotations

from numbers import Integral
from typing import Any


STATE_WIDTH = 4
STATE_COUNT = 10**STATE_WIDTH
DIGIT_SYMBOLS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _validate_state(n: int) -> int:
    """Validate and return a four-digit state as a built-in integer."""

    if isinstance(n, bool) or not isinstance(n, Integral):
        raise TypeError("a Kaprekar state must be an integer")
    value = int(n)
    if not 0 <= value < STATE_COUNT:
        raise ValueError("a four-digit state must be between 0 and 9999")
    return value


def normalize_state(n: int) -> str:
    """Return *n* as an exactly four-character, zero-padded string."""

    return f"{_validate_state(n):04d}"


def kaprekar_components(n: int) -> tuple[int, int, int]:
    """Return the descending number, ascending number, and their difference."""

    digits = normalize_state(n)
    ascending_digits = "".join(sorted(digits))
    descending_digits = ascending_digits[::-1]
    descending = int(descending_digits)
    ascending = int(ascending_digits)
    return descending, ascending, descending - ascending


def kaprekar_step(n: int) -> int:
    """Perform exactly one four-digit Kaprekar transformation."""

    return kaprekar_components(n)[2]


def is_repdigit(n: int) -> bool:
    """Return whether all four digits of *n*, including leading zeros, agree."""

    return len(set(normalize_state(n))) == 1


def digit_multiset(n: int) -> str:
    """Return the state's digits sorted in ascending order."""

    return "".join(sorted(normalize_state(n)))


def digit_differences(n: int) -> tuple[int, int]:
    """Return ``(x, y) = (a-d, b-c)`` for digits ``a >= b >= c >= d``."""

    d, c, b, a = (int(digit) for digit in digit_multiset(n))
    return a - d, b - c


def trajectory(n: int, max_iterations: int = 100) -> list[int]:
    """Return states from *n* through the first repeated state.

    The initial state and the repeated closing state are both included.  Thus
    ``trajectory(6174)`` is ``[6174, 6174]``.  A ``RuntimeError`` is raised if
    no repeat is found within *max_iterations* transformations.
    """

    current = _validate_state(n)
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral):
        raise TypeError("max_iterations must be an integer")
    if max_iterations < 1:
        raise ValueError("max_iterations must be positive")

    seen = {current}
    states = [current]
    for _ in range(int(max_iterations)):
        current = kaprekar_step(current)
        states.append(current)
        if current in seen:
            return states
        seen.add(current)
    raise RuntimeError(f"no repeated state found within {max_iterations} transformations")


def canonical_cycle(cycle: tuple[int, ...]) -> tuple[int, ...]:
    """Rotate a directed cycle to its lexicographically smallest form."""

    if not cycle:
        raise ValueError("a cycle cannot be empty")
    rotations = tuple(cycle[index:] + cycle[:index] for index in range(len(cycle)))
    return min(rotations)


def analyze_start(n: int, max_iterations: int = 100) -> dict[str, Any]:
    """Classify a starting state and return its trajectory statistics.

    No attractor is presumed.  The first repeated state is located, the cycle
    is extracted from the trajectory, and distances are derived from positions
    in that trajectory.
    """

    start = _validate_state(n)
    states = trajectory(start, max_iterations=max_iterations)
    repeated = states[-1]
    cycle_start = states[:-1].index(repeated)
    cycle = canonical_cycle(tuple(states[cycle_start:-1]))
    descending, ascending, first_result = kaprekar_components(start)
    state_text = normalize_state(start)
    distinct_digits = len(set(state_text))

    unique_path = states[:-1]
    distance_to_6174 = unique_path.index(6174) if 6174 in unique_path else None
    distance_to_0000 = unique_path.index(0) if 0 in unique_path else None
    reaches_6174 = distance_to_6174 is not None
    reaches_0000 = distance_to_0000 is not None

    return {
        "start": start,
        "start_state": state_text,
        "is_repdigit": distinct_digits == 1,
        "distinct_digits": distinct_digits,
        "digit_multiset": digit_multiset(start),
        "x": digit_differences(start)[0],
        "y": digit_differences(start)[1],
        "descending": descending,
        "ascending": ascending,
        "first_result": first_result,
        "first_result_state": normalize_state(first_result),
        "trajectory": states,
        "final_attractor": cycle[0],
        "iterations_to_attractor": cycle_start,
        "distance_to_6174": distance_to_6174,
        "distance_to_0000": distance_to_0000,
        "reaches_6174": reaches_6174,
        "reaches_0000": reaches_0000,
        "other_cycle_detected": not reaches_6174 and not reaches_0000,
        "cycle_length": len(cycle),
        "cycle": cycle,
    }


def _validate_general_parameters(n: int, digits: int, base: int) -> tuple[int, int, int]:
    """Validate a generalized state specification."""

    for name, value in (("n", n), ("digits", digits), ("base", base)):
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError(f"{name} must be an integer")
    n, digits, base = int(n), int(digits), int(base)
    if digits < 1:
        raise ValueError("digits must be positive")
    if not 2 <= base <= 36:
        raise ValueError("base must be between 2 and 36")
    if not 0 <= n < base**digits:
        raise ValueError(f"n must be in 0..{base**digits - 1} for this system")
    return n, digits, base


def digits_in_base(n: int, digits: int, base: int = 10) -> tuple[int, ...]:
    """Return the fixed-width, most-significant-first digits of *n* in *base*."""

    n, digits, base = _validate_general_parameters(n, digits, base)
    remaining = n
    values = [0] * digits
    for index in range(digits - 1, -1, -1):
        remaining, values[index] = divmod(remaining, base)
    return tuple(values)


def integer_from_digits(values: tuple[int, ...] | list[int], base: int = 10) -> int:
    """Convert a non-empty sequence of base-*base* digits to an integer."""

    if isinstance(base, bool) or not isinstance(base, Integral):
        raise TypeError("base must be an integer")
    base = int(base)
    if not 2 <= base <= 36:
        raise ValueError("base must be between 2 and 36")
    if not values:
        raise ValueError("values must contain at least one digit")
    result = 0
    for digit in values:
        if isinstance(digit, bool) or not isinstance(digit, Integral):
            raise TypeError("each digit must be an integer")
        digit = int(digit)
        if not 0 <= digit < base:
            raise ValueError(f"each digit must be in 0..{base - 1}")
        result = result * base + digit
    return result


def format_state_general(n: int, digits: int, base: int = 10) -> str:
    """Format *n* as an exactly *digits*-character state in bases 2 through 36."""

    return "".join(DIGIT_SYMBOLS[digit] for digit in digits_in_base(n, digits, base))


def kaprekar_components_general(n: int, digits: int, base: int = 10) -> tuple[int, int, int]:
    """Return descending, ascending, and difference for a fixed-width system."""

    values = sorted(digits_in_base(n, digits, base))
    ascending = integer_from_digits(values, base)
    descending = integer_from_digits(list(reversed(values)), base)
    return descending, ascending, descending - ascending


def kaprekar_step_general(n: int, digits: int, base: int = 10) -> int:
    """Perform one Kaprekar step for an arbitrary width and numeric base."""

    return kaprekar_components_general(n, digits, base)[2]


def trajectory_general(
    n: int,
    digits: int,
    base: int = 10,
    max_iterations: int | None = None,
) -> list[int]:
    """Return a generalized trajectory through its first repeated state."""

    n, digits, base = _validate_general_parameters(n, digits, base)
    if max_iterations is None:
        max_iterations = base**digits + 1
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral):
        raise TypeError("max_iterations must be an integer or None")
    if max_iterations < 1:
        raise ValueError("max_iterations must be positive")

    current = n
    seen = {current}
    states = [current]
    for _ in range(int(max_iterations)):
        current = kaprekar_step_general(current, digits, base)
        states.append(current)
        if current in seen:
            return states
        seen.add(current)
    raise RuntimeError(f"no repeated state found within {max_iterations} transformations")


def analyze_general_state(n: int, digits: int, base: int = 10) -> dict[str, Any]:
    """Return a generic cycle and transient analysis for one generalized state."""

    n, digits, base = _validate_general_parameters(n, digits, base)
    states = trajectory_general(n, digits, base)
    repeated = states[-1]
    cycle_start = states[:-1].index(repeated)
    cycle = canonical_cycle(tuple(states[cycle_start:-1]))
    values = digits_in_base(n, digits, base)
    descending, ascending, successor = kaprekar_components_general(n, digits, base)
    return {
        "start": n,
        "start_state": format_state_general(n, digits, base),
        "base": base,
        "digits": digits,
        "is_repdigit": len(set(values)) == 1,
        "descending": descending,
        "ascending": ascending,
        "successor": successor,
        "trajectory": states,
        "transient_depth": cycle_start,
        "cycle": cycle,
        "cycle_length": len(cycle),
    }
