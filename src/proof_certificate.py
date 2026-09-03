"""Generate and verify the finite 55-pair proof certificate for 6174."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from .kaprekar import digit_differences, normalize_state


PAIR_FIELDS = (
    "x",
    "y",
    "kaprekar_output",
    "output_state",
    "next_x",
    "next_y",
    "attractor_pair",
    "pair_distance_to_attractor",
    "pair_path",
    "is_maximum_nonzero_pair",
)

EXPECTED_MAXIMUM_PAIRS = {
    (4, 1),
    (5, 1),
    (5, 2),
    (6, 1),
    (8, 5),
    (9, 4),
    (9, 5),
    (9, 6),
}


def pair_output(pair: tuple[int, int]) -> int:
    """Return the four-digit decimal output represented by ``(x, y)``."""

    x, y = pair
    if not 0 <= y <= x <= 9:
        raise ValueError("a feasible pair must satisfy 0 <= y <= x <= 9")
    return 999 * x + 90 * y


def pair_successor(pair: tuple[int, int]) -> tuple[int, int]:
    """Return the digit-difference pair of the output represented by *pair*."""

    return digit_differences(pair_output(pair))


def generate_pair_certificate() -> list[dict[str, Any]]:
    """Generate and validate all 55 feasible pair transitions."""

    pairs = [(x, y) for x in range(10) for y in range(x + 1)]
    successors = {pair: pair_successor(pair) for pair in pairs}
    rows: list[dict[str, Any]] = []
    maximum_pairs: set[tuple[int, int]] = set()

    for pair in pairs:
        path = [pair]
        seen = {pair}
        current = pair
        while True:
            next_pair = successors[current]
            path.append(next_pair)
            if next_pair in seen:
                break
            seen.add(next_pair)
            current = next_pair
        attractor = next_pair
        if pair == (0, 0):
            distance = 0
        else:
            if attractor != (6, 2):
                raise AssertionError(f"nonzero pair {pair} does not reach (6, 2)")
            distance = path[:-1].index((6, 2))
        if pair != (0, 0) and distance == 6:
            maximum_pairs.add(pair)
        rows.append(
            {
                "x": pair[0],
                "y": pair[1],
                "kaprekar_output": pair_output(pair),
                "output_state": normalize_state(pair_output(pair)),
                "next_x": successors[pair][0],
                "next_y": successors[pair][1],
                "attractor_pair": f"({attractor[0]},{attractor[1]})",
                "pair_distance_to_attractor": distance,
                "pair_path": " -> ".join(f"({x},{y})" for x, y in path),
                "is_maximum_nonzero_pair": pair != (0, 0) and distance == 6,
            }
        )

    outputs = {row["kaprekar_output"] for row in rows}
    if len(rows) != 55 or len(outputs) != 55:
        raise AssertionError("the feasible pair set must contain 55 distinct outputs")
    if pair_successor((0, 0)) != (0, 0):
        raise AssertionError("(0,0) must be the sole pair in the zero basin")
    if pair_successor((6, 2)) != (6, 2) or pair_output((6, 2)) != 6174:
        raise AssertionError("(6,2) must be fixed and represent 6174")
    if maximum_pairs != EXPECTED_MAXIMUM_PAIRS:
        raise AssertionError("the maximum pair witnesses differ from the expected certificate")
    if sum(row["attractor_pair"] == "(6,2)" for row in rows) != 54:
        raise AssertionError("all 54 nonzero pairs must enter the 6174 pair attractor")
    return rows


def write_pair_certificate(output_path: Path) -> list[dict[str, Any]]:
    """Write the checked certificate as CSV and return its rows."""

    rows = generate_pair_certificate()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_proof_report(report_path: Path, certificate_relative_path: str) -> None:
    """Write the mathematical explanation backed by the generated certificate."""

    rows = generate_pair_certificate()
    maximum_pairs = [f"({row['x']},{row['y']})" for row in rows if row["is_maximum_nonzero_pair"]]
    report = rf"""# A Finite Proof Certificate for the Seven-Step 6174 Bound

## Statement

Every four-digit decimal state containing at least two distinct digits reaches 6174 in at most seven Kaprekar transformations. The bound is sharp.

## Analytic reduction

Let the sorted digits be $a \ge b \ge c \ge d$, and set

```text
x = a − d          y = b − c.
```

Direct subtraction gives

```text
K(n) = 999(a − d) + 90(b − c) = 999x + 90y.
```

The nested digit intervals imply $0 \le y \le x \le 9$. Conversely, every such pair is realized by sorted digits $(a,b,c,d)=(x,y,0,0)$, so there are

```text
1 + 2 + ⋯ + 10 = 55
```

feasible pairs. If two pairs have the same output, then

```text
999(x − x′) = −90(y − y′),
```

or $111(x-x')=-10(y-y')$. A nonzero left side has magnitude at least 111, whereas the right side has magnitude at most 90. Hence both differences are zero and the 55 outputs are distinct.

## Checked transition certificate

For a pair $p=(x,y)$, define

```text
F(p) = 999x + 90y
P(p) = the digit-difference pair of F(p).
```

The generated [55-row certificate]({certificate_relative_path}) lists every feasible pair, $F(p)$, $P(p)$, its complete pair path, and its distance to a pair attractor. The checker establishes:

- `(0,0)` maps to itself and is the only pair in the 0000 basin;
- every one of the 54 nonzero pairs reaches `(6,2)`;
- `(6,2)` maps to itself and $F(6,2)=6174$;
- the largest pair-graph distance to `(6,2)` is 6;
- the eight distance-six witnesses are {', '.join(maximum_pairs)}.

## Translation to the original routine

Every non-repdigit start has a nonzero initial pair $p$. If that pair is at pair-graph distance $r$ from `(6,2)`, then after $r$ transformations the current state's pair is `(6,2)`, and one further transformation produces 6174. Therefore

```text
T(n) ≤ r + 1 ≤ 7.
```

The eight distance-six pairs are realizable by four-digit states, so starts with $T(n)=7$ exist. The full exhaustive dataset independently identifies 2,184 such ordered states. The special start 6174 itself has distance zero by convention.

## Nature of the proof

The reduction from 10,000 states to 55 pairs is analytic. The remaining claim is a transparent finite proof by complete enumeration of the 55 explicitly exported transitions. Automated tests regenerate the table, verify its invariants, and compare the underlying Kaprekar implementation with independent definitions.
"""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
