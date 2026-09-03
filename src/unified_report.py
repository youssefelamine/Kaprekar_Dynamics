"""Assemble the project reports into one canonical reading document."""

from __future__ import annotations

import html
import re
from pathlib import Path


def _between(text: str, start: str, end: str | None = None) -> str:
    """Return report content between two top-level section markers."""

    content = text.split(start, 1)[1]
    if end is not None:
        content = content.split(end, 1)[0]
    return content.strip()


def _lower_headings(text: str) -> str:
    """Nest an existing report beneath a part heading."""

    return re.sub(r"^(#{2,})", lambda match: "#" + match.group(1), text, flags=re.MULTILINE)


def _repair_classical_math(text: str) -> str:
    """Repair escape sequences damaged in the legacy generated Markdown."""

    replacements = {
        "(S={0,1,ldots,9999})": r"\(S=\{0,1,\ldots,9999\}\)",
        "(nin S)": r"\(n\in S\)",
        "(D(n))": r"\(D(n)\)",
        "(A(n))": r"\(A(n)\)",
        "(K:S\nightarrow S)": r"$K:S\rightarrow S$",
        "(S) is finite": "$S$ is finite",
        "(d(n))": "$d(n)$",
        "(T(n))": "`T(n)`",
        "(tgeq0)": "`t ≥ 0`",
        "(K^t(n)=6174)": "`K^t(n) = 6174`",
        "(T(3524)=3)": "`T(3524) = 3`",
        "(T(6174)=0)": "`T(6174) = 0`",
        "ageq bgeq cgeq d": r"a\geq b\geq c\geq d",
        "((x,y)mapsto K(n))": "`(x,y) → K(n)`",
        "(0leq yleq xleq9)": "`0 ≤ y ≤ x ≤ 9`",
        "(999x+90y)": "`999x + 90y`",
    }
    for damaged, repaired in replacements.items():
        text = text.replace(damaged, repaired)
    return text


def _size_readme_figures(text: str, width: int = 560) -> str:
    """Render plots at a restrained width while linking to full resolution."""

    image_pattern = re.compile(r"^!\[([^\]]*)\]\(([^)]+)\)$", re.MULTILINE)

    def replacement(match: re.Match[str]) -> str:
        alt = html.escape(match.group(1), quote=True)
        source = html.escape(match.group(2), quote=True)
        return (
            '<p align="center">\n'
            f'  <a href="{source}"><img src="{source}" alt="{alt}" width="{width}"></a>\n'
            '</p>'
        )

    return image_pattern.sub(replacement, text)


def write_unified_report(project_root: Path) -> Path:
    """Build the repository's single Markdown document: ``README.md``."""

    report_dir = project_root / "report"
    classical = (project_root / "src" / "classical_report_source.txt").read_text(encoding="utf-8")
    generalized = (report_dir / "generalized_kaprekar_report.md").read_text(encoding="utf-8")
    proof = (report_dir / "6174_finite_proof.md").read_text(encoding="utf-8")

    classical_body = _repair_classical_math(_lower_headings(
        "## 1. Introduction\n\n" + _between(classical, "## 1. Introduction", "## Data and artifact index")
    )).replace("../", "")
    generalized_body = _lower_headings(
        "## 1. Research design\n\n"
        + _between(generalized, "## 1. Research design", "## References")
    ).replace(
        "The companion [finite proof chapter](6174_finite_proof.md)",
        "Part III below",
    )
    generalized_body = generalized_body.replace("../", "")
    proof_body = _lower_headings(_between(proof, "## Statement")).replace("../", "")
    references = _between(classical, "## References")

    report = f"""# Kaprekar Dynamics

### An exhaustive study of 6174 and generalized Kaprekar systems

An exhaustive computational study of Kaprekar's routine, from the classical four-digit 6174 problem to generalized fixed-width systems in bases 2–16. The project combines complete state-space enumeration, graph analysis, symmetry reduction, and a finite proof of the classical seven-step bound.

## Key findings

- **10,000** classical decimal states were exhaustively analyzed.
- **9,990 out of 9,990** non-repdigit states converge to `6174`.
- **7 transformations** is the exact maximum; 2,184 states attain it.
- The classical state space contracts from **10,000 ordered states to 715 digit multisets to 55 first-step outputs**.
- **75 generalized systems** were analyzed exactly, revealing **199 attractors**, transient depths up to **31**, and cycles up to length **14**.

## Key insight

For descending digits $a \\ge b \\ge c \\ge d$, define $x=a-d$ and $y=b-c$. Then

```text
K(n) = 999(a − d) + 90(b − c) = 999x + 90y

10,000 ordered states  →  715 digit multisets  →  55 first-step outputs
```

This many-to-one contraction is the central structural reason the classical routine converges so quickly.

## Why this is interesting

The classical 6174 phenomenon is often presented as a numerical curiosity: repeatedly rearranging and subtracting four decimal digits somehow leads to the same constant. The exhaustive state-space analysis shows that the behavior is less mysterious when viewed as a finite dynamical system.

Digit permutation immediately removes most positional information, and the four-digit transformation can be expressed using only two digit differences. The generalized census also shows that this behavior is not universal: other bases and widths can have multiple attractors and non-trivial cycles. The familiar four-digit decimal case is therefore a particularly simple member of a much richer family of finite dynamical systems.

> **Scope of contribution.** Convergence to 6174 and the seven-step upper bound are established mathematical results. This repository independently reproduces them exhaustively, characterizes the complete state space and maximum-depth cases, supplies a computer-checkable 55-pair certificate, and extends the analysis to 75 base-and-width systems.

## Contents

- [Quick start](#quick-start)
- [Repository navigation](#repository-navigation)
- [Part I — The classical four-digit decimal system](#part-i--the-classical-four-digit-decimal-system)
- [Part II — Generalized Kaprekar systems](#part-ii--generalized-kaprekar-systems)
- [Part III — Finite proof certificate](#part-iii--finite-proof-certificate-for-the-seven-step-6174-bound)
- [Reproducibility](#reproducibility)
- [Artifact index](#artifact-index)
- [References](#references)

## Quick start

The exact tested environment is Python 3.13.12 with the package versions pinned in `requirements.txt`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 -m src.pipeline
python3 -m src.generalized_pipeline
python3 -m unittest discover -s tests -v
```

## Repository navigation

- [`src/`](src/): transformation logic, exhaustive analyses, visualization, proof, and pipeline code
- [`tests/`](tests/): independent reference checks and generalized weighted-census tests
- [`data/`](data/): complete classical states and generalized weighted-class records
- [`tables/`](tables/): summaries, basins, depth distributions, and proof certificate
- [`figures/`](figures/): all plots in PNG and PDF formats
- [`requirements.txt`](requirements.txt): exact tested Python package versions

## Part I — The classical four-digit decimal system

{classical_body}

## Part II — Generalized Kaprekar systems

{generalized_body}

## Part III — Finite proof certificate for the seven-step 6174 bound

{proof_body}

## Reproducibility

The pipelines are deterministic. `src.pipeline` regenerates the complete 10,000-state decimal study. `src.generalized_pipeline` regenerates the 75-system census, proof certificate, figures, and this README. Custom inclusive grids are supported, for example:

```bash
python3 -m src.generalized_pipeline --bases 2:16 --digits 2:6
```

### Project structure

```text
src/                 Core maps, analysis, pipelines, figures, and report generation
tests/               Classical and generalized verification suites
data/                Full classical dataset and generalized weighted-class census
tables/              Classical and generalized machine-readable result tables
figures/             Classical and generalized PNG/PDF visualizations
requirements.txt     Python dependencies
.python-version      Exact Python version used for the published run
README.md            Complete report and repository documentation
```

### Validation strategy

The 24 tests compare the classical transformation with an independent literal implementation on all 10,000 states, verify the identity `K(n) = 999x + 90y` exhaustively, and compare the specialized and generalized decimal implementations. Generalized tests compare weighted symmetry results with brute-force ordered enumeration, verify multinomial totals and cycle-class depth splitting, test bases through 16, and regenerate every invariant in the 55-pair proof certificate.

## Artifact index

The principal machine-readable artifacts are:

- Classical 10,000-state dataset: [`data/kaprekar_results.csv`](data/kaprekar_results.csv)
- Classical numerical summary: [`data/analysis_summary.json`](data/analysis_summary.json)
- All 715 classical permutation classes: [`tables/permutation_classes.csv`](tables/permutation_classes.csv)
- All 55 classical difference pairs: [`tables/xy_reduced_states.csv`](tables/xy_reduced_states.csv)
- All 2,184 seven-step states: [`tables/maximum_distance_states.csv`](tables/maximum_distance_states.csv)
- Generalized 75-system summary: [`tables/generalized/system_summary.csv`](tables/generalized/system_summary.csv)
- Every generalized attractor and basin: [`tables/generalized/attractors_and_basins.csv`](tables/generalized/attractors_and_basins.csv)
- Exact weighted depth distributions: [`tables/generalized/depth_distributions.csv`](tables/generalized/depth_distributions.csv)
- All 244,999 weighted classes: [`data/generalized/weighted_classes.csv`](data/generalized/weighted_classes.csv)
- Generalized JSON summary: [`data/generalized/generalized_summary.json`](data/generalized/generalized_summary.json)
- The 55-row proof certificate: [`tables/generalized/kaprekar_6174_pair_certificate.csv`](tables/generalized/kaprekar_6174_pair_certificate.csv)

## References

{references}
"""
    report = _size_readme_figures(report)
    output = project_root / "README.md"
    output.write_text(report, encoding="utf-8")
    return output
