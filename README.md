# Kelvin–Voigt RMI numerical solver

This source release computes candidate temporal eigenvalues of the two-medium
Kelvin–Voigt dispersion relation, verifies them against the original equation,
and applies the semi-infinite spatial-decay conditions. It includes the frozen
symbolic polynomial, executable tests, and small reproducible examples.

**Release candidate:** authorship, license, and paper/DOI metadata are pending
the project owner's decision. This folder is a local review candidate; no license has been granted yet.
See [publication steps](release/README.md) and `release/metadata.json`.

## Setup

Use Python 3.11. The current release candidate was validated with Python
3.11.15 on Linux x86_64. Other environments require their own validation.
Run these commands from the repository root:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m unittest discover -s tests -v
```

`requirements-lock.txt` records the exact tested Python dependencies, including
plotting and symbolic-comparison dependencies; it does not pin the OS, BLAS,
or Python interpreter. `requirements.txt` contains only solver dependencies.
`requirements-dev.txt` provides broader dependency ranges for development;
these ranges are not a claim that every version combination was tested.
Mathematica is not needed to solve or run the Python tests. Re-running the
original symbolic proof scripts requires a Wolfram installation.

## Solve a point

Use the arbitrary-precision reference path for difficult material contrasts and
closely spaced roots. All inputs below are passed as decimal strings:

```bash
python scripts/analyze_mp_parameter_point.py \
  --Arho -0.9974 --Amu -0.999998 --AG -1 --Ek 0 --precision 80
```

This produces JSON with polynomial/genuine/admissible counts, accepted roots,
both spatial decay factors, and normalized source residuals. The air–hydrogel
example has the admissible conjugate pair approximately
`-0.456786491471073 +/- 0.839594002196374j`. Acceptance is evaluated at the
requested precision; JSON root coordinates are subsequently rounded to
binary64. A high-precision residual describes the internal root, not a promise
that the rounded JSON root achieves that residual.

The Python entry point is:

```python
from src.mp_reference import solve_reference_point

result = solve_reference_point(
    "-0.9974", "-0.999998", "-1", "0", decimal_precision=80
)
for s, q_plus, q_minus, relative_residual in result.accepted:
    print(s, q_plus, q_minus, relative_residual)
```

The fast float64 workflow provides more detailed rejection classifications:

```bash
python scripts/analyze_parameter_point.py \
  --Arho 0.2 --Amu -0.3 --AG 0.4 --Lambda 1.7 \
  --enable-high-precision-refinement
```

Its high-precision option refines ambiguous candidates only; it does not
replace float64 candidate generation with a full arbitrary-precision solve.

## Reproduce a small reference sweep

```bash
python scripts/run_mp_reference_sweep.py \
  --Arho -0.9974 --Amu -0.999998 --samples 3 --workers 1 \
  --precision 80 --output-dir results/reference_smoke
```

This evaluates the Cartesian grid in `(AG, Ek)` and writes `grid_data.npz`,
`admissible_roots.csv`, and `metadata.json`. The full 501 by 501 manuscript
campaign is expensive and is not part of the quick-start validation. Existing
manuscript data and figures are not included in this solver archive.

## Numerical methods and validation

- [Solver methods and limitations](docs/SOLVER_METHODS.md)
- [Detailed root verification workflow](root_verification_workflow.md)
- [Nondimensionalization](SPEC.md)
- [Symbolic provenance](symbolic/README.md)
- [Release validation](docs/RELEASE_VALIDATION.md)

## Repository layout

- `src/`: numerical solver and workflow modules.
- `scripts/`: point solves, sweeps, summaries, and archive building.
- `tests/`: numerical and workflow regression tests.
- `symbolic/`, `derivation/`: frozen coefficients and mathematical provenance.
- `examples/`: a small reference output with reproduction instructions.
- `docs/`: numerical methods and validation evidence.
- `release/`: archive selection, citation template, and publication metadata.

See [source provenance](PROVENANCE.md) for the imported candidate. To build a
source archive, run `python scripts/build_solver_release.py --output
dist/rmi-solver-candidate.zip`. `SOURCE_MANIFEST.json` inside the archive
records the SHA-256 of each exported file. Local results and archives are
ignored by Git. Automated checks are in `.github/workflows/tests.yml`.

The recorded 31 by 31 asymmetric benchmark can be reproduced using the
instructions in [examples/README.md](examples/README.md).
