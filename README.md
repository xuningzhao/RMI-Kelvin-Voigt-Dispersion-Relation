# Kelvin–Voigt RMI dispersion relation

Numerical solver for the two-medium Kelvin–Voigt dispersion relation. The
solver generates degree-14 polynomial candidates, verifies them against the
original dispersion relation, and applies the semi-infinite spatial-decay
conditions. Both float64 and full arbitrary-precision workflows are included.

## Install and test

Use Python 3.11. Run these commands from the repository root:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

`requirements.txt` records the exact tested solver, plotting, and test
dependencies. Validation used Python 3.11.15 on Linux x86_64; the dependency
file does not pin the operating system, BLAS, or Python itself. Mathematica
is needed only to rerun the Wolfram proof scripts, not to run the solver or
Python tests.

## Solve one point

For difficult contrasts and closely spaced roots, use the full
arbitrary-precision reference path with decimal-string inputs:

```bash
python scripts/analyze_mp_parameter_point.py \
  --Arho -0.9974 --Amu -0.999998 --AG -1 --Ek 0 --precision 80
```

The JSON output contains candidate/genuine/admissible counts, accepted roots,
spatial decay factors, and normalized source residuals. This air–hydrogel
example has two admissible roots, approximately
`-0.456786491471073 +/- 0.839594002196374j`. The recorded output is in
[examples/air_hydrogel_reference.json](examples/air_hydrogel_reference.json).
Acceptance uses the requested precision; exported root coordinates are
rounded to binary64. Reported high-precision residuals describe the internal
roots, not their rounded JSON representations.

Python API:

```python
from src.mp_reference import solve_reference_point

result = solve_reference_point(
    "-0.9974", "-0.999998", "-1", "0", decimal_precision=80
)
for s, q_plus, q_minus, residual in result.accepted:
    print(s, q_plus, q_minus, residual)
```

The faster float64 workflow provides detailed rejection classifications:

```bash
python scripts/analyze_parameter_point.py \
  --Arho 0.2 --Amu -0.3 --AG 0.4 --Lambda 1.7 \
  --enable-high-precision-refinement
```

Its refinement option revisits ambiguous candidates only. It cannot recover
all candidates potentially lost during float64 generation or earlier rejection.

## Parameter sweeps

Small full-MP reference sweep:

```bash
python scripts/run_mp_reference_sweep.py \
  --Arho -0.9974 --Amu -0.999998 --samples 3 --workers 1 \
  --precision 80 --output-dir results/reference_smoke
```

This writes `grid_data.npz`, `admissible_roots.csv`, and `metadata.json` on a
3 by 3 `(AG, Ek)` grid. Generated results are ignored by Git. The expensive
501 by 501 manuscript datasets are not distributed with the source.

<details>
<summary>Reproduce the 31 by 31 asymmetric float64 benchmark</summary>


`benchmarks/verify_Arho_m0p95_Amu_m0p8.yaml` records this 31 by 31 study. The
CLI does not read YAML; these Bash commands construct the corresponding grid:

```bash
AG_VALUES=$(python -c 'import numpy as np; print(" ".join(f"{x:.17g}" for x in np.linspace(-0.8, 0.8, 31)))')
EK_VALUES=$(python -c 'import numpy as np; print(" ".join(f"{x:.17g}" for x in np.linspace(-0.8, 0.8, 31)))')
python scripts/run_parameter_sweep.py \
  --Arho -0.95 --Amu -0.8 --sweep AG $AG_VALUES --sweep Ek $EK_VALUES \
  --output results/verify_Arho_m0p95_Amu_m0p8/sweep_output.jsonl \
  --enable-high-precision-refinement
python scripts/summarize_parameter_sweep.py \
  results/verify_Arho_m0p95_Amu_m0p8/sweep_output.jsonl \
  --output-dir results/verify_Arho_m0p95_Amu_m0p8
python scripts/plot_parameter_sweep.py \
  results/verify_Arho_m0p95_Amu_m0p8 \
  --output-dir results/verify_Arho_m0p95_Amu_m0p8/figures --all --color-by Ek
```

This benchmark is separate from the tiny MP smoke test; it was not rerun as
part of standalone-folder preparation. It is included as a reproduction recipe.

</details>

<details>
<summary>Manifest-driven production campaigns</summary>

This optional workflow uses float64 with configurable targeted refinement.
It is separate from the full-MP reference sweep above.

`config/production_sweep_v1.json` describes Group A (343 cases with 501 Ek
points) and three Group B families (49 maps per family, 501 by 501 points per
map). Group A uses atomic case payloads; Group B uses checksummed tiles.

### Build manifests

```bash
python scripts/build_production_manifests.py \
  --config config/production_sweep_v1.json \
  --campaign-root results/production_v1
```

This writes manifests under `results/production_v1/manifests/` and records
configuration and environment provenance. It does not launch solver jobs.

Use `python scripts/run_production_bundle.py --help` to select a manifest,
case-worker count, cases per worker, and optional work limit. The number of
manifest records must equal workers times cases per worker. Before running a
full campaign, select resources and a configuration appropriate to your
machine. Cluster-specific submission scripts are not part of this release.

### Validate completed output

```bash
python scripts/validate_production_results.py \
  --config config/production_sweep_v1.json \
  --manifest results/production_v1/manifests/a.jsonl \
  --campaign-root results/production_v1
python scripts/update_production_index.py results/production_v1
python scripts/build_retry_manifest.py \
  --manifest results/production_v1/manifests/a.jsonl \
  --campaign-root results/production_v1 \
  --output results/production_v1/manifests/a_retry.jsonl
```

Validation fails while selected cases are incomplete. Restart skips only
payloads whose configuration hash, checksum, schema, and internal array
structure validate. Failed or damaged output units are recomputed. Keep the
original configuration and manifest with each campaign.

</details>

## Input domain and precision

The reference path requires finite real inputs with `|Arho| < 1`,
`|Amu| < 1`, `|AG| <= 1`, and `-1 < Ek < 1`. The float64 command accepts
`Lambda`, related to `Ek` by `Lambda = ((1-Ek)/(1+Ek))^2`.

Check closely spaced roots at increased precision. Accepted binary64 outputs
can merge at the deduplication tolerance. Use separate processes rather than
threads for concurrent mpmath solves.

## Validation and contributions

The 118-test suite covers dimensional equivalence, frozen symbolic coefficient
comparisons, domain exclusions, branch choices, reduced degree, conjugacy,
medium-exchange symmetry, 60/80-digit reference convergence, and workflow
storage/restart behavior. Tests and the single-point examples passed from an
independently extracted source archive; the 3 by 3 MP sweep also completed.
Full manuscript sweeps and Wolfram proofs were not rerun in these checks.
The GitHub workflow in `.github/workflows/tests.yml` runs tests and examples
on pushes and pull requests; its current result is available in the repository's
Actions tab.

For a numerical change, add a regression test and check source residuals,
spatial admissibility, and increased-precision agreement. Update symbolic
sources and their validation when changing coefficient formulas. Include
parameters, backend, precision, and dependency versions in bug reports.

## Source layout and release

- `src/`, `scripts/`, `tests/`: solver, commands, and tests.
- `symbolic/`, `derivation/`: executable symbolic sources and LaTeX derivation.
- `config/`, `benchmarks/`, `examples/`: reproducibility inputs and reference output.

The current version is **0.1.0rc1**. Author-approved authorship, license, and
citation metadata remain pending; no software license has been granted yet.
Add `LICENSE` and `CITATION.cff` when those decisions are made. The archive
builder includes both automatically if present:

```bash
python scripts/build_solver_release.py --output dist/rmi-solver-candidate.zip
```

The source selection is in `config/source_manifest.json`. The archive includes
a `SOURCE_MANIFEST.json` of per-file SHA-256 hashes; an adjacent `.zip.sha256`
checksums the archive. Existing archives are not overwritten. Before a tagged
release, run the tests from an extracted copy and record the chosen version
and citation metadata. This is a source-code release, not a PyPI package.

This repository starts a separate Git history from the research project.
The imported candidate was `rmi-solver-candidate-20260923.zip`, SHA-256
`4ae078f7b5d1022254dd6c1cede5a68e1f9459b054ec0c14fca1dbd63d6b3f73`.
Earlier preparation records and imported file hashes remain in Git history.
