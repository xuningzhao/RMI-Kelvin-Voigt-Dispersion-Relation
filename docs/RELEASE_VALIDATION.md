# Solver release validation

Validation environment: Linux x86_64, CPython 3.11.15, NumPy 2.4.6,
mpmath 1.3.0, Matplotlib 3.11.1, SymPy 1.14.0. Exact dependency versions
are in `requirements-lock.txt`.

Executed during release preparation:

- `python -m unittest discover -s tests`: 118 tests passed. This includes
  the two critical-transition tests previously missed by unittest discovery.
- Arbitrary-precision single-point JSON example at `Arho=-0.9974`,
  `Amu=-0.999998`, `AG=-1`, `Ek=0`, 80 digits: 14 polynomial candidates,
  four genuine candidates, two physically admissible roots. The accepted
  internal normalized residuals were approximately `3.45e-73`.
- Tests compare the air–hydrogel accepted spectrum at 60 and 80 digits,
  verify medium-exchange symmetry, and check that the public MP entry point
  restores the caller's precision and rejects out-of-scope inputs.
- The documented 3 by 3 MP sweep completed with one worker at 80 digits,
  producing all three expected artifacts and a 3 by 3 count array.
- A candidate archive was extracted to a separate directory. All 118 tests
  and both documented single-point examples passed from the extracted tree,
  using the recorded environment. Individual exported file hashes were
  verified against its manifest. A fresh dependency installation and other
  operating systems were not tested.

The test suite checks the frozen Wolfram exports through SymPy. It does not
re-run the Mathematica proofs or regenerate the expensive manuscript sweeps.
Passing these tests is numerical regression evidence, not a certified
completeness or error-bound proof.

## Build the review archive

```bash
python scripts/build_solver_release.py \
  --output dist/rmi-solver-candidate.zip
```

The archive builder uses an explicit file selection and refuses to overwrite
an existing archive. It does not stage, commit, tag, push, or publish anything.
Generated data, raw logs, private environment files, historical experiments,
and cluster submission scripts are excluded. The source snapshot contains
the present solver changes, even if they have not yet been committed.

`SOURCE_MANIFEST.json` records individual content hashes; the adjacent
`.zip.sha256` records the complete archive hash. Extract into a fresh directory
and rerun the tests from that directory before release.

## Remaining publication decisions

`release/metadata.json` deliberately leaves authors, license, paper citation,
and DOI unspecified. A final public release needs owner-approved metadata
and an appropriate license file. The candidate should also receive a stable
version/tag after review. No author identity, license grant, publication DOI,
or clean commit provenance is inferred from the local environment.

## Standalone repository preparation

The public-folder candidate also includes a GitHub Actions workflow and a
small reference JSON example. GitHub-hosted CI has not run yet. Publication
metadata is still pending in `release/metadata.json`; the citation file is a
template rather than an asserted author list.

Standalone-folder validation: all 118 tests passed in the recorded Python
environment, the JSON example reproduced the two expected admissible roots,
and all local Markdown links resolved. SHA-256 comparisons confirmed that
37 numerical and symbolic source files matched the imported candidate.
