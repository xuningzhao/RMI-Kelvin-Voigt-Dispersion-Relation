# Contributing

Use the setup commands in README.md, make a branch, and submit a pull request.
Before proposing a change, run:

```bash
python -m unittest discover -s tests -v
```

For a numerical bug, include the four input parameters, solver backend,
working precision, dependency versions, and expected versus observed results.
For a solver change, include a regression test and examine source residuals,
spatial admissibility, and convergence with increased precision.

The polynomial coefficients are mechanically derived from the frozen Wolfram
exports. Algebra changes must update and validate the symbolic source too.
Preserve the distinction between polynomial candidates, genuine roots, and
physically admissible roots.

Keep generated sweep data, environments, and local logs out of commits. Update
`release/source_manifest.json` if a new file belongs in the source archive.
Run the archive builder into `dist/` to inspect the distribution contents.
Authorship and license decisions remain pending for this candidate.
