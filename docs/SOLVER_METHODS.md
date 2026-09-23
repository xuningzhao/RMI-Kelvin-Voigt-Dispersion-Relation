# Numerical solver methods

## Scope and coordinates

The viscous time scale is `t_v = rho_T / (mu_T k^2)`, and the dimensionless
temporal eigenvalue is `s = gamma t_v`. The controls are `Arho`, `Amu`, `AG`,
and `Lambda`; the bounded rate coordinate satisfies
`Lambda = ((1 - Ek)/(1 + Ek))^2`. Figure notation `Ev` denotes this `Ek`
coordinate; it is not a second independent parameter. Dimensionless time in
the amplitude figures is `t' = t/t_v`.

The implemented symbolic theorem applies for finite real controls with
`|Arho| < 1`, `|Amu| < 1`, `|AG| <= 1`, and `Lambda > 0`.
The reference entry point uses `-1 < Ek < 1`. It rejects endpoints and
out-of-scope material contrasts. Candidate-level exclusions are `s = 0`,
`E_plus = 0`, `E_minus = 0`, `B_plus = 0`, and `B_minus = 0`.

## Candidate generation

`src/polynomial_coefficients.py` is the numerical transcription of the frozen
Wolfram coefficient export. Coefficients are stored in ascending powers of
`s`; the generic degree is 14, with exact reductions at special contrasts.
Do not modify these formulas independently of their symbolic source.

The float64 path (`src/polynomial.py`) uses effective-degree handling and a
scaled polynomial solve. Targeted high-precision refinement only revisits
mathematically ambiguous candidates. It cannot guarantee recovery of roots
that float64 candidate generation or earlier rejection has already lost.

The reference path (`src/mp_reference.py`) evaluates the same coefficient
expressions with mpmath, computing Lambda directly from decimal Ek. It
removes exactly zero leading coefficients, deflates exact zero roots while
retaining their multiplicity, balances the variable, and calls `mp.polyroots`
with 200 extra digits and at most 10,000 iterations. Failure to converge
raises an exception; it is not interpreted as absence of physical roots.

## Verification and physical acceptance

For each reference candidate, define

```text
E_plus  = s (1 + Amu) + Lambda (1 + AG)
E_minus = s (1 - Amu) + Lambda (1 - AG)
q_plus  = sqrt(1 + (1 + Arho) s^2 / E_plus)
q_minus = sqrt(1 + (1 - Arho) s^2 / E_minus)
B_minus = E_plus + E_minus q_minus
B_plus  = E_minus + E_plus q_plus
D = s^2/B_minus + s^2/B_plus + 2
```

Both square roots use the principal branch. At working precision `p`, the
reference implementation excludes denominator/domain magnitudes at or below
`10^(-floor(p/2))`. It accepts the mathematical source check when

```text
abs(D) / max(1, abs(s^2/B_minus) + abs(s^2/B_plus) + 2)
    <= 10^(-floor(p/2))
```

Physical acceptance additionally requires, separately for each q,

```text
Re(q) > 1e-12 + 1e-10 max(1, abs(q))
```

The sign of `Re(s)` describes temporal growth or decay only. It is not the
spatial admissibility criterion. The float64 workflow has its own documented
relative/absolute domain and source thresholds in `src/workflow.py`; do not
assume that both backends use identical mathematical acceptance thresholds.

## Counts, precision, and limitations

Polynomial and genuine counts include multiplicity. Accepted reference roots
are converted to binary64 and deduplicated with distance threshold
`1e-12 max(1, abs(s))`, as in the manuscript sweep implementation. Extremely
close distinct roots can therefore merge in the exported accepted list.
Root coordinates and q values in JSON/CSV are not arbitrary-precision output.

The public reference entry point restores the caller's mpmath precision after
a solve. The underlying mpmath context is shared: use separate processes,
not concurrent threads, for parallel reference solves. Lower-level MP helper
functions still set the global precision themselves.

The symbolic no-root-loss theorem does not certify the completeness of a
finite-precision numerical solve. Small residuals alone do not certify root
accuracy near repeated or ill-conditioned roots. For publication-sensitive
points, compare spectra at increased precision and report unresolved or
marginal cases. The current implementation is not an interval-certified
root-isolation solver.

Critical-Ek map searches use float64 coarse scans and MP bracket refinement.
A finite coarse scan can miss narrow transitions; transition tolerance is
not a bound on all possible unresolved crossings. Likewise, the least-stable
mode amplitude figures illustrate a selected modal response, not a complete
initial-value solution with all modal coefficients determined.

## Validation evidence

Tests compare Python polynomial coefficients to frozen Mathematica exports,
check dimensional equivalence, principal branches, domain exclusions,
reduced degree, refinement failure handling, conjugacy, medium exchange, and
reference precision convergence. Production checkpoint and summary tests
cover reproducibility plumbing separately from numerical correctness.
See `docs/RELEASE_VALIDATION.md` for the environment and executed checks.
