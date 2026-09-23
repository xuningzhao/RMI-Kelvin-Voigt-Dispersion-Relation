# Phase I symbolic framework review

This note reviews the Phase I Mathematica notebook as a future publication
artifact. It does not begin Phase II and does not derive the polynomial
representation.

## Overall assessment

The Phase I framework now verifies the nondimensionalization

$$
D(\gamma)\longrightarrow D^*_\Lambda(s)
$$

in a reproducible way and records the nonzero scaling factor

$$
D^*_\Lambda=\frac{\rho_T}{2k^2}D(\gamma).
$$

The strongest feature is that the notebook and the batch audit share the same
symbolic definitions through `verification/Phase1Verification.wl`, reducing the
risk of notebook/script drift.

## Issues found in the first exported PDF

The first PDF was mathematically correct in its final check, but it read too
much like an implementation transcript:

- The first audit cell printed the entire audit association, overwhelming the
  derivation.
- Mathematica emitted harmless `CreateDirectory::eexist` messages when
  `logs/` and `exports/` already existed.
- The final equivalence check was present, but the intermediate algebraic
  transformations were not sufficiently decomposed for publication-level audit.
- Domain information was present but not yet separated into original expression
  domain, nondimensional simplification assumptions, and future algebraic
  manipulation assumptions.

## Improvements made

The Phase I notebook and audit framework were revised so that:

- The audit load is silent in the notebook and explicit summaries are displayed
  afterward.
- Existing log/export directories no longer generate warning messages.
- The notebook now checks the transformation component by component:
  - Kelvin--Voigt effective-viscosity scaling;
  - radicand mapping;
  - principal square-root mapping;
  - reciprocal-denominator scaling;
  - constant-term scaling;
  - reciprocal-term scaling;
  - whole-expression equivalence.
- The original expression domain is recorded separately from the
  nondimensional assumptions.
- The final PASS/FAIL summary now requires all component checks to pass, not
  only the final expression-level simplification.

## Recommendations before Phase II

## Theorem-scope clarification

The Phase I symbolic theorem is established under the explicit assumption
$\Lambda>0$. The limiting case $\Lambda=0$ (equivalently $E_k=1$) corresponds,
within the chosen nondimensional parameterization, to both shear moduli
vanishing simultaneously. This limiting point is outside the scope of the
current theorem and should be analyzed separately if it is needed. This is a
scope statement, not a change to the verified nondimensionalization.

Before deriving the polynomial representation, the framework should keep the
following conventions:

1. Do not clear denominators without recording every excluded denominator.
2. Do not square equations without recording the implication direction and the
   possible spurious-root mechanism.
3. Track radical signs and branch selection before eliminating square roots.
4. Preserve a table of assumptions accumulated by each transformation.
5. Export a compact proof summary after each phase, separate from raw symbolic
   expressions.

With the current Phase I revisions, the framework is ready to be re-executed
locally and re-exported as a cleaner PDF before Phase II begins.
