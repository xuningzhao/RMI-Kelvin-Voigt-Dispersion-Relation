# Polynomial-candidate root verification workflow

This document defines the mathematical terminology and implementation contract
for the production Python workflow. It starts from the verified symbolic result

$$
D^*_\Lambda(s)=0
\Longrightarrow
P_{14}(s)=0
$$

under the recorded definition and nonzero assumptions. The reverse implication
is not guaranteed. Therefore roots of $P_{14}$ are candidates only and must be
filtered by direct substitution into the original nondimensional dispersion
relation.

## Source equations and branch convention

The source equation remains the nondimensional dispersion relation

$$
D^*_\Lambda(s;A_\rho,A_\mu,A_G,\Lambda)=0.
$$

The polynomial candidate equation is

$$
P_{14}(s;A_\rho,A_\mu,A_G,\Lambda)=0.
$$

The polynomial is used only to generate a finite candidate set. The original
dispersion relation is the authority for deciding whether a candidate is a
genuine root.

All square roots use the principal complex square-root branch, matching both:

- Mathematica `Sqrt` in the symbolic derivation; and
- NumPy `sqrt` in `src/dispersion.py`.

## Compact quantities

Define

$$
E_+=s(1+A_\mu)+\Lambda(1+A_G),
$$

$$
E_-=s(1-A_\mu)+\Lambda(1-A_G).
$$

The radicands are

$$
R_-=
1+\frac{(1-A_\rho)s^2}{E_-},
$$

$$
R_+=
1+\frac{(1+A_\rho)s^2}{E_+}.
$$

The principal square-root quantities are

$$
q_-=\sqrt{R_-},
\qquad
q_+=\sqrt{R_+}.
$$

The reciprocal denominators used in the polynomial derivation are

$$
B_- = E_+ + E_-q_-,
$$

$$
B_+ = E_- + E_+q_+.
$$

## Root terminology

### Polynomial candidate root

A polynomial candidate root is any finite complex root of the structurally
effective polynomial equation

$$
P_{\mathrm{eff}}(s)=0,
$$

where $P_{\mathrm{eff}}$ is obtained from $P_{14}$ after structural
effective-degree detection. Leading coefficients are not discarded merely
because they are small relative to lower-order coefficients.

A polynomial candidate root is not automatically a root of the original
dispersion relation.

### Mathematically genuine root

A mathematically genuine root is a polynomial candidate root that:

1. satisfies the recorded mathematical definition and nonzero assumptions;
2. evaluates the original nondimensional dispersion relation without poles,
   invalid values, or branch inconsistency; and
3. satisfies

   $$
   D^*_\Lambda(s)=0
   $$

   within the documented numerical residual tolerance.

Thus every retained mathematically genuine root is checked in the source
equation, not merely in the polynomial equation.

### Physically admissible root

A physically admissible root is a mathematically genuine root that also
satisfies the far-field spatial decay condition of the semi-infinite
normal-mode solution.

The spatial fields contain $e^{q_-y}$ in the lower domain and $e^{-q_+y}$ in
the upper domain. Decay away from the interface therefore requires

$$
\operatorname{Re}(q_-)>0,
\qquad
\operatorname{Re}(q_+)>0.
$$

The quantities $q_-$ and $q_+$ must be computed from the same principal
square-root convention used in the symbolic derivation and the source
dispersion evaluator.

## Temporal behavior versus spatial admissibility

Temporal behavior is determined by the real part of the response rate:

$$
\operatorname{Re}(s).
$$

For a physically admissible root:

- $\operatorname{Re}(s)>0$ indicates temporal growth;
- $\operatorname{Re}(s)<0$ indicates temporal decay;
- $\operatorname{Re}(s)$ near zero indicates temporal marginality.

Spatial admissibility is determined by

$$
\operatorname{Re}(q_-),
\qquad
\operatorname{Re}(q_+).
$$

These are different concepts and must never be conflated. In particular,
$\operatorname{Re}(s)$ must not be used to accept or reject a root as
physically admissible. It only classifies temporal stability after the root has
passed mathematical genuineness and spatial decay checks.

## Complete workflow

For a prescribed parameter tuple

$$
(A_\rho,A_\mu,A_G,\Lambda),
$$

the implementation should perform the following steps.

### 1. Check theorem-scope parameters

The verified symbolic theorem assumes

$$
\Lambda>0,
$$

$$
-1<A_\rho<1,
\qquad
-1<A_\mu<1,
\qquad
-1\le A_G\le1.
$$

The limiting case $\Lambda=0$, equivalently $E_k=1$, lies outside the current
symbolic theorem and should be classified separately unless a dedicated limit
analysis is added.

This is the parameter-level theorem scope. Candidate-level theorem conditions
such as $s\ne0$, $E_\pm\ne0$, and $B_\pm\ne0$ are checked separately for each
polynomial candidate by the definition-domain filter.

### 2. Evaluate polynomial coefficients

Evaluate the verified coefficient list

$$
P_{14}(s)=\sum_{j=0}^{14}c_j s^j
$$

at the given parameter tuple.

The exported Mathematica coefficient list is ordered low-to-high:

$$
\{c_0,c_1,\ldots,c_{14}\}.
$$

Most numerical root routines expect high-to-low ordering, so the implementation
must reverse the list before calling such routines.

### 3. Detect effective degree

The generic degree is

$$
\deg_s P_{14}=14,
$$

with leading coefficient

$$
c_{14}=(A_\rho+A_\mu)^2.
$$

The degree may drop when

$$
A_\rho+A_\mu=0.
$$

The implementation must preserve degree 14 whenever the leading coefficient is
structurally nonzero. In particular, the degree must not be reduced merely
because lower-order coefficients are much larger in high-$\Lambda$ cases.
Reduced degree is not an error by itself, but it must be reported and should
correspond to structural cancellation, such as $A_\rho+A_\mu=0$.

For numerical conditioning, the Python implementation solves a scaled-variable
polynomial. It chooses a positive scale $\alpha$, solves in $z$ with

$$
s=\alpha z,
$$

and maps roots back to $s$. This improves conditioning without deleting
mathematically essential coefficients.

### 4. Compute finite polynomial roots

Compute all finite roots of the effective polynomial. If the effective degree
is zero or the coefficient vector is numerically unusable, classify the
parameter point as degenerate or unresolved rather than forcing roots.

Duplicate roots produced by numerical roundoff should be clustered using a
documented distance tolerance.

### 5. Apply definition and nonzero checks

For each polynomial candidate root, check the mathematical definition domain
and no-root-loss denominator conditions. At minimum, reject or mark as
domain-failed candidates for which any of the following are zero or
numerically too close to zero:

$$
s,
\qquad
E_+,
\qquad
E_-,
\qquad
B_-,
\qquad
B_+.
$$

These checks are numerical versions of the recorded symbolic assumptions.

### 6. Verify against the source dispersion relation

Evaluate

$$
D^*_\Lambda(s;A_\rho,A_\mu,A_G,\Lambda)
$$

directly using the source implementation. A candidate becomes a mathematically
genuine root only if the source residual passes the documented tolerance.

Polynomial residual alone is insufficient for accepting a root.

### 6a. Targeted high-precision refinement for ambiguous candidates

The default workflow is double precision. If a candidate lands in the
mathematical residual ambiguity band, the implementation may optionally invoke
a targeted high-precision fallback. This fallback is applied only to candidates
whose double-precision verification status is `numerically_ambiguous`.

The fallback:

1. recomputes the P14 coefficients at arbitrary precision;
2. locally refines the nearby polynomial root;
3. evaluates the original dispersion relation at the same high precision using
   the principal square-root branch; and
4. resolves the candidate as genuine, spurious, or still unresolved.

Confidently accepted or rejected double-precision candidates are not altered.
This refinement improves root accuracy for clustered polynomial roots without
relaxing the global residual, domain, or physical-admissibility tolerances.

### 7. Compute spatial decay quantities

For each mathematically genuine root, compute

$$
q_-=\sqrt{1+\frac{(1-A_\rho)s^2}{E_-}},
$$

$$
q_+=\sqrt{1+\frac{(1+A_\rho)s^2}{E_+}},
$$

using the principal square-root branch.

### 8. Apply physical admissibility

Classify a mathematically genuine root as physically admissible if

$$
\operatorname{Re}(q_-)>0,
\qquad
\operatorname{Re}(q_+)>0.
$$

If either real part lies within numerical tolerance of zero, the root should be
classified as a spatial-decay boundary or marginal-admissibility case rather
than forced into admissible or inadmissible.

## Recommended parameter-point classification

Each parameter point should receive one primary classification, with detailed
counts stored separately.

Implemented categories:

- `outside_theorem_scope`: parameters violate the current symbolic theorem
  assumptions, for example $\Lambda\le0$.
- `identically_zero_polynomial`: after effective-degree detection, the
  polynomial is numerically indistinguishable from zero.
- `constant_nonzero_polynomial`: after effective-degree detection, the
  polynomial is nonzero constant and therefore has no finite polynomial roots.
- `only_domain_invalid_candidates`: all polynomial candidates fail definition or
  nonzero checks.
- `only_spurious_candidates`: candidates pass domain checks but fail the source
  residual test for $D^*_\Lambda$.
- `genuine_roots_nonadmissible`: at least one mathematically genuine root
  exists, but none satisfy the spatial decay conditions.
- `marginal_spatial_decay`: at least one mathematically genuine root lies
  within tolerance of $\operatorname{Re}(q_-)=0$ or
  $\operatorname{Re}(q_+)=0$, and no confidently admissible root is present.
- `physically_admissible_roots`: at least one mathematically genuine root
  satisfies both spatial decay inequalities.
- `numerically_unresolved`: numerical failures or tolerance conflicts prevent a
  reliable classification.

Suggested precedence:

1. `outside_theorem_scope`
2. `identically_zero_polynomial`
3. `constant_nonzero_polynomial`
4. `numerically_unresolved` for ambiguous mathematical residuals, evaluation
   failures, missing domain intermediates, or incomplete polynomial roots
5. `only_domain_invalid_candidates`
6. `only_spurious_candidates`
7. `physically_admissible_roots`
8. `marginal_spatial_decay`
9. `genuine_roots_nonadmissible`

The implementation stores the underlying counts:

- polynomial candidate count;
- domain-passing candidate count;
- mathematically genuine root count;
- physically admissible root count;
- marginal spatial-decay count;
- rejected spurious candidate count.

## Numerical tolerance policy

The Python implementation defines separate, scale-aware tolerances for:

- structural effective-degree detection;
- candidate duplicate clustering;
- domain/nonzero checks for $s$, $E_\pm$, and $B_\pm$;
- polynomial residual diagnostics;
- source residual acceptance for $D^*_\Lambda$;
- spatial decay boundary checks for $\operatorname{Re}(q_\pm)$.

The implemented defaults use separate named tolerances rather than one global
value:

```text
polynomial_relative_tolerance = 1e-12
polynomial_absolute_tolerance = 0
domain_relative_tolerance = 1e-10
domain_absolute_tolerance = 1e-14
residual_relative_tolerance = 1e-8
residual_absolute_tolerance = 1e-10
q_relative_tolerance = 1e-10
q_absolute_tolerance = 1e-12
enable_high_precision_refinement = False
high_precision_initial_decimal_precision = 80
high_precision_maximum_decimal_precision = 120
high_precision_precision_step = 20
```

For spatial admissibility, the implemented scale-aware boundary tolerances are

$$
\varepsilon_{q_\pm}
=
\varepsilon_{\mathrm{abs}}
+
\tau_q\max(1,|q_\pm|),
$$

with $\tau_q=10^{-10}$ and
$\varepsilon_{\mathrm{abs}}=10^{-12}$ by default. Then:

- admissible if
  $\operatorname{Re}(q_+)>\varepsilon_{q_+}$ and
  $\operatorname{Re}(q_-)>\varepsilon_{q_-}$;
- marginal if either real part lies within its corresponding tolerance of
  zero;
- otherwise nonadmissible.
