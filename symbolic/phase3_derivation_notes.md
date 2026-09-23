# Phase III polynomial-completeness audit

Phase III starts after Phase II is frozen. Its purpose is to package and audit
the final polynomial-completeness theorem: the degree-14 polynomial candidate
equation does not lose any root of the original nondimensional dispersion
relation.

Phase II is the authoritative source for the radical elimination,
denominator-clearing steps, coefficient verification, and forward-preserving
symbolic identities. Phase III reuses those verified Phase II checks rather
than repeating the full derivation under a second name.

The central theorem is the one-way implication

$$
\Omega_{\mathrm{definition}}
\land
D^*_\Lambda(s)=0
\Longrightarrow
P_{14}(s)=0.
$$

Equivalently, every root of the original dispersion relation is represented
among the roots of $P_{14}$ under the stated definition-domain and derivation
assumptions.

Phase III does **not** claim the reverse implication

$$
P_{14}(s)=0
\Longrightarrow
D^*_\Lambda(s)=0.
$$

Extra polynomial roots are expected because the derivation clears denominators
and squares equations. Those candidates are filtered later by direct
substitution into the original dispersion relation.

## Definition domain

The original nondimensional expression requires:

$$
s\ne0,
$$

$$
E_+\ne0,\qquad E_-\ne0,
$$

and the two reciprocal denominators to be nonzero.

These conditions define $\Omega_{\mathrm{definition}}$, the domain on which
$D^*_\Lambda(s)$ is mathematically well-defined.

## No-root-loss conditions

The no-root-loss audit records only the assumptions required to preserve the
forward implication from the original equation to the polynomial equation.

Important distinctions:

- multiplying by an expression cannot lose a root, although it may add
  extraneous roots;
- squaring cannot lose a root satisfying the pre-squared equation, although it
  may add extraneous roots;
- taking the numerator of a rational expression requires the rational
  denominator to be nonzero at the root;
- dividing by the nonzero constant normalization factor $4$ is harmless.

Branch or reverse-equivalence conditions are not part of the main Phase III
theorem.

## Forward-step audit

The audited Phase II chain is

$$
D^*_\Lambda(s)=0
\Longrightarrow
F_0=0
\Longrightarrow
F_1=0
\Longrightarrow
F_2=0
\Longrightarrow
F_3=0
\Longrightarrow
F_4=0
\Longrightarrow
P_{14}(s)=0.
$$

Each step is classified in Phase II by whether it preserves the forward
implication and whether it can introduce extraneous polynomial candidates.
Phase III imports those classifications into the final theorem summary.

## Runtime root-selection workflow

The later numerical workflow is:

1. solve $P_{14}(s)=0$;
2. obtain polynomial candidate roots;
3. reject candidates outside $\Omega_{\mathrm{definition}}$;
4. substitute each remaining candidate into $D^*_\Lambda(s)$;
5. retain candidates satisfying the original dispersion relation within
   tolerance;
6. apply physical admissibility criteria later.

Physical admissibility is intentionally outside this polynomial-completeness
proof.

## Stop boundary

Phase III does not implement the production polynomial-root filtering
algorithm, does not construct physical spectra, and does not claim reverse
equivalence.
