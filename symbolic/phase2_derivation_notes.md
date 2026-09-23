# Phase II derivation notes

Phase II starts from the frozen Phase I theorem:

$$
D^*_\Lambda(s)=\frac{\rho_T}{2k^2}D(\gamma).
$$

The goal is not to prove equivalence between the polynomial and the original
dispersion relation. The Phase II theorem is only the one-way implication

$$
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0
$$

on the admissible domain.

## Compact notation

Define

$$
E_+=s(1+A_\mu)+\Lambda(1+A_G),
$$

$$
E_-=s(1-A_\mu)+\Lambda(1-A_G).
$$

Then

$$
\eta_+=\frac{E_+}{s},\qquad
\eta_-=\frac{E_-}{s}.
$$

The two radicands are

$$
R_- = 1+\frac{(1-A_\rho)s^2}{E_-},
\qquad
R_+ = 1+\frac{(1+A_\rho)s^2}{E_+},
$$

with principal square roots

$$
q_-=\sqrt{R_-},\qquad q_+=\sqrt{R_+}.
$$

The reciprocal denominators in the $s^2$ clock form are

$$
B_- = E_+ + E_-q_-,
\qquad
B_+ = E_- + E_+q_+.
$$

Thus

$$
F_0
=s^2\left(\frac1{B_-}+\frac1{B_+}\right)+2.
$$

## Admissible-domain records

The original nondimensional expression requires

$$
s\ne0,\qquad E_+\ne0,\qquad E_-\ne0,
$$

and

$$
B_-\ne0,\qquad B_+\ne0.
$$

The physical Phase II assumptions are inherited from Phase I:

$$
k>0,\quad \rho_T>0,\quad \mu_T>0,\quad \Lambda>0,
$$

$$
-1<A_\rho<1,\qquad -1<A_\mu<1,\qquad -1\le A_G\le1,
$$

$$
s\in\mathbb C,\qquad s\ne0.
$$

## Step chain

Clearing the reciprocal denominators gives

$$
F_1=s^2(B_-+B_+)+2B_-B_+=0.
$$

This multiplication is valid in the forward direction on the admissible domain.
It may introduce spurious roots if used backward.

Collecting radicals gives

$$
F_1=\alpha+\beta q_-+\gamma q_+ +\delta q_-q_+=0,
$$

where

$$
\alpha=s^2(E_++E_-)+2E_+E_-,
$$

$$
\beta=E_-(s^2+2E_-),
\qquad
\gamma=E_+(s^2+2E_+),
$$

$$
\delta=2E_+E_-.
$$

Isolate $q_+$ and square:

$$
F_2=(\alpha+\beta q_-)^2
-R_+(\gamma+\delta q_-)^2=0.
$$

This preserves forward implication but may introduce spurious roots.

Write

$$
F_2=u_0+u_1q_-+u_2q_-^2.
$$

Using $q_-^2=R_-$ gives

$$
F_3=u_0+u_1q_-+u_2R_-=0.
$$

Then isolate $q_-$ and square:

$$
F_4=(u_0+u_2R_-)^2-u_1^2R_-=0.
$$

This again preserves forward implication but may introduce spurious roots.

Finally clear rational denominators. The factor

$$
E_+^2E_-^2
$$

is a denominator-clearing/admissibility factor and is not divided out of the
final polynomial. Only the harmless nonzero constant factor $4$ is removed from
the cleared numerator. The resulting polynomial is the candidate

$$
P_{14}(s;A_\rho,A_\mu,A_G,\Lambda)=0.
$$

## Generic degree

The Phase II symbolic audit is designed to confirm:

$$
\deg_s P_{14}=14,
$$

with leading coefficient

$$
c_{14}=(A_\rho+A_\mu)^2.
$$

Therefore the degree may drop on the parameter hypersurface

$$
A_\rho+A_\mu=0.
$$

Additional lower-degree cancellations on special parameter strata should be
treated separately in a later admissible-domain analysis.

## Stop boundary

Phase II does not claim

$$
P_{14}(s)=0\Longrightarrow D^*_\Lambda(s)=0.
$$

It also does not implement a polynomial-root filtering algorithm. Those tasks
belong to later phases.
