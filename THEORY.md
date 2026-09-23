# Theory and numerical methods

This reference combines the numerical-method documentation, derivations, and
frozen symbolic verification records. The first section describes the current
solver. Later review and audit records retain their original scientific scope;
installation and release instructions are maintained in [README.md](README.md).

- [Numerical methods and limitations](#numerical-methods)
- [Nondimensional formulation](#specification)
- [Detailed root verification](#root-verification)
- [Algebraic derivation](#derivation)
- [Buckingham analysis](#buckingham)
- [Nondimensional dispersion relation](#dispersion)
- [Symbolic verification framework](#symbolic-verification)
- [Phase I review record](#phase1-review)
- [Phase II derivation record](#phase2-derivation)
- [Phase III derivation record](#phase3-derivation)
- [Canonical exports](#canonical-exports)
- [Frozen Phase II report](#phase2-report)
- [Frozen Phase III report](#phase3-report)


<a id="numerical-methods"></a>

## Numerical methods and limitations

### Numerical solver methods

#### Scope and coordinates

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

#### Candidate generation

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

#### Verification and physical acceptance

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

#### Counts, precision, and limitations

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

#### Validation evidence

Tests compare Python polynomial coefficients to frozen Mathematica exports,
check dimensional equivalence, principal branches, domain exclusions,
reduced degree, refinement failure handling, conjugacy, medium exchange, and
reference precision convergence. Production checkpoint and summary tests
cover reproducibility plumbing separately from numerical correctness.
See `README.md` for the environment and executed checks.


<a id="specification"></a>

## Nondimensional formulation

### Specification: canonical nondimensionalization

#### Physical problem

Two Kelvin–Voigt media meet at a perturbed interface. For medium $`i`$, let
$`\rho_i`$ be the solid density, $`\mu_i`$ the dynamic viscosity, and $`G_i`$
the shear modulus. Let $`k`$ be the interfacial wavenumber and $`\gamma`$ the
possibly complex temporal growth rate. Define

```math
\rho_T=\rho_1+\rho_2,\qquad
\mu_T=\mu_1+\mu_2,\qquad
G_T=G_1+G_2.
```

The derivation assumes $`k>0`$, $`\rho_i>0`$, and, for the chosen canonical
clock, $`\mu_T>0`$. An individual modulus may vanish. The normalized
$`g_i`$ and $`A_G`$ coordinates assume $`G_T>0`$; the purely viscous point
$`G_1=G_2=0`$ is their continuous $`\Lambda\to0`$ boundary, where $`A_G`$ is
irrelevant and not identifiable. Whenever division by an individual
property is used, that property is additionally assumed positive. Square-root
branches are the same branches as in the dimensional dispersion relation;
nondimensionalization does not select or alter them.

#### Original dimensional dispersion relation

An independent transcription of the supplied equation is

```math
\mathcal D(\gamma,k)={}\gamma\left[
\frac{1}{
\mu_1+G_1/\gamma+
(\mu_2+G_2/\gamma)
\sqrt{1+\frac{\rho_2\gamma}{(\mu_2+G_2/\gamma)k^2}}}
+\frac{1}{
\mu_2+G_2/\gamma+
(\mu_1+G_1/\gamma)
\sqrt{1+\frac{\rho_1\gamma}{(\mu_1+G_1/\gamma)k^2}}}
\right]
+\frac{4k^2}{\rho_T}=0.
```

This is identical to the supplied relation: superscripted medium labels have
only been changed to subscripts, and the second denominator is the exact
$`1\leftrightarrow2`$ image of the first. No term, sign, factor, or radical has
changed.

The eight dimensional quantities and their dimensions are:

| Quantity | Meaning | Dimensions |
|---|---|---|
| $`\gamma`$ | response/growth rate | $`T^{-1}`$ |
| $`k`$ | wavenumber | $`L^{-1}`$ |
| $`\rho_1,\rho_2`$ | densities | $`M L^{-3}`$ |
| $`\mu_1,\mu_2`$ | dynamic viscosities | $`M L^{-1}T^{-1}`$ |
| $`G_1,G_2`$ | shear moduli | $`M L^{-1}T^{-2}`$ |

The two addends in $`\mathcal D`$ both have dimension $`L/M`$. In particular,
$`[\gamma/\mu]=L/M`$ and $`[k^2/\rho]=L/M`$.

#### Buckingham $`\Pi`$ analysis

There are $`n=8`$ dimensional quantities. Their exponent vectors in the
$`(M,L,T)`$ basis form the matrix (columns are ordered as
$`\gamma,k,\rho_1,\rho_2,\mu_1,\mu_2,G_1,G_2`$)

```math
\mathbf D=
\begin{pmatrix}
0&0&1&1&1&1&1&1\\
0&-1&-3&-3&-1&-1&-1&-1\\
-1&0&0&0&-1&-1&-2&-2
\end{pmatrix}.
```

For example, the columns
$`[\rho]=(1,-3,0)`$, $`[\mu]=(1,-1,-1)`$, and
$`[k]=(0,-1,0)`$ form a $`3\times3`$ matrix with determinant $`-1`$.
Thus $`\mathrm{rank}\,\mathbf D\ge3`$; since there are only three base
dimensions, $`\mathrm{rank}\,\mathbf D=3`$. Buckingham's theorem gives
$`n-r=8-3=5`$ independent groups. One is the dimensionless response, leaving
four independent control parameters.

A strict monomial basis, using $`\rho_2,\mu_2,k`$ as repeating variables, is

```math
\frac{\rho_2\gamma}{\mu_2k^2},\quad
\frac{\rho_1}{\rho_2},\quad
\frac{\mu_1}{\mu_2},\quad
\frac{\rho_2G_1}{\mu_2^2k^2},\quad
\frac{\rho_2G_2}{\mu_2^2k^2}.
```

It is valid but privileges medium 2. A symmetric and physically clearer set of
coordinates is

```math
s_v=\frac{\rho_T\gamma}{\mu_Tk^2},\qquad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2},\qquad
A_\rho,\ A_\mu,\ A_G,
```

where

```math
A_x=\frac{x_1-x_2}{x_1+x_2},\qquad
x_1/x_T=\frac{1+A_x}{2},\quad
x_2/x_T=\frac{1-A_x}{2}.
```

These five coordinates are invertibly related to the monomial basis in the
positive-property interior, so they are independent. Because sums and bounded
contrasts are not monomials in the original variables, this is best called a
natural symmetric $`\Pi`$-coordinate system rather than *the* unique
Buckingham basis.

Other equally complete bases replace $`(s_v,\Lambda)`$ by

```math
(s_e,\chi),\qquad (s_R,\Lambda),\qquad
\chi=\frac{r_e}{r_v}=\sqrt{\Lambda},
```

with the three contrasts unchanged.

#### Characteristic time scales

The total-property rates and times are defined exactly by

```math
r_v=\frac{\mu_Tk^2}{\rho_T},\quad t_v=r_v^{-1},\qquad
r_e=k\sqrt{\frac{G_T}{\rho_T}},\quad t_e=r_e^{-1},
```

and the Kelvin–Voigt crossover or relaxation time is

```math
t_R=\frac{\mu_T}{G_T}.
```

They are not three independent clocks:

```math
t_R=\frac{t_e^2}{t_v},\qquad
\chi=\frac{r_e}{r_v}=\frac{t_v}{t_e},\qquad
\Lambda=\chi^2.
```

The proposed bounded rate contrast naturally accompanies the rate-sum clock,

```math
E_k=\frac{r_v-r_e}{r_v+r_e}=\frac{1-\chi}{1+\chi},\qquad
t_+=\frac{1}{r_v+r_e}.
```

For nonnegative material properties, $`\chi\in[0,\infty)`$ maps to
$`E_k\in(-1,1]`$; the lower endpoint is approached as $`\chi\to\infty`$.

Constituent clocks also occur:

```math
t_{v,i}=\frac{\rho_i}{\mu_i k^2},\qquad
t_{e,i}=\frac1k\sqrt{\frac{\rho_i}{G_i}},\qquad
t_{R,i}=\frac{\mu_i}{G_i}.
```

They introduce no new dimensional freedom: their ratios to total clocks are
functions of $`\Lambda`$ and the contrasts. No capillary, gravitational,
acoustic, or imposed shock time appears because the corresponding dimensional
quantity is absent from the given dispersion relation. Infinitely many clocks
can be manufactured as $`t_v f(\Lambda,A_\rho,A_\mu,A_G)`$; they are
reparameterizations, not new balances.

#### Candidate nondimensionalizations

Set

```math
r_i=\rho_i/\rho_T,\qquad m_i=\mu_i/\mu_T,\qquad
g_i=G_i/G_T,
```

so each pair sums to one. For any clock $`t_c`$, define

```math
s=\gamma t_c,\quad p=\frac{t_v}{t_c},\quad
q=\frac{G_Tt_c}{\mu_T},\quad H_i=m_i+\frac{qg_i}{s}.
```

Direct substitution and multiplication by $`\rho_T/k^2`$ gives the universal
dimensionless relation

```math
p s\left[
\frac1{H_1+H_2\sqrt{1+p r_2s/H_2}}+
\frac1{H_2+H_1\sqrt{1+p r_1s/H_1}}
\right]+4=0,                                      \tag{1}
```

with $`pq=\Lambda`$. Equivalently, for $`s\ne0`$, let
$`a_i=m_is+qg_i`$. Then

```math
p s^2\left[
\frac1{a_1+a_2\sqrt{1+p r_2s^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+p r_1s^2/a_1}}
\right]+4=0.                                      \tag{2}
```

The un-cleared form (1) is the direct identity; (2) is convenient away from
the original singular representation at $`\gamma=0`$.

The choices are:

| Clock | $`s`$ | $`p`$ | $`q`$ |
|---|---:|---:|---:|
| viscous $`t_v`$ | $`s_v`$ | $`1`$ | $`\Lambda`$ |
| elastic $`t_e`$ | $`s_e`$ | $`\chi`$ | $`\chi`$ |
| relaxation $`t_R`$ | $`s_R`$ | $`\Lambda`$ | $`1`$ |
| rate sum $`t_+`$ | $`s_+`$ | $`1+\chi`$ | $`\chi^2/(1+\chi)`$ |

Thus (2) yields, respectively (the sums below contain the two ordered pairs
$`(i,j)=(1,2),(2,1)`$),

```math
s_v^2\sum_{i\ne j}\frac1{a_i+a_j\sqrt{1+r_js_v^2/a_j}}+4=0,
\quad a_i=m_is_v+\Lambda g_i,                     \tag{V}
```

```math
\chi s_e^2\sum_{i\ne j}\frac1{a_i+a_j\sqrt{1+\chi r_js_e^2/a_j}}+4=0,
\quad a_i=m_is_e+\chi g_i,                       \tag{E}
```

```math
\Lambda s_R^2\sum_{i\ne j}\frac1{a_i+a_j\sqrt{1+\Lambda r_js_R^2/a_j}}+4=0,
\quad a_i=m_is_R+g_i.                            \tag{R}
```

For the rate-sum form, substitute its $`p,q`$ from the table into (2). In
terms of $`E_k`$,

```math
p=\frac{2}{1+E_k},\qquad
q=\frac{(1-E_k)^2}{2(1+E_k)}.                    \tag{+}
```

Writing $`p_+=1+\chi`$, $`q_+=\chi^2/(1+\chi)`$, and
$`a_i=m_is_++q_+g_i`$, its equation is explicitly

```math
p_+s_+^2\left[
\frac1{a_1+a_2\sqrt{1+p_+r_2s_+^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+p_+r_1s_+^2/a_1}}
\right]+4=0.                                    \tag{S}
```

For a constituent clock, no fresh derivation is needed because (1)–(2) are
general. The relevant substitutions are

| $`t_c`$ | $`p`$ | $`q`$ |
| --- | --- | --- |
| $`t_{v,j}`$ | $`m_j/r_j`$ | $`\Lambda r_j/m_j`$ |
| $`t_{e,j}`$ | $`\chi\sqrt{g_j/r_j}`$ | $`\chi\sqrt{r_j/g_j}`$ |
| $`t_{R,j}`$ | $`\Lambda g_j/m_j`$ | $`m_j/g_j`$ |

The response variables obey

```math
s_e=\frac{s_v}{\chi},\qquad
s_R=\frac{s_v}{\Lambda},\qquad
s_+=\frac{s_v}{1+\chi}.
```

#### Comparison

| Form | Algebra | Symmetry | Controls | Interpretation and limits |
|---|---|---|---:|---|
| viscous | simplest outer and radical coefficients | exchange-symmetric | 4 | direct diffusion clock; regular as $`G_T\to0`$ |
| elastic | balanced $`p=q=\chi`$ | exchange-symmetric | 4 | direct wave clock; undefined at $`G_T=0`$ |
| relaxation | normalizes elastic/viscous crossover | exchange-symmetric | 4 | useful for constitutive relaxation; undefined at $`G_T=0`$, poor long-wave scaling |
| rate sum | two coefficient functions | exchange-symmetric | 4 | bounded $`E_k`$, but no primitive balance and more algebra |
| constituent | extra contrast factors | privileges one medium | 4 | useful only when one medium is a reference; singular if its selected property vanishes |

All total-property formulations have the same number of independent
parameters. A change of clock cannot reduce the Buckingham count; it merely
moves powers of $`\chi`$ or $`\Lambda`$ between the response and coefficients.
All are invariant under medium exchange when
$`(r_1,m_1,g_1)\leftrightarrow(r_2,m_2,g_2)`$, equivalently when all three
contrasts change sign.

#### Status of $`A_\rho,A_\mu,A_G,E_k`$

The four quantities are a valid, independent set of **control coordinates**
for positive properties. Adding one response such as $`s_v`$ completes the
five groups required by Buckingham's theorem. They are not a unique or strict
monomial $`\Pi`$ basis:

1. Buckingham bases are nonunique.
2. Each $`A_x`$ is a bounded transform of the ratio $`x_1/x_2`$.
3. $`E_k`$ is a bounded transform of $`\chi=\sqrt\Lambda`$:
   $`\chi=(1-E_k)/(1+E_k)`$.
4. Choosing $`E_k`$ does not itself specify the dimensionless response clock.

The contrasts arise naturally from exchange symmetry and total-property
normalization, not from dimensional analysis alone. $`E_k`$ is convenient for
bounded plots, while $`\Lambda`$ is more canonical algebraically because it is
the direct dimensionless coefficient generated by the viscous scaling.

#### Recommended formulation

Adopt (V) with

```math
\boxed{
s=\frac{\rho_T\gamma}{\mu_Tk^2},\quad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2},\quad
(A_\rho,A_\mu,A_G)
}
```

and reconstruct fractions with $`(1\pm A_x)/2`$. This formulation is
dimensionally complete, exchange-symmetric, minimal in parameter count,
algebraically shortest, and continuous into the purely viscous limit. At the
exact point $`G_T=0`$, use $`\Lambda=0`$; $`A_G`$ can be omitted because it
has no physical effect there. Store
or compute $`\chi=\sqrt\Lambda`$ and $`E_k=(1-\chi)/(1+\chi)`$ only as derived
coordinates. If a bounded sweep coordinate is desirable, sweeping $`E_k`$ is
perfectly legitimate, but the governing specification should remain in
$`\Lambda`$.

The elastic form (E) is an equally valid secondary presentation when wave
physics is the focus. It is not recommended as the repository's defining form
because it loses its clock at zero total modulus.

### Phase II: Parameter-space stability analysis

> **Current implementation note.** This section records the original numerical
> exploration plan written before the symbolic polynomial-completeness theorem
> was finalized. The production workflow is now the polynomial-candidate method
> defined in [root-verification reference](THEORY.md#root-verification):
> solve the verified $`P_{14}`$ polynomial, filter candidates by the original
> nondimensional dispersion relation, and then apply the $`q_\pm`$ spatial-decay
> admissibility conditions. Finite-budget direct searches remain useful only as
> exploratory diagnostics.

#### Scope and scientific objective

Phase II begins from the final nondimensional dispersion relation documented in
[derivation/final_nondimensional_dispersion.md](THEORY.md#dispersion).
The Phase I theory and its recommended viscous-clock nondimensionalization
remain the mathematical foundation of the project.

The primary objective is to determine how stability is organized over the
nondimensional material-parameter space. The canonical coordinates are

```math
s,\qquad A_\rho,\qquad A_\mu,\qquad A_G,\qquad \Lambda,
```

and the bounded sweep-coordinate representation is

```math
s,\qquad A_\rho,\qquad A_\mu,\qquad A_G,\qquad E_k,
```

with

```math
\Lambda=\left(\frac{1-E_k}{1+E_k}\right)^2.
```

The ultimate scientific objective is not merely to compute roots. It is to
understand how density contrast, viscosity contrast, elastic contrast, and
elastic–viscous rate competition organize the RMI stability behavior.

#### Numerical problem

For each prescribed parameter set, solve

```math
\mathcal D_\Lambda^*(s;A_\rho,A_\mu,A_G,\Lambda)=0
```

or, equivalently,

```math
\mathcal D_{E_k}^*(s;A_\rho,A_\mu,A_G,E_k)=0,
```

where $`s\in\mathbb C`$. Within each declared search domain, identify all
verified roots when possible. If the verified roots are
$`\{s_i\}_{i=1}^N`$, identify the dominant root by

```math
s_{\mathrm{dom}}
\in\underset{s_i}{\mathrm{arg\,max}}\ \mathrm{Re}(s_i).
```

If several roots share the same largest real part within the stated numerical
tolerance, report the complete tied set rather than silently selecting one.
Classify stability from the dominant real part:

```math
\begin{array}{ll}
\mathrm{Re}(s_{\mathrm{dom}})>0
&\text{unstable},\\
\mathrm{Re}(s_{\mathrm{dom}})<0
&\text{stable or decaying},\\
\mathrm{Re}(s_{\mathrm{dom}})=0
&\text{neutral or marginal}.
\end{array}
```

Numerically, these comparisons must use a documented stability tolerance.
Cases lying within that tolerance of zero must be labeled marginal rather than
assigned a sign that is not resolved by the computation.

#### Required outputs

Each parameter study must produce, where applicable:

- all verified roots within the declared search domain;
- the dominant root $`s_{\mathrm{dom}}`$;
- the dominant growth rate $`\mathrm{Re}(s_{\mathrm{dom}})`$;
- the dominant frequency $`\mathrm{Im}(s_{\mathrm{dom}})`$;
- root trajectories in the complex $`s`$-plane under parameter continuation;
- stability maps over selected nondimensional parameter planes; and
- diagnostic reports for failed or unresolved cases.

Every reported result must retain the parameter values, search domain,
square-root branch convention, solver settings, verification tolerances, and
classification status needed to reproduce it.

#### Root-finding outcome categories

Every attempted parameter set must receive one of the following numerical
outcome categories:

- **root_found**: at least one candidate root passes the required residual and
  consistency checks;
- **likely_no_root**: a documented and sufficiently broad search found no
  verified root in the declared search domain, and the diagnostics support
  that conclusion; or
- **unresolved**: the search did not establish either a verified root or a
  sufficiently supported no-root result.

These are numerical classifications, not automatic physical conclusions.
Failure of a root-finding method is not evidence that the dispersion relation
has no roots. In particular, nonconvergence, branch inconsistency, inadequate
initial guesses, insufficient search extent, and loss of numerical precision
must lead to explicit diagnostics and may require the **unresolved** category.
The category **likely_no_root** is always relative to the documented search
domain and verification procedure.

#### Numerical principles

The Phase II implementation must follow these principles:

1. The final nondimensional dispersion relation is the source of truth.
2. Every candidate root must be verified by evaluating the residual
   $`|\mathcal D^*(s)|`$ in the source equation. The absolute and relative
   residual tolerances must be stated.
3. Transformed, cleared, rationalized, squared, or polynomialized equations may
   be used only as search aids or diagnostics. Any candidate they produce must
   be checked in the source equation because such transformations may add,
   remove, or obscure roots.
4. Numerical failures and incomplete searches must be reported explicitly.
5. No root may be forced, inferred from solver termination alone, or retained
   merely because an iteration returned a finite value.
6. Square-root branches must be selected and tracked consistently. Branch
   changes, branch-cut encounters, and continuation discontinuities must be
   diagnosed rather than silently absorbed.
7. Duplicate roots found from multiple initial guesses must be clustered using
   a documented tolerance, while genuinely distinct nearby roots must be
   preserved.
8. Search completeness must always be stated relative to a specified region
   of the complex $`s`$-plane; it must not be claimed globally without a
   supporting mathematical argument.

#### Recommended first numerical implementation

The first numerical implementation should:

1. implement `D_star_Lambda(s, Arho, Amu, AG, Lambda)` directly from
   $`\mathcal D_\Lambda^*`$;
2. implement `D_star_Ek(s, Arho, Amu, AG, Ek)` directly from
   $`\mathcal D_{E_k}^*`$;
3. verify over representative admissible parameters and complex values of $`s`$
   that

   ```math
   \mathcal{D}_{E_k}^{*}(s,A_\rho,A_\mu,A_G,E_k)
   =\mathcal{D}_{\Lambda}^{*}\left(
   s,A_\rho,A_\mu,A_G,
   \left(\frac{1-E_k}{1+E_k}\right)^2
   \right);
   ```

4. keep compact helper variables internal to the implementation, using them
   only to evaluate the final mathematical definition rather than replacing
   it; and
5. test medium-exchange symmetry, limiting cases, branch consistency, and
   residual evaluation before introducing parameter sweeps.

This section specifies the future implementation only. No root solver is
introduced in this phase of the documentation update.

#### Initial parameter-exploration plan

Begin with one-dimensional sweeps in $`E_k`$ at fixed
$`A_\rho`$, $`A_\mu`$, and $`A_G`$. These sweeps should establish reliable root
seeding, branch tracking, dominant-root switching behavior, and diagnostic
reporting before moving to higher-dimensional studies.

After the one-dimensional workflow is verified, construct two-dimensional
maps including:

- $`E_k`$ versus $`A_\rho`$;
- $`E_k`$ versus $`A_\mu`$;
- $`E_k`$ versus $`A_G`$; and
- $`A_\rho`$ versus $`A_G`$ at fixed $`E_k`$.

For every sweep or map, record the fixed parameters, sampled ranges,
continuation direction, search domain, root classification, and unresolved
regions. Stability boundaries should be inferred only from verified dominant
roots and should retain uncertainty information where the numerical
classification changes or remains unresolved.


<a id="root-verification"></a>

## Detailed root verification

### Polynomial-candidate root verification workflow

This document defines the mathematical terminology and implementation contract
for the production Python workflow. It starts from the verified symbolic result

```math
D^*_\Lambda(s)=0
\Longrightarrow
P_{14}(s)=0
```

under the recorded definition and nonzero assumptions. The reverse implication
is not guaranteed. Therefore roots of $`P_{14}`$ are candidates only and must be
filtered by direct substitution into the original nondimensional dispersion
relation.

#### Source equations and branch convention

The source equation remains the nondimensional dispersion relation

```math
D^*_\Lambda(s;A_\rho,A_\mu,A_G,\Lambda)=0.
```

The polynomial candidate equation is

```math
P_{14}(s;A_\rho,A_\mu,A_G,\Lambda)=0.
```

The polynomial is used only to generate a finite candidate set. The original
dispersion relation is the authority for deciding whether a candidate is a
genuine root.

All square roots use the principal complex square-root branch, matching both:

- Mathematica `Sqrt` in the symbolic derivation; and
- NumPy `sqrt` in `src/dispersion.py`.

#### Compact quantities

Define

```math
E_+=s(1+A_\mu)+\Lambda(1+A_G),
```

```math
E_-=s(1-A_\mu)+\Lambda(1-A_G).
```

The radicands are

```math
R_-=
1+\frac{(1-A_\rho)s^2}{E_-},
```

```math
R_+=
1+\frac{(1+A_\rho)s^2}{E_+}.
```

The principal square-root quantities are

```math
q_-=\sqrt{R_-},
\qquad
q_+=\sqrt{R_+}.
```

The reciprocal denominators used in the polynomial derivation are

```math
B_- = E_+ + E_-q_-,
```

```math
B_+ = E_- + E_+q_+.
```

#### Root terminology

##### Polynomial candidate root

A polynomial candidate root is any finite complex root of the structurally
effective polynomial equation

```math
P_{\mathrm{eff}}(s)=0,
```

where $`P_{\mathrm{eff}}`$ is obtained from $`P_{14}`$ after structural
effective-degree detection. Leading coefficients are not discarded merely
because they are small relative to lower-order coefficients.

A polynomial candidate root is not automatically a root of the original
dispersion relation.

##### Mathematically genuine root

A mathematically genuine root is a polynomial candidate root that:

1. satisfies the recorded mathematical definition and nonzero assumptions;
2. evaluates the original nondimensional dispersion relation without poles,
   invalid values, or branch inconsistency; and
3. satisfies

   ```math
   D^*_\Lambda(s)=0
   ```

   within the documented numerical residual tolerance.

Thus every retained mathematically genuine root is checked in the source
equation, not merely in the polynomial equation.

##### Physically admissible root

A physically admissible root is a mathematically genuine root that also
satisfies the far-field spatial decay condition of the semi-infinite
normal-mode solution.

The spatial fields contain $`e^{q_-y}`$ in the lower domain and $`e^{-q_+y}`$ in
the upper domain. Decay away from the interface therefore requires

```math
\mathrm{Re}(q_-)>0,
\qquad
\mathrm{Re}(q_+)>0.
```

The quantities $`q_-`$ and $`q_+`$ must be computed from the same principal
square-root convention used in the symbolic derivation and the source
dispersion evaluator.

#### Temporal behavior versus spatial admissibility

Temporal behavior is determined by the real part of the response rate:

```math
\mathrm{Re}(s).
```

For a physically admissible root:

- $`\mathrm{Re}(s)>0`$ indicates temporal growth;
- $`\mathrm{Re}(s)<0`$ indicates temporal decay;
- $`\mathrm{Re}(s)`$ near zero indicates temporal marginality.

Spatial admissibility is determined by

```math
\mathrm{Re}(q_-),
\qquad
\mathrm{Re}(q_+).
```

These are different concepts and must never be conflated. In particular,
$`\mathrm{Re}(s)`$ must not be used to accept or reject a root as
physically admissible. It only classifies temporal stability after the root has
passed mathematical genuineness and spatial decay checks.

#### Complete workflow

For a prescribed parameter tuple

```math
(A_\rho,A_\mu,A_G,\Lambda),
```

the implementation should perform the following steps.

##### 1. Check theorem-scope parameters

The verified symbolic theorem assumes

```math
\Lambda>0,
```

```math
-1<A_\rho<1,
\qquad
-1<A_\mu<1,
\qquad
-1\le A_G\le1.
```

The limiting case $`\Lambda=0`$, equivalently $`E_k=1`$, lies outside the current
symbolic theorem and should be classified separately unless a dedicated limit
analysis is added.

This is the parameter-level theorem scope. Candidate-level theorem conditions
such as $`s\ne0`$, $`E_\pm\ne0`$, and $`B_\pm\ne0`$ are checked separately for each
polynomial candidate by the definition-domain filter.

##### 2. Evaluate polynomial coefficients

Evaluate the verified coefficient list

```math
P_{14}(s)=\sum_{j=0}^{14}c_j s^j
```

at the given parameter tuple.

The exported Mathematica coefficient list is ordered low-to-high:

```math
\{c_0,c_1,\ldots,c_{14}\}.
```

Most numerical root routines expect high-to-low ordering, so the implementation
must reverse the list before calling such routines.

##### 3. Detect effective degree

The generic degree is

```math
\deg_s P_{14}=14,
```

with leading coefficient

```math
c_{14}=(A_\rho+A_\mu)^2.
```

The degree may drop when

```math
A_\rho+A_\mu=0.
```

The implementation must preserve degree 14 whenever the leading coefficient is
structurally nonzero. In particular, the degree must not be reduced merely
because lower-order coefficients are much larger in high-$`\Lambda`$ cases.
Reduced degree is not an error by itself, but it must be reported and should
correspond to structural cancellation, such as $`A_\rho+A_\mu=0`$.

For numerical conditioning, the Python implementation solves a scaled-variable
polynomial. It chooses a positive scale $`\alpha`$, solves in $`z`$ with

```math
s=\alpha z,
```

and maps roots back to $`s`$. This improves conditioning without deleting
mathematically essential coefficients.

##### 4. Compute finite polynomial roots

Compute all finite roots of the effective polynomial. If the effective degree
is zero or the coefficient vector is numerically unusable, classify the
parameter point as degenerate or unresolved rather than forcing roots.

Duplicate roots produced by numerical roundoff should be clustered using a
documented distance tolerance.

##### 5. Apply definition and nonzero checks

For each polynomial candidate root, check the mathematical definition domain
and no-root-loss denominator conditions. At minimum, reject or mark as
domain-failed candidates for which any of the following are zero or
numerically too close to zero:

```math
s,
\qquad
E_+,
\qquad
E_-,
\qquad
B_-,
\qquad
B_+.
```

These checks are numerical versions of the recorded symbolic assumptions.

##### 6. Verify against the source dispersion relation

Evaluate

```math
D^*_\Lambda(s;A_\rho,A_\mu,A_G,\Lambda)
```

directly using the source implementation. A candidate becomes a mathematically
genuine root only if the source residual passes the documented tolerance.

Polynomial residual alone is insufficient for accepting a root.

##### 6a. Targeted high-precision refinement for ambiguous candidates

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

##### 7. Compute spatial decay quantities

For each mathematically genuine root, compute

```math
q_-=\sqrt{1+\frac{(1-A_\rho)s^2}{E_-}},
```

```math
q_+=\sqrt{1+\frac{(1+A_\rho)s^2}{E_+}},
```

using the principal square-root branch.

##### 8. Apply physical admissibility

Classify a mathematically genuine root as physically admissible if

```math
\mathrm{Re}(q_-)>0,
\qquad
\mathrm{Re}(q_+)>0.
```

If either real part lies within numerical tolerance of zero, the root should be
classified as a spatial-decay boundary or marginal-admissibility case rather
than forced into admissible or inadmissible.

#### Recommended parameter-point classification

Each parameter point should receive one primary classification, with detailed
counts stored separately.

Implemented categories:

- `outside_theorem_scope`: parameters violate the current symbolic theorem
  assumptions, for example $`\Lambda\le0`$.
- `identically_zero_polynomial`: after effective-degree detection, the
  polynomial is numerically indistinguishable from zero.
- `constant_nonzero_polynomial`: after effective-degree detection, the
  polynomial is nonzero constant and therefore has no finite polynomial roots.
- `only_domain_invalid_candidates`: all polynomial candidates fail definition or
  nonzero checks.
- `only_spurious_candidates`: candidates pass domain checks but fail the source
  residual test for $`D^*_\Lambda`$.
- `genuine_roots_nonadmissible`: at least one mathematically genuine root
  exists, but none satisfy the spatial decay conditions.
- `marginal_spatial_decay`: at least one mathematically genuine root lies
  within tolerance of $`\mathrm{Re}(q_-)=0`$ or
  $`\mathrm{Re}(q_+)=0`$, and no confidently admissible root is present.
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

#### Numerical tolerance policy

The Python implementation defines separate, scale-aware tolerances for:

- structural effective-degree detection;
- candidate duplicate clustering;
- domain/nonzero checks for $`s`$, $`E_\pm`$, and $`B_\pm`$;
- polynomial residual diagnostics;
- source residual acceptance for $`D^*_\Lambda`$;
- spatial decay boundary checks for $`\mathrm{Re}(q_\pm)`$.

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

```math
\varepsilon_{q_\pm}
=
\varepsilon_{\mathrm{abs}}
+
\tau_q\max(1,|q_\pm|),
```

with $`\tau_q=10^{-10}`$ and
$`\varepsilon_{\mathrm{abs}}=10^{-12}`$ by default. Then:

- admissible if
  $`\mathrm{Re}(q_+)>\varepsilon_{q_+}`$ and
  $`\mathrm{Re}(q_-)>\varepsilon_{q_-}`$;
- marginal if either real part lies within its corresponding tolerance of
  zero;
- otherwise nonadmissible.


<a id="derivation"></a>

## Algebraic derivation

### Detailed derivation

#### 1. Verified starting equation

Write the Kelvin–Voigt effective viscosity as

```math
\eta_i(\gamma)=\mu_i+\frac{G_i}{\gamma}.
```

The supplied equation is exactly

```math
\gamma\left[
\frac1{\eta_1+\eta_2\sqrt{1+\rho_2\gamma/(\eta_2k^2)}}+
\frac1{\eta_2+\eta_1\sqrt{1+\rho_1\gamma/(\eta_1k^2)}}
\right]+\frac{4k^2}{\rho_T}=0.                  \tag{D}
```

Expanding each $`\eta_i`$ reproduces the original equation term by term.

#### 2. General change of variables

For an arbitrary nonzero characteristic time $`t_c`$, introduce

```math
\gamma=\frac{s}{t_c},\quad
\rho_i=\rho_T r_i,\quad
\mu_i=\mu_Tm_i,\quad
G_i=G_Tg_i.
```

Define

```math
p=\frac{\rho_T}{\mu_Tk^2t_c}=\frac{t_v}{t_c},
\qquad q=\frac{G_Tt_c}{\mu_T}.
```

Then, without omitting an algebraic step,

```math
\begin{aligned}
\eta_i
&=\mu_Tm_i+\frac{G_Tg_i}{s/t_c}\\
&=\mu_T\left(m_i+\frac{G_Tt_c}{\mu_T}\frac{g_i}{s}\right)\\
&=\mu_TH_i,
\qquad H_i=m_i+\frac{qg_i}{s},
\end{aligned}
```

and

```math
\begin{aligned}
1+\frac{\rho_i\gamma}{\eta_i k^2}
&=1+\frac{\rho_T r_i(s/t_c)}{\mu_TH_i k^2}\\
&=1+p\frac{r_is}{H_i}.
\end{aligned}
```

The first reciprocal denominator in (D) becomes

```math
\frac1{\eta_1+\eta_2\sqrt{1+\rho_2\gamma/(\eta_2k^2)}}
=\frac1{\mu_T}
\frac1{H_1+H_2\sqrt{1+pr_2s/H_2}}.
```

The other term follows by exchanging 1 and 2. Finally multiply (D) by
$`\rho_T/k^2`$. Since

```math
\frac{\rho_T}{k^2}\frac{\gamma}{\mu_T}
=\frac{\rho_Ts}{\mu_Tk^2t_c}=ps,
```

the exact general dimensionless equation is

```math
ps\left[
\frac1{H_1+H_2\sqrt{1+pr_2s/H_2}}+
\frac1{H_2+H_1\sqrt{1+pr_1s/H_1}}
\right]+4=0.                                    \tag{G1}
```

Also,

```math
pq=\frac{\rho_TG_T}{\mu_T^2k^2}=\Lambda,
```

so only one independent rate-ratio parameter is present.

For a form with no explicit $`1/s`$, define

```math
a_i=sH_i=m_is+qg_i.
```

For $`s\ne0`$,

```math
H_i=\frac{a_i}{s},\qquad
1+\frac{pr_is}{H_i}=1+\frac{pr_is^2}{a_i},
```

and

```math
\frac1{H_i+H_j\sqrt{1+pr_js/H_j}}
=\frac{s}{a_i+a_j\sqrt{1+pr_js^2/a_j}}.
```

Thus

```math
ps^2\left[
\frac1{a_1+a_2\sqrt{1+pr_2s^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+pr_1s^2/a_1}}
\right]+4=0.                                    \tag{G2}
```

#### 3. Viscous clock

Take

```math
t_c=t_v=\frac{\rho_T}{\mu_Tk^2}.
```

Then

```math
p=\frac{t_v}{t_v}=1,qquad
q=\frac{G_T}{\mu_T}\frac{\rho_T}{\mu_Tk^2}
=\Lambda.
```

In (G1),

```math
H_i=m_i+\frac{\Lambda g_i}{s_v},
```

giving

```math
s_v\left[
\frac1{H_1+H_2\sqrt{1+r_2s_v/H_2}}+
\frac1{H_2+H_1\sqrt{1+r_1s_v/H_1}}
\right]+4=0.
```

Equivalently, with $`a_i=m_is_v+\Lambda g_i`$,

```math
s_v^2\left[
\frac1{a_1+a_2\sqrt{1+r_2s_v^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+r_1s_v^2/a_1}}
\right]+4=0.                                    \tag{V}
```

#### 4. Elastic-wave clock

Take

```math
t_c=t_e=\frac1k\sqrt{\frac{\rho_T}{G_T}},\qquad
\chi=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}=\sqrt\Lambda.
```

Direct calculation gives

```math
p=\frac{\rho_T}{\mu_Tk^2t_e}
=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}=\chi,
```

and

```math
q=\frac{G_Tt_e}{\mu_T}
=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}=\chi.
```

Consequently $`a_i=m_is_e+\chi g_i`$, and (G2) becomes

```math
\chi s_e^2\left[
\frac1{a_1+a_2\sqrt{1+\chi r_2s_e^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+\chi r_1s_e^2/a_1}}
\right]+4=0.                                    \tag{E}
```

#### 5. Kelvin–Voigt relaxation clock

Take $`t_c=t_R=\mu_T/G_T`$. Then

```math
p=\frac{\rho_T}{\mu_Tk^2}\frac{G_T}{\mu_T}
=\Lambda,
\qquad q=\frac{G_T}{\mu_T}\frac{\mu_T}{G_T}=1.
```

Thus $`a_i=m_is_R+g_i`$, and

```math
\Lambda s_R^2\left[
\frac1{a_1+a_2\sqrt{1+\Lambda r_2s_R^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+\Lambda r_1s_R^2/a_1}}
\right]+4=0.                                    \tag{R}
```

The identity $`t_R=t_e^2/t_v`$ follows immediately from the definitions.

#### 6. Rate-sum clock and $`E_k`$

Let

```math
t_+=\frac1{r_v+r_e},\quad \chi=\frac{r_e}{r_v},\quad
E_k=\frac{1-\chi}{1+\chi}.
```

Because $`t_+=t_v/(1+\chi)`$,

```math
p=\frac{t_v}{t_+}=1+\chi,qquad
q=\frac{\Lambda}{p}=\frac{\chi^2}{1+\chi}.
```

Solving the definition of $`E_k`$ gives

```math
\chi=\frac{1-E_k}{1+E_k},quad
1+\chi=\frac2{1+E_k},quad
\frac{\chi^2}{1+\chi}=\frac{(1-E_k)^2}{2(1+E_k)}.
```

Substitution of these $`p,q`$ and
$`a_i=m_is_++qg_i`$ into (G2) gives explicitly

```math
(1+\chi)s_+^2\left[
\frac1{a_1+a_2\sqrt{1+(1+\chi)r_2s_+^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+(1+\chi)r_1s_+^2/a_1}}
\right]+4=0,
```

where $`a_i=m_is_++\chi^2g_i/(1+\chi)`$. It has no fewer independent
controls than (V); it only bounds the rate ratio.

#### 7. Constituent clocks

From $`r_i=\rho_i/\rho_T`$, $`m_i=\mu_i/\mu_T`$, and $`g_i=G_i/G_T`$,

```math
t_{v,j}=t_v\frac{r_j}{m_j},\qquad
t_{e,j}=t_e\sqrt{\frac{r_j}{g_j}},\qquad
t_{R,j}=t_R\frac{m_j}{g_j}.
```

Using $`p=t_v/t_c`$, $`q=\Lambda/p`$, these yield

| $`t_c`$ | $`p`$ | $`q`$ |
| --- | --- | --- |
| $`t_{v,j}`$ | $`m_j/r_j`$ | $`\Lambda r_j/m_j`$ |
| $`t_{e,j}`$ | $`\chi\sqrt{g_j/r_j}`$ | $`\chi\sqrt{r_j/g_j}`$ |
| $`t_{R,j}`$ | $`\Lambda g_j/m_j`$ | $`m_j/g_j`$ |

Putting any row into (G1) or (G2) gives its full dimensionless relation.
These clocks privilege medium $`j`$ and add contrast factors without reducing
the number of independent controls.

#### 8. Contrast reconstruction and symmetry

For $`x\in\{\rho,\mu,G\}`$,

```math
x_1=\frac{x_T}{2}(1+A_x),\qquad
x_2=\frac{x_T}{2}(1-A_x).
```

Therefore

```math
r_{1,2}=\frac{1\pm A_\rho}{2},\quad
m_{1,2}=\frac{1\pm A_\mu}{2},\quad
g_{1,2}=\frac{1\pm A_G}{2}.
```

Exchanging the media sends all three contrasts to their negatives and swaps
the two reciprocal terms, leaving every total-property equation invariant.


<a id="buckingham"></a>

## Buckingham analysis

### Mathematical structure of the Buckingham $`\Pi`$ space

#### Scope and relation to the existing derivation

This document studies the linear-algebraic content of Buckingham's theorem for
the dimensional variables

```math
(\gamma,k,\rho_1,\rho_2,\mu_1,\mu_2,G_1,G_2).
```

It does not change the dispersion relation or the recommended
nondimensionalization in [nondimensional formulation](THEORY.md#specification). The detailed substitution
into the dispersion relation remains in [DERIVATION.md](THEORY.md#derivation). No
repeating variables are used to obtain the null space below; repeating
variables are introduced only afterward to explain the classical construction.

For the geometric interpretation, the dimensional variables are first viewed
on a positive coordinate chart so that arbitrary real powers and logarithms
are well defined. The exponent algebra itself is formal and remains valid for
a nonzero complex response $`\gamma`$; the dispersion relation supplies that
analytic continuation. Zero material-property values can be included afterward
by taking appropriate limits.

#### 1. Dimensional matrix

Use base-dimension order $`(M,L,T)`$ and variable order

```math
\boldsymbol{x}=(\gamma,k,\rho_1,\rho_2,\mu_1,\mu_2,G_1,G_2).
```

Every column is the exponent vector of one variable:

| Variable | Physical dimension | Column of $`D`$ |
| --- | --- | --- |
| $`\gamma`$ | $`T^{-1}`$ | $`(0,0,-1)^\mathsf T`$ |
| $`k`$ | $`L^{-1}`$ | $`(0,-1,0)^\mathsf T`$ |
| $`\rho_1`$ | $`ML^{-3}`$ | $`(1,-3,0)^\mathsf T`$ |
| $`\rho_2`$ | $`ML^{-3}`$ | $`(1,-3,0)^\mathsf T`$ |
| $`\mu_1`$ | $`ML^{-1}T^{-1}`$ | $`(1,-1,-1)^\mathsf T`$ |
| $`\mu_2`$ | $`ML^{-1}T^{-1}`$ | $`(1,-1,-1)^\mathsf T`$ |
| $`G_1`$ | $`ML^{-1}T^{-2}`$ | $`(1,-1,-2)^\mathsf T`$ |
| $`G_2`$ | $`ML^{-1}T^{-2}`$ | $`(1,-1,-2)^\mathsf T`$ |

Therefore

```math
D=
\begin{pmatrix}
0&0&1&1&1&1&1&1\\
0&-1&-3&-3&-1&-1&-1&-1\\
-1&0&0&0&-1&-1&-2&-2
\end{pmatrix}.                                      \tag{1}
```

For an exponent vector

```math
\boldsymbol a=(a_\gamma,a_k,a_{\rho_1},a_{\rho_2},
a_{\mu_1},a_{\mu_2},a_{G_1},a_{G_2})^\mathsf T,
```

the monomial

```math
\Pi_{\boldsymbol a}
=\gamma^{a_\gamma}k^{a_k}
\rho_1^{a_{\rho_1}}\rho_2^{a_{\rho_2}}
\mu_1^{a_{\mu_1}}\mu_2^{a_{\mu_2}}
G_1^{a_{G_1}}G_2^{a_{G_2}}                         \tag{2}
```

is dimensionless exactly when $`D\boldsymbol a=0`$.

#### 2. Direct symbolic solution of $`D\boldsymbol a=0`$

Writing out the three rows of (1) gives

```math
\begin{aligned}
a_{\rho_1}+a_{\rho_2}+a_{\mu_1}+a_{\mu_2}+a_{G_1}+a_{G_2}&=0,\\
-a_k-3a_{\rho_1}-3a_{\rho_2}-a_{\mu_1}-a_{\mu_2}-a_{G_1}-a_{G_2}&=0,\\
-a_\gamma-a_{\mu_1}-a_{\mu_2}-2a_{G_1}-2a_{G_2}&=0.
\end{aligned}                                      \tag{3}
```

Row reduction, performed without selecting repeating variables, gives

```math
\mathrm{rref}(D)=
\begin{pmatrix}
1&0&0&0&1&1&2&2\\
0&1&0&0&-2&-2&-2&-2\\
0&0&1&1&1&1&1&1
\end{pmatrix}.                                      \tag{4}
```

The pivot variables are
$`a_\gamma,a_k,a_{\rho_1}`$. Introduce arbitrary free parameters

```math
c_1=a_{\rho_2},\quad c_2=a_{\mu_1},\quad
c_3=a_{\mu_2},\quad c_4=a_{G_1},\quad c_5=a_{G_2}.
```

Equation (4) then gives

```math
\begin{aligned}
a_\gamma&=-c_2-c_3-2c_4-2c_5,\\
a_k&=2c_2+2c_3+2c_4+2c_5,\\
a_{\rho_1}&=-c_1-c_2-c_3-c_4-c_5,\\
a_{\rho_2}&=c_1,\quad a_{\mu_1}=c_2,\quad
a_{\mu_2}=c_3,\quad a_{G_1}=c_4,\quad a_{G_2}=c_5.
\end{aligned}                                      \tag{5}
```

Equivalently,

```math
\boldsymbol a=c_1\boldsymbol n_1+c_2\boldsymbol n_2
+c_3\boldsymbol n_3+c_4\boldsymbol n_4+c_5\boldsymbol n_5, \tag{6}
```

where

```math
\begin{aligned}
\boldsymbol n_1&=(0,0,-1,1,0,0,0,0)^\mathsf T,\\
\boldsymbol n_2&=(-1,2,-1,0,1,0,0,0)^\mathsf T,\\
\boldsymbol n_3&=(-1,2,-1,0,0,1,0,0)^\mathsf T,\\
\boldsymbol n_4&=(-2,2,-1,0,0,0,1,0)^\mathsf T,\\
\boldsymbol n_5&=(-2,2,-1,0,0,0,0,1)^\mathsf T.
\end{aligned}                                      \tag{7}
```

Direct multiplication verifies $`D\boldsymbol n_j=0`$ for all five vectors.
Because each vector has a unique unit entry in one of the five free-variable
positions, they are linearly independent.

#### 3. What the five-dimensional null space means

The columns for $`\rho,\mu,k`$, for example, contain the nonsingular minor

```math
\det
\begin{pmatrix}
1&1&0\\
-3&-1&-1\\
0&-1&0
\end{pmatrix}=-1.
```

Hence $`\mathrm{rank}\,D\ge3`$. There are only three base dimensions, so
$`\mathrm{rank}\,D=3`$, and rank-nullity gives

```math
\dim\ker D=8-3=5.                                  \tag{8}
```

The five dimensions admit a particularly useful physical decomposition.
Since each material pair has identical columns,

```math
\begin{aligned}
\boldsymbol d_\rho&=(0,0,1,-1,0,0,0,0)^\mathsf T,\\
\boldsymbol d_\mu&=(0,0,0,0,1,-1,0,0)^\mathsf T,\\
\boldsymbol d_G&=(0,0,0,0,0,0,1,-1)^\mathsf T
\end{aligned}                                      \tag{9}
```

are automatically in $`\ker D`$. They generate the dimensionless ratios
$`\rho_1/\rho_2`$, $`\mu_1/\mu_2`$, and $`G_1/G_2`$. Thus three of the five
directions encode how each property is partitioned between the media. The
remaining two directions may be chosen as a dimensionless temporal response
and one elastic-to-viscous balance.

This explains why a property pair does not disappear under
nondimensionalization. Changing units multiplies both members of a pair by the
same conversion factor. It can change their numerical magnitudes, but it
cannot change their ratio. Algebraically, equal columns imply that their
column difference is zero; geometrically, unit rescaling acts along the common
direction and leaves the relative direction untouched.

More formally, let

```math
K_-=\mathrm{span}\,\{\boldsymbol d_\rho,
\boldsymbol d_\mu,\boldsymbol d_G\}.
```

Then $`K_-\subset\ker D`$, $`\dim K_-=3`$, and the quotient
$`\ker D/K_-`$ is two-dimensional. The assignment “three contrasts, one
response, one rate balance” is therefore a natural decomposition, although
the particular coordinates chosen within each part remain nonunique.

##### Vector space versus invariant coordinates

The precise vector space produced by linear algebra is the exponent space
$`\ker D`$. Its addition law corresponds to multiplication of monomials:

```math
\Pi_{\boldsymbol a+\boldsymbol b}
=\Pi_{\boldsymbol a}\Pi_{\boldsymbol b},\qquad
\Pi_{c\boldsymbol a}=\Pi_{\boldsymbol a}^{c}.       \tag{10}
```

Once any five independent monomial groups are known, every dimensionless
quantity can be written locally as a function of them. Such functions do not
themselves form the same linear exponent space. It is therefore useful to
distinguish:

1. the five-dimensional linear space $`\ker D`$ of monomial exponents;
2. a chosen vector-space basis of $`\ker D`$;
3. the five-dimensional quotient or invariant space, on which arbitrary
   invertible nonlinear coordinates may be used.

Buckingham's theorem fixes the dimension of this space, not a preferred basis
or preferred coordinates.

#### 4. A basis obtained directly from row reduction

Using (7) in (2) gives the direct null-space basis

```math
\boxed{
N_1=\frac{\rho_2}{\rho_1},\quad
N_2=\frac{\mu_1k^2}{\gamma\rho_1},\quad
N_3=\frac{\mu_2k^2}{\gamma\rho_1},\quad
N_4=\frac{G_1k^2}{\gamma^2\rho_1},\quad
N_5=\frac{G_2k^2}{\gamma^2\rho_1}.}                \tag{11}
```

This basis is a mechanically correct output of a particular row-reduction
convention. It is not physically canonical. A different column order changes
the pivot columns and therefore changes the basis returned by the same
algorithm. Moreover, any nonsingular matrix $`C\in GL(5,\mathbb R)`$ produces
another basis

```math
B=NC,                                               \tag{12}
```

where $`N`$ is the $`8\times5`$ matrix with columns $`\boldsymbol n_j`$.
At the group level,

```math
\widetilde\Pi_j=\prod_{i=1}^5N_i^{C_{ij}}.          \tag{13}
```

Thus nonuniqueness is the ordinary nonuniqueness of a basis in a
five-dimensional vector space.

#### 5. Repeating variables as a basis selection

The classical construction chooses three dimensionally independent repeating
variables. Choose $`(\rho_2,\mu_2,k)`$. Their dimension columns form an
invertible $`3\times3`$ matrix. For each remaining variable, solving for the
three repeating-variable exponents constructs one kernel vector. The result is

```math
\begin{aligned}
P&=\frac{\rho_2\gamma}{\mu_2k^2},&
R_\rho&=\frac{\rho_1}{\rho_2},&
R_\mu&=\frac{\mu_1}{\mu_2},\\
Q_1&=\frac{\rho_2G_1}{\mu_2^2k^2},&
Q_2&=\frac{\rho_2G_2}{\mu_2^2k^2}.&&
\end{aligned}\tag{14}
```

This construction has not found a different invariant space. Indeed, the
groups in (14) are explicit products of the direct groups in (11):

```math
P=\frac{N_1}{N_3},\qquad
R_\rho=N_1^{-1},\qquad
R_\mu=\frac{N_2}{N_3},\qquad
Q_1=\frac{N_1N_4}{N_3^2},\qquad
Q_2=\frac{N_1N_5}{N_3^2}.                          \tag{15}
```

The corresponding change-of-basis matrix is

```math
C=
\begin{pmatrix}
1&-1&0&1&1\\
0&0&1&0&0\\
-1&0&-1&-2&-2\\
0&0&0&1&0\\
0&0&0&0&1
\end{pmatrix},\qquad \det C=1.                     \tag{16}
```

Because $`C`$ is invertible, the exponent vectors in (14) and (11) span
exactly the same kernel. In general, choosing repeating variables amounts to
choosing an invertible $`3\times3`$ column minor of $`D`$, then using the
remaining five variables as free coordinates. It is a convenient algorithm
for selecting one kernel basis, not an additional theorem and not a uniqueness
principle.

#### 6. Several valid bases and coordinate systems

##### 6.1 Direct row-reduction basis

Equation (11) is a strict monomial basis. Its advantages are that it follows
immediately from linear algebra and makes no preliminary physical choice. Its
disadvantages are that the result depends on column ordering, privileges
$`\rho_1`$, mixes $`\gamma`$ into four groups, and obscures exchange symmetry.

##### 6.2 Classical repeating-variable basis

Equation (14) is also a strict monomial basis. It has integer exponents and
isolates one response group. It is convenient when medium 2 is a genuine
reference material. It privileges medium 2, becomes unsuitable if a selected
repeating property vanishes, and does not display exchange symmetry.

##### 6.3 Reference-scale ratio basis

Replacing $`Q_1`$ in (14) by $`R_G=G_1/G_2=Q_1/Q_2`$ gives

```math
\boxed{P,\ R_\rho,\ R_\mu,\ R_G,\ Q_2}.            \tag{17}
```

This is a strict monomial basis because the replacement is an invertible basis
change. It cleanly exposes one degree of freedom for every material pair,
plus one response and one elastic-strength group. It still uses medium 2 as
the dimensional reference, and all three ratios are unbounded.

##### 6.4 Exchange-adapted geometric-mean basis

Define geometric means

```math
\rho_g=\sqrt{\rho_1\rho_2},\qquad
\mu_g=\sqrt{\mu_1\mu_2},\qquad
G_g=\sqrt{G_1G_2}.
```

Then

```math
\boxed{
P_g=\frac{\gamma\rho_g}{\mu_gk^2},\quad
L_g=\frac{\rho_gG_g}{\mu_g^2k^2},\quad
R_\rho,\quad R_\mu,\quad R_G}                     \tag{18}
```

is a strict monomial basis, allowing half-integer exponents. Direct
substitution of its five exponent vectors into (1) gives zero, and their
$`8\times5`$ exponent matrix has rank five.

Under the exchange operator $`1\leftrightarrow2`$,

```math
P_g\mapsto P_g,\quad L_g\mapsto L_g,\quad
R_x\mapsto R_x^{-1}.                               \tag{19}
```

Equivalently, $`\log R_x\mapsto-\log R_x`$. This basis diagonalizes exchange
symmetry into a two-dimensional even sector and a three-dimensional odd
sector in logarithmic coordinates. It is the most exchange-adapted *monomial*
basis. Its costs are fractional powers, unbounded ratios, and loss of a useful
coordinate when a property vanishes.

##### 6.5 Total-property and contrast coordinates

Let

```math
\rho_T=\rho_1+\rho_2,\quad
\mu_T=\mu_1+\mu_2,\quad
G_T=G_1+G_2,
```

and define

```math
s=\frac{\rho_T\gamma}{\mu_Tk^2},\qquad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2},\qquad
A_x=\frac{x_1-x_2}{x_1+x_2}.                      \tag{20}
```

The five quantities

```math
\boxed{s,\ \Lambda,\ A_\rho,\ A_\mu,\ A_G}       \tag{21}
```

are dimensionless and independent on the positive-property interior. They are
not monomials in the eight original variables because they contain sums.
Consequently, they do not correspond to five exponent vectors and are not a
vector-space basis of $`\ker D`$. They are an invertible nonlinear coordinate
system on the same five-dimensional invariant space. They preserve exchange
symmetry in the especially transparent form

```math
s\mapsto s,\quad\Lambda\mapsto\Lambda,\quad
(A_\rho,A_\mu,A_G)\mapsto(-A_\rho,-A_\mu,-A_G).    \tag{22}
```

These coordinates are adapted to the sums that actually occur in the
dispersion relation. The contrasts are bounded for positive properties. At
the exactly purely viscous point $`G_T=0`$, $`A_G`$ is unidentifiable but also
physically irrelevant, as explained in the existing specification.

#### 7. Ratios and contrasts are invertible coordinates

For any positive pair $`(x_1,x_2)`$, let

```math
R_x=\frac{x_1}{x_2},\qquad
A_x=\frac{x_1-x_2}{x_1+x_2}.
```

Dividing the numerator and denominator of $`A_x`$ by $`x_2`$ gives

```math
A_x=\frac{R_x-1}{R_x+1}.                           \tag{23}
```

Solving for the ratio gives the inverse

```math
R_x=\frac{1+A_x}{1-A_x}.                           \tag{24}
```

Thus (23) is a bijection from $`R_x\in(0,\infty)`$ to
$`A_x\in(-1,1)`$. It extends to zero/infinite ratios by including the endpoints
$`-1`$ and $`+1`$. In particular,

```math
(R_\rho,R_\mu,R_G)
\longleftrightarrow(A_\rho,A_\mu,A_G)
```

is an invertible componentwise coordinate transformation. The contrasts do
not add new dimensionless information; they give bounded, exchange-odd
coordinates for the same three kernel directions.

#### 8. Precise status of the preferred variables

The relation to the classical basis (14) makes the answer exact. From
$`R_G=Q_1/Q_2`$,

```math
\begin{aligned}
s&=P\frac{1+R_\rho}{1+R_\mu},\\
\Lambda&=\frac{(1+R_\rho)(Q_1+Q_2)}{(1+R_\mu)^2},\\
A_\rho&=\frac{R_\rho-1}{R_\rho+1},\qquad
A_\mu=\frac{R_\mu-1}{R_\mu+1},\\
A_G&=\frac{Q_1-Q_2}{Q_1+Q_2}.
\end{aligned}                                      \tag{25}
```

Conversely, from the five preferred coordinates,

```math
R_\rho=\frac{1+A_\rho}{1-A_\rho},\quad
R_\mu=\frac{1+A_\mu}{1-A_\mu},\quad
R_G=\frac{1+A_G}{1-A_G},                           \tag{26}
```

```math
P=s\frac{1+R_\mu}{1+R_\rho},\qquad
Q_\Sigma=Q_1+Q_2
=\Lambda\frac{(1+R_\mu)^2}{1+R_\rho},            \tag{27}
```

and

```math
Q_1=\frac{1+A_G}{2}Q_\Sigma,\qquad
Q_2=\frac{1-A_G}{2}Q_\Sigma.                      \tag{28}
```

Equations (25)–(28) prove invertibility on the positive-property interior.
Therefore:

- In the strict linear-algebraic sense, (21) is **not a basis of
  $`\ker D`$**, because its members are not monomials and have no single
  exponent vectors.
- In standard applied Buckingham terminology, it may be called a
  **transformed $`\Pi`$ basis**, because it is an invertible transformation of
  any monomial basis.
- Most precisely, it is a **convenient nonlinear coordinate system on the
  five-dimensional positive $`\Pi`$-space**.

These statements are compatible rather than contradictory; they refer to
different meanings of the word “basis.”

#### 9. Is one basis mathematically preferable?

There is no basis that is optimal under every proposed criterion. The criteria
select different structures.

| Choice | Exchange symmetry | Algebraic simplicity | Bounded contrast coordinates | Physical interpretation |
|---|---|---|---|---|
| direct RREF (11) | hidden | simple to compute, poor for the equation | no | weak |
| repeating variables (14) | privileges medium 2 | integer monomials | no | good with a reference medium |
| ratio/reference (17) | ratios invert | simple monomials | no | exposes three pair freedoms |
| geometric means (18) | exact and diagonal in log coordinates | compact but fractional powers | no | strong for positive two-medium symmetry |
| totals/contrasts (21) | exact; contrasts change sign | simplest for this dispersion relation | yes, for all three contrasts | strongest for this problem |

The conclusions behind the table are mathematical:

1. **Exchange symmetry.** The exchange map is a linear involution on exponent
   space. Its $`-1`$ eigenspace is the three-dimensional pair-difference space
   (9), while its $`+1`$ sector inside $`\ker D`$ has dimension two. The
   geometric basis respects this decomposition with
   $`(P_g,L_g)`$ even and $`(\log R_\rho,\log R_\mu,\log R_G)`$ odd. The
   total/contrast coordinates express the same decomposition without
   logarithms.

2. **Algebraic simplicity.** The dimensional dispersion relation contains
   $`\rho_1+\rho_2`$, and the canonical derivation naturally normalizes by
   $`\rho_T,\mu_T,G_T`$. Substitution therefore yields fractions
   $`(1\pm A_x)/2`$ and a single rate parameter $`\Lambda`$. Reference or
   geometric bases require repeated rational conversions to these sums. For
   this equation—not by dimensional analysis alone—the total coordinates are
   algebraically preferable.

3. **Boundedness.** A nonconstant positive monomial is unbounded on the full
   positive invariant space: along a suitable logarithmic kernel direction it
   is $`\exp(ct)`$, which approaches either zero or infinity. Hence bounded
   contrasts cannot be obtained by merely choosing another linear basis of
   $`\ker D`$; a nonlinear transformation such as (23) is necessary. The
   contrast variables achieve boundedness without losing invertibility.

4. **Physical interpretation.** The quotient decomposition gives exactly
   three material-partition coordinates plus two common-scale coordinates.
   Choosing the latter as $`s`$ and $`\Lambda=(r_e/r_v)^2`$ identifies the
   response and the elastic-to-viscous rate competition directly. This is more
   closely tied to the governing equation than the raw RREF groups.

Accordingly, the geometric-mean construction (18) is the mathematically
cleanest strict monomial basis when exchange symmetry is the sole priority.
The total/contrast variables (21) are the preferable coordinates when all four
criteria are considered together. This null-space analysis therefore supports
the existing recommendation; it does not reveal an inconsistency or a reason
to change the canonical nondimensional formulation.


<a id="dispersion"></a>

## Nondimensional dispersion relation

### Final explicit nondimensional dispersion relations

#### Purpose

This document expands the recommended viscous-clock nondimensional dispersion
relation entirely in the variables

```math
s,\quad A_\rho,\quad A_\mu,\quad A_G,\quad \Lambda
```

and then replaces $`\Lambda`$ by the bounded rate coordinate $`E_k`$. The physics
and the viscous-clock scaling are unchanged. No numerical root finding is
performed.

The total properties and nondimensional response are

```math
\rho_T=\rho_1+\rho_2,\qquad
\mu_T=\mu_1+\mu_2,\qquad
G_T=G_1+G_2,\qquad
s=\frac{\rho_T\gamma}{\mu_Tk^2}.
```

The material contrasts and elastic–viscous parameter are

```math
A_\rho=\frac{\rho_1-\rho_2}{\rho_T},\qquad
A_\mu=\frac{\mu_1-\mu_2}{\mu_T},\qquad
A_G=\frac{G_1-G_2}{G_T},\qquad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2}.
```

#### 1. Starting viscous-clock equation

The recommended viscous-clock equation derived in the existing specification
can be written without its former denominator shorthand as

```math
\begin{aligned}
0={}&s^2\Bigg[
\frac{1}{
m_1s+\Lambda g_1
+(m_2s+\Lambda g_2)
\sqrt{1+\dfrac{r_2s^2}{m_2s+\Lambda g_2}}
}
\\
&\qquad\qquad+
\frac{1}{
m_2s+\Lambda g_2
+(m_1s+\Lambda g_1)
\sqrt{1+\dfrac{r_1s^2}{m_1s+\Lambda g_1}}
}
\Bigg]+4,
\end{aligned}                                      \tag{1}
```

where, only for the purpose of identifying the starting equation,

```math
r_i=\frac{\rho_i}{\rho_T},\qquad
m_i=\frac{\mu_i}{\mu_T},\qquad
g_i=\frac{G_i}{G_T}.
```

Equation (1) is the cleared version of the direct viscous-clock equation. The
fully explicit final equations below contain none of these normalized-fraction
symbols.

#### 2. Substitution of the contrast reconstructions

Use

```math
r_1=\frac{1+A_\rho}{2},\qquad
r_2=\frac{1-A_\rho}{2},
```

```math
m_1=\frac{1+A_\mu}{2},\qquad
m_2=\frac{1-A_\mu}{2},
```

and

```math
g_1=\frac{1+A_G}{2},\qquad
g_2=\frac{1-A_G}{2}.
```

The two distinct material combinations occurring in (1) become

```math
m_1s+\Lambda g_1
=\frac{s(1+A_\mu)+\Lambda(1+A_G)}{2},
```

```math
m_2s+\Lambda g_2
=\frac{s(1-A_\mu)+\Lambda(1-A_G)}{2},             \tag{2}
```

and the radical arguments become

```math
1+\frac{r_2s^2}{m_2s+\Lambda g_2}
=1+\frac{(1-A_\rho)s^2}
{s(1-A_\mu)+\Lambda(1-A_G)},
```

```math
1+\frac{r_1s^2}{m_1s+\Lambda g_1}
=1+\frac{(1+A_\rho)s^2}
{s(1+A_\mu)+\Lambda(1+A_G)}.                     \tag{3}
```

Each complete reciprocal denominator in (1) therefore has an overall factor
$`1/2`$. Removing that factor from both reciprocals transforms (1) into

```math
\begin{aligned}
0={}&2s^2\Bigg[
\frac{1}{
s(1+A_\mu)+\Lambda(1+A_G)
+\left[s(1-A_\mu)+\Lambda(1-A_G)\right]
\sqrt{1+\dfrac{(1-A_\rho)s^2}
{s(1-A_\mu)+\Lambda(1-A_G)}}
}
\\
&\qquad+
\frac{1}{
s(1-A_\mu)+\Lambda(1-A_G)
+\left[s(1+A_\mu)+\Lambda(1+A_G)\right]
\sqrt{1+\dfrac{(1+A_\rho)s^2}
{s(1+A_\mu)+\Lambda(1+A_G)}}
}
\Bigg]+4.
\end{aligned}                                      \tag{4}
```

Dividing (4) by two gives the equally valid cleared normalization

```math
\begin{aligned}
0={}&s^2\Bigg[
\frac{1}{
s(1+A_\mu)+\Lambda(1+A_G)
+\left[s(1-A_\mu)+\Lambda(1-A_G)\right]
\sqrt{1+\dfrac{(1-A_\rho)s^2}
{s(1-A_\mu)+\Lambda(1-A_G)}}
}
\\
&\qquad+
\frac{1}{
s(1-A_\mu)+\Lambda(1-A_G)
+\left[s(1+A_\mu)+\Lambda(1+A_G)\right]
\sqrt{1+\dfrac{(1+A_\rho)s^2}
{s(1+A_\mu)+\Lambda(1+A_G)}}
}
\Bigg]+2.
\end{aligned}                                      \tag{5}
```

For $`s\ne0`$, factoring one $`s`$ out of every denominator in (5) gives the
direct-viscosity form used for the final equations. This step changes the
outer factor from $`s^2`$ to $`s`$ but changes no radical argument.

#### 3. Final explicit $`\Lambda`$ form

The completely expanded dispersion function in the algebraically canonical
parameter $`\Lambda`$ is

```math
\boxed{
\begin{aligned}
\mathcal D_\Lambda^*
(s;A_\rho,A_\mu,A_G,\Lambda)
={}&s\Bigg[
\frac{1}{
(1+A_\mu)+\dfrac{\Lambda}{s}(1+A_G)
+\left[(1-A_\mu)+\dfrac{\Lambda}{s}(1-A_G)\right]
\sqrt{
1+\dfrac{(1-A_\rho)s}
{(1-A_\mu)+\dfrac{\Lambda}{s}(1-A_G)}
}
}
\\
&\quad+
\frac{1}{
(1-A_\mu)+\dfrac{\Lambda}{s}(1-A_G)
+\left[(1+A_\mu)+\dfrac{\Lambda}{s}(1+A_G)\right]
\sqrt{
1+\dfrac{(1+A_\rho)s}
{(1+A_\mu)+\dfrac{\Lambda}{s}(1+A_G)}
}
}
\Bigg]+2=0.
\end{aligned}}
                                                               \tag{6}
```

Equation (6) contains only $`s`$, $`A_\rho`$, $`A_\mu`$, $`A_G`$, and $`\Lambda`$.
The first reciprocal contains the inertia of medium 2 in its radical, hence
the factors $`1-A_\rho`$, $`1-A_\mu`$, and $`1-A_G`$ there. The second reciprocal is
its exact $`1\leftrightarrow2`$ counterpart. This provides a direct sign check.

#### 4. From $`\Lambda`$ to $`E_k`$

The viscous and elastic response rates are

```math
r_v=\frac{\mu_Tk^2}{\rho_T},\qquad
r_e=k\sqrt{\frac{G_T}{\rho_T}},
```

so

```math
\frac{r_e}{r_v}
=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}
=\sqrt{\Lambda}.
```

By definition,

```math
E_k=\frac{r_v-r_e}{r_v+r_e}
=\frac{1-\sqrt{\Lambda}}{1+\sqrt{\Lambda}}.
```

Solving for the rate ratio gives

```math
\sqrt{\Lambda}=\frac{1-E_k}{1+E_k},\qquad
\Lambda=\left(\frac{1-E_k}{1+E_k}\right)^2.       \tag{7}
```

Thus the $`E_k`$ equation is obtained from (6) by making only the substitution
(7). The response variable $`s`$ remains the viscous-clock response; changing
from $`\Lambda`$ to $`E_k`$ does not change the clock.

#### 5. Final explicit $`E_k`$ form

Substitution of (7) into every occurrence of $`\Lambda`$ in (6) gives

```math
\boxed{
\begin{aligned}
\mathcal D_{E_k}^*
(s;A_\rho,A_\mu,A_G,E_k)
={}&s\Bigg[
\frac{1}{
(1+A_\mu)+\dfrac{(1-E_k)^2}{s(1+E_k)^2}(1+A_G)
+\left[(1-A_\mu)+\dfrac{(1-E_k)^2}{s(1+E_k)^2}(1-A_G)\right]
\sqrt{
1+\dfrac{(1-A_\rho)s}
{(1-A_\mu)+\dfrac{(1-E_k)^2}{s(1+E_k)^2}(1-A_G)}
}
}
\\
&\quad+
\frac{1}{
(1-A_\mu)+\dfrac{(1-E_k)^2}{s(1+E_k)^2}(1-A_G)
+\left[(1+A_\mu)+\dfrac{(1-E_k)^2}{s(1+E_k)^2}(1+A_G)\right]
\sqrt{
1+\dfrac{(1+A_\rho)s}
{(1+A_\mu)+\dfrac{(1-E_k)^2}{s(1+E_k)^2}(1+A_G)}
}
}
\Bigg]+2=0.
\end{aligned}}
                                                               \tag{8}
```

There is no sign or branch change between (6) and (8). In particular,

```math
\mathcal D_{E_k}^*(s;A_\rho,A_\mu,A_G,E_k)
=\mathcal D_\Lambda^*\!\left(
s;A_\rho,A_\mu,A_G,
\left(\frac{1-E_k}{1+E_k}\right)^2
\right).                                                        \tag{9}
```

#### 6. Constant-prefactor and normalization audit

The constant can be traced without relying on the cleared equation. Starting
from the dimensional relation, multiply by $`\rho_T/k^2`$. Under the viscous
clock

```math
s=\frac{\rho_T\gamma}{\mu_Tk^2},
```

the dimensional constant becomes exactly

```math
\frac{\rho_T}{k^2}\frac{4k^2}{\rho_T}=4.          \tag{10}
```

Before contrast substitution, the direct dimensionless equation consequently
has outer factor $`s`$ and constant $`4`$. Each normalized material contribution
contains a factor $`1/2`$ after the contrast reconstructions. Therefore each
complete reciprocal is multiplied by $`2`$, and the equation is

```math
2s[\text{two fully expanded reciprocal terms}]+4=0.             \tag{11}
```

Dividing the entire equation by the nonzero constant $`2`$ gives

```math
s[\text{the same two fully expanded reciprocal terms}]+2=0,     \tag{12}
```

which is the normalization used in the boxed equations (6) and (8). Hence
$`2s[\cdots]+4=0`$ and $`s[\cdots]+2=0`$ are exactly equivalent; neither represents
a change in physics. The cleared normalizations (4) and (5) are likewise
equivalent for $`s\ne0`$. The original Kelvin–Voigt representation already
contains $`G_i/\gamma`$, so $`s=0`$ must be understood through the appropriate
limit rather than by treating the uncleared formula as an ordinary value.

#### 7. Compact implementation form

Only after establishing the fully explicit equations, it is convenient for an
implementation to define

```math
\eta_1=(1+A_\mu)+\frac{\Lambda}{s}(1+A_G),\qquad
\eta_2=(1-A_\mu)+\frac{\Lambda}{s}(1-A_G),
```

```math
q_1=\sqrt{1+\frac{(1+A_\rho)s}{\eta_1}},\qquad
q_2=\sqrt{1+\frac{(1-A_\rho)s}{\eta_2}}.
```

Then the $`\Lambda`$ form can be evaluated as

```math
\mathcal D_\Lambda^*
=s\left(\frac{1}{\eta_1+\eta_2q_2}
+\frac{1}{\eta_2+\eta_1q_1}\right)+2.             \tag{13}
```

For the $`E_k`$ form, use the same implementation after assigning

```math
\Lambda=\left(\frac{1-E_k}{1+E_k}\right)^2.
```

The square-root branch convention must remain the same as in the dimensional
dispersion relation.

#### 8. Comparison of the two forms

The two equations contain the same four independent controls and the same
dimensionless response $`s`$.

- $`\Lambda=\rho_TG_T/(\mu_T^2k^2)`$ is algebraically canonical for the viscous
  clock. It appears directly when dimensional factors are collected, and it
  keeps the equation rational in the material parameters outside the square
  roots.
- $`E_k=(r_v-r_e)/(r_v+r_e)`$ is a bounded reparameterization of the positive
  rate ratio. For nonnegative material properties, $`-1<E_k\leq1`$, with
  $`E_k=1`$ at $`\Lambda=0`$ and $`E_k\to-1`$ as $`\Lambda\to\infty`$.
- Equations (6) and (8) are mathematically equivalent by the invertible map
  (7) on the finite positive-parameter interior. Choosing one or the other
  changes only the coordinate used for parameter studies, not the dispersion
  relation or its solutions.


<a id="symbolic-verification"></a>

## Symbolic verification framework

### Symbolic dispersion-relation verification

This directory contains the Mathematica/Wolfram Language symbolic audit
framework for the Kelvin--Voigt dispersion relation.

The current symbolic workflow is organized into three maintained phases:

1. **Phase I** verifies the dimensional-to-nondimensional transformation
   $`D(\gamma)\mapsto D^*_\Lambda(s)`$.
2. **Phase II** is the authoritative radical-elimination and polynomial
   derivation phase. It constructs $`P_{14}(s)`$ and verifies the forward
   no-root-loss chain.
3. **Phase III** is a thin polynomial-completeness audit built from the
   verified Phase II identities. It does not repeat the full derivation.

The Phase I target is

```math
D(\gamma)\longrightarrow D^*(s)
```

The Phase II/III target is the one-way theorem

```math
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0.
```

The reverse implication is not claimed.

#### Directory layout

- `DispersionRelationVerification.nb` is the primary Mathematica notebook
  interface. It is organized as a derivation with explanatory sections and
  auditable symbolic checks.
- `audit.wl` is the batch verification script for reproducible command-line
  Phase I auditing.
- `audit_phase2.wl` is the batch verification script for the Phase II
  polynomial-candidate derivation audit.
- `audit_phase3.wl` is the batch verification script for the Phase III
  polynomial-completeness/no-root-loss audit.
- `verification/Phase1Verification.wl` contains the reusable symbolic
  definitions used by both the notebook and the Phase I audit script.
- `verification/Phase2PolynomialDerivation.wl` contains the Phase II symbolic
  derivation from $`D^*_\Lambda(s)=0`$ to the polynomial candidate
  $`P_{14}(s)=0`$.
- `verification/Phase3PolynomialCompleteness.wl` contains the maintained
  Phase III completeness audit. It packages the final no-root-loss theorem
  using the authoritative Phase II checks.
- Legacy runtime-filtering infrastructure is retained in the upstream research
  archive and is not included here or loaded by `audit_phase3.wl`.
- `phase1_publication_review.md` records publication-readiness review notes for
  the Phase I framework.
- `phase2_derivation_notes.md` records the Phase II derivation sequence and
  admissible-domain logic.
- `phase3_derivation_notes.md` records the Phase III
  polynomial-completeness/no-root-loss audit and the later runtime filtering
  workflow.
- `THEORY.md#root-verification` defines the project-wide terminology and
  implementation contract for polynomial candidate roots, mathematically
  genuine roots, and physically admissible roots.
- `logs/` stores newly generated audit logs.
- `exports/canonical/` stores the latest frozen implementation inputs.
- older timestamped logs and exports remain in the upstream research archive;
  they are not required by this standalone solver release.

#### Phase I scope

Implemented sections:

1. Section 0 — Definitions
2. Section 1 — Original dimensional dispersion relation
3. Section 2 — Nondimensionalization

Not yet implemented:

1. Phase IV — Production polynomial-root filtering
2. Phase V — Parameter-region and branch-structure analysis

#### Verification target

The framework verifies the exact normalization used by the numerical module:

```math
D^*_\Lambda(s;A_\rho,A_\mu,A_G,\Lambda)
=
\frac{\rho_T}{2k^2}
D(\gamma)
```

after applying the substitutions

```math
\gamma=\frac{s\mu_Tk^2}{\rho_T},
```

```math
\rho_1=\frac{\rho_T(1+A_\rho)}{2},\qquad
\rho_2=\frac{\rho_T(1-A_\rho)}{2},
```

```math
\mu_1=\frac{\mu_T(1+A_\mu)}{2},\qquad
\mu_2=\frac{\mu_T(1-A_\mu)}{2},
```

```math
G_1=\frac{\Lambda\mu_T^2k^2(1+A_G)}{2\rho_T},\qquad
G_2=\frac{\Lambda\mu_T^2k^2(1-A_G)}{2\rho_T}.
```

Equivalently,

```math
D^*_\Lambda-\frac{\rho_T}{2k^2}D=0
```

under the Phase I assumptions.

The present symbolic verification is established under the explicit theorem
assumption $`\Lambda>0`$. The limiting case $`\Lambda=0`$ corresponds, within this
chosen parameterization, to both shear moduli vanishing simultaneously
($`E_k=1`$). That boundary case lies outside the scope of the current Phase I
theorem and may be analyzed separately if needed. Individual elastic moduli may
still vanish within the theorem through the boundary values $`A_G=\pm1`$ when
$`\Lambda>0`$.

The audit also verifies the nondimensionalization component by component:

1. Kelvin--Voigt effective-viscosity scaling
2. Radicand mapping
3. Principal square-root mapping
4. Reciprocal-denominator scaling
5. Constant-term scaling
6. Reciprocal-term scaling
7. Whole-expression equivalence

#### Running the audit

From the repository root, run the frozen Phase I audit with:

```bash
wolframscript -script symbolic/audit.wl
```

Run the Phase II polynomial-candidate audit with:

```bash
wolframscript -script symbolic/audit_phase2.wl
```

Run the Phase III polynomial-completeness audit with:

```bash
wolframscript -script symbolic/audit_phase3.wl
```

The script writes timestamped logs to `symbolic/logs/` and exported symbolic
expressions to `symbolic/exports/`.

If `wolframscript` is not installed, open
`symbolic/DispersionRelationVerification.nb` in Mathematica and evaluate the
notebook cells manually.

From a fresh Mathematica notebook, the same audit can be run as:

```wolfram
projectRoot = "/absolute/path/to/RMI-dispersion-relation";
Get[FileNameJoin[{projectRoot, "symbolic", "audit.wl"}]]
```

For Phase II:

```wolfram
projectRoot = "/absolute/path/to/RMI-dispersion-relation";
Get[FileNameJoin[{projectRoot, "symbolic", "audit_phase2.wl"}]]
```

For Phase III:

```wolfram
projectRoot = "/absolute/path/to/RMI-dispersion-relation";
Get[FileNameJoin[{projectRoot, "symbolic", "audit_phase3.wl"}]]
```

When loaded through `Get[...]` in a notebook, `audit.wl` returns an Association
with fields such as `"EquivalenceVerified"`, `"PassFailSummary"`, and the paths
of the generated log/export files. It does not call `Exit[]` in a front-end
notebook session.

#### Phase II theorem scope

Phase II establishes only the forward implication

```math
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0
```

on the recorded admissible domain. It does not claim the reverse implication.
Spurious roots may enter when reciprocal denominators are cleared and when the
two radical-elimination squarings are performed.

#### Phase III theorem scope

Phase III packages the final polynomial-completeness theorem:

```math
\Omega_{\mathrm{definition}}
\land
\Omega_{\mathrm{no\ root\ loss}}
\land
D^*_\Lambda(s)=0
\Longrightarrow
P_{14}(s)=0.
```

Phase II is authoritative for the symbolic identities and algebraic
transformations; Phase III reuses those results. It does not prove or claim
the reverse implication.

#### Runtime candidate filtering

Polynomial roots are only candidates. The full mathematical contract for the
later runtime workflow is documented in
[root-verification reference](THEORY.md#root-verification).

In summary, the workflow is:

1. solve $`P_{14}(s)=0`$;
2. reject candidates outside the mathematical definition domain;
3. substitute each remaining candidate into $`D^*_\Lambda(s)`$;
4. retain only true roots of the original dispersion relation;
5. compute $`q_-`$ and $`q_+`$ on the principal square-root branch; and
6. apply the spatial decay conditions

   ```math
   \mathrm{Re}(q_-)>0,
   \qquad
   \mathrm{Re}(q_+)>0.
   ```

This runtime filtering workflow is separate from the symbolic no-root-loss
theorem.

Temporal behavior is classified by $`\mathrm{Re}(s)`$. Physical spatial
admissibility is classified by $`\mathrm{Re}(q_\pm)`$. These are different
concepts and must not be conflated.

#### Mathematica GUI validation checklist

If `wolframscript` is unavailable, validate from Mathematica GUI:

1. Quit the kernel.
2. Open `DispersionRelationVerification.nb`.
3. Evaluate Phase I and confirm every component check reports PASS.
4. Evaluate Phase II and confirm:
   - `P14IsPolynomialInS -> True`;
   - degree in `s` is 14 generically;
   - coefficient count is 15;
   - leading coefficient is `(Arho + Amu)^2`;
   - every `ForwardPreservingChecks` entry has `"Status" -> "PASS"`;
   - no `"UNRESOLVED"` check is treated as PASS.
5. Evaluate Phase III and confirm the polynomial-completeness audit reports
   PASS.
6. Regenerate logs and exports.
7. Quit the kernel and rerun once to check for hidden state dependence.

#### Branch convention

The symbolic expressions use Mathematica's principal square-root branch,
matching the NumPy principal branch used by `src/dispersion.py`. Phase I checks
verify algebraic equivalence of the dimensional and nondimensional expressions
on the same branch convention. Later admissible-domain work will make branch
selection and possible spurious roots explicit.


<a id="phase1-review"></a>

## Phase I review record

### Phase I symbolic framework review

This note reviews the Phase I Mathematica notebook as a future publication
artifact. It does not begin Phase II and does not derive the polynomial
representation.

#### Overall assessment

The Phase I framework now verifies the nondimensionalization

```math
D(\gamma)\longrightarrow D^*_\Lambda(s)
```

in a reproducible way and records the nonzero scaling factor

```math
D^*_\Lambda=\frac{\rho_T}{2k^2}D(\gamma).
```

The strongest feature is that the notebook and the batch audit share the same
symbolic definitions through `verification/Phase1Verification.wl`, reducing the
risk of notebook/script drift.

#### Issues found in the first exported PDF

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

#### Improvements made

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

#### Recommendations before Phase II

#### Theorem-scope clarification

The Phase I symbolic theorem is established under the explicit assumption
$`\Lambda>0`$. The limiting case $`\Lambda=0`$ (equivalently $`E_k=1`$) corresponds,
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


<a id="phase2-derivation"></a>

## Phase II derivation record

### Phase II derivation notes

Phase II starts from the frozen Phase I theorem:

```math
D^*_\Lambda(s)=\frac{\rho_T}{2k^2}D(\gamma).
```

The goal is not to prove equivalence between the polynomial and the original
dispersion relation. The Phase II theorem is only the one-way implication

```math
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0
```

on the admissible domain.

#### Compact notation

Define

```math
E_+=s(1+A_\mu)+\Lambda(1+A_G),
```

```math
E_-=s(1-A_\mu)+\Lambda(1-A_G).
```

Then

```math
\eta_+=\frac{E_+}{s},\qquad
\eta_-=\frac{E_-}{s}.
```

The two radicands are

```math
R_- = 1+\frac{(1-A_\rho)s^2}{E_-},
\qquad
R_+ = 1+\frac{(1+A_\rho)s^2}{E_+},
```

with principal square roots

```math
q_-=\sqrt{R_-},\qquad q_+=\sqrt{R_+}.
```

The reciprocal denominators in the $`s^2`$ clock form are

```math
B_- = E_+ + E_-q_-,
\qquad
B_+ = E_- + E_+q_+.
```

Thus

```math
F_0
=s^2\left(\frac1{B_-}+\frac1{B_+}\right)+2.
```

#### Admissible-domain records

The original nondimensional expression requires

```math
s\ne0,\qquad E_+\ne0,\qquad E_-\ne0,
```

and

```math
B_-\ne0,\qquad B_+\ne0.
```

The physical Phase II assumptions are inherited from Phase I:

```math
k>0,\quad \rho_T>0,\quad \mu_T>0,\quad \Lambda>0,
```

```math
-1<A_\rho<1,\qquad -1<A_\mu<1,\qquad -1\le A_G\le1,
```

```math
s\in\mathbb C,\qquad s\ne0.
```

#### Step chain

Clearing the reciprocal denominators gives

```math
F_1=s^2(B_-+B_+)+2B_-B_+=0.
```

This multiplication is valid in the forward direction on the admissible domain.
It may introduce spurious roots if used backward.

Collecting radicals gives

```math
F_1=\alpha+\beta q_-+\gamma q_+ +\delta q_-q_+=0,
```

where

```math
\alpha=s^2(E_++E_-)+2E_+E_-,
```

```math
\beta=E_-(s^2+2E_-),
\qquad
\gamma=E_+(s^2+2E_+),
```

```math
\delta=2E_+E_-.
```

Isolate $`q_+`$ and square:

```math
F_2=(\alpha+\beta q_-)^2
-R_+(\gamma+\delta q_-)^2=0.
```

This preserves forward implication but may introduce spurious roots.

Write

```math
F_2=u_0+u_1q_-+u_2q_-^2.
```

Using $`q_-^2=R_-`$ gives

```math
F_3=u_0+u_1q_-+u_2R_-=0.
```

Then isolate $`q_-`$ and square:

```math
F_4=(u_0+u_2R_-)^2-u_1^2R_-=0.
```

This again preserves forward implication but may introduce spurious roots.

Finally clear rational denominators. The factor

```math
E_+^2E_-^2
```

is a denominator-clearing/admissibility factor and is not divided out of the
final polynomial. Only the harmless nonzero constant factor $`4`$ is removed from
the cleared numerator. The resulting polynomial is the candidate

```math
P_{14}(s;A_\rho,A_\mu,A_G,\Lambda)=0.
```

#### Generic degree

The Phase II symbolic audit is designed to confirm:

```math
\deg_s P_{14}=14,
```

with leading coefficient

```math
c_{14}=(A_\rho+A_\mu)^2.
```

Therefore the degree may drop on the parameter hypersurface

```math
A_\rho+A_\mu=0.
```

Additional lower-degree cancellations on special parameter strata should be
treated separately in a later admissible-domain analysis.

#### Stop boundary

Phase II does not claim

```math
P_{14}(s)=0\Longrightarrow D^*_\Lambda(s)=0.
```

It also does not implement a polynomial-root filtering algorithm. Those tasks
belong to later phases.


<a id="phase3-derivation"></a>

## Phase III derivation record

### Phase III polynomial-completeness audit

Phase III starts after Phase II is frozen. Its purpose is to package and audit
the final polynomial-completeness theorem: the degree-14 polynomial candidate
equation does not lose any root of the original nondimensional dispersion
relation.

Phase II is the authoritative source for the radical elimination,
denominator-clearing steps, coefficient verification, and forward-preserving
symbolic identities. Phase III reuses those verified Phase II checks rather
than repeating the full derivation under a second name.

The central theorem is the one-way implication

```math
\Omega_{\mathrm{definition}}
\land
D^*_\Lambda(s)=0
\Longrightarrow
P_{14}(s)=0.
```

Equivalently, every root of the original dispersion relation is represented
among the roots of $`P_{14}`$ under the stated definition-domain and derivation
assumptions.

Phase III does **not** claim the reverse implication

```math
P_{14}(s)=0
\Longrightarrow
D^*_\Lambda(s)=0.
```

Extra polynomial roots are expected because the derivation clears denominators
and squares equations. Those candidates are filtered later by direct
substitution into the original dispersion relation.

#### Definition domain

The original nondimensional expression requires:

```math
s\ne0,
```

```math
E_+\ne0,\qquad E_-\ne0,
```

and the two reciprocal denominators to be nonzero.

These conditions define $`\Omega_{\mathrm{definition}}`$, the domain on which
$`D^*_\Lambda(s)`$ is mathematically well-defined.

#### No-root-loss conditions

The no-root-loss audit records only the assumptions required to preserve the
forward implication from the original equation to the polynomial equation.

Important distinctions:

- multiplying by an expression cannot lose a root, although it may add
  extraneous roots;
- squaring cannot lose a root satisfying the pre-squared equation, although it
  may add extraneous roots;
- taking the numerator of a rational expression requires the rational
  denominator to be nonzero at the root;
- dividing by the nonzero constant normalization factor $`4`$ is harmless.

Branch or reverse-equivalence conditions are not part of the main Phase III
theorem.

#### Forward-step audit

The audited Phase II chain is

```math
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
```

Each step is classified in Phase II by whether it preserves the forward
implication and whether it can introduce extraneous polynomial candidates.
Phase III imports those classifications into the final theorem summary.

#### Runtime root-selection workflow

The later numerical workflow is:

1. solve $`P_{14}(s)=0`$;
2. obtain polynomial candidate roots;
3. reject candidates outside $`\Omega_{\mathrm{definition}}`$;
4. substitute each remaining candidate into $`D^*_\Lambda(s)`$;
5. retain candidates satisfying the original dispersion relation within
   tolerance;
6. apply physical admissibility criteria later.

Physical admissibility is intentionally outside this polynomial-completeness
proof.

#### Stop boundary

Phase III does not implement the production polynomial-root filtering
algorithm, does not construct physical spectra, and does not claim reverse
equivalence.


<a id="canonical-exports"></a>

## Canonical exports

### Canonical symbolic exports

This directory contains the latest frozen symbolic exports used as
implementation inputs for the production Python workflow.

The files here are normalized, stable copies of the latest verified timestamped
exports:

- `P14_coefficients.wl`
- `P14_polynomial.wl`
- `P14_polynomial.tex`
- `NoRootLossConditions.wl`
- `Phase2ForwardCompletenessChecks.wl`
- `DStarLambda_phase1.tex`
- `ScaledDimensionalD_phase1.tex`
- `phase2_report.md`
- `phase3_forward_numerical_validation.wl`
- `phase3_polynomial_completeness_report.md`

Older timestamped symbolic artifacts remain in the upstream research archive
and are not included in this standalone release.


<a id="phase2-report"></a>

## Frozen Phase II report

### Phase II polynomial derivation audit

The established theorem is one-way:

```math
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0.
```

The reverse implication is not claimed.

- PASS: True
- Forward theorem PASS: True
- Numerical spot-check PASS: True
- Degree in $`s`$: 14
- Coefficient count: 15
- Leading coefficient: `(Amu + Arho)^2`
- Polynomial normalization factor: `4`
- Denominator-clearing eta factor retained in admissibility conditions: `((1 - AG)*Lambda + (1 - Amu)*s)^2*((1 + AG)*Lambda + (1 + Amu)*s)^2`

#### Forward-preserving checks

```wolfram
{<|"Name" -> "DStarLambda equals F0ClockForm", "Status" -> "PASS", "Statement" -> HoldForm[2 + s*((1 + Amu + ((1 + AG)*Lambda)/s + (1 - Amu + ((1 - AG)*Lambda)/s)*Sqrt[1 + ((1 - Arho)*s)/(1 - Amu + ((1 - AG)*Lambda)/s)])^(-1) + (1 - Amu + ((1 - AG)*Lambda)/s + (1 + Amu + ((1 + AG)*Lambda)/s)*Sqrt[1 + ((1 + Arho)*s)/(1 + Amu + ((1 + AG)*Lambda)/s)])^(-1)) == 2 + s^2*(((1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])^(-1) + ((1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)])^(-1))], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Exact rewriting of the nondimensional dispersion relation into the s^2 clock form."|>, <|"Name" -> "First denominator-clearing step", "Status" -> "PASS", "Statement" -> HoldForm[((1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])*((1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)])*(2 + s^2*(((1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])^(-1) + ((1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)])^(-1))) == 2*Lambda^2 - 2*AG^2*Lambda^2 + 4*Lambda*s - 4*AG*Amu*Lambda*s + 2*s^2 - 2*Amu^2*s^2 + 2*Lambda*s^2 + 2*s^3 + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - AG*Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - Amu*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + AG*Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Amu*s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)]], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Multiplication by reciprocal denominators. This preserves roots in the forward direction on the original-expression domain."|>, <|"Name" -> "Collected F1 radical form", "Status" -> "PASS", "Statement" -> HoldForm[2*Lambda^2 - 2*AG^2*Lambda^2 + 4*Lambda*s - 4*AG*Amu*Lambda*s + 2*s^2 - 2*Amu^2*s^2 + 2*Lambda*s^2 + 2*s^3 + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - AG*Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - Amu*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + AG*Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Amu*s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] == 2*((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s) + s^2*((1 - AG)*Lambda + (1 + AG)*Lambda + (1 - Amu)*s + (1 + Amu)*s) + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*(s^2 + 2*((1 - AG)*Lambda + (1 - Amu)*s)) + 2*((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)]*(s^2 + 2*((1 + AG)*Lambda + (1 + Amu)*s))], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Exact collection into alpha0 + betaMinus qMinus + gammaPlus qPlus + delta qMinus qPlus."|>, <|"Name" -> "First radical elimination preserves forward implication", "Status" -> "PASS", "Statement" -> HoldForm[Implies[F1LinearRadicals[] == 0, F2AfterFirstSquaring[] == 0]], "Assumptions" -> True, "Details" -> "The step squares an equation satisfied by any original root. Squaring may add roots but cannot lose roots satisfying the pre-squared equation."|>, <|"Name" -> "F2 single-radical rewrite", "Status" -> "PASS", "Statement" -> HoldForm[-32*AG*Lambda^4 - 32*AG^3*Lambda^4 - 96*AG*Lambda^3*s - 32*AG^3*Lambda^3*s - 32*Amu*Lambda^3*s - 96*AG^2*Amu*Lambda^3*s - 96*AG*Lambda^2*s^2 - 96*Amu*Lambda^2*s^2 - 96*AG^2*Amu*Lambda^2*s^2 - 96*AG*Amu^2*Lambda^2*s^2 + 8*Lambda^3*s^2 - 24*AG*Lambda^3*s^2 - 8*AG^2*Lambda^3*s^2 - 8*AG^3*Lambda^3*s^2 - 32*AG*Lambda*s^3 - 96*Amu*Lambda*s^3 - 96*AG*Amu^2*Lambda*s^3 - 32*Amu^3*Lambda*s^3 + 24*Lambda^2*s^3 - 48*AG*Lambda^2*s^3 - 8*AG^2*Lambda^2*s^3 - 24*Amu*Lambda^2*s^3 - 16*AG*Amu*Lambda^2*s^3 - 24*AG^2*Amu*Lambda^2*s^3 - 32*Amu*s^4 - 32*Amu^3*s^4 + 24*Lambda*s^4 - 24*AG*Lambda*s^4 - 48*Amu*Lambda*s^4 - 16*AG*Amu*Lambda*s^4 - 8*Amu^2*Lambda*s^4 - 24*AG*Amu^2*Lambda*s^4 + 4*Lambda^2*s^4 - 4*AG*Lambda^2*s^4 + 8*s^5 - 24*Amu*s^5 - 8*Amu^2*s^5 - 8*Amu^3*s^5 + 8*Lambda*s^5 - 4*AG*Lambda*s^5 - 4*Amu*Lambda*s^5 + 4*s^6 - 4*Amu*s^6 - (16*AG*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG^2*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG^2*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (128*AG*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (128*AG*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (AG^2*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Amu^3*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Amu^3*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + s^8/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (Amu^2*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Amu^2*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (8*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - s^8/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*AG^4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*AG^4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^3*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^2*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^3*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (32*AG*Amu*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*AG^2*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (32*AG*Amu*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*AG^2*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Amu^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu^3*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Amu^2*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu^3*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*Amu^4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Amu^4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - 32*AG*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^3*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*AG*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*AG*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG^2*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^3*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - (8*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) == -32*AG*Lambda^4 - 32*AG^3*Lambda^4 - 96*AG*Lambda^3*s - 32*AG^3*Lambda^3*s - 32*Amu*Lambda^3*s - 96*AG^2*Amu*Lambda^3*s - 96*AG*Lambda^2*s^2 - 96*Amu*Lambda^2*s^2 - 96*AG^2*Amu*Lambda^2*s^2 - 96*AG*Amu^2*Lambda^2*s^2 + 8*Lambda^3*s^2 - 24*AG*Lambda^3*s^2 - 8*AG^2*Lambda^3*s^2 - 8*AG^3*Lambda^3*s^2 - 32*AG*Lambda*s^3 - 96*Amu*Lambda*s^3 - 96*AG*Amu^2*Lambda*s^3 - 32*Amu^3*Lambda*s^3 + 24*Lambda^2*s^3 - 48*AG*Lambda^2*s^3 - 8*AG^2*Lambda^2*s^3 - 24*Amu*Lambda^2*s^3 - 16*AG*Amu*Lambda^2*s^3 - 24*AG^2*Amu*Lambda^2*s^3 - 32*Amu*s^4 - 32*Amu^3*s^4 + 24*Lambda*s^4 - 24*AG*Lambda*s^4 - 48*Amu*Lambda*s^4 - 16*AG*Amu*Lambda*s^4 - 8*Amu^2*Lambda*s^4 - 24*AG*Amu^2*Lambda*s^4 + 4*Lambda^2*s^4 - 4*AG*Lambda^2*s^4 + 8*s^5 - 24*Amu*s^5 - 8*Amu^2*s^5 - 8*Amu^3*s^5 + 8*Lambda*s^5 - 4*AG*Lambda*s^5 - 4*Amu*Lambda*s^5 + 4*s^6 - 4*Amu*s^6 - (16*AG*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG^2*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG^2*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (128*AG*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (128*AG*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (AG^2*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Amu^3*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Amu^3*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + s^8/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (Amu^2*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Amu^2*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (8*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - s^8/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*AG^4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*AG^4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^3*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^2*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^3*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (32*AG*Amu*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*AG^2*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (32*AG*Amu*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*AG^2*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Amu^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu^3*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Amu^2*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu^3*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*Amu^4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Amu^4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - 32*AG*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^3*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*AG*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*AG*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG^2*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^3*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - (8*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s)], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Uses qMinus^2 = radicandMinus with Mathematica's principal Sqrt convention; no PowerExpand or branch replacement is used."|>, <|"Name" -> "Second radical elimination preserves forward implication", "Status" -> "PASS", "Statement" -> HoldForm[Implies[F3SingleRadical[] == 0, F4AfterSecondSquaring[] == 0]], "Assumptions" -> True, "Details" -> "The second squaring may add roots but cannot lose roots satisfying F3 == 0."|>, <|"Name" -> "Rational numerator extraction preserves forward implication", "Status" -> "PASS", "Statement" -> HoldForm[Implies[F4AfterSecondSquaring[] == 0 && Denominator[Together[F4AfterSecondSquaring[]]] != 0, rawPolynomialNumerator[] == 0]], "Assumptions" -> HoldForm[Denominator[Together[F4AfterSecondSquaring[]]] != 0], "Details" -> "If a defined rational expression is zero and its denominator is nonzero, then its numerator is zero."|>, <|"Name" -> "Final expression equals P14 up to constant normalization", "Status" -> "PASS", "Statement" -> HoldForm[4*(1024*AG^4*Lambda^8 + 4096*AG^4*Lambda^7*s + 4096*AG^3*Amu*Lambda^7*s + 6144*AG^4*Lambda^6*s^2 + 16384*AG^3*Amu*Lambda^6*s^2 + 6144*AG^2*Amu^2*Lambda^6*s^2 + 2048*AG^4*Lambda^7*s^2 + 1024*AG^3*Arho*Lambda^7*s^2 + 4096*AG^4*Lambda^5*s^3 + 24576*AG^3*Amu*Lambda^5*s^3 + 24576*AG^2*Amu^2*Lambda^5*s^3 + 4096*AG*Amu^3*Lambda^5*s^3 + 6144*AG^4*Lambda^6*s^3 + 8192*AG^3*Amu*Lambda^6*s^3 + 4096*AG^3*Arho*Lambda^6*s^3 + 3072*AG^2*Amu*Arho*Lambda^6*s^3 + 1024*AG^4*Lambda^4*s^4 + 16384*AG^3*Amu*Lambda^4*s^4 + 36864*AG^2*Amu^2*Lambda^4*s^4 + 16384*AG*Amu^3*Lambda^4*s^4 + 1024*Amu^4*Lambda^4*s^4 + 6144*AG^4*Lambda^5*s^4 + 24576*AG^3*Amu*Lambda^5*s^4 + 12288*AG^2*Amu^2*Lambda^5*s^4 + 6144*AG^3*Arho*Lambda^5*s^4 + 12288*AG^2*Amu*Arho*Lambda^5*s^4 + 3072*AG*Amu^2*Arho*Lambda^5*s^4 + 64*AG^2*Lambda^6*s^4 + 1472*AG^4*Lambda^6*s^4 + 1536*AG^3*Arho*Lambda^6*s^4 + 384*AG^2*Arho^2*Lambda^6*s^4 - 128*AG^4*Arho^2*Lambda^6*s^4 + 4096*AG^3*Amu*Lambda^3*s^5 + 24576*AG^2*Amu^2*Lambda^3*s^5 + 24576*AG*Amu^3*Lambda^3*s^5 + 4096*Amu^4*Lambda^3*s^5 + 2048*AG^4*Lambda^4*s^5 + 24576*AG^3*Amu*Lambda^4*s^5 + 36864*AG^2*Amu^2*Lambda^4*s^5 + 8192*AG*Amu^3*Lambda^4*s^5 + 4096*AG^3*Arho*Lambda^4*s^5 + 18432*AG^2*Amu*Arho*Lambda^4*s^5 + 12288*AG*Amu^2*Arho*Lambda^4*s^5 + 1024*Amu^3*Arho*Lambda^4*s^5 + 256*AG^2*Lambda^5*s^5 + 2944*AG^4*Lambda^5*s^5 + 128*AG*Amu*Lambda^5*s^5 + 5888*AG^3*Amu*Lambda^5*s^5 + 4608*AG^3*Arho*Lambda^5*s^5 + 4608*AG^2*Amu*Arho*Lambda^5*s^5 + 1536*AG^2*Arho^2*Lambda^5*s^5 - 256*AG^4*Arho^2*Lambda^5*s^5 + 768*AG*Amu*Arho^2*Lambda^5*s^5 - 512*AG^3*Amu*Arho^2*Lambda^5*s^5 + 6144*AG^2*Amu^2*Lambda^2*s^6 + 16384*AG*Amu^3*Lambda^2*s^6 + 6144*Amu^4*Lambda^2*s^6 + 8192*AG^3*Amu*Lambda^3*s^6 + 36864*AG^2*Amu^2*Lambda^3*s^6 + 24576*AG*Amu^3*Lambda^3*s^6 + 2048*Amu^4*Lambda^3*s^6 + 1024*AG^3*Arho*Lambda^3*s^6 + 12288*AG^2*Amu*Arho*Lambda^3*s^6 + 18432*AG*Amu^2*Arho*Lambda^3*s^6 + 4096*Amu^3*Arho*Lambda^3*s^6 + 384*AG^2*Lambda^4*s^6 + 1472*AG^4*Lambda^4*s^6 + 512*AG*Amu*Lambda^4*s^6 + 11776*AG^3*Amu*Lambda^4*s^6 + 64*Amu^2*Lambda^4*s^6 + 8832*AG^2*Amu^2*Lambda^4*s^6 + 4608*AG^3*Arho*Lambda^4*s^6 + 13824*AG^2*Amu*Arho*Lambda^4*s^6 + 4608*AG*Amu^2*Arho*Lambda^4*s^6 + 2304*AG^2*Arho^2*Lambda^4*s^6 - 128*AG^4*Arho^2*Lambda^4*s^6 + 3072*AG*Amu*Arho^2*Lambda^4*s^6 - 1024*AG^3*Amu*Arho^2*Lambda^4*s^6 + 384*Amu^2*Arho^2*Lambda^4*s^6 - 768*AG^2*Amu^2*Arho^2*Lambda^4*s^6 + 128*AG^2*Lambda^5*s^6 + 448*AG^4*Lambda^5*s^6 + 32*AG*Arho*Lambda^5*s^6 + 800*AG^3*Arho*Lambda^5*s^6 + 384*AG^2*Arho^2*Lambda^5*s^6 - 128*AG^4*Arho^2*Lambda^5*s^6 + 64*AG*Arho^3*Lambda^5*s^6 - 64*AG^3*Arho^3*Lambda^5*s^6 + 4096*AG*Amu^3*Lambda*s^7 + 4096*Amu^4*Lambda*s^7 + 12288*AG^2*Amu^2*Lambda^2*s^7 + 24576*AG*Amu^3*Lambda^2*s^7 + 6144*Amu^4*Lambda^2*s^7 + 3072*AG^2*Amu*Arho*Lambda^2*s^7 + 12288*AG*Amu^2*Arho*Lambda^2*s^7 + 6144*Amu^3*Arho*Lambda^2*s^7 + 256*AG^2*Lambda^3*s^7 + 768*AG*Amu*Lambda^3*s^7 + 5888*AG^3*Amu*Lambda^3*s^7 + 256*Amu^2*Lambda^3*s^7 + 17664*AG^2*Amu^2*Lambda^3*s^7 + 5888*AG*Amu^3*Lambda^3*s^7 + 1536*AG^3*Arho*Lambda^3*s^7 + 13824*AG^2*Amu*Arho*Lambda^3*s^7 + 13824*AG*Amu^2*Arho*Lambda^3*s^7 + 1536*Amu^3*Arho*Lambda^3*s^7 + 1536*AG^2*Arho^2*Lambda^3*s^7 + 4608*AG*Amu*Arho^2*Lambda^3*s^7 - 512*AG^3*Amu*Arho^2*Lambda^3*s^7 + 1536*Amu^2*Arho^2*Lambda^3*s^7 - 1536*AG^2*Amu^2*Arho^2*Lambda^3*s^7 - 512*AG*Amu^3*Arho^2*Lambda^3*s^7 + 384*AG^2*Lambda^4*s^7 + 448*AG^4*Lambda^4*s^7 + 256*AG*Amu*Lambda^4*s^7 + 1792*AG^3*Amu*Lambda^4*s^7 + 128*AG*Arho*Lambda^4*s^7 + 1600*AG^3*Arho*Lambda^4*s^7 + 32*Amu*Arho*Lambda^4*s^7 + 2400*AG^2*Amu*Arho*Lambda^4*s^7 + 1152*AG^2*Arho^2*Lambda^4*s^7 - 128*AG^4*Arho^2*Lambda^4*s^7 + 768*AG*Amu*Arho^2*Lambda^4*s^7 - 512*AG^3*Amu*Arho^2*Lambda^4*s^7 + 256*AG*Arho^3*Lambda^4*s^7 - 128*AG^3*Arho^3*Lambda^4*s^7 + 64*Amu*Arho^3*Lambda^4*s^7 - 192*AG^2*Amu*Arho^3*Lambda^4*s^7 + 1024*Amu^4*s^8 + 8192*AG*Amu^3*Lambda*s^8 + 6144*Amu^4*Lambda*s^8 + 3072*AG*Amu^2*Arho*Lambda*s^8 + 4096*Amu^3*Arho*Lambda*s^8 + 64*AG^2*Lambda^2*s^8 + 512*AG*Amu*Lambda^2*s^8 + 384*Amu^2*Lambda^2*s^8 + 8832*AG^2*Amu^2*Lambda^2*s^8 + 11776*AG*Amu^3*Lambda^2*s^8 + 1472*Amu^4*Lambda^2*s^8 + 4608*AG^2*Amu*Arho*Lambda^2*s^8 + 13824*AG*Amu^2*Arho*Lambda^2*s^8 + 4608*Amu^3*Arho*Lambda^2*s^8 + 384*AG^2*Arho^2*Lambda^2*s^8 + 3072*AG*Amu*Arho^2*Lambda^2*s^8 + 2304*Amu^2*Arho^2*Lambda^2*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^8 - 1024*AG*Amu^3*Arho^2*Lambda^2*s^8 - 128*Amu^4*Arho^2*Lambda^2*s^8 + 384*AG^2*Lambda^3*s^8 + 768*AG*Amu*Lambda^3*s^8 + 1792*AG^3*Amu*Lambda^3*s^8 + 128*Amu^2*Lambda^3*s^8 + 2688*AG^2*Amu^2*Lambda^3*s^8 + 192*AG*Arho*Lambda^3*s^8 + 800*AG^3*Arho*Lambda^3*s^8 + 128*Amu*Arho*Lambda^3*s^8 + 4800*AG^2*Amu*Arho*Lambda^3*s^8 + 2400*AG*Amu^2*Arho*Lambda^3*s^8 + 1152*AG^2*Arho^2*Lambda^3*s^8 + 2304*AG*Amu*Arho^2*Lambda^3*s^8 - 512*AG^3*Amu*Arho^2*Lambda^3*s^8 + 384*Amu^2*Arho^2*Lambda^3*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^3*s^8 + 384*AG*Arho^3*Lambda^3*s^8 - 64*AG^3*Arho^3*Lambda^3*s^8 + 256*Amu*Arho^3*Lambda^3*s^8 - 384*AG^2*Amu*Arho^3*Lambda^3*s^8 - 192*AG*Amu^2*Arho^3*Lambda^3*s^8 + 80*AG^2*Lambda^4*s^8 + 48*AG^4*Lambda^4*s^8 + 48*AG*Arho*Lambda^4*s^8 + 176*AG^3*Arho*Lambda^4*s^8 + 4*Arho^2*Lambda^4*s^8 + 120*AG^2*Arho^2*Lambda^4*s^8 - 28*AG^4*Arho^2*Lambda^4*s^8 + 32*AG*Arho^3*Lambda^4*s^8 - 32*AG^3*Arho^3*Lambda^4*s^8 + 4*Arho^4*Lambda^4*s^8 - 8*AG^2*Arho^4*Lambda^4*s^8 + 4*AG^4*Arho^4*Lambda^4*s^8 + 2048*Amu^4*s^9 + 1024*Amu^3*Arho*s^9 + 128*AG*Amu*Lambda*s^9 + 256*Amu^2*Lambda*s^9 + 5888*AG*Amu^3*Lambda*s^9 + 2944*Amu^4*Lambda*s^9 + 4608*AG*Amu^2*Arho*Lambda*s^9 + 4608*Amu^3*Arho*Lambda*s^9 + 768*AG*Amu*Arho^2*Lambda*s^9 + 1536*Amu^2*Arho^2*Lambda*s^9 - 512*AG*Amu^3*Arho^2*Lambda*s^9 - 256*Amu^4*Arho^2*Lambda*s^9 + 128*AG^2*Lambda^2*s^9 + 768*AG*Amu*Lambda^2*s^9 + 384*Amu^2*Lambda^2*s^9 + 2688*AG^2*Amu^2*Lambda^2*s^9 + 1792*AG*Amu^3*Lambda^2*s^9 + 128*AG*Arho*Lambda^2*s^9 + 192*Amu*Arho*Lambda^2*s^9 + 2400*AG^2*Amu*Arho*Lambda^2*s^9 + 4800*AG*Amu^2*Arho*Lambda^2*s^9 + 800*Amu^3*Arho*Lambda^2*s^9 + 384*AG^2*Arho^2*Lambda^2*s^9 + 2304*AG*Amu*Arho^2*Lambda^2*s^9 + 1152*Amu^2*Arho^2*Lambda^2*s^9 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^9 - 512*AG*Amu^3*Arho^2*Lambda^2*s^9 + 256*AG*Arho^3*Lambda^2*s^9 + 384*Amu*Arho^3*Lambda^2*s^9 - 192*AG^2*Amu*Arho^3*Lambda^2*s^9 - 384*AG*Amu^2*Arho^3*Lambda^2*s^9 - 64*Amu^3*Arho^3*Lambda^2*s^9 + 160*AG^2*Lambda^3*s^9 + 160*AG*Amu*Lambda^3*s^9 + 192*AG^3*Amu*Lambda^3*s^9 + 144*AG*Arho*Lambda^3*s^9 + 176*AG^3*Arho*Lambda^3*s^9 + 48*Amu*Arho*Lambda^3*s^9 + 528*AG^2*Amu*Arho*Lambda^3*s^9 + 16*Arho^2*Lambda^3*s^9 + 240*AG^2*Arho^2*Lambda^3*s^9 + 240*AG*Amu*Arho^2*Lambda^3*s^9 - 112*AG^3*Amu*Arho^2*Lambda^3*s^9 + 96*AG*Arho^3*Lambda^3*s^9 - 32*AG^3*Arho^3*Lambda^3*s^9 + 32*Amu*Arho^3*Lambda^3*s^9 - 96*AG^2*Amu*Arho^3*Lambda^3*s^9 + 16*Arho^4*Lambda^3*s^9 - 16*AG^2*Arho^4*Lambda^3*s^9 - 16*AG*Amu*Arho^4*Lambda^3*s^9 + 16*AG^3*Amu*Arho^4*Lambda^3*s^9 + 64*Amu^2*s^10 + 1472*Amu^4*s^10 + 1536*Amu^3*Arho*s^10 + 384*Amu^2*Arho^2*s^10 - 128*Amu^4*Arho^2*s^10 + 256*AG*Amu*Lambda*s^10 + 384*Amu^2*Lambda*s^10 + 1792*AG*Amu^3*Lambda*s^10 + 448*Amu^4*Lambda*s^10 + 32*AG*Arho*Lambda*s^10 + 128*Amu*Arho*Lambda*s^10 + 2400*AG*Amu^2*Arho*Lambda*s^10 + 1600*Amu^3*Arho*Lambda*s^10 + 768*AG*Amu*Arho^2*Lambda*s^10 + 1152*Amu^2*Arho^2*Lambda*s^10 - 512*AG*Amu^3*Arho^2*Lambda*s^10 - 128*Amu^4*Arho^2*Lambda*s^10 + 64*AG*Arho^3*Lambda*s^10 + 256*Amu*Arho^3*Lambda*s^10 - 192*AG*Amu^2*Arho^3*Lambda*s^10 - 128*Amu^3*Arho^3*Lambda*s^10 + 80*AG^2*Lambda^2*s^10 + 320*AG*Amu*Lambda^2*s^10 + 80*Amu^2*Lambda^2*s^10 + 288*AG^2*Amu^2*Lambda^2*s^10 + 144*AG*Arho*Lambda^2*s^10 + 144*Amu*Arho*Lambda^2*s^10 + 528*AG^2*Amu*Arho*Lambda^2*s^10 + 528*AG*Amu^2*Arho*Lambda^2*s^10 + 24*Arho^2*Lambda^2*s^10 + 120*AG^2*Arho^2*Lambda^2*s^10 + 480*AG*Amu*Arho^2*Lambda^2*s^10 + 120*Amu^2*Arho^2*Lambda^2*s^10 - 168*AG^2*Amu^2*Arho^2*Lambda^2*s^10 + 96*AG*Arho^3*Lambda^2*s^10 + 96*Amu*Arho^3*Lambda^2*s^10 - 96*AG^2*Amu*Arho^3*Lambda^2*s^10 - 96*AG*Amu^2*Arho^3*Lambda^2*s^10 + 24*Arho^4*Lambda^2*s^10 - 8*AG^2*Arho^4*Lambda^2*s^10 - 32*AG*Amu*Arho^4*Lambda^2*s^10 - 8*Amu^2*Arho^4*Lambda^2*s^10 + 24*AG^2*Amu^2*Arho^4*Lambda^2*s^10 + 16*AG^2*Lambda^3*s^10 + 20*AG*Arho*Lambda^3*s^10 + 12*AG^3*Arho*Lambda^3*s^10 + 4*Arho^2*Lambda^3*s^10 + 12*AG^2*Arho^2*Lambda^3*s^10 + 4*AG*Arho^3*Lambda^3*s^10 - 4*AG^3*Arho^3*Lambda^3*s^10 + 128*Amu^2*s^11 + 448*Amu^4*s^11 + 32*Amu*Arho*s^11 + 800*Amu^3*Arho*s^11 + 384*Amu^2*Arho^2*s^11 - 128*Amu^4*Arho^2*s^11 + 64*Amu*Arho^3*s^11 - 64*Amu^3*Arho^3*s^11 + 160*AG*Amu*Lambda*s^11 + 160*Amu^2*Lambda*s^11 + 192*AG*Amu^3*Lambda*s^11 + 48*AG*Arho*Lambda*s^11 + 144*Amu*Arho*Lambda*s^11 + 528*AG*Amu^2*Arho*Lambda*s^11 + 176*Amu^3*Arho*Lambda*s^11 + 16*Arho^2*Lambda*s^11 + 240*AG*Amu*Arho^2*Lambda*s^11 + 240*Amu^2*Arho^2*Lambda*s^11 - 112*AG*Amu^3*Arho^2*Lambda*s^11 + 32*AG*Arho^3*Lambda*s^11 + 96*Amu*Arho^3*Lambda*s^11 - 96*AG*Amu^2*Arho^3*Lambda*s^11 - 32*Amu^3*Arho^3*Lambda*s^11 + 16*Arho^4*Lambda*s^11 - 16*AG*Amu*Arho^4*Lambda*s^11 - 16*Amu^2*Arho^4*Lambda*s^11 + 16*AG*Amu^3*Arho^4*Lambda*s^11 + 16*AG^2*Lambda^2*s^11 + 32*AG*Amu*Lambda^2*s^11 + 40*AG*Arho*Lambda^2*s^11 + 20*Amu*Arho*Lambda^2*s^11 + 36*AG^2*Amu*Arho*Lambda^2*s^11 + 12*Arho^2*Lambda^2*s^11 + 12*AG^2*Arho^2*Lambda^2*s^11 + 24*AG*Amu*Arho^2*Lambda^2*s^11 + 8*AG*Arho^3*Lambda^2*s^11 + 4*Amu*Arho^3*Lambda^2*s^11 - 12*AG^2*Amu*Arho^3*Lambda^2*s^11 + 80*Amu^2*s^12 + 48*Amu^4*s^12 + 48*Amu*Arho*s^12 + 176*Amu^3*Arho*s^12 + 4*Arho^2*s^12 + 120*Amu^2*Arho^2*s^12 - 28*Amu^4*Arho^2*s^12 + 32*Amu*Arho^3*s^12 - 32*Amu^3*Arho^3*s^12 + 4*Arho^4*s^12 - 8*Amu^2*Arho^4*s^12 + 4*Amu^4*Arho^4*s^12 + 32*AG*Amu*Lambda*s^12 + 16*Amu^2*Lambda*s^12 + 20*AG*Arho*Lambda*s^12 + 40*Amu*Arho*Lambda*s^12 + 36*AG*Amu^2*Arho*Lambda*s^12 + 12*Arho^2*Lambda*s^12 + 24*AG*Amu*Arho^2*Lambda*s^12 + 12*Amu^2*Arho^2*Lambda*s^12 + 4*AG*Arho^3*Lambda*s^12 + 8*Amu*Arho^3*Lambda*s^12 - 12*AG*Amu^2*Arho^3*Lambda*s^12 + AG^2*Lambda^2*s^12 + 2*AG*Arho*Lambda^2*s^12 + Arho^2*Lambda^2*s^12 + 16*Amu^2*s^13 + 20*Amu*Arho*s^13 + 12*Amu^3*Arho*s^13 + 4*Arho^2*s^13 + 12*Amu^2*Arho^2*s^13 + 4*Amu*Arho^3*s^13 - 4*Amu^3*Arho^3*s^13 + 2*AG*Amu*Lambda*s^13 + 2*AG*Arho*Lambda*s^13 + 2*Amu*Arho*Lambda*s^13 + 2*Arho^2*Lambda*s^13 + Amu^2*s^14 + 2*Amu*Arho*s^14 + Arho^2*s^14) == 4*(1024*AG^4*Lambda^8 + 4096*AG^4*Lambda^7*s + 4096*AG^3*Amu*Lambda^7*s + 6144*AG^4*Lambda^6*s^2 + 16384*AG^3*Amu*Lambda^6*s^2 + 6144*AG^2*Amu^2*Lambda^6*s^2 + 2048*AG^4*Lambda^7*s^2 + 1024*AG^3*Arho*Lambda^7*s^2 + 4096*AG^4*Lambda^5*s^3 + 24576*AG^3*Amu*Lambda^5*s^3 + 24576*AG^2*Amu^2*Lambda^5*s^3 + 4096*AG*Amu^3*Lambda^5*s^3 + 6144*AG^4*Lambda^6*s^3 + 8192*AG^3*Amu*Lambda^6*s^3 + 4096*AG^3*Arho*Lambda^6*s^3 + 3072*AG^2*Amu*Arho*Lambda^6*s^3 + 1024*AG^4*Lambda^4*s^4 + 16384*AG^3*Amu*Lambda^4*s^4 + 36864*AG^2*Amu^2*Lambda^4*s^4 + 16384*AG*Amu^3*Lambda^4*s^4 + 1024*Amu^4*Lambda^4*s^4 + 6144*AG^4*Lambda^5*s^4 + 24576*AG^3*Amu*Lambda^5*s^4 + 12288*AG^2*Amu^2*Lambda^5*s^4 + 6144*AG^3*Arho*Lambda^5*s^4 + 12288*AG^2*Amu*Arho*Lambda^5*s^4 + 3072*AG*Amu^2*Arho*Lambda^5*s^4 + 64*AG^2*Lambda^6*s^4 + 1472*AG^4*Lambda^6*s^4 + 1536*AG^3*Arho*Lambda^6*s^4 + 384*AG^2*Arho^2*Lambda^6*s^4 - 128*AG^4*Arho^2*Lambda^6*s^4 + 4096*AG^3*Amu*Lambda^3*s^5 + 24576*AG^2*Amu^2*Lambda^3*s^5 + 24576*AG*Amu^3*Lambda^3*s^5 + 4096*Amu^4*Lambda^3*s^5 + 2048*AG^4*Lambda^4*s^5 + 24576*AG^3*Amu*Lambda^4*s^5 + 36864*AG^2*Amu^2*Lambda^4*s^5 + 8192*AG*Amu^3*Lambda^4*s^5 + 4096*AG^3*Arho*Lambda^4*s^5 + 18432*AG^2*Amu*Arho*Lambda^4*s^5 + 12288*AG*Amu^2*Arho*Lambda^4*s^5 + 1024*Amu^3*Arho*Lambda^4*s^5 + 256*AG^2*Lambda^5*s^5 + 2944*AG^4*Lambda^5*s^5 + 128*AG*Amu*Lambda^5*s^5 + 5888*AG^3*Amu*Lambda^5*s^5 + 4608*AG^3*Arho*Lambda^5*s^5 + 4608*AG^2*Amu*Arho*Lambda^5*s^5 + 1536*AG^2*Arho^2*Lambda^5*s^5 - 256*AG^4*Arho^2*Lambda^5*s^5 + 768*AG*Amu*Arho^2*Lambda^5*s^5 - 512*AG^3*Amu*Arho^2*Lambda^5*s^5 + 6144*AG^2*Amu^2*Lambda^2*s^6 + 16384*AG*Amu^3*Lambda^2*s^6 + 6144*Amu^4*Lambda^2*s^6 + 8192*AG^3*Amu*Lambda^3*s^6 + 36864*AG^2*Amu^2*Lambda^3*s^6 + 24576*AG*Amu^3*Lambda^3*s^6 + 2048*Amu^4*Lambda^3*s^6 + 1024*AG^3*Arho*Lambda^3*s^6 + 12288*AG^2*Amu*Arho*Lambda^3*s^6 + 18432*AG*Amu^2*Arho*Lambda^3*s^6 + 4096*Amu^3*Arho*Lambda^3*s^6 + 384*AG^2*Lambda^4*s^6 + 1472*AG^4*Lambda^4*s^6 + 512*AG*Amu*Lambda^4*s^6 + 11776*AG^3*Amu*Lambda^4*s^6 + 64*Amu^2*Lambda^4*s^6 + 8832*AG^2*Amu^2*Lambda^4*s^6 + 4608*AG^3*Arho*Lambda^4*s^6 + 13824*AG^2*Amu*Arho*Lambda^4*s^6 + 4608*AG*Amu^2*Arho*Lambda^4*s^6 + 2304*AG^2*Arho^2*Lambda^4*s^6 - 128*AG^4*Arho^2*Lambda^4*s^6 + 3072*AG*Amu*Arho^2*Lambda^4*s^6 - 1024*AG^3*Amu*Arho^2*Lambda^4*s^6 + 384*Amu^2*Arho^2*Lambda^4*s^6 - 768*AG^2*Amu^2*Arho^2*Lambda^4*s^6 + 128*AG^2*Lambda^5*s^6 + 448*AG^4*Lambda^5*s^6 + 32*AG*Arho*Lambda^5*s^6 + 800*AG^3*Arho*Lambda^5*s^6 + 384*AG^2*Arho^2*Lambda^5*s^6 - 128*AG^4*Arho^2*Lambda^5*s^6 + 64*AG*Arho^3*Lambda^5*s^6 - 64*AG^3*Arho^3*Lambda^5*s^6 + 4096*AG*Amu^3*Lambda*s^7 + 4096*Amu^4*Lambda*s^7 + 12288*AG^2*Amu^2*Lambda^2*s^7 + 24576*AG*Amu^3*Lambda^2*s^7 + 6144*Amu^4*Lambda^2*s^7 + 3072*AG^2*Amu*Arho*Lambda^2*s^7 + 12288*AG*Amu^2*Arho*Lambda^2*s^7 + 6144*Amu^3*Arho*Lambda^2*s^7 + 256*AG^2*Lambda^3*s^7 + 768*AG*Amu*Lambda^3*s^7 + 5888*AG^3*Amu*Lambda^3*s^7 + 256*Amu^2*Lambda^3*s^7 + 17664*AG^2*Amu^2*Lambda^3*s^7 + 5888*AG*Amu^3*Lambda^3*s^7 + 1536*AG^3*Arho*Lambda^3*s^7 + 13824*AG^2*Amu*Arho*Lambda^3*s^7 + 13824*AG*Amu^2*Arho*Lambda^3*s^7 + 1536*Amu^3*Arho*Lambda^3*s^7 + 1536*AG^2*Arho^2*Lambda^3*s^7 + 4608*AG*Amu*Arho^2*Lambda^3*s^7 - 512*AG^3*Amu*Arho^2*Lambda^3*s^7 + 1536*Amu^2*Arho^2*Lambda^3*s^7 - 1536*AG^2*Amu^2*Arho^2*Lambda^3*s^7 - 512*AG*Amu^3*Arho^2*Lambda^3*s^7 + 384*AG^2*Lambda^4*s^7 + 448*AG^4*Lambda^4*s^7 + 256*AG*Amu*Lambda^4*s^7 + 1792*AG^3*Amu*Lambda^4*s^7 + 128*AG*Arho*Lambda^4*s^7 + 1600*AG^3*Arho*Lambda^4*s^7 + 32*Amu*Arho*Lambda^4*s^7 + 2400*AG^2*Amu*Arho*Lambda^4*s^7 + 1152*AG^2*Arho^2*Lambda^4*s^7 - 128*AG^4*Arho^2*Lambda^4*s^7 + 768*AG*Amu*Arho^2*Lambda^4*s^7 - 512*AG^3*Amu*Arho^2*Lambda^4*s^7 + 256*AG*Arho^3*Lambda^4*s^7 - 128*AG^3*Arho^3*Lambda^4*s^7 + 64*Amu*Arho^3*Lambda^4*s^7 - 192*AG^2*Amu*Arho^3*Lambda^4*s^7 + 1024*Amu^4*s^8 + 8192*AG*Amu^3*Lambda*s^8 + 6144*Amu^4*Lambda*s^8 + 3072*AG*Amu^2*Arho*Lambda*s^8 + 4096*Amu^3*Arho*Lambda*s^8 + 64*AG^2*Lambda^2*s^8 + 512*AG*Amu*Lambda^2*s^8 + 384*Amu^2*Lambda^2*s^8 + 8832*AG^2*Amu^2*Lambda^2*s^8 + 11776*AG*Amu^3*Lambda^2*s^8 + 1472*Amu^4*Lambda^2*s^8 + 4608*AG^2*Amu*Arho*Lambda^2*s^8 + 13824*AG*Amu^2*Arho*Lambda^2*s^8 + 4608*Amu^3*Arho*Lambda^2*s^8 + 384*AG^2*Arho^2*Lambda^2*s^8 + 3072*AG*Amu*Arho^2*Lambda^2*s^8 + 2304*Amu^2*Arho^2*Lambda^2*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^8 - 1024*AG*Amu^3*Arho^2*Lambda^2*s^8 - 128*Amu^4*Arho^2*Lambda^2*s^8 + 384*AG^2*Lambda^3*s^8 + 768*AG*Amu*Lambda^3*s^8 + 1792*AG^3*Amu*Lambda^3*s^8 + 128*Amu^2*Lambda^3*s^8 + 2688*AG^2*Amu^2*Lambda^3*s^8 + 192*AG*Arho*Lambda^3*s^8 + 800*AG^3*Arho*Lambda^3*s^8 + 128*Amu*Arho*Lambda^3*s^8 + 4800*AG^2*Amu*Arho*Lambda^3*s^8 + 2400*AG*Amu^2*Arho*Lambda^3*s^8 + 1152*AG^2*Arho^2*Lambda^3*s^8 + 2304*AG*Amu*Arho^2*Lambda^3*s^8 - 512*AG^3*Amu*Arho^2*Lambda^3*s^8 + 384*Amu^2*Arho^2*Lambda^3*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^3*s^8 + 384*AG*Arho^3*Lambda^3*s^8 - 64*AG^3*Arho^3*Lambda^3*s^8 + 256*Amu*Arho^3*Lambda^3*s^8 - 384*AG^2*Amu*Arho^3*Lambda^3*s^8 - 192*AG*Amu^2*Arho^3*Lambda^3*s^8 + 80*AG^2*Lambda^4*s^8 + 48*AG^4*Lambda^4*s^8 + 48*AG*Arho*Lambda^4*s^8 + 176*AG^3*Arho*Lambda^4*s^8 + 4*Arho^2*Lambda^4*s^8 + 120*AG^2*Arho^2*Lambda^4*s^8 - 28*AG^4*Arho^2*Lambda^4*s^8 + 32*AG*Arho^3*Lambda^4*s^8 - 32*AG^3*Arho^3*Lambda^4*s^8 + 4*Arho^4*Lambda^4*s^8 - 8*AG^2*Arho^4*Lambda^4*s^8 + 4*AG^4*Arho^4*Lambda^4*s^8 + 2048*Amu^4*s^9 + 1024*Amu^3*Arho*s^9 + 128*AG*Amu*Lambda*s^9 + 256*Amu^2*Lambda*s^9 + 5888*AG*Amu^3*Lambda*s^9 + 2944*Amu^4*Lambda*s^9 + 4608*AG*Amu^2*Arho*Lambda*s^9 + 4608*Amu^3*Arho*Lambda*s^9 + 768*AG*Amu*Arho^2*Lambda*s^9 + 1536*Amu^2*Arho^2*Lambda*s^9 - 512*AG*Amu^3*Arho^2*Lambda*s^9 - 256*Amu^4*Arho^2*Lambda*s^9 + 128*AG^2*Lambda^2*s^9 + 768*AG*Amu*Lambda^2*s^9 + 384*Amu^2*Lambda^2*s^9 + 2688*AG^2*Amu^2*Lambda^2*s^9 + 1792*AG*Amu^3*Lambda^2*s^9 + 128*AG*Arho*Lambda^2*s^9 + 192*Amu*Arho*Lambda^2*s^9 + 2400*AG^2*Amu*Arho*Lambda^2*s^9 + 4800*AG*Amu^2*Arho*Lambda^2*s^9 + 800*Amu^3*Arho*Lambda^2*s^9 + 384*AG^2*Arho^2*Lambda^2*s^9 + 2304*AG*Amu*Arho^2*Lambda^2*s^9 + 1152*Amu^2*Arho^2*Lambda^2*s^9 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^9 - 512*AG*Amu^3*Arho^2*Lambda^2*s^9 + 256*AG*Arho^3*Lambda^2*s^9 + 384*Amu*Arho^3*Lambda^2*s^9 - 192*AG^2*Amu*Arho^3*Lambda^2*s^9 - 384*AG*Amu^2*Arho^3*Lambda^2*s^9 - 64*Amu^3*Arho^3*Lambda^2*s^9 + 160*AG^2*Lambda^3*s^9 + 160*AG*Amu*Lambda^3*s^9 + 192*AG^3*Amu*Lambda^3*s^9 + 144*AG*Arho*Lambda^3*s^9 + 176*AG^3*Arho*Lambda^3*s^9 + 48*Amu*Arho*Lambda^3*s^9 + 528*AG^2*Amu*Arho*Lambda^3*s^9 + 16*Arho^2*Lambda^3*s^9 + 240*AG^2*Arho^2*Lambda^3*s^9 + 240*AG*Amu*Arho^2*Lambda^3*s^9 - 112*AG^3*Amu*Arho^2*Lambda^3*s^9 + 96*AG*Arho^3*Lambda^3*s^9 - 32*AG^3*Arho^3*Lambda^3*s^9 + 32*Amu*Arho^3*Lambda^3*s^9 - 96*AG^2*Amu*Arho^3*Lambda^3*s^9 + 16*Arho^4*Lambda^3*s^9 - 16*AG^2*Arho^4*Lambda^3*s^9 - 16*AG*Amu*Arho^4*Lambda^3*s^9 + 16*AG^3*Amu*Arho^4*Lambda^3*s^9 + 64*Amu^2*s^10 + 1472*Amu^4*s^10 + 1536*Amu^3*Arho*s^10 + 384*Amu^2*Arho^2*s^10 - 128*Amu^4*Arho^2*s^10 + 256*AG*Amu*Lambda*s^10 + 384*Amu^2*Lambda*s^10 + 1792*AG*Amu^3*Lambda*s^10 + 448*Amu^4*Lambda*s^10 + 32*AG*Arho*Lambda*s^10 + 128*Amu*Arho*Lambda*s^10 + 2400*AG*Amu^2*Arho*Lambda*s^10 + 1600*Amu^3*Arho*Lambda*s^10 + 768*AG*Amu*Arho^2*Lambda*s^10 + 1152*Amu^2*Arho^2*Lambda*s^10 - 512*AG*Amu^3*Arho^2*Lambda*s^10 - 128*Amu^4*Arho^2*Lambda*s^10 + 64*AG*Arho^3*Lambda*s^10 + 256*Amu*Arho^3*Lambda*s^10 - 192*AG*Amu^2*Arho^3*Lambda*s^10 - 128*Amu^3*Arho^3*Lambda*s^10 + 80*AG^2*Lambda^2*s^10 + 320*AG*Amu*Lambda^2*s^10 + 80*Amu^2*Lambda^2*s^10 + 288*AG^2*Amu^2*Lambda^2*s^10 + 144*AG*Arho*Lambda^2*s^10 + 144*Amu*Arho*Lambda^2*s^10 + 528*AG^2*Amu*Arho*Lambda^2*s^10 + 528*AG*Amu^2*Arho*Lambda^2*s^10 + 24*Arho^2*Lambda^2*s^10 + 120*AG^2*Arho^2*Lambda^2*s^10 + 480*AG*Amu*Arho^2*Lambda^2*s^10 + 120*Amu^2*Arho^2*Lambda^2*s^10 - 168*AG^2*Amu^2*Arho^2*Lambda^2*s^10 + 96*AG*Arho^3*Lambda^2*s^10 + 96*Amu*Arho^3*Lambda^2*s^10 - 96*AG^2*Amu*Arho^3*Lambda^2*s^10 - 96*AG*Amu^2*Arho^3*Lambda^2*s^10 + 24*Arho^4*Lambda^2*s^10 - 8*AG^2*Arho^4*Lambda^2*s^10 - 32*AG*Amu*Arho^4*Lambda^2*s^10 - 8*Amu^2*Arho^4*Lambda^2*s^10 + 24*AG^2*Amu^2*Arho^4*Lambda^2*s^10 + 16*AG^2*Lambda^3*s^10 + 20*AG*Arho*Lambda^3*s^10 + 12*AG^3*Arho*Lambda^3*s^10 + 4*Arho^2*Lambda^3*s^10 + 12*AG^2*Arho^2*Lambda^3*s^10 + 4*AG*Arho^3*Lambda^3*s^10 - 4*AG^3*Arho^3*Lambda^3*s^10 + 128*Amu^2*s^11 + 448*Amu^4*s^11 + 32*Amu*Arho*s^11 + 800*Amu^3*Arho*s^11 + 384*Amu^2*Arho^2*s^11 - 128*Amu^4*Arho^2*s^11 + 64*Amu*Arho^3*s^11 - 64*Amu^3*Arho^3*s^11 + 160*AG*Amu*Lambda*s^11 + 160*Amu^2*Lambda*s^11 + 192*AG*Amu^3*Lambda*s^11 + 48*AG*Arho*Lambda*s^11 + 144*Amu*Arho*Lambda*s^11 + 528*AG*Amu^2*Arho*Lambda*s^11 + 176*Amu^3*Arho*Lambda*s^11 + 16*Arho^2*Lambda*s^11 + 240*AG*Amu*Arho^2*Lambda*s^11 + 240*Amu^2*Arho^2*Lambda*s^11 - 112*AG*Amu^3*Arho^2*Lambda*s^11 + 32*AG*Arho^3*Lambda*s^11 + 96*Amu*Arho^3*Lambda*s^11 - 96*AG*Amu^2*Arho^3*Lambda*s^11 - 32*Amu^3*Arho^3*Lambda*s^11 + 16*Arho^4*Lambda*s^11 - 16*AG*Amu*Arho^4*Lambda*s^11 - 16*Amu^2*Arho^4*Lambda*s^11 + 16*AG*Amu^3*Arho^4*Lambda*s^11 + 16*AG^2*Lambda^2*s^11 + 32*AG*Amu*Lambda^2*s^11 + 40*AG*Arho*Lambda^2*s^11 + 20*Amu*Arho*Lambda^2*s^11 + 36*AG^2*Amu*Arho*Lambda^2*s^11 + 12*Arho^2*Lambda^2*s^11 + 12*AG^2*Arho^2*Lambda^2*s^11 + 24*AG*Amu*Arho^2*Lambda^2*s^11 + 8*AG*Arho^3*Lambda^2*s^11 + 4*Amu*Arho^3*Lambda^2*s^11 - 12*AG^2*Amu*Arho^3*Lambda^2*s^11 + 80*Amu^2*s^12 + 48*Amu^4*s^12 + 48*Amu*Arho*s^12 + 176*Amu^3*Arho*s^12 + 4*Arho^2*s^12 + 120*Amu^2*Arho^2*s^12 - 28*Amu^4*Arho^2*s^12 + 32*Amu*Arho^3*s^12 - 32*Amu^3*Arho^3*s^12 + 4*Arho^4*s^12 - 8*Amu^2*Arho^4*s^12 + 4*Amu^4*Arho^4*s^12 + 32*AG*Amu*Lambda*s^12 + 16*Amu^2*Lambda*s^12 + 20*AG*Arho*Lambda*s^12 + 40*Amu*Arho*Lambda*s^12 + 36*AG*Amu^2*Arho*Lambda*s^12 + 12*Arho^2*Lambda*s^12 + 24*AG*Amu*Arho^2*Lambda*s^12 + 12*Amu^2*Arho^2*Lambda*s^12 + 4*AG*Arho^3*Lambda*s^12 + 8*Amu*Arho^3*Lambda*s^12 - 12*AG*Amu^2*Arho^3*Lambda*s^12 + AG^2*Lambda^2*s^12 + 2*AG*Arho*Lambda^2*s^12 + Arho^2*Lambda^2*s^12 + 16*Amu^2*s^13 + 20*Amu*Arho*s^13 + 12*Amu^3*Arho*s^13 + 4*Arho^2*s^13 + 12*Amu^2*Arho^2*s^13 + 4*Amu*Arho^3*s^13 - 4*Amu^3*Arho^3*s^13 + 2*AG*Amu*Lambda*s^13 + 2*AG*Arho*Lambda*s^13 + 2*Amu*Arho*Lambda*s^13 + 2*Arho^2*Lambda*s^13 + Amu^2*s^14 + 2*Amu*Arho*s^14 + Arho^2*s^14)], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0, "Details" -> "Only the nonzero constant factor 4 is divided out. E_+^2 E_-^2 is not divided out of P14."|>}
```

Spurious roots may be introduced by denominator clearing and by the two squaring operations.

The full polynomial and coefficient list are exported as Wolfram expressions in `symbolic/exports/`.


<a id="phase3-report"></a>

## Frozen Phase III report

### Phase III polynomial-completeness audit

Phase III packages the final no-root-loss theorem built from the Phase II derivation.

The certified theorem is:

```math
\Omega_{\mathrm{definition}}\land\Omega_{\mathrm{no\ root\ loss}}\land D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0.
```

Equivalently, every root of the original nondimensional dispersion relation is represented among the roots of the polynomial candidate.

The reverse implication is not claimed. Extra polynomial roots are expected.

#### Authoritative source of checks

Phase II is authoritative for radical elimination, denominator clearing, polynomial construction, coefficient verification, and forward-preserving symbolic identities.

Phase III reuses the Phase II result rather than duplicating the algebra.

#### OmegaDefinition

```wolfram
s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0
```

#### NoRootLossConditions

```wolfram
(1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0 && ((1 - AG)*Lambda + (1 - Amu)*s)^2*((1 + AG)*Lambda + (1 + Amu)*s)^2 != 0
```

#### Theorem

```wolfram
HoldForm[Implies[ForwardImplicationConditions[] && DStarLambda[] == 0, P14Polynomial[] == 0]]
```

#### Phase II forward-preserving checks

```wolfram
{<|"Name" -> "DStarLambda equals F0ClockForm", "Status" -> "PASS", "Statement" -> HoldForm[2 + s*((1 + Amu + ((1 + AG)*Lambda)/s + (1 - Amu + ((1 - AG)*Lambda)/s)*Sqrt[1 + ((1 - Arho)*s)/(1 - Amu + ((1 - AG)*Lambda)/s)])^(-1) + (1 - Amu + ((1 - AG)*Lambda)/s + (1 + Amu + ((1 + AG)*Lambda)/s)*Sqrt[1 + ((1 + Arho)*s)/(1 + Amu + ((1 + AG)*Lambda)/s)])^(-1)) == 2 + s^2*(((1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])^(-1) + ((1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)])^(-1))], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Exact rewriting of the nondimensional dispersion relation into the s^2 clock form."|>, <|"Name" -> "First denominator-clearing step", "Status" -> "PASS", "Statement" -> HoldForm[((1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])*((1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)])*(2 + s^2*(((1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])^(-1) + ((1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)])^(-1))) == 2*Lambda^2 - 2*AG^2*Lambda^2 + 4*Lambda*s - 4*AG*Amu*Lambda*s + 2*s^2 - 2*Amu^2*s^2 + 2*Lambda*s^2 + 2*s^3 + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - AG*Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - Amu*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + AG*Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Amu*s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)]], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Multiplication by reciprocal denominators. This preserves roots in the forward direction on the original-expression domain."|>, <|"Name" -> "Collected F1 radical form", "Status" -> "PASS", "Statement" -> HoldForm[2*Lambda^2 - 2*AG^2*Lambda^2 + 4*Lambda*s - 4*AG*Amu*Lambda*s + 2*s^2 - 2*Amu^2*s^2 + 2*Lambda*s^2 + 2*s^3 + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - AG*Lambda*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - Amu*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*AG^2*Lambda^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Amu*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Amu^2*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + AG*Lambda*s^2*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + Amu*s^3*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*AG^2*Lambda^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 4*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 4*AG*Amu*Lambda*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + 2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] - 2*Amu^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] == 2*((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s) + s^2*((1 - AG)*Lambda + (1 + AG)*Lambda + (1 - Amu)*s + (1 + Amu)*s) + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*(s^2 + 2*((1 - AG)*Lambda + (1 - Amu)*s)) + 2*((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)]*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)]*(s^2 + 2*((1 + AG)*Lambda + (1 + Amu)*s))], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Exact collection into alpha0 + betaMinus qMinus + gammaPlus qPlus + delta qMinus qPlus."|>, <|"Name" -> "First radical elimination preserves forward implication", "Status" -> "PASS", "Statement" -> HoldForm[Implies[F1LinearRadicals[] == 0, F2AfterFirstSquaring[] == 0]], "Assumptions" -> True, "Details" -> "The step squares an equation satisfied by any original root. Squaring may add roots but cannot lose roots satisfying the pre-squared equation."|>, <|"Name" -> "F2 single-radical rewrite", "Status" -> "PASS", "Statement" -> HoldForm[-32*AG*Lambda^4 - 32*AG^3*Lambda^4 - 96*AG*Lambda^3*s - 32*AG^3*Lambda^3*s - 32*Amu*Lambda^3*s - 96*AG^2*Amu*Lambda^3*s - 96*AG*Lambda^2*s^2 - 96*Amu*Lambda^2*s^2 - 96*AG^2*Amu*Lambda^2*s^2 - 96*AG*Amu^2*Lambda^2*s^2 + 8*Lambda^3*s^2 - 24*AG*Lambda^3*s^2 - 8*AG^2*Lambda^3*s^2 - 8*AG^3*Lambda^3*s^2 - 32*AG*Lambda*s^3 - 96*Amu*Lambda*s^3 - 96*AG*Amu^2*Lambda*s^3 - 32*Amu^3*Lambda*s^3 + 24*Lambda^2*s^3 - 48*AG*Lambda^2*s^3 - 8*AG^2*Lambda^2*s^3 - 24*Amu*Lambda^2*s^3 - 16*AG*Amu*Lambda^2*s^3 - 24*AG^2*Amu*Lambda^2*s^3 - 32*Amu*s^4 - 32*Amu^3*s^4 + 24*Lambda*s^4 - 24*AG*Lambda*s^4 - 48*Amu*Lambda*s^4 - 16*AG*Amu*Lambda*s^4 - 8*Amu^2*Lambda*s^4 - 24*AG*Amu^2*Lambda*s^4 + 4*Lambda^2*s^4 - 4*AG*Lambda^2*s^4 + 8*s^5 - 24*Amu*s^5 - 8*Amu^2*s^5 - 8*Amu^3*s^5 + 8*Lambda*s^5 - 4*AG*Lambda*s^5 - 4*Amu*Lambda*s^5 + 4*s^6 - 4*Amu*s^6 - (16*AG*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG^2*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG^2*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (128*AG*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (128*AG*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (AG^2*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Amu^3*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Amu^3*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + s^8/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (Amu^2*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Amu^2*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (8*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - s^8/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*AG^4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*AG^4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^3*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^2*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^3*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (32*AG*Amu*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*AG^2*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (32*AG*Amu*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*AG^2*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Amu^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu^3*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Amu^2*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu^3*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*Amu^4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Amu^4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - 32*AG*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^3*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*AG*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*AG*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG^2*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^3*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - (8*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) == -32*AG*Lambda^4 - 32*AG^3*Lambda^4 - 96*AG*Lambda^3*s - 32*AG^3*Lambda^3*s - 32*Amu*Lambda^3*s - 96*AG^2*Amu*Lambda^3*s - 96*AG*Lambda^2*s^2 - 96*Amu*Lambda^2*s^2 - 96*AG^2*Amu*Lambda^2*s^2 - 96*AG*Amu^2*Lambda^2*s^2 + 8*Lambda^3*s^2 - 24*AG*Lambda^3*s^2 - 8*AG^2*Lambda^3*s^2 - 8*AG^3*Lambda^3*s^2 - 32*AG*Lambda*s^3 - 96*Amu*Lambda*s^3 - 96*AG*Amu^2*Lambda*s^3 - 32*Amu^3*Lambda*s^3 + 24*Lambda^2*s^3 - 48*AG*Lambda^2*s^3 - 8*AG^2*Lambda^2*s^3 - 24*Amu*Lambda^2*s^3 - 16*AG*Amu*Lambda^2*s^3 - 24*AG^2*Amu*Lambda^2*s^3 - 32*Amu*s^4 - 32*Amu^3*s^4 + 24*Lambda*s^4 - 24*AG*Lambda*s^4 - 48*Amu*Lambda*s^4 - 16*AG*Amu*Lambda*s^4 - 8*Amu^2*Lambda*s^4 - 24*AG*Amu^2*Lambda*s^4 + 4*Lambda^2*s^4 - 4*AG*Lambda^2*s^4 + 8*s^5 - 24*Amu*s^5 - 8*Amu^2*s^5 - 8*Amu^3*s^5 + 8*Lambda*s^5 - 4*AG*Lambda*s^5 - 4*Amu*Lambda*s^5 + 4*s^6 - 4*Amu*s^6 - (16*AG*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG^2*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG^2*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*AG^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (128*AG*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*AG^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (128*AG*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*AG*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*AG*Amu*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (64*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*AG*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*AG*Amu*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (64*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*Lambda*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (32*Amu^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (16*Amu^3*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (32*Amu^2*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (16*Amu^3*Arho*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*AG*Amu*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (24*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (AG^2*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu^2*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Amu^3*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (4*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (12*Amu*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (12*Amu^2*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (4*Amu^3*Arho*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Amu*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*AG*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 - AG)*Lambda + (1 - Amu)*s) + s^8/((1 - AG)*Lambda + (1 - Amu)*s) - (2*Amu*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (Amu^2*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) + (2*Amu*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (Amu^2*Arho*s^8)/((1 - AG)*Lambda + (1 - Amu)*s) - (8*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG^4*Arho*Lambda^4*s^2)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^2*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^3*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG^3*Amu*Arho*Lambda^3*s^3)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (64*AG*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG^2*Amu^2*Arho*Lambda^2*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG^3*Arho*Lambda^3*s^4)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Amu^2*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (32*AG*Amu^3*Arho*Lambda*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG^2*Amu*Arho*Lambda^2*s^5)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^2*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu^3*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu^4*Arho*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (24*AG*Amu*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*AG*Amu^2*Arho*Lambda*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (AG^2*Arho*Lambda^2*s^6)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Amu^2*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu^3*Arho*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*AG*Amu*Arho*Lambda*s^7)/((1 + AG)*Lambda + (1 + Amu)*s) - s^8/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (2*Amu*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (Amu^2*Arho*s^8)/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*AG^4*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*AG^4*Arho^2*Lambda^4*s^4)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^3*Amu*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG^2*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG^3*Amu*Arho^2*Lambda^3*s^5)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*AG^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (32*AG*Amu*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (24*AG^2*Amu^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*AG^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (32*AG*Amu*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (24*AG^2*Amu^2*Arho^2*Lambda^2*s^6)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Amu^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu^3*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*AG*Amu*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (16*Amu^2*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (16*AG*Amu^3*Arho^2*Lambda*s^7)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (8*Amu^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (4*Amu^4*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - (8*Amu^2*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) + (4*Amu^4*Arho^2*s^8)/(((1 - AG)*Lambda + (1 - Amu)*s)*((1 + AG)*Lambda + (1 + Amu)*s)) - 32*AG*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*AG^3*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^3*s*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*AG*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG^2*Amu*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda^2*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^3*Lambda^3*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*AG*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 96*Amu*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 96*AG*Amu^2*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*Lambda*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*AG*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*AG^2*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG^2*Amu*Lambda^2*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 32*Amu*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 32*Amu^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*AG*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 48*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 16*AG*Amu*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 24*AG*Amu^2*Lambda*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 24*Amu*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Amu^3*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 8*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] + 4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - 4*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] - (8*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG^4*Arho*Lambda^4*s^2*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*AG^3*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG^3*Amu*Arho*Lambda^3*s^3*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*AG*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG^2*Amu^2*Arho*Lambda^2*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^3*Arho*Lambda^3*s^4*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (32*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*AG*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (48*Amu*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (48*AG*Amu^2*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (32*AG*Amu^3*Arho*Lambda*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*AG*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*AG^2*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG^2*Amu*Arho*Lambda^2*s^5*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (16*Amu*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (16*Amu^3*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*Amu^4*Arho*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (12*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*AG*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (8*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (8*AG*Amu*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (12*AG*Amu^2*Arho*Lambda*s^6*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) - (4*Amu*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^2*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s) + (4*Amu^3*Arho*s^7*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)])/((1 + AG)*Lambda + (1 + Amu)*s)], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0 && s != 0 && (1 + AG)*Lambda + (1 + Amu)*s != 0 && (1 - AG)*Lambda + (1 - Amu)*s != 0 && (1 + AG)*Lambda + (1 + Amu)*s + ((1 - AG)*Lambda + (1 - Amu)*s)*Sqrt[1 + ((1 - Arho)*s^2)/((1 - AG)*Lambda + (1 - Amu)*s)] != 0 && (1 - AG)*Lambda + (1 - Amu)*s + ((1 + AG)*Lambda + (1 + Amu)*s)*Sqrt[1 + ((1 + Arho)*s^2)/((1 + AG)*Lambda + (1 + Amu)*s)] != 0, "Details" -> "Uses qMinus^2 = radicandMinus with Mathematica's principal Sqrt convention; no PowerExpand or branch replacement is used."|>, <|"Name" -> "Second radical elimination preserves forward implication", "Status" -> "PASS", "Statement" -> HoldForm[Implies[F3SingleRadical[] == 0, F4AfterSecondSquaring[] == 0]], "Assumptions" -> True, "Details" -> "The second squaring may add roots but cannot lose roots satisfying F3 == 0."|>, <|"Name" -> "Rational numerator extraction preserves forward implication", "Status" -> "PASS", "Statement" -> HoldForm[Implies[F4AfterSecondSquaring[] == 0 && Denominator[Together[F4AfterSecondSquaring[]]] != 0, rawPolynomialNumerator[] == 0]], "Assumptions" -> HoldForm[Denominator[Together[F4AfterSecondSquaring[]]] != 0], "Details" -> "If a defined rational expression is zero and its denominator is nonzero, then its numerator is zero."|>, <|"Name" -> "Final expression equals P14 up to constant normalization", "Status" -> "PASS", "Statement" -> HoldForm[4*(1024*AG^4*Lambda^8 + 4096*AG^4*Lambda^7*s + 4096*AG^3*Amu*Lambda^7*s + 6144*AG^4*Lambda^6*s^2 + 16384*AG^3*Amu*Lambda^6*s^2 + 6144*AG^2*Amu^2*Lambda^6*s^2 + 2048*AG^4*Lambda^7*s^2 + 1024*AG^3*Arho*Lambda^7*s^2 + 4096*AG^4*Lambda^5*s^3 + 24576*AG^3*Amu*Lambda^5*s^3 + 24576*AG^2*Amu^2*Lambda^5*s^3 + 4096*AG*Amu^3*Lambda^5*s^3 + 6144*AG^4*Lambda^6*s^3 + 8192*AG^3*Amu*Lambda^6*s^3 + 4096*AG^3*Arho*Lambda^6*s^3 + 3072*AG^2*Amu*Arho*Lambda^6*s^3 + 1024*AG^4*Lambda^4*s^4 + 16384*AG^3*Amu*Lambda^4*s^4 + 36864*AG^2*Amu^2*Lambda^4*s^4 + 16384*AG*Amu^3*Lambda^4*s^4 + 1024*Amu^4*Lambda^4*s^4 + 6144*AG^4*Lambda^5*s^4 + 24576*AG^3*Amu*Lambda^5*s^4 + 12288*AG^2*Amu^2*Lambda^5*s^4 + 6144*AG^3*Arho*Lambda^5*s^4 + 12288*AG^2*Amu*Arho*Lambda^5*s^4 + 3072*AG*Amu^2*Arho*Lambda^5*s^4 + 64*AG^2*Lambda^6*s^4 + 1472*AG^4*Lambda^6*s^4 + 1536*AG^3*Arho*Lambda^6*s^4 + 384*AG^2*Arho^2*Lambda^6*s^4 - 128*AG^4*Arho^2*Lambda^6*s^4 + 4096*AG^3*Amu*Lambda^3*s^5 + 24576*AG^2*Amu^2*Lambda^3*s^5 + 24576*AG*Amu^3*Lambda^3*s^5 + 4096*Amu^4*Lambda^3*s^5 + 2048*AG^4*Lambda^4*s^5 + 24576*AG^3*Amu*Lambda^4*s^5 + 36864*AG^2*Amu^2*Lambda^4*s^5 + 8192*AG*Amu^3*Lambda^4*s^5 + 4096*AG^3*Arho*Lambda^4*s^5 + 18432*AG^2*Amu*Arho*Lambda^4*s^5 + 12288*AG*Amu^2*Arho*Lambda^4*s^5 + 1024*Amu^3*Arho*Lambda^4*s^5 + 256*AG^2*Lambda^5*s^5 + 2944*AG^4*Lambda^5*s^5 + 128*AG*Amu*Lambda^5*s^5 + 5888*AG^3*Amu*Lambda^5*s^5 + 4608*AG^3*Arho*Lambda^5*s^5 + 4608*AG^2*Amu*Arho*Lambda^5*s^5 + 1536*AG^2*Arho^2*Lambda^5*s^5 - 256*AG^4*Arho^2*Lambda^5*s^5 + 768*AG*Amu*Arho^2*Lambda^5*s^5 - 512*AG^3*Amu*Arho^2*Lambda^5*s^5 + 6144*AG^2*Amu^2*Lambda^2*s^6 + 16384*AG*Amu^3*Lambda^2*s^6 + 6144*Amu^4*Lambda^2*s^6 + 8192*AG^3*Amu*Lambda^3*s^6 + 36864*AG^2*Amu^2*Lambda^3*s^6 + 24576*AG*Amu^3*Lambda^3*s^6 + 2048*Amu^4*Lambda^3*s^6 + 1024*AG^3*Arho*Lambda^3*s^6 + 12288*AG^2*Amu*Arho*Lambda^3*s^6 + 18432*AG*Amu^2*Arho*Lambda^3*s^6 + 4096*Amu^3*Arho*Lambda^3*s^6 + 384*AG^2*Lambda^4*s^6 + 1472*AG^4*Lambda^4*s^6 + 512*AG*Amu*Lambda^4*s^6 + 11776*AG^3*Amu*Lambda^4*s^6 + 64*Amu^2*Lambda^4*s^6 + 8832*AG^2*Amu^2*Lambda^4*s^6 + 4608*AG^3*Arho*Lambda^4*s^6 + 13824*AG^2*Amu*Arho*Lambda^4*s^6 + 4608*AG*Amu^2*Arho*Lambda^4*s^6 + 2304*AG^2*Arho^2*Lambda^4*s^6 - 128*AG^4*Arho^2*Lambda^4*s^6 + 3072*AG*Amu*Arho^2*Lambda^4*s^6 - 1024*AG^3*Amu*Arho^2*Lambda^4*s^6 + 384*Amu^2*Arho^2*Lambda^4*s^6 - 768*AG^2*Amu^2*Arho^2*Lambda^4*s^6 + 128*AG^2*Lambda^5*s^6 + 448*AG^4*Lambda^5*s^6 + 32*AG*Arho*Lambda^5*s^6 + 800*AG^3*Arho*Lambda^5*s^6 + 384*AG^2*Arho^2*Lambda^5*s^6 - 128*AG^4*Arho^2*Lambda^5*s^6 + 64*AG*Arho^3*Lambda^5*s^6 - 64*AG^3*Arho^3*Lambda^5*s^6 + 4096*AG*Amu^3*Lambda*s^7 + 4096*Amu^4*Lambda*s^7 + 12288*AG^2*Amu^2*Lambda^2*s^7 + 24576*AG*Amu^3*Lambda^2*s^7 + 6144*Amu^4*Lambda^2*s^7 + 3072*AG^2*Amu*Arho*Lambda^2*s^7 + 12288*AG*Amu^2*Arho*Lambda^2*s^7 + 6144*Amu^3*Arho*Lambda^2*s^7 + 256*AG^2*Lambda^3*s^7 + 768*AG*Amu*Lambda^3*s^7 + 5888*AG^3*Amu*Lambda^3*s^7 + 256*Amu^2*Lambda^3*s^7 + 17664*AG^2*Amu^2*Lambda^3*s^7 + 5888*AG*Amu^3*Lambda^3*s^7 + 1536*AG^3*Arho*Lambda^3*s^7 + 13824*AG^2*Amu*Arho*Lambda^3*s^7 + 13824*AG*Amu^2*Arho*Lambda^3*s^7 + 1536*Amu^3*Arho*Lambda^3*s^7 + 1536*AG^2*Arho^2*Lambda^3*s^7 + 4608*AG*Amu*Arho^2*Lambda^3*s^7 - 512*AG^3*Amu*Arho^2*Lambda^3*s^7 + 1536*Amu^2*Arho^2*Lambda^3*s^7 - 1536*AG^2*Amu^2*Arho^2*Lambda^3*s^7 - 512*AG*Amu^3*Arho^2*Lambda^3*s^7 + 384*AG^2*Lambda^4*s^7 + 448*AG^4*Lambda^4*s^7 + 256*AG*Amu*Lambda^4*s^7 + 1792*AG^3*Amu*Lambda^4*s^7 + 128*AG*Arho*Lambda^4*s^7 + 1600*AG^3*Arho*Lambda^4*s^7 + 32*Amu*Arho*Lambda^4*s^7 + 2400*AG^2*Amu*Arho*Lambda^4*s^7 + 1152*AG^2*Arho^2*Lambda^4*s^7 - 128*AG^4*Arho^2*Lambda^4*s^7 + 768*AG*Amu*Arho^2*Lambda^4*s^7 - 512*AG^3*Amu*Arho^2*Lambda^4*s^7 + 256*AG*Arho^3*Lambda^4*s^7 - 128*AG^3*Arho^3*Lambda^4*s^7 + 64*Amu*Arho^3*Lambda^4*s^7 - 192*AG^2*Amu*Arho^3*Lambda^4*s^7 + 1024*Amu^4*s^8 + 8192*AG*Amu^3*Lambda*s^8 + 6144*Amu^4*Lambda*s^8 + 3072*AG*Amu^2*Arho*Lambda*s^8 + 4096*Amu^3*Arho*Lambda*s^8 + 64*AG^2*Lambda^2*s^8 + 512*AG*Amu*Lambda^2*s^8 + 384*Amu^2*Lambda^2*s^8 + 8832*AG^2*Amu^2*Lambda^2*s^8 + 11776*AG*Amu^3*Lambda^2*s^8 + 1472*Amu^4*Lambda^2*s^8 + 4608*AG^2*Amu*Arho*Lambda^2*s^8 + 13824*AG*Amu^2*Arho*Lambda^2*s^8 + 4608*Amu^3*Arho*Lambda^2*s^8 + 384*AG^2*Arho^2*Lambda^2*s^8 + 3072*AG*Amu*Arho^2*Lambda^2*s^8 + 2304*Amu^2*Arho^2*Lambda^2*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^8 - 1024*AG*Amu^3*Arho^2*Lambda^2*s^8 - 128*Amu^4*Arho^2*Lambda^2*s^8 + 384*AG^2*Lambda^3*s^8 + 768*AG*Amu*Lambda^3*s^8 + 1792*AG^3*Amu*Lambda^3*s^8 + 128*Amu^2*Lambda^3*s^8 + 2688*AG^2*Amu^2*Lambda^3*s^8 + 192*AG*Arho*Lambda^3*s^8 + 800*AG^3*Arho*Lambda^3*s^8 + 128*Amu*Arho*Lambda^3*s^8 + 4800*AG^2*Amu*Arho*Lambda^3*s^8 + 2400*AG*Amu^2*Arho*Lambda^3*s^8 + 1152*AG^2*Arho^2*Lambda^3*s^8 + 2304*AG*Amu*Arho^2*Lambda^3*s^8 - 512*AG^3*Amu*Arho^2*Lambda^3*s^8 + 384*Amu^2*Arho^2*Lambda^3*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^3*s^8 + 384*AG*Arho^3*Lambda^3*s^8 - 64*AG^3*Arho^3*Lambda^3*s^8 + 256*Amu*Arho^3*Lambda^3*s^8 - 384*AG^2*Amu*Arho^3*Lambda^3*s^8 - 192*AG*Amu^2*Arho^3*Lambda^3*s^8 + 80*AG^2*Lambda^4*s^8 + 48*AG^4*Lambda^4*s^8 + 48*AG*Arho*Lambda^4*s^8 + 176*AG^3*Arho*Lambda^4*s^8 + 4*Arho^2*Lambda^4*s^8 + 120*AG^2*Arho^2*Lambda^4*s^8 - 28*AG^4*Arho^2*Lambda^4*s^8 + 32*AG*Arho^3*Lambda^4*s^8 - 32*AG^3*Arho^3*Lambda^4*s^8 + 4*Arho^4*Lambda^4*s^8 - 8*AG^2*Arho^4*Lambda^4*s^8 + 4*AG^4*Arho^4*Lambda^4*s^8 + 2048*Amu^4*s^9 + 1024*Amu^3*Arho*s^9 + 128*AG*Amu*Lambda*s^9 + 256*Amu^2*Lambda*s^9 + 5888*AG*Amu^3*Lambda*s^9 + 2944*Amu^4*Lambda*s^9 + 4608*AG*Amu^2*Arho*Lambda*s^9 + 4608*Amu^3*Arho*Lambda*s^9 + 768*AG*Amu*Arho^2*Lambda*s^9 + 1536*Amu^2*Arho^2*Lambda*s^9 - 512*AG*Amu^3*Arho^2*Lambda*s^9 - 256*Amu^4*Arho^2*Lambda*s^9 + 128*AG^2*Lambda^2*s^9 + 768*AG*Amu*Lambda^2*s^9 + 384*Amu^2*Lambda^2*s^9 + 2688*AG^2*Amu^2*Lambda^2*s^9 + 1792*AG*Amu^3*Lambda^2*s^9 + 128*AG*Arho*Lambda^2*s^9 + 192*Amu*Arho*Lambda^2*s^9 + 2400*AG^2*Amu*Arho*Lambda^2*s^9 + 4800*AG*Amu^2*Arho*Lambda^2*s^9 + 800*Amu^3*Arho*Lambda^2*s^9 + 384*AG^2*Arho^2*Lambda^2*s^9 + 2304*AG*Amu*Arho^2*Lambda^2*s^9 + 1152*Amu^2*Arho^2*Lambda^2*s^9 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^9 - 512*AG*Amu^3*Arho^2*Lambda^2*s^9 + 256*AG*Arho^3*Lambda^2*s^9 + 384*Amu*Arho^3*Lambda^2*s^9 - 192*AG^2*Amu*Arho^3*Lambda^2*s^9 - 384*AG*Amu^2*Arho^3*Lambda^2*s^9 - 64*Amu^3*Arho^3*Lambda^2*s^9 + 160*AG^2*Lambda^3*s^9 + 160*AG*Amu*Lambda^3*s^9 + 192*AG^3*Amu*Lambda^3*s^9 + 144*AG*Arho*Lambda^3*s^9 + 176*AG^3*Arho*Lambda^3*s^9 + 48*Amu*Arho*Lambda^3*s^9 + 528*AG^2*Amu*Arho*Lambda^3*s^9 + 16*Arho^2*Lambda^3*s^9 + 240*AG^2*Arho^2*Lambda^3*s^9 + 240*AG*Amu*Arho^2*Lambda^3*s^9 - 112*AG^3*Amu*Arho^2*Lambda^3*s^9 + 96*AG*Arho^3*Lambda^3*s^9 - 32*AG^3*Arho^3*Lambda^3*s^9 + 32*Amu*Arho^3*Lambda^3*s^9 - 96*AG^2*Amu*Arho^3*Lambda^3*s^9 + 16*Arho^4*Lambda^3*s^9 - 16*AG^2*Arho^4*Lambda^3*s^9 - 16*AG*Amu*Arho^4*Lambda^3*s^9 + 16*AG^3*Amu*Arho^4*Lambda^3*s^9 + 64*Amu^2*s^10 + 1472*Amu^4*s^10 + 1536*Amu^3*Arho*s^10 + 384*Amu^2*Arho^2*s^10 - 128*Amu^4*Arho^2*s^10 + 256*AG*Amu*Lambda*s^10 + 384*Amu^2*Lambda*s^10 + 1792*AG*Amu^3*Lambda*s^10 + 448*Amu^4*Lambda*s^10 + 32*AG*Arho*Lambda*s^10 + 128*Amu*Arho*Lambda*s^10 + 2400*AG*Amu^2*Arho*Lambda*s^10 + 1600*Amu^3*Arho*Lambda*s^10 + 768*AG*Amu*Arho^2*Lambda*s^10 + 1152*Amu^2*Arho^2*Lambda*s^10 - 512*AG*Amu^3*Arho^2*Lambda*s^10 - 128*Amu^4*Arho^2*Lambda*s^10 + 64*AG*Arho^3*Lambda*s^10 + 256*Amu*Arho^3*Lambda*s^10 - 192*AG*Amu^2*Arho^3*Lambda*s^10 - 128*Amu^3*Arho^3*Lambda*s^10 + 80*AG^2*Lambda^2*s^10 + 320*AG*Amu*Lambda^2*s^10 + 80*Amu^2*Lambda^2*s^10 + 288*AG^2*Amu^2*Lambda^2*s^10 + 144*AG*Arho*Lambda^2*s^10 + 144*Amu*Arho*Lambda^2*s^10 + 528*AG^2*Amu*Arho*Lambda^2*s^10 + 528*AG*Amu^2*Arho*Lambda^2*s^10 + 24*Arho^2*Lambda^2*s^10 + 120*AG^2*Arho^2*Lambda^2*s^10 + 480*AG*Amu*Arho^2*Lambda^2*s^10 + 120*Amu^2*Arho^2*Lambda^2*s^10 - 168*AG^2*Amu^2*Arho^2*Lambda^2*s^10 + 96*AG*Arho^3*Lambda^2*s^10 + 96*Amu*Arho^3*Lambda^2*s^10 - 96*AG^2*Amu*Arho^3*Lambda^2*s^10 - 96*AG*Amu^2*Arho^3*Lambda^2*s^10 + 24*Arho^4*Lambda^2*s^10 - 8*AG^2*Arho^4*Lambda^2*s^10 - 32*AG*Amu*Arho^4*Lambda^2*s^10 - 8*Amu^2*Arho^4*Lambda^2*s^10 + 24*AG^2*Amu^2*Arho^4*Lambda^2*s^10 + 16*AG^2*Lambda^3*s^10 + 20*AG*Arho*Lambda^3*s^10 + 12*AG^3*Arho*Lambda^3*s^10 + 4*Arho^2*Lambda^3*s^10 + 12*AG^2*Arho^2*Lambda^3*s^10 + 4*AG*Arho^3*Lambda^3*s^10 - 4*AG^3*Arho^3*Lambda^3*s^10 + 128*Amu^2*s^11 + 448*Amu^4*s^11 + 32*Amu*Arho*s^11 + 800*Amu^3*Arho*s^11 + 384*Amu^2*Arho^2*s^11 - 128*Amu^4*Arho^2*s^11 + 64*Amu*Arho^3*s^11 - 64*Amu^3*Arho^3*s^11 + 160*AG*Amu*Lambda*s^11 + 160*Amu^2*Lambda*s^11 + 192*AG*Amu^3*Lambda*s^11 + 48*AG*Arho*Lambda*s^11 + 144*Amu*Arho*Lambda*s^11 + 528*AG*Amu^2*Arho*Lambda*s^11 + 176*Amu^3*Arho*Lambda*s^11 + 16*Arho^2*Lambda*s^11 + 240*AG*Amu*Arho^2*Lambda*s^11 + 240*Amu^2*Arho^2*Lambda*s^11 - 112*AG*Amu^3*Arho^2*Lambda*s^11 + 32*AG*Arho^3*Lambda*s^11 + 96*Amu*Arho^3*Lambda*s^11 - 96*AG*Amu^2*Arho^3*Lambda*s^11 - 32*Amu^3*Arho^3*Lambda*s^11 + 16*Arho^4*Lambda*s^11 - 16*AG*Amu*Arho^4*Lambda*s^11 - 16*Amu^2*Arho^4*Lambda*s^11 + 16*AG*Amu^3*Arho^4*Lambda*s^11 + 16*AG^2*Lambda^2*s^11 + 32*AG*Amu*Lambda^2*s^11 + 40*AG*Arho*Lambda^2*s^11 + 20*Amu*Arho*Lambda^2*s^11 + 36*AG^2*Amu*Arho*Lambda^2*s^11 + 12*Arho^2*Lambda^2*s^11 + 12*AG^2*Arho^2*Lambda^2*s^11 + 24*AG*Amu*Arho^2*Lambda^2*s^11 + 8*AG*Arho^3*Lambda^2*s^11 + 4*Amu*Arho^3*Lambda^2*s^11 - 12*AG^2*Amu*Arho^3*Lambda^2*s^11 + 80*Amu^2*s^12 + 48*Amu^4*s^12 + 48*Amu*Arho*s^12 + 176*Amu^3*Arho*s^12 + 4*Arho^2*s^12 + 120*Amu^2*Arho^2*s^12 - 28*Amu^4*Arho^2*s^12 + 32*Amu*Arho^3*s^12 - 32*Amu^3*Arho^3*s^12 + 4*Arho^4*s^12 - 8*Amu^2*Arho^4*s^12 + 4*Amu^4*Arho^4*s^12 + 32*AG*Amu*Lambda*s^12 + 16*Amu^2*Lambda*s^12 + 20*AG*Arho*Lambda*s^12 + 40*Amu*Arho*Lambda*s^12 + 36*AG*Amu^2*Arho*Lambda*s^12 + 12*Arho^2*Lambda*s^12 + 24*AG*Amu*Arho^2*Lambda*s^12 + 12*Amu^2*Arho^2*Lambda*s^12 + 4*AG*Arho^3*Lambda*s^12 + 8*Amu*Arho^3*Lambda*s^12 - 12*AG*Amu^2*Arho^3*Lambda*s^12 + AG^2*Lambda^2*s^12 + 2*AG*Arho*Lambda^2*s^12 + Arho^2*Lambda^2*s^12 + 16*Amu^2*s^13 + 20*Amu*Arho*s^13 + 12*Amu^3*Arho*s^13 + 4*Arho^2*s^13 + 12*Amu^2*Arho^2*s^13 + 4*Amu*Arho^3*s^13 - 4*Amu^3*Arho^3*s^13 + 2*AG*Amu*Lambda*s^13 + 2*AG*Arho*Lambda*s^13 + 2*Amu*Arho*Lambda*s^13 + 2*Arho^2*Lambda*s^13 + Amu^2*s^14 + 2*Amu*Arho*s^14 + Arho^2*s^14) == 4*(1024*AG^4*Lambda^8 + 4096*AG^4*Lambda^7*s + 4096*AG^3*Amu*Lambda^7*s + 6144*AG^4*Lambda^6*s^2 + 16384*AG^3*Amu*Lambda^6*s^2 + 6144*AG^2*Amu^2*Lambda^6*s^2 + 2048*AG^4*Lambda^7*s^2 + 1024*AG^3*Arho*Lambda^7*s^2 + 4096*AG^4*Lambda^5*s^3 + 24576*AG^3*Amu*Lambda^5*s^3 + 24576*AG^2*Amu^2*Lambda^5*s^3 + 4096*AG*Amu^3*Lambda^5*s^3 + 6144*AG^4*Lambda^6*s^3 + 8192*AG^3*Amu*Lambda^6*s^3 + 4096*AG^3*Arho*Lambda^6*s^3 + 3072*AG^2*Amu*Arho*Lambda^6*s^3 + 1024*AG^4*Lambda^4*s^4 + 16384*AG^3*Amu*Lambda^4*s^4 + 36864*AG^2*Amu^2*Lambda^4*s^4 + 16384*AG*Amu^3*Lambda^4*s^4 + 1024*Amu^4*Lambda^4*s^4 + 6144*AG^4*Lambda^5*s^4 + 24576*AG^3*Amu*Lambda^5*s^4 + 12288*AG^2*Amu^2*Lambda^5*s^4 + 6144*AG^3*Arho*Lambda^5*s^4 + 12288*AG^2*Amu*Arho*Lambda^5*s^4 + 3072*AG*Amu^2*Arho*Lambda^5*s^4 + 64*AG^2*Lambda^6*s^4 + 1472*AG^4*Lambda^6*s^4 + 1536*AG^3*Arho*Lambda^6*s^4 + 384*AG^2*Arho^2*Lambda^6*s^4 - 128*AG^4*Arho^2*Lambda^6*s^4 + 4096*AG^3*Amu*Lambda^3*s^5 + 24576*AG^2*Amu^2*Lambda^3*s^5 + 24576*AG*Amu^3*Lambda^3*s^5 + 4096*Amu^4*Lambda^3*s^5 + 2048*AG^4*Lambda^4*s^5 + 24576*AG^3*Amu*Lambda^4*s^5 + 36864*AG^2*Amu^2*Lambda^4*s^5 + 8192*AG*Amu^3*Lambda^4*s^5 + 4096*AG^3*Arho*Lambda^4*s^5 + 18432*AG^2*Amu*Arho*Lambda^4*s^5 + 12288*AG*Amu^2*Arho*Lambda^4*s^5 + 1024*Amu^3*Arho*Lambda^4*s^5 + 256*AG^2*Lambda^5*s^5 + 2944*AG^4*Lambda^5*s^5 + 128*AG*Amu*Lambda^5*s^5 + 5888*AG^3*Amu*Lambda^5*s^5 + 4608*AG^3*Arho*Lambda^5*s^5 + 4608*AG^2*Amu*Arho*Lambda^5*s^5 + 1536*AG^2*Arho^2*Lambda^5*s^5 - 256*AG^4*Arho^2*Lambda^5*s^5 + 768*AG*Amu*Arho^2*Lambda^5*s^5 - 512*AG^3*Amu*Arho^2*Lambda^5*s^5 + 6144*AG^2*Amu^2*Lambda^2*s^6 + 16384*AG*Amu^3*Lambda^2*s^6 + 6144*Amu^4*Lambda^2*s^6 + 8192*AG^3*Amu*Lambda^3*s^6 + 36864*AG^2*Amu^2*Lambda^3*s^6 + 24576*AG*Amu^3*Lambda^3*s^6 + 2048*Amu^4*Lambda^3*s^6 + 1024*AG^3*Arho*Lambda^3*s^6 + 12288*AG^2*Amu*Arho*Lambda^3*s^6 + 18432*AG*Amu^2*Arho*Lambda^3*s^6 + 4096*Amu^3*Arho*Lambda^3*s^6 + 384*AG^2*Lambda^4*s^6 + 1472*AG^4*Lambda^4*s^6 + 512*AG*Amu*Lambda^4*s^6 + 11776*AG^3*Amu*Lambda^4*s^6 + 64*Amu^2*Lambda^4*s^6 + 8832*AG^2*Amu^2*Lambda^4*s^6 + 4608*AG^3*Arho*Lambda^4*s^6 + 13824*AG^2*Amu*Arho*Lambda^4*s^6 + 4608*AG*Amu^2*Arho*Lambda^4*s^6 + 2304*AG^2*Arho^2*Lambda^4*s^6 - 128*AG^4*Arho^2*Lambda^4*s^6 + 3072*AG*Amu*Arho^2*Lambda^4*s^6 - 1024*AG^3*Amu*Arho^2*Lambda^4*s^6 + 384*Amu^2*Arho^2*Lambda^4*s^6 - 768*AG^2*Amu^2*Arho^2*Lambda^4*s^6 + 128*AG^2*Lambda^5*s^6 + 448*AG^4*Lambda^5*s^6 + 32*AG*Arho*Lambda^5*s^6 + 800*AG^3*Arho*Lambda^5*s^6 + 384*AG^2*Arho^2*Lambda^5*s^6 - 128*AG^4*Arho^2*Lambda^5*s^6 + 64*AG*Arho^3*Lambda^5*s^6 - 64*AG^3*Arho^3*Lambda^5*s^6 + 4096*AG*Amu^3*Lambda*s^7 + 4096*Amu^4*Lambda*s^7 + 12288*AG^2*Amu^2*Lambda^2*s^7 + 24576*AG*Amu^3*Lambda^2*s^7 + 6144*Amu^4*Lambda^2*s^7 + 3072*AG^2*Amu*Arho*Lambda^2*s^7 + 12288*AG*Amu^2*Arho*Lambda^2*s^7 + 6144*Amu^3*Arho*Lambda^2*s^7 + 256*AG^2*Lambda^3*s^7 + 768*AG*Amu*Lambda^3*s^7 + 5888*AG^3*Amu*Lambda^3*s^7 + 256*Amu^2*Lambda^3*s^7 + 17664*AG^2*Amu^2*Lambda^3*s^7 + 5888*AG*Amu^3*Lambda^3*s^7 + 1536*AG^3*Arho*Lambda^3*s^7 + 13824*AG^2*Amu*Arho*Lambda^3*s^7 + 13824*AG*Amu^2*Arho*Lambda^3*s^7 + 1536*Amu^3*Arho*Lambda^3*s^7 + 1536*AG^2*Arho^2*Lambda^3*s^7 + 4608*AG*Amu*Arho^2*Lambda^3*s^7 - 512*AG^3*Amu*Arho^2*Lambda^3*s^7 + 1536*Amu^2*Arho^2*Lambda^3*s^7 - 1536*AG^2*Amu^2*Arho^2*Lambda^3*s^7 - 512*AG*Amu^3*Arho^2*Lambda^3*s^7 + 384*AG^2*Lambda^4*s^7 + 448*AG^4*Lambda^4*s^7 + 256*AG*Amu*Lambda^4*s^7 + 1792*AG^3*Amu*Lambda^4*s^7 + 128*AG*Arho*Lambda^4*s^7 + 1600*AG^3*Arho*Lambda^4*s^7 + 32*Amu*Arho*Lambda^4*s^7 + 2400*AG^2*Amu*Arho*Lambda^4*s^7 + 1152*AG^2*Arho^2*Lambda^4*s^7 - 128*AG^4*Arho^2*Lambda^4*s^7 + 768*AG*Amu*Arho^2*Lambda^4*s^7 - 512*AG^3*Amu*Arho^2*Lambda^4*s^7 + 256*AG*Arho^3*Lambda^4*s^7 - 128*AG^3*Arho^3*Lambda^4*s^7 + 64*Amu*Arho^3*Lambda^4*s^7 - 192*AG^2*Amu*Arho^3*Lambda^4*s^7 + 1024*Amu^4*s^8 + 8192*AG*Amu^3*Lambda*s^8 + 6144*Amu^4*Lambda*s^8 + 3072*AG*Amu^2*Arho*Lambda*s^8 + 4096*Amu^3*Arho*Lambda*s^8 + 64*AG^2*Lambda^2*s^8 + 512*AG*Amu*Lambda^2*s^8 + 384*Amu^2*Lambda^2*s^8 + 8832*AG^2*Amu^2*Lambda^2*s^8 + 11776*AG*Amu^3*Lambda^2*s^8 + 1472*Amu^4*Lambda^2*s^8 + 4608*AG^2*Amu*Arho*Lambda^2*s^8 + 13824*AG*Amu^2*Arho*Lambda^2*s^8 + 4608*Amu^3*Arho*Lambda^2*s^8 + 384*AG^2*Arho^2*Lambda^2*s^8 + 3072*AG*Amu*Arho^2*Lambda^2*s^8 + 2304*Amu^2*Arho^2*Lambda^2*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^8 - 1024*AG*Amu^3*Arho^2*Lambda^2*s^8 - 128*Amu^4*Arho^2*Lambda^2*s^8 + 384*AG^2*Lambda^3*s^8 + 768*AG*Amu*Lambda^3*s^8 + 1792*AG^3*Amu*Lambda^3*s^8 + 128*Amu^2*Lambda^3*s^8 + 2688*AG^2*Amu^2*Lambda^3*s^8 + 192*AG*Arho*Lambda^3*s^8 + 800*AG^3*Arho*Lambda^3*s^8 + 128*Amu*Arho*Lambda^3*s^8 + 4800*AG^2*Amu*Arho*Lambda^3*s^8 + 2400*AG*Amu^2*Arho*Lambda^3*s^8 + 1152*AG^2*Arho^2*Lambda^3*s^8 + 2304*AG*Amu*Arho^2*Lambda^3*s^8 - 512*AG^3*Amu*Arho^2*Lambda^3*s^8 + 384*Amu^2*Arho^2*Lambda^3*s^8 - 768*AG^2*Amu^2*Arho^2*Lambda^3*s^8 + 384*AG*Arho^3*Lambda^3*s^8 - 64*AG^3*Arho^3*Lambda^3*s^8 + 256*Amu*Arho^3*Lambda^3*s^8 - 384*AG^2*Amu*Arho^3*Lambda^3*s^8 - 192*AG*Amu^2*Arho^3*Lambda^3*s^8 + 80*AG^2*Lambda^4*s^8 + 48*AG^4*Lambda^4*s^8 + 48*AG*Arho*Lambda^4*s^8 + 176*AG^3*Arho*Lambda^4*s^8 + 4*Arho^2*Lambda^4*s^8 + 120*AG^2*Arho^2*Lambda^4*s^8 - 28*AG^4*Arho^2*Lambda^4*s^8 + 32*AG*Arho^3*Lambda^4*s^8 - 32*AG^3*Arho^3*Lambda^4*s^8 + 4*Arho^4*Lambda^4*s^8 - 8*AG^2*Arho^4*Lambda^4*s^8 + 4*AG^4*Arho^4*Lambda^4*s^8 + 2048*Amu^4*s^9 + 1024*Amu^3*Arho*s^9 + 128*AG*Amu*Lambda*s^9 + 256*Amu^2*Lambda*s^9 + 5888*AG*Amu^3*Lambda*s^9 + 2944*Amu^4*Lambda*s^9 + 4608*AG*Amu^2*Arho*Lambda*s^9 + 4608*Amu^3*Arho*Lambda*s^9 + 768*AG*Amu*Arho^2*Lambda*s^9 + 1536*Amu^2*Arho^2*Lambda*s^9 - 512*AG*Amu^3*Arho^2*Lambda*s^9 - 256*Amu^4*Arho^2*Lambda*s^9 + 128*AG^2*Lambda^2*s^9 + 768*AG*Amu*Lambda^2*s^9 + 384*Amu^2*Lambda^2*s^9 + 2688*AG^2*Amu^2*Lambda^2*s^9 + 1792*AG*Amu^3*Lambda^2*s^9 + 128*AG*Arho*Lambda^2*s^9 + 192*Amu*Arho*Lambda^2*s^9 + 2400*AG^2*Amu*Arho*Lambda^2*s^9 + 4800*AG*Amu^2*Arho*Lambda^2*s^9 + 800*Amu^3*Arho*Lambda^2*s^9 + 384*AG^2*Arho^2*Lambda^2*s^9 + 2304*AG*Amu*Arho^2*Lambda^2*s^9 + 1152*Amu^2*Arho^2*Lambda^2*s^9 - 768*AG^2*Amu^2*Arho^2*Lambda^2*s^9 - 512*AG*Amu^3*Arho^2*Lambda^2*s^9 + 256*AG*Arho^3*Lambda^2*s^9 + 384*Amu*Arho^3*Lambda^2*s^9 - 192*AG^2*Amu*Arho^3*Lambda^2*s^9 - 384*AG*Amu^2*Arho^3*Lambda^2*s^9 - 64*Amu^3*Arho^3*Lambda^2*s^9 + 160*AG^2*Lambda^3*s^9 + 160*AG*Amu*Lambda^3*s^9 + 192*AG^3*Amu*Lambda^3*s^9 + 144*AG*Arho*Lambda^3*s^9 + 176*AG^3*Arho*Lambda^3*s^9 + 48*Amu*Arho*Lambda^3*s^9 + 528*AG^2*Amu*Arho*Lambda^3*s^9 + 16*Arho^2*Lambda^3*s^9 + 240*AG^2*Arho^2*Lambda^3*s^9 + 240*AG*Amu*Arho^2*Lambda^3*s^9 - 112*AG^3*Amu*Arho^2*Lambda^3*s^9 + 96*AG*Arho^3*Lambda^3*s^9 - 32*AG^3*Arho^3*Lambda^3*s^9 + 32*Amu*Arho^3*Lambda^3*s^9 - 96*AG^2*Amu*Arho^3*Lambda^3*s^9 + 16*Arho^4*Lambda^3*s^9 - 16*AG^2*Arho^4*Lambda^3*s^9 - 16*AG*Amu*Arho^4*Lambda^3*s^9 + 16*AG^3*Amu*Arho^4*Lambda^3*s^9 + 64*Amu^2*s^10 + 1472*Amu^4*s^10 + 1536*Amu^3*Arho*s^10 + 384*Amu^2*Arho^2*s^10 - 128*Amu^4*Arho^2*s^10 + 256*AG*Amu*Lambda*s^10 + 384*Amu^2*Lambda*s^10 + 1792*AG*Amu^3*Lambda*s^10 + 448*Amu^4*Lambda*s^10 + 32*AG*Arho*Lambda*s^10 + 128*Amu*Arho*Lambda*s^10 + 2400*AG*Amu^2*Arho*Lambda*s^10 + 1600*Amu^3*Arho*Lambda*s^10 + 768*AG*Amu*Arho^2*Lambda*s^10 + 1152*Amu^2*Arho^2*Lambda*s^10 - 512*AG*Amu^3*Arho^2*Lambda*s^10 - 128*Amu^4*Arho^2*Lambda*s^10 + 64*AG*Arho^3*Lambda*s^10 + 256*Amu*Arho^3*Lambda*s^10 - 192*AG*Amu^2*Arho^3*Lambda*s^10 - 128*Amu^3*Arho^3*Lambda*s^10 + 80*AG^2*Lambda^2*s^10 + 320*AG*Amu*Lambda^2*s^10 + 80*Amu^2*Lambda^2*s^10 + 288*AG^2*Amu^2*Lambda^2*s^10 + 144*AG*Arho*Lambda^2*s^10 + 144*Amu*Arho*Lambda^2*s^10 + 528*AG^2*Amu*Arho*Lambda^2*s^10 + 528*AG*Amu^2*Arho*Lambda^2*s^10 + 24*Arho^2*Lambda^2*s^10 + 120*AG^2*Arho^2*Lambda^2*s^10 + 480*AG*Amu*Arho^2*Lambda^2*s^10 + 120*Amu^2*Arho^2*Lambda^2*s^10 - 168*AG^2*Amu^2*Arho^2*Lambda^2*s^10 + 96*AG*Arho^3*Lambda^2*s^10 + 96*Amu*Arho^3*Lambda^2*s^10 - 96*AG^2*Amu*Arho^3*Lambda^2*s^10 - 96*AG*Amu^2*Arho^3*Lambda^2*s^10 + 24*Arho^4*Lambda^2*s^10 - 8*AG^2*Arho^4*Lambda^2*s^10 - 32*AG*Amu*Arho^4*Lambda^2*s^10 - 8*Amu^2*Arho^4*Lambda^2*s^10 + 24*AG^2*Amu^2*Arho^4*Lambda^2*s^10 + 16*AG^2*Lambda^3*s^10 + 20*AG*Arho*Lambda^3*s^10 + 12*AG^3*Arho*Lambda^3*s^10 + 4*Arho^2*Lambda^3*s^10 + 12*AG^2*Arho^2*Lambda^3*s^10 + 4*AG*Arho^3*Lambda^3*s^10 - 4*AG^3*Arho^3*Lambda^3*s^10 + 128*Amu^2*s^11 + 448*Amu^4*s^11 + 32*Amu*Arho*s^11 + 800*Amu^3*Arho*s^11 + 384*Amu^2*Arho^2*s^11 - 128*Amu^4*Arho^2*s^11 + 64*Amu*Arho^3*s^11 - 64*Amu^3*Arho^3*s^11 + 160*AG*Amu*Lambda*s^11 + 160*Amu^2*Lambda*s^11 + 192*AG*Amu^3*Lambda*s^11 + 48*AG*Arho*Lambda*s^11 + 144*Amu*Arho*Lambda*s^11 + 528*AG*Amu^2*Arho*Lambda*s^11 + 176*Amu^3*Arho*Lambda*s^11 + 16*Arho^2*Lambda*s^11 + 240*AG*Amu*Arho^2*Lambda*s^11 + 240*Amu^2*Arho^2*Lambda*s^11 - 112*AG*Amu^3*Arho^2*Lambda*s^11 + 32*AG*Arho^3*Lambda*s^11 + 96*Amu*Arho^3*Lambda*s^11 - 96*AG*Amu^2*Arho^3*Lambda*s^11 - 32*Amu^3*Arho^3*Lambda*s^11 + 16*Arho^4*Lambda*s^11 - 16*AG*Amu*Arho^4*Lambda*s^11 - 16*Amu^2*Arho^4*Lambda*s^11 + 16*AG*Amu^3*Arho^4*Lambda*s^11 + 16*AG^2*Lambda^2*s^11 + 32*AG*Amu*Lambda^2*s^11 + 40*AG*Arho*Lambda^2*s^11 + 20*Amu*Arho*Lambda^2*s^11 + 36*AG^2*Amu*Arho*Lambda^2*s^11 + 12*Arho^2*Lambda^2*s^11 + 12*AG^2*Arho^2*Lambda^2*s^11 + 24*AG*Amu*Arho^2*Lambda^2*s^11 + 8*AG*Arho^3*Lambda^2*s^11 + 4*Amu*Arho^3*Lambda^2*s^11 - 12*AG^2*Amu*Arho^3*Lambda^2*s^11 + 80*Amu^2*s^12 + 48*Amu^4*s^12 + 48*Amu*Arho*s^12 + 176*Amu^3*Arho*s^12 + 4*Arho^2*s^12 + 120*Amu^2*Arho^2*s^12 - 28*Amu^4*Arho^2*s^12 + 32*Amu*Arho^3*s^12 - 32*Amu^3*Arho^3*s^12 + 4*Arho^4*s^12 - 8*Amu^2*Arho^4*s^12 + 4*Amu^4*Arho^4*s^12 + 32*AG*Amu*Lambda*s^12 + 16*Amu^2*Lambda*s^12 + 20*AG*Arho*Lambda*s^12 + 40*Amu*Arho*Lambda*s^12 + 36*AG*Amu^2*Arho*Lambda*s^12 + 12*Arho^2*Lambda*s^12 + 24*AG*Amu*Arho^2*Lambda*s^12 + 12*Amu^2*Arho^2*Lambda*s^12 + 4*AG*Arho^3*Lambda*s^12 + 8*Amu*Arho^3*Lambda*s^12 - 12*AG*Amu^2*Arho^3*Lambda*s^12 + AG^2*Lambda^2*s^12 + 2*AG*Arho*Lambda^2*s^12 + Arho^2*Lambda^2*s^12 + 16*Amu^2*s^13 + 20*Amu*Arho*s^13 + 12*Amu^3*Arho*s^13 + 4*Arho^2*s^13 + 12*Amu^2*Arho^2*s^13 + 4*Amu*Arho^3*s^13 - 4*Amu^3*Arho^3*s^13 + 2*AG*Amu*Lambda*s^13 + 2*AG*Arho*Lambda*s^13 + 2*Amu*Arho*Lambda*s^13 + 2*Arho^2*Lambda*s^13 + Amu^2*s^14 + 2*Amu*Arho*s^14 + Arho^2*s^14)], "Assumptions" -> Element[k | rhoT | muT | Arho | Amu | AG | Lambda, Reals] && Element[s, Complexes] && k > 0 && rhoT > 0 && muT > 0 && Lambda > 0 && -1 < Arho < 1 && -1 < Amu < 1 && -1 <= AG <= 1 && s != 0, "Details" -> "Only the nonzero constant factor 4 is divided out. E_+^2 E_-^2 is not divided out of P14."|>}
```

#### Runtime filtering workflow

1. Solve P14Polynomial[] == 0 to obtain polynomial candidate roots.
2. Reject candidates outside OmegaDefinition[].
3. Substitute each remaining candidate into DStarLambda[].
4. Retain only candidates satisfying the original dispersion relation within tolerance.
5. Apply physical admissibility criteria afterward.

Physical admissibility is intentionally outside this symbolic no-root-loss proof.
