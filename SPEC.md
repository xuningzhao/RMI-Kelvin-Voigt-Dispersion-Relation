# Specification: canonical nondimensionalization

## Physical problem

Two Kelvin–Voigt media meet at a perturbed interface. For medium $i$, let
$\rho_i$ be the solid density, $\mu_i$ the dynamic viscosity, and $G_i$
the shear modulus. Let $k$ be the interfacial wavenumber and $\gamma$ the
possibly complex temporal growth rate. Define

$$
\rho_T=\rho_1+\rho_2,\qquad
\mu_T=\mu_1+\mu_2,\qquad
G_T=G_1+G_2.
$$

The derivation assumes $k>0$, $\rho_i>0$, and, for the chosen canonical
clock, $\mu_T>0$. An individual modulus may vanish. The normalized
$g_i$ and $A_G$ coordinates assume $G_T>0$; the purely viscous point
$G_1=G_2=0$ is their continuous $\Lambda\to0$ boundary, where $A_G$ is
irrelevant and not identifiable. Whenever division by an individual
property is used, that property is additionally assumed positive. Square-root
branches are the same branches as in the dimensional dispersion relation;
nondimensionalization does not select or alter them.

## Original dimensional dispersion relation

An independent transcription of the supplied equation is

$$
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
$$

This is identical to the supplied relation: superscripted medium labels have
only been changed to subscripts, and the second denominator is the exact
$1\leftrightarrow2$ image of the first. No term, sign, factor, or radical has
changed.

The eight dimensional quantities and their dimensions are:

| Quantity | Meaning | Dimensions |
|---|---|---|
| $\gamma$ | response/growth rate | $T^{-1}$ |
| $k$ | wavenumber | $L^{-1}$ |
| $\rho_1,\rho_2$ | densities | $M L^{-3}$ |
| $\mu_1,\mu_2$ | dynamic viscosities | $M L^{-1}T^{-1}$ |
| $G_1,G_2$ | shear moduli | $M L^{-1}T^{-2}$ |

The two addends in $\mathcal D$ both have dimension $L/M$. In particular,
$[\gamma/\mu]=L/M$ and $[k^2/\rho]=L/M$.

## Buckingham $\Pi$ analysis

There are $n=8$ dimensional quantities. Their exponent vectors in the
$(M,L,T)$ basis form the matrix (columns are ordered as
$\gamma,k,\rho_1,\rho_2,\mu_1,\mu_2,G_1,G_2$)

$$
\mathbf D=
\begin{pmatrix}
0&0&1&1&1&1&1&1\\
0&-1&-3&-3&-1&-1&-1&-1\\
-1&0&0&0&-1&-1&-2&-2
\end{pmatrix}.
$$

For example, the columns
$[\rho]=(1,-3,0)$, $[\mu]=(1,-1,-1)$, and
$[k]=(0,-1,0)$ form a $3\times3$ matrix with determinant $-1$.
Thus $\operatorname{rank}\mathbf D\ge3$; since there are only three base
dimensions, $\operatorname{rank}\mathbf D=3$. Buckingham's theorem gives
$n-r=8-3=5$ independent groups. One is the dimensionless response, leaving
four independent control parameters.

A strict monomial basis, using $\rho_2,\mu_2,k$ as repeating variables, is

$$
\frac{\rho_2\gamma}{\mu_2k^2},\quad
\frac{\rho_1}{\rho_2},\quad
\frac{\mu_1}{\mu_2},\quad
\frac{\rho_2G_1}{\mu_2^2k^2},\quad
\frac{\rho_2G_2}{\mu_2^2k^2}.
$$

It is valid but privileges medium 2. A symmetric and physically clearer set of
coordinates is

$$
s_v=\frac{\rho_T\gamma}{\mu_Tk^2},\qquad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2},\qquad
A_\rho,\ A_\mu,\ A_G,
$$

where

$$
A_x=\frac{x_1-x_2}{x_1+x_2},\qquad
x_1/x_T=\frac{1+A_x}{2},\quad
x_2/x_T=\frac{1-A_x}{2}.
$$

These five coordinates are invertibly related to the monomial basis in the
positive-property interior, so they are independent. Because sums and bounded
contrasts are not monomials in the original variables, this is best called a
natural symmetric $\Pi$-coordinate system rather than *the* unique
Buckingham basis.

Other equally complete bases replace $(s_v,\Lambda)$ by

$$
(s_e,\chi),\qquad (s_R,\Lambda),\qquad
\chi=\frac{r_e}{r_v}=\sqrt{\Lambda},
$$

with the three contrasts unchanged.

## Characteristic time scales

The total-property rates and times are defined exactly by

$$
r_v=\frac{\mu_Tk^2}{\rho_T},\quad t_v=r_v^{-1},\qquad
r_e=k\sqrt{\frac{G_T}{\rho_T}},\quad t_e=r_e^{-1},
$$

and the Kelvin–Voigt crossover or relaxation time is

$$
t_R=\frac{\mu_T}{G_T}.
$$

They are not three independent clocks:

$$
t_R=\frac{t_e^2}{t_v},\qquad
\chi=\frac{r_e}{r_v}=\frac{t_v}{t_e},\qquad
\Lambda=\chi^2.
$$

The proposed bounded rate contrast naturally accompanies the rate-sum clock,

$$
E_k=\frac{r_v-r_e}{r_v+r_e}=\frac{1-\chi}{1+\chi},\qquad
t_+=\frac{1}{r_v+r_e}.
$$

For nonnegative material properties, $\chi\in[0,\infty)$ maps to
$E_k\in(-1,1]$; the lower endpoint is approached as $\chi\to\infty$.

Constituent clocks also occur:

$$
t_{v,i}=\frac{\rho_i}{\mu_i k^2},\qquad
t_{e,i}=\frac1k\sqrt{\frac{\rho_i}{G_i}},\qquad
t_{R,i}=\frac{\mu_i}{G_i}.
$$

They introduce no new dimensional freedom: their ratios to total clocks are
functions of $\Lambda$ and the contrasts. No capillary, gravitational,
acoustic, or imposed shock time appears because the corresponding dimensional
quantity is absent from the given dispersion relation. Infinitely many clocks
can be manufactured as $t_v f(\Lambda,A_\rho,A_\mu,A_G)$; they are
reparameterizations, not new balances.

## Candidate nondimensionalizations

Set

$$
r_i=\rho_i/\rho_T,\qquad m_i=\mu_i/\mu_T,\qquad
g_i=G_i/G_T,
$$

so each pair sums to one. For any clock $t_c$, define

$$
s=\gamma t_c,\quad p=\frac{t_v}{t_c},\quad
q=\frac{G_Tt_c}{\mu_T},\quad H_i=m_i+\frac{qg_i}{s}.
$$

Direct substitution and multiplication by $\rho_T/k^2$ gives the universal
dimensionless relation

$$
p s\left[
\frac1{H_1+H_2\sqrt{1+p r_2s/H_2}}+
\frac1{H_2+H_1\sqrt{1+p r_1s/H_1}}
\right]+4=0,                                      \tag{1}
$$

with $pq=\Lambda$. Equivalently, for $s\ne0$, let
$a_i=m_is+qg_i$. Then

$$
p s^2\left[
\frac1{a_1+a_2\sqrt{1+p r_2s^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+p r_1s^2/a_1}}
\right]+4=0.                                      \tag{2}
$$

The un-cleared form (1) is the direct identity; (2) is convenient away from
the original singular representation at $\gamma=0$.

The choices are:

| Clock | $s$ | $p$ | $q$ |
|---|---:|---:|---:|
| viscous $t_v$ | $s_v$ | $1$ | $\Lambda$ |
| elastic $t_e$ | $s_e$ | $\chi$ | $\chi$ |
| relaxation $t_R$ | $s_R$ | $\Lambda$ | $1$ |
| rate sum $t_+$ | $s_+$ | $1+\chi$ | $\chi^2/(1+\chi)$ |

Thus (2) yields, respectively (the sums below contain the two ordered pairs
$(i,j)=(1,2),(2,1)$),

$$
s_v^2\sum_{i\ne j}\frac1{a_i+a_j\sqrt{1+r_js_v^2/a_j}}+4=0,
\quad a_i=m_is_v+\Lambda g_i,                     \tag{V}
$$

$$
\chi s_e^2\sum_{i\ne j}\frac1{a_i+a_j\sqrt{1+\chi r_js_e^2/a_j}}+4=0,
\quad a_i=m_is_e+\chi g_i,                       \tag{E}
$$

$$
\Lambda s_R^2\sum_{i\ne j}\frac1{a_i+a_j\sqrt{1+\Lambda r_js_R^2/a_j}}+4=0,
\quad a_i=m_is_R+g_i.                            \tag{R}
$$

For the rate-sum form, substitute its $p,q$ from the table into (2). In
terms of $E_k$,

$$
p=\frac{2}{1+E_k},\qquad
q=\frac{(1-E_k)^2}{2(1+E_k)}.                    \tag{+}
$$

Writing $p_+=1+\chi$, $q_+=\chi^2/(1+\chi)$, and
$a_i=m_is_++q_+g_i$, its equation is explicitly

$$
p_+s_+^2\left[
\frac1{a_1+a_2\sqrt{1+p_+r_2s_+^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+p_+r_1s_+^2/a_1}}
\right]+4=0.                                    \tag{S}
$$

For a constituent clock, no fresh derivation is needed because (1)–(2) are
general. The relevant substitutions are

$$
\begin{array}{c|cc}
t_c&p&q\\ \hline
t_{v,j}&m_j/r_j&\Lambda r_j/m_j\\
t_{e,j}&\chi\sqrt{g_j/r_j}&\chi\sqrt{r_j/g_j}\\
t_{R,j}&\Lambda g_j/m_j&m_j/g_j.
\end{array}
$$

The response variables obey

$$
s_e=\frac{s_v}{\chi},\qquad
s_R=\frac{s_v}{\Lambda},\qquad
s_+=\frac{s_v}{1+\chi}.
$$

## Comparison

| Form | Algebra | Symmetry | Controls | Interpretation and limits |
|---|---|---|---:|---|
| viscous | simplest outer and radical coefficients | exchange-symmetric | 4 | direct diffusion clock; regular as $G_T\to0$ |
| elastic | balanced $p=q=\chi$ | exchange-symmetric | 4 | direct wave clock; undefined at $G_T=0$ |
| relaxation | normalizes elastic/viscous crossover | exchange-symmetric | 4 | useful for constitutive relaxation; undefined at $G_T=0$, poor long-wave scaling |
| rate sum | two coefficient functions | exchange-symmetric | 4 | bounded $E_k$, but no primitive balance and more algebra |
| constituent | extra contrast factors | privileges one medium | 4 | useful only when one medium is a reference; singular if its selected property vanishes |

All total-property formulations have the same number of independent
parameters. A change of clock cannot reduce the Buckingham count; it merely
moves powers of $\chi$ or $\Lambda$ between the response and coefficients.
All are invariant under medium exchange when
$(r_1,m_1,g_1)\leftrightarrow(r_2,m_2,g_2)$, equivalently when all three
contrasts change sign.

## Status of $A_\rho,A_\mu,A_G,E_k$

The four quantities are a valid, independent set of **control coordinates**
for positive properties. Adding one response such as $s_v$ completes the
five groups required by Buckingham's theorem. They are not a unique or strict
monomial $\Pi$ basis:

1. Buckingham bases are nonunique.
2. Each $A_x$ is a bounded transform of the ratio $x_1/x_2$.
3. $E_k$ is a bounded transform of $\chi=\sqrt\Lambda$:
   $\chi=(1-E_k)/(1+E_k)$.
4. Choosing $E_k$ does not itself specify the dimensionless response clock.

The contrasts arise naturally from exchange symmetry and total-property
normalization, not from dimensional analysis alone. $E_k$ is convenient for
bounded plots, while $\Lambda$ is more canonical algebraically because it is
the direct dimensionless coefficient generated by the viscous scaling.

## Recommended formulation

Adopt (V) with

$$
\boxed{
s=\frac{\rho_T\gamma}{\mu_Tk^2},\quad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2},\quad
(A_\rho,A_\mu,A_G)
}
$$

and reconstruct fractions with $(1\pm A_x)/2$. This formulation is
dimensionally complete, exchange-symmetric, minimal in parameter count,
algebraically shortest, and continuous into the purely viscous limit. At the
exact point $G_T=0$, use $\Lambda=0$; $A_G$ can be omitted because it
has no physical effect there. Store
or compute $\chi=\sqrt\Lambda$ and $E_k=(1-\chi)/(1+\chi)$ only as derived
coordinates. If a bounded sweep coordinate is desirable, sweeping $E_k$ is
perfectly legitimate, but the governing specification should remain in
$\Lambda$.

The elastic form (E) is an equally valid secondary presentation when wave
physics is the focus. It is not recommended as the repository's defining form
because it loses its clock at zero total modulus.

# Phase II: Parameter-space stability analysis

> **Current implementation note.** This section records the original numerical
> exploration plan written before the symbolic polynomial-completeness theorem
> was finalized. The production workflow is now the polynomial-candidate method
> defined in [root_verification_workflow.md](root_verification_workflow.md):
> solve the verified $P_{14}$ polynomial, filter candidates by the original
> nondimensional dispersion relation, and then apply the $q_\pm$ spatial-decay
> admissibility conditions. Finite-budget direct searches remain useful only as
> exploratory diagnostics.

## Scope and scientific objective

Phase II begins from the final nondimensional dispersion relation documented in
[derivation/final_nondimensional_dispersion.md](derivation/final_nondimensional_dispersion.md).
The Phase I theory and its recommended viscous-clock nondimensionalization
remain the mathematical foundation of the project.

The primary objective is to determine how stability is organized over the
nondimensional material-parameter space. The canonical coordinates are

$$
s,\qquad A_\rho,\qquad A_\mu,\qquad A_G,\qquad \Lambda,
$$

and the bounded sweep-coordinate representation is

$$
s,\qquad A_\rho,\qquad A_\mu,\qquad A_G,\qquad E_k,
$$

with

$$
\Lambda=\left(\frac{1-E_k}{1+E_k}\right)^2.
$$

The ultimate scientific objective is not merely to compute roots. It is to
understand how density contrast, viscosity contrast, elastic contrast, and
elastic–viscous rate competition organize the RMI stability behavior.

## Numerical problem

For each prescribed parameter set, solve

$$
\mathcal D_\Lambda^*(s;A_\rho,A_\mu,A_G,\Lambda)=0
$$

or, equivalently,

$$
\mathcal D_{E_k}^*(s;A_\rho,A_\mu,A_G,E_k)=0,
$$

where $s\in\mathbb C$. Within each declared search domain, identify all
verified roots when possible. If the verified roots are
$\{s_i\}_{i=1}^N$, identify the dominant root by

$$
s_{\mathrm{dom}}
\in\underset{s_i}{\operatorname{arg\,max}}\ \operatorname{Re}(s_i).
$$

If several roots share the same largest real part within the stated numerical
tolerance, report the complete tied set rather than silently selecting one.
Classify stability from the dominant real part:

$$
\begin{array}{ll}
\operatorname{Re}(s_{\mathrm{dom}})>0
&\text{unstable},\\
\operatorname{Re}(s_{\mathrm{dom}})<0
&\text{stable or decaying},\\
\operatorname{Re}(s_{\mathrm{dom}})=0
&\text{neutral or marginal}.
\end{array}
$$

Numerically, these comparisons must use a documented stability tolerance.
Cases lying within that tolerance of zero must be labeled marginal rather than
assigned a sign that is not resolved by the computation.

## Required outputs

Each parameter study must produce, where applicable:

- all verified roots within the declared search domain;
- the dominant root $s_{\mathrm{dom}}$;
- the dominant growth rate $\operatorname{Re}(s_{\mathrm{dom}})$;
- the dominant frequency $\operatorname{Im}(s_{\mathrm{dom}})$;
- root trajectories in the complex $s$-plane under parameter continuation;
- stability maps over selected nondimensional parameter planes; and
- diagnostic reports for failed or unresolved cases.

Every reported result must retain the parameter values, search domain,
square-root branch convention, solver settings, verification tolerances, and
classification status needed to reproduce it.

## Root-finding outcome categories

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

## Numerical principles

The Phase II implementation must follow these principles:

1. The final nondimensional dispersion relation is the source of truth.
2. Every candidate root must be verified by evaluating the residual
   $|\mathcal D^*(s)|$ in the source equation. The absolute and relative
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
   of the complex $s$-plane; it must not be claimed globally without a
   supporting mathematical argument.

## Recommended first numerical implementation

The first numerical implementation should:

1. implement D_star_Lambda(s, Arho, Amu, AG, Lambda) directly from
   $\mathcal D_\Lambda^*$;
2. implement D_star_Ek(s, Arho, Amu, AG, Ek) directly from
   $\mathcal D_{E_k}^*$;
3. verify over representative admissible parameters and complex values of $s$
   that

   $$
   \mathrm{D\_star\_Ek}(s,A_\rho,A_\mu,A_G,E_k)
   =
   \mathrm{D\_star\_Lambda}\left(
   s,A_\rho,A_\mu,A_G,
   \left(\frac{1-E_k}{1+E_k}\right)^2
   \right);
   $$

4. keep compact helper variables internal to the implementation, using them
   only to evaluate the final mathematical definition rather than replacing
   it; and
5. test medium-exchange symmetry, limiting cases, branch consistency, and
   residual evaluation before introducing parameter sweeps.

This section specifies the future implementation only. No root solver is
introduced in this phase of the documentation update.

## Initial parameter-exploration plan

Begin with one-dimensional sweeps in $E_k$ at fixed
$A_\rho$, $A_\mu$, and $A_G$. These sweeps should establish reliable root
seeding, branch tracking, dominant-root switching behavior, and diagnostic
reporting before moving to higher-dimensional studies.

After the one-dimensional workflow is verified, construct two-dimensional
maps including:

- $E_k$ versus $A_\rho$;
- $E_k$ versus $A_\mu$;
- $E_k$ versus $A_G$; and
- $A_\rho$ versus $A_G$ at fixed $E_k$.

For every sweep or map, record the fixed parameters, sampled ranges,
continuation direction, search domain, root classification, and unresolved
regions. Stability boundaries should be inferred only from verified dominant
roots and should retain uncertainty information where the numerical
classification changes or remains unresolved.
