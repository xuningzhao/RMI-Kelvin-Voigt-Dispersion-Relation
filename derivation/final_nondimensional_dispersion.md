# Final explicit nondimensional dispersion relations

## Purpose

This document expands the recommended viscous-clock nondimensional dispersion
relation entirely in the variables

$$
s,\quad A_\rho,\quad A_\mu,\quad A_G,\quad \Lambda
$$

and then replaces $\Lambda$ by the bounded rate coordinate $E_k$. The physics
and the viscous-clock scaling are unchanged. No numerical root finding is
performed.

The total properties and nondimensional response are

$$
\rho_T=\rho_1+\rho_2,\qquad
\mu_T=\mu_1+\mu_2,\qquad
G_T=G_1+G_2,\qquad
s=\frac{\rho_T\gamma}{\mu_Tk^2}.
$$

The material contrasts and elastic–viscous parameter are

$$
A_\rho=\frac{\rho_1-\rho_2}{\rho_T},\qquad
A_\mu=\frac{\mu_1-\mu_2}{\mu_T},\qquad
A_G=\frac{G_1-G_2}{G_T},\qquad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2}.
$$

## 1. Starting viscous-clock equation

The recommended viscous-clock equation derived in the existing specification
can be written without its former denominator shorthand as

$$
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
$$

where, only for the purpose of identifying the starting equation,

$$
r_i=\frac{\rho_i}{\rho_T},\qquad
m_i=\frac{\mu_i}{\mu_T},\qquad
g_i=\frac{G_i}{G_T}.
$$

Equation (1) is the cleared version of the direct viscous-clock equation. The
fully explicit final equations below contain none of these normalized-fraction
symbols.

## 2. Substitution of the contrast reconstructions

Use

$$
r_1=\frac{1+A_\rho}{2},\qquad
r_2=\frac{1-A_\rho}{2},
$$

$$
m_1=\frac{1+A_\mu}{2},\qquad
m_2=\frac{1-A_\mu}{2},
$$

and

$$
g_1=\frac{1+A_G}{2},\qquad
g_2=\frac{1-A_G}{2}.
$$

The two distinct material combinations occurring in (1) become

$$
m_1s+\Lambda g_1
=\frac{s(1+A_\mu)+\Lambda(1+A_G)}{2},
$$

$$
m_2s+\Lambda g_2
=\frac{s(1-A_\mu)+\Lambda(1-A_G)}{2},             \tag{2}
$$

and the radical arguments become

$$
1+\frac{r_2s^2}{m_2s+\Lambda g_2}
=1+\frac{(1-A_\rho)s^2}
{s(1-A_\mu)+\Lambda(1-A_G)},
$$

$$
1+\frac{r_1s^2}{m_1s+\Lambda g_1}
=1+\frac{(1+A_\rho)s^2}
{s(1+A_\mu)+\Lambda(1+A_G)}.                     \tag{3}
$$

Each complete reciprocal denominator in (1) therefore has an overall factor
$1/2$. Removing that factor from both reciprocals transforms (1) into

$$
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
$$

Dividing (4) by two gives the equally valid cleared normalization

$$
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
$$

For $s\ne0$, factoring one $s$ out of every denominator in (5) gives the
direct-viscosity form used for the final equations. This step changes the
outer factor from $s^2$ to $s$ but changes no radical argument.

## 3. Final explicit $\Lambda$ form

The completely expanded dispersion function in the algebraically canonical
parameter $\Lambda$ is

$$
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
$$

Equation (6) contains only $s$, $A_\rho$, $A_\mu$, $A_G$, and $\Lambda$.
The first reciprocal contains the inertia of medium 2 in its radical, hence
the factors $1-A_\rho$, $1-A_\mu$, and $1-A_G$ there. The second reciprocal is
its exact $1\leftrightarrow2$ counterpart. This provides a direct sign check.

## 4. From $\Lambda$ to $E_k$

The viscous and elastic response rates are

$$
r_v=\frac{\mu_Tk^2}{\rho_T},\qquad
r_e=k\sqrt{\frac{G_T}{\rho_T}},
$$

so

$$
\frac{r_e}{r_v}
=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}
=\sqrt{\Lambda}.
$$

By definition,

$$
E_k=\frac{r_v-r_e}{r_v+r_e}
=\frac{1-\sqrt{\Lambda}}{1+\sqrt{\Lambda}}.
$$

Solving for the rate ratio gives

$$
\sqrt{\Lambda}=\frac{1-E_k}{1+E_k},\qquad
\Lambda=\left(\frac{1-E_k}{1+E_k}\right)^2.       \tag{7}
$$

Thus the $E_k$ equation is obtained from (6) by making only the substitution
(7). The response variable $s$ remains the viscous-clock response; changing
from $\Lambda$ to $E_k$ does not change the clock.

## 5. Final explicit $E_k$ form

Substitution of (7) into every occurrence of $\Lambda$ in (6) gives

$$
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
$$

There is no sign or branch change between (6) and (8). In particular,

$$
\mathcal D_{E_k}^*(s;A_\rho,A_\mu,A_G,E_k)
=\mathcal D_\Lambda^*\!\left(
s;A_\rho,A_\mu,A_G,
\left(\frac{1-E_k}{1+E_k}\right)^2
\right).                                                        \tag{9}
$$

## 6. Constant-prefactor and normalization audit

The constant can be traced without relying on the cleared equation. Starting
from the dimensional relation, multiply by $\rho_T/k^2$. Under the viscous
clock

$$
s=\frac{\rho_T\gamma}{\mu_Tk^2},
$$

the dimensional constant becomes exactly

$$
\frac{\rho_T}{k^2}\frac{4k^2}{\rho_T}=4.          \tag{10}
$$

Before contrast substitution, the direct dimensionless equation consequently
has outer factor $s$ and constant $4$. Each normalized material contribution
contains a factor $1/2$ after the contrast reconstructions. Therefore each
complete reciprocal is multiplied by $2$, and the equation is

$$
2s[\text{two fully expanded reciprocal terms}]+4=0.             \tag{11}
$$

Dividing the entire equation by the nonzero constant $2$ gives

$$
s[\text{the same two fully expanded reciprocal terms}]+2=0,     \tag{12}
$$

which is the normalization used in the boxed equations (6) and (8). Hence
$2s[\cdots]+4=0$ and $s[\cdots]+2=0$ are exactly equivalent; neither represents
a change in physics. The cleared normalizations (4) and (5) are likewise
equivalent for $s\ne0$. The original Kelvin–Voigt representation already
contains $G_i/\gamma$, so $s=0$ must be understood through the appropriate
limit rather than by treating the uncleared formula as an ordinary value.

## 7. Compact implementation form

Only after establishing the fully explicit equations, it is convenient for an
implementation to define

$$
\eta_1=(1+A_\mu)+\frac{\Lambda}{s}(1+A_G),\qquad
\eta_2=(1-A_\mu)+\frac{\Lambda}{s}(1-A_G),
$$

$$
q_1=\sqrt{1+\frac{(1+A_\rho)s}{\eta_1}},\qquad
q_2=\sqrt{1+\frac{(1-A_\rho)s}{\eta_2}}.
$$

Then the $\Lambda$ form can be evaluated as

$$
\mathcal D_\Lambda^*
=s\left(\frac{1}{\eta_1+\eta_2q_2}
+\frac{1}{\eta_2+\eta_1q_1}\right)+2.             \tag{13}
$$

For the $E_k$ form, use the same implementation after assigning

$$
\Lambda=\left(\frac{1-E_k}{1+E_k}\right)^2.
$$

The square-root branch convention must remain the same as in the dimensional
dispersion relation.

## 8. Comparison of the two forms

The two equations contain the same four independent controls and the same
dimensionless response $s$.

- $\Lambda=\rho_TG_T/(\mu_T^2k^2)$ is algebraically canonical for the viscous
  clock. It appears directly when dimensional factors are collected, and it
  keeps the equation rational in the material parameters outside the square
  roots.
- $E_k=(r_v-r_e)/(r_v+r_e)$ is a bounded reparameterization of the positive
  rate ratio. For nonnegative material properties, $-1<E_k\leq1$, with
  $E_k=1$ at $\Lambda=0$ and $E_k\to-1$ as $\Lambda\to\infty$.
- Equations (6) and (8) are mathematically equivalent by the invertible map
  (7) on the finite positive-parameter interior. Choosing one or the other
  changes only the coordinate used for parameter studies, not the dispersion
  relation or its solutions.
