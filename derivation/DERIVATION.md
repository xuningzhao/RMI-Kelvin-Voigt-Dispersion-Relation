# Detailed derivation

## 1. Verified starting equation

Write the Kelvin–Voigt effective viscosity as

$$
\eta_i(\gamma)=\mu_i+\frac{G_i}{\gamma}.
$$

The supplied equation is exactly

$$
\gamma\left[
\frac1{\eta_1+\eta_2\sqrt{1+\rho_2\gamma/(\eta_2k^2)}}+
\frac1{\eta_2+\eta_1\sqrt{1+\rho_1\gamma/(\eta_1k^2)}}
\right]+\frac{4k^2}{\rho_T}=0.                  \tag{D}
$$

Expanding each $\eta_i$ reproduces the original equation term by term.

## 2. General change of variables

For an arbitrary nonzero characteristic time $t_c$, introduce

$$
\gamma=\frac{s}{t_c},\quad
\rho_i=\rho_T r_i,\quad
\mu_i=\mu_Tm_i,\quad
G_i=G_Tg_i.
$$

Define

$$
p=\frac{\rho_T}{\mu_Tk^2t_c}=\frac{t_v}{t_c},
\qquad q=\frac{G_Tt_c}{\mu_T}.
$$

Then, without omitting an algebraic step,

$$
\begin{aligned}
\eta_i
&=\mu_Tm_i+\frac{G_Tg_i}{s/t_c}\\
&=\mu_T\left(m_i+\frac{G_Tt_c}{\mu_T}\frac{g_i}{s}\right)\\
&=\mu_TH_i,
\qquad H_i=m_i+\frac{qg_i}{s},
\end{aligned}
$$

and

$$
\begin{aligned}
1+\frac{\rho_i\gamma}{\eta_i k^2}
&=1+\frac{\rho_T r_i(s/t_c)}{\mu_TH_i k^2}\\
&=1+p\frac{r_is}{H_i}.
\end{aligned}
$$

The first reciprocal denominator in (D) becomes

$$
\frac1{\eta_1+\eta_2\sqrt{1+\rho_2\gamma/(\eta_2k^2)}}
=\frac1{\mu_T}
\frac1{H_1+H_2\sqrt{1+pr_2s/H_2}}.
$$

The other term follows by exchanging 1 and 2. Finally multiply (D) by
$\rho_T/k^2$. Since

$$
\frac{\rho_T}{k^2}\frac{\gamma}{\mu_T}
=\frac{\rho_Ts}{\mu_Tk^2t_c}=ps,
$$

the exact general dimensionless equation is

$$
ps\left[
\frac1{H_1+H_2\sqrt{1+pr_2s/H_2}}+
\frac1{H_2+H_1\sqrt{1+pr_1s/H_1}}
\right]+4=0.                                    \tag{G1}
$$

Also,

$$
pq=\frac{\rho_TG_T}{\mu_T^2k^2}=\Lambda,
$$

so only one independent rate-ratio parameter is present.

For a form with no explicit $1/s$, define

$$
a_i=sH_i=m_is+qg_i.
$$

For $s\ne0$,

$$
H_i=\frac{a_i}{s},\qquad
1+\frac{pr_is}{H_i}=1+\frac{pr_is^2}{a_i},
$$

and

$$
\frac1{H_i+H_j\sqrt{1+pr_js/H_j}}
=\frac{s}{a_i+a_j\sqrt{1+pr_js^2/a_j}}.
$$

Thus

$$
ps^2\left[
\frac1{a_1+a_2\sqrt{1+pr_2s^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+pr_1s^2/a_1}}
\right]+4=0.                                    \tag{G2}
$$

## 3. Viscous clock

Take

$$
t_c=t_v=\frac{\rho_T}{\mu_Tk^2}.
$$

Then

$$
p=\frac{t_v}{t_v}=1,qquad
q=\frac{G_T}{\mu_T}\frac{\rho_T}{\mu_Tk^2}
=\Lambda.
$$

In (G1),

$$
H_i=m_i+\frac{\Lambda g_i}{s_v},
$$

giving

$$
s_v\left[
\frac1{H_1+H_2\sqrt{1+r_2s_v/H_2}}+
\frac1{H_2+H_1\sqrt{1+r_1s_v/H_1}}
\right]+4=0.
$$

Equivalently, with $a_i=m_is_v+\Lambda g_i$,

$$
s_v^2\left[
\frac1{a_1+a_2\sqrt{1+r_2s_v^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+r_1s_v^2/a_1}}
\right]+4=0.                                    \tag{V}
$$

## 4. Elastic-wave clock

Take

$$
t_c=t_e=\frac1k\sqrt{\frac{\rho_T}{G_T}},\qquad
\chi=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}=\sqrt\Lambda.
$$

Direct calculation gives

$$
p=\frac{\rho_T}{\mu_Tk^2t_e}
=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}=\chi,
$$

and

$$
q=\frac{G_Tt_e}{\mu_T}
=\frac{\sqrt{\rho_TG_T}}{\mu_Tk}=\chi.
$$

Consequently $a_i=m_is_e+\chi g_i$, and (G2) becomes

$$
\chi s_e^2\left[
\frac1{a_1+a_2\sqrt{1+\chi r_2s_e^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+\chi r_1s_e^2/a_1}}
\right]+4=0.                                    \tag{E}
$$

## 5. Kelvin–Voigt relaxation clock

Take $t_c=t_R=\mu_T/G_T$. Then

$$
p=\frac{\rho_T}{\mu_Tk^2}\frac{G_T}{\mu_T}
=\Lambda,
\qquad q=\frac{G_T}{\mu_T}\frac{\mu_T}{G_T}=1.
$$

Thus $a_i=m_is_R+g_i$, and

$$
\Lambda s_R^2\left[
\frac1{a_1+a_2\sqrt{1+\Lambda r_2s_R^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+\Lambda r_1s_R^2/a_1}}
\right]+4=0.                                    \tag{R}
$$

The identity $t_R=t_e^2/t_v$ follows immediately from the definitions.

## 6. Rate-sum clock and $E_k$

Let

$$
t_+=\frac1{r_v+r_e},\quad \chi=\frac{r_e}{r_v},\quad
E_k=\frac{1-\chi}{1+\chi}.
$$

Because $t_+=t_v/(1+\chi)$,

$$
p=\frac{t_v}{t_+}=1+\chi,qquad
q=\frac{\Lambda}{p}=\frac{\chi^2}{1+\chi}.
$$

Solving the definition of $E_k$ gives

$$
\chi=\frac{1-E_k}{1+E_k},quad
1+\chi=\frac2{1+E_k},quad
\frac{\chi^2}{1+\chi}=\frac{(1-E_k)^2}{2(1+E_k)}.
$$

Substitution of these $p,q$ and
$a_i=m_is_++qg_i$ into (G2) gives explicitly

$$
(1+\chi)s_+^2\left[
\frac1{a_1+a_2\sqrt{1+(1+\chi)r_2s_+^2/a_2}}+
\frac1{a_2+a_1\sqrt{1+(1+\chi)r_1s_+^2/a_1}}
\right]+4=0,
$$

where $a_i=m_is_++\chi^2g_i/(1+\chi)$. It has no fewer independent
controls than (V); it only bounds the rate ratio.

## 7. Constituent clocks

From $r_i=\rho_i/\rho_T$, $m_i=\mu_i/\mu_T$, and $g_i=G_i/G_T$,

$$
t_{v,j}=t_v\frac{r_j}{m_j},\qquad
t_{e,j}=t_e\sqrt{\frac{r_j}{g_j}},\qquad
t_{R,j}=t_R\frac{m_j}{g_j}.
$$

Using $p=t_v/t_c$, $q=\Lambda/p$, these yield

$$
\begin{array}{c|cc}
t_c&p&q\\ \hline
t_{v,j}&m_j/r_j&\Lambda r_j/m_j\\
t_{e,j}&\chi\sqrt{g_j/r_j}&\chi\sqrt{r_j/g_j}\\
t_{R,j}&\Lambda g_j/m_j&m_j/g_j.
\end{array}
$$

Putting any row into (G1) or (G2) gives its full dimensionless relation.
These clocks privilege medium $j$ and add contrast factors without reducing
the number of independent controls.

## 8. Contrast reconstruction and symmetry

For $x\in\{\rho,\mu,G\}$,

$$
x_1=\frac{x_T}{2}(1+A_x),\qquad
x_2=\frac{x_T}{2}(1-A_x).
$$

Therefore

$$
r_{1,2}=\frac{1\pm A_\rho}{2},\quad
m_{1,2}=\frac{1\pm A_\mu}{2},\quad
g_{1,2}=\frac{1\pm A_G}{2}.
$$

Exchanging the media sends all three contrasts to their negatives and swaps
the two reciprocal terms, leaving every total-property equation invariant.
