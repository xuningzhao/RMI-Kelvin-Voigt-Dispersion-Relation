# Mathematical structure of the Buckingham $\Pi$ space

## Scope and relation to the existing derivation

This document studies the linear-algebraic content of Buckingham's theorem for
the dimensional variables

$$
(\gamma,k,\rho_1,\rho_2,\mu_1,\mu_2,G_1,G_2).
$$

It does not change the dispersion relation or the recommended
nondimensionalization in [../SPEC.md](../SPEC.md). The detailed substitution
into the dispersion relation remains in [DERIVATION.md](DERIVATION.md). No
repeating variables are used to obtain the null space below; repeating
variables are introduced only afterward to explain the classical construction.

For the geometric interpretation, the dimensional variables are first viewed
on a positive coordinate chart so that arbitrary real powers and logarithms
are well defined. The exponent algebra itself is formal and remains valid for
a nonzero complex response $\gamma$; the dispersion relation supplies that
analytic continuation. Zero material-property values can be included afterward
by taking appropriate limits.

## 1. Dimensional matrix

Use base-dimension order $(M,L,T)$ and variable order

$$
\boldsymbol{x}=(\gamma,k,\rho_1,\rho_2,\mu_1,\mu_2,G_1,G_2).
$$

Every column is the exponent vector of one variable:

$$
\begin{array}{c|c|c}
\text{variable}&\text{physical dimension}&\text{column of }D\\ \hline
\gamma&T^{-1}&(0,0,-1)^\mathsf T\\
k&L^{-1}&(0,-1,0)^\mathsf T\\
\rho_1&ML^{-3}&(1,-3,0)^\mathsf T\\
\rho_2&ML^{-3}&(1,-3,0)^\mathsf T\\
\mu_1&ML^{-1}T^{-1}&(1,-1,-1)^\mathsf T\\
\mu_2&ML^{-1}T^{-1}&(1,-1,-1)^\mathsf T\\
G_1&ML^{-1}T^{-2}&(1,-1,-2)^\mathsf T\\
G_2&ML^{-1}T^{-2}&(1,-1,-2)^\mathsf T.
\end{array}
$$

Therefore

$$
D=
\begin{pmatrix}
0&0&1&1&1&1&1&1\\
0&-1&-3&-3&-1&-1&-1&-1\\
-1&0&0&0&-1&-1&-2&-2
\end{pmatrix}.                                      \tag{1}
$$

For an exponent vector

$$
\boldsymbol a=(a_\gamma,a_k,a_{\rho_1},a_{\rho_2},
a_{\mu_1},a_{\mu_2},a_{G_1},a_{G_2})^\mathsf T,
$$

the monomial

$$
\Pi_{\boldsymbol a}
=\gamma^{a_\gamma}k^{a_k}
\rho_1^{a_{\rho_1}}\rho_2^{a_{\rho_2}}
\mu_1^{a_{\mu_1}}\mu_2^{a_{\mu_2}}
G_1^{a_{G_1}}G_2^{a_{G_2}}                         \tag{2}
$$

is dimensionless exactly when $D\boldsymbol a=0$.

## 2. Direct symbolic solution of $D\boldsymbol a=0$

Writing out the three rows of (1) gives

$$
\begin{aligned}
a_{\rho_1}+a_{\rho_2}+a_{\mu_1}+a_{\mu_2}+a_{G_1}+a_{G_2}&=0,\\
-a_k-3a_{\rho_1}-3a_{\rho_2}-a_{\mu_1}-a_{\mu_2}-a_{G_1}-a_{G_2}&=0,\\
-a_\gamma-a_{\mu_1}-a_{\mu_2}-2a_{G_1}-2a_{G_2}&=0.
\end{aligned}                                      \tag{3}
$$

Row reduction, performed without selecting repeating variables, gives

$$
\operatorname{rref}(D)=
\begin{pmatrix}
1&0&0&0&1&1&2&2\\
0&1&0&0&-2&-2&-2&-2\\
0&0&1&1&1&1&1&1
\end{pmatrix}.                                      \tag{4}
$$

The pivot variables are
$a_\gamma,a_k,a_{\rho_1}$. Introduce arbitrary free parameters

$$
c_1=a_{\rho_2},\quad c_2=a_{\mu_1},\quad
c_3=a_{\mu_2},\quad c_4=a_{G_1},\quad c_5=a_{G_2}.
$$

Equation (4) then gives

$$
\begin{aligned}
a_\gamma&=-c_2-c_3-2c_4-2c_5,\\
a_k&=2c_2+2c_3+2c_4+2c_5,\\
a_{\rho_1}&=-c_1-c_2-c_3-c_4-c_5,\\
a_{\rho_2}&=c_1,\quad a_{\mu_1}=c_2,\quad
a_{\mu_2}=c_3,\quad a_{G_1}=c_4,\quad a_{G_2}=c_5.
\end{aligned}                                      \tag{5}
$$

Equivalently,

$$
\boldsymbol a=c_1\boldsymbol n_1+c_2\boldsymbol n_2
+c_3\boldsymbol n_3+c_4\boldsymbol n_4+c_5\boldsymbol n_5, \tag{6}
$$

where

$$
\begin{aligned}
\boldsymbol n_1&=(0,0,-1,1,0,0,0,0)^\mathsf T,\\
\boldsymbol n_2&=(-1,2,-1,0,1,0,0,0)^\mathsf T,\\
\boldsymbol n_3&=(-1,2,-1,0,0,1,0,0)^\mathsf T,\\
\boldsymbol n_4&=(-2,2,-1,0,0,0,1,0)^\mathsf T,\\
\boldsymbol n_5&=(-2,2,-1,0,0,0,0,1)^\mathsf T.
\end{aligned}                                      \tag{7}
$$

Direct multiplication verifies $D\boldsymbol n_j=0$ for all five vectors.
Because each vector has a unique unit entry in one of the five free-variable
positions, they are linearly independent.

## 3. What the five-dimensional null space means

The columns for $\rho,\mu,k$, for example, contain the nonsingular minor

$$
\det
\begin{pmatrix}
1&1&0\\
-3&-1&-1\\
0&-1&0
\end{pmatrix}=-1.
$$

Hence $\operatorname{rank}D\ge3$. There are only three base dimensions, so
$\operatorname{rank}D=3$, and rank-nullity gives

$$
\dim\ker D=8-3=5.                                  \tag{8}
$$

The five dimensions admit a particularly useful physical decomposition.
Since each material pair has identical columns,

$$
\begin{aligned}
\boldsymbol d_\rho&=(0,0,1,-1,0,0,0,0)^\mathsf T,\\
\boldsymbol d_\mu&=(0,0,0,0,1,-1,0,0)^\mathsf T,\\
\boldsymbol d_G&=(0,0,0,0,0,0,1,-1)^\mathsf T
\end{aligned}                                      \tag{9}
$$

are automatically in $\ker D$. They generate the dimensionless ratios
$\rho_1/\rho_2$, $\mu_1/\mu_2$, and $G_1/G_2$. Thus three of the five
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

$$
K_-=\operatorname{span}\{\boldsymbol d_\rho,
\boldsymbol d_\mu,\boldsymbol d_G\}.
$$

Then $K_-\subset\ker D$, $\dim K_-=3$, and the quotient
$\ker D/K_-$ is two-dimensional. The assignment “three contrasts, one
response, one rate balance” is therefore a natural decomposition, although
the particular coordinates chosen within each part remain nonunique.

### Vector space versus invariant coordinates

The precise vector space produced by linear algebra is the exponent space
$\ker D$. Its addition law corresponds to multiplication of monomials:

$$
\Pi_{\boldsymbol a+\boldsymbol b}
=\Pi_{\boldsymbol a}\Pi_{\boldsymbol b},\qquad
\Pi_{c\boldsymbol a}=\Pi_{\boldsymbol a}^{c}.       \tag{10}
$$

Once any five independent monomial groups are known, every dimensionless
quantity can be written locally as a function of them. Such functions do not
themselves form the same linear exponent space. It is therefore useful to
distinguish:

1. the five-dimensional linear space $\ker D$ of monomial exponents;
2. a chosen vector-space basis of $\ker D$;
3. the five-dimensional quotient or invariant space, on which arbitrary
   invertible nonlinear coordinates may be used.

Buckingham's theorem fixes the dimension of this space, not a preferred basis
or preferred coordinates.

## 4. A basis obtained directly from row reduction

Using (7) in (2) gives the direct null-space basis

$$
\boxed{
N_1=\frac{\rho_2}{\rho_1},\quad
N_2=\frac{\mu_1k^2}{\gamma\rho_1},\quad
N_3=\frac{\mu_2k^2}{\gamma\rho_1},\quad
N_4=\frac{G_1k^2}{\gamma^2\rho_1},\quad
N_5=\frac{G_2k^2}{\gamma^2\rho_1}.}                \tag{11}
$$

This basis is a mechanically correct output of a particular row-reduction
convention. It is not physically canonical. A different column order changes
the pivot columns and therefore changes the basis returned by the same
algorithm. Moreover, any nonsingular matrix $C\in GL(5,\mathbb R)$ produces
another basis

$$
B=NC,                                               \tag{12}
$$

where $N$ is the $8\times5$ matrix with columns $\boldsymbol n_j$.
At the group level,

$$
\widetilde\Pi_j=\prod_{i=1}^5N_i^{C_{ij}}.          \tag{13}
$$

Thus nonuniqueness is the ordinary nonuniqueness of a basis in a
five-dimensional vector space.

## 5. Repeating variables as a basis selection

The classical construction chooses three dimensionally independent repeating
variables. Choose $(\rho_2,\mu_2,k)$. Their dimension columns form an
invertible $3\times3$ matrix. For each remaining variable, solving for the
three repeating-variable exponents constructs one kernel vector. The result is

$$
\begin{aligned}
P&=\frac{\rho_2\gamma}{\mu_2k^2},&
R_\rho&=\frac{\rho_1}{\rho_2},&
R_\mu&=\frac{\mu_1}{\mu_2},\\
Q_1&=\frac{\rho_2G_1}{\mu_2^2k^2},&
Q_2&=\frac{\rho_2G_2}{\mu_2^2k^2}.&&              \tag{14}
\end{aligned}
$$

This construction has not found a different invariant space. Indeed, the
groups in (14) are explicit products of the direct groups in (11):

$$
P=\frac{N_1}{N_3},\qquad
R_\rho=N_1^{-1},\qquad
R_\mu=\frac{N_2}{N_3},\qquad
Q_1=\frac{N_1N_4}{N_3^2},\qquad
Q_2=\frac{N_1N_5}{N_3^2}.                          \tag{15}
$$

The corresponding change-of-basis matrix is

$$
C=
\begin{pmatrix}
1&-1&0&1&1\\
0&0&1&0&0\\
-1&0&-1&-2&-2\\
0&0&0&1&0\\
0&0&0&0&1
\end{pmatrix},qquad \det C=1.                     \tag{16}
$$

Because $C$ is invertible, the exponent vectors in (14) and (11) span
exactly the same kernel. In general, choosing repeating variables amounts to
choosing an invertible $3\times3$ column minor of $D$, then using the
remaining five variables as free coordinates. It is a convenient algorithm
for selecting one kernel basis, not an additional theorem and not a uniqueness
principle.

## 6. Several valid bases and coordinate systems

### 6.1 Direct row-reduction basis

Equation (11) is a strict monomial basis. Its advantages are that it follows
immediately from linear algebra and makes no preliminary physical choice. Its
disadvantages are that the result depends on column ordering, privileges
$\rho_1$, mixes $\gamma$ into four groups, and obscures exchange symmetry.

### 6.2 Classical repeating-variable basis

Equation (14) is also a strict monomial basis. It has integer exponents and
isolates one response group. It is convenient when medium 2 is a genuine
reference material. It privileges medium 2, becomes unsuitable if a selected
repeating property vanishes, and does not display exchange symmetry.

### 6.3 Reference-scale ratio basis

Replacing $Q_1$ in (14) by $R_G=G_1/G_2=Q_1/Q_2$ gives

$$
\boxed{P,\ R_\rho,\ R_\mu,\ R_G,\ Q_2}.            \tag{17}
$$

This is a strict monomial basis because the replacement is an invertible basis
change. It cleanly exposes one degree of freedom for every material pair,
plus one response and one elastic-strength group. It still uses medium 2 as
the dimensional reference, and all three ratios are unbounded.

### 6.4 Exchange-adapted geometric-mean basis

Define geometric means

$$
\rho_g=\sqrt{\rho_1\rho_2},\qquad
\mu_g=\sqrt{\mu_1\mu_2},\qquad
G_g=\sqrt{G_1G_2}.
$$

Then

$$
\boxed{
P_g=\frac{\gamma\rho_g}{\mu_gk^2},\quad
L_g=\frac{\rho_gG_g}{\mu_g^2k^2},\quad
R_\rho,\quad R_\mu,\quad R_G}                     \tag{18}
$$

is a strict monomial basis, allowing half-integer exponents. Direct
substitution of its five exponent vectors into (1) gives zero, and their
$8\times5$ exponent matrix has rank five.

Under the exchange operator $1\leftrightarrow2$,

$$
P_g\mapsto P_g,\quad L_g\mapsto L_g,\quad
R_x\mapsto R_x^{-1}.                               \tag{19}
$$

Equivalently, $\log R_x\mapsto-\log R_x$. This basis diagonalizes exchange
symmetry into a two-dimensional even sector and a three-dimensional odd
sector in logarithmic coordinates. It is the most exchange-adapted *monomial*
basis. Its costs are fractional powers, unbounded ratios, and loss of a useful
coordinate when a property vanishes.

### 6.5 Total-property and contrast coordinates

Let

$$
\rho_T=\rho_1+\rho_2,\quad
\mu_T=\mu_1+\mu_2,\quad
G_T=G_1+G_2,
$$

and define

$$
s=\frac{\rho_T\gamma}{\mu_Tk^2},\qquad
\Lambda=\frac{\rho_TG_T}{\mu_T^2k^2},\qquad
A_x=\frac{x_1-x_2}{x_1+x_2}.                      \tag{20}
$$

The five quantities

$$
\boxed{s,\ \Lambda,\ A_\rho,\ A_\mu,\ A_G}       \tag{21}
$$

are dimensionless and independent on the positive-property interior. They are
not monomials in the eight original variables because they contain sums.
Consequently, they do not correspond to five exponent vectors and are not a
vector-space basis of $\ker D$. They are an invertible nonlinear coordinate
system on the same five-dimensional invariant space. They preserve exchange
symmetry in the especially transparent form

$$
s\mapsto s,\quad\Lambda\mapsto\Lambda,\quad
(A_\rho,A_\mu,A_G)\mapsto(-A_\rho,-A_\mu,-A_G).    \tag{22}
$$

These coordinates are adapted to the sums that actually occur in the
dispersion relation. The contrasts are bounded for positive properties. At
the exactly purely viscous point $G_T=0$, $A_G$ is unidentifiable but also
physically irrelevant, as explained in the existing specification.

## 7. Ratios and contrasts are invertible coordinates

For any positive pair $(x_1,x_2)$, let

$$
R_x=\frac{x_1}{x_2},\qquad
A_x=\frac{x_1-x_2}{x_1+x_2}.
$$

Dividing the numerator and denominator of $A_x$ by $x_2$ gives

$$
A_x=\frac{R_x-1}{R_x+1}.                           \tag{23}
$$

Solving for the ratio gives the inverse

$$
R_x=\frac{1+A_x}{1-A_x}.                           \tag{24}
$$

Thus (23) is a bijection from $R_x\in(0,\infty)$ to
$A_x\in(-1,1)$. It extends to zero/infinite ratios by including the endpoints
$-1$ and $+1$. In particular,

$$
(R_\rho,R_\mu,R_G)
\longleftrightarrow(A_\rho,A_\mu,A_G)
$$

is an invertible componentwise coordinate transformation. The contrasts do
not add new dimensionless information; they give bounded, exchange-odd
coordinates for the same three kernel directions.

## 8. Precise status of the preferred variables

The relation to the classical basis (14) makes the answer exact. From
$R_G=Q_1/Q_2$,

$$
\begin{aligned}
s&=P\frac{1+R_\rho}{1+R_\mu},\\
\Lambda&=\frac{(1+R_\rho)(Q_1+Q_2)}{(1+R_\mu)^2},\\
A_\rho&=\frac{R_\rho-1}{R_\rho+1},\qquad
A_\mu=\frac{R_\mu-1}{R_\mu+1},\\
A_G&=\frac{Q_1-Q_2}{Q_1+Q_2}.
\end{aligned}                                      \tag{25}
$$

Conversely, from the five preferred coordinates,

$$
R_\rho=\frac{1+A_\rho}{1-A_\rho},\quad
R_\mu=\frac{1+A_\mu}{1-A_\mu},\quad
R_G=\frac{1+A_G}{1-A_G},                           \tag{26}
$$

$$
P=s\frac{1+R_\mu}{1+R_\rho},\qquad
Q_\Sigma=Q_1+Q_2
=\Lambda\frac{(1+R_\mu)^2}{1+R_\rho},            \tag{27}
$$

and

$$
Q_1=\frac{1+A_G}{2}Q_\Sigma,\qquad
Q_2=\frac{1-A_G}{2}Q_\Sigma.                      \tag{28}
$$

Equations (25)–(28) prove invertibility on the positive-property interior.
Therefore:

- In the strict linear-algebraic sense, (21) is **not a basis of
  $\ker D$**, because its members are not monomials and have no single
  exponent vectors.
- In standard applied Buckingham terminology, it may be called a
  **transformed $\Pi$ basis**, because it is an invertible transformation of
  any monomial basis.
- Most precisely, it is a **convenient nonlinear coordinate system on the
  five-dimensional positive $\Pi$-space**.

These statements are compatible rather than contradictory; they refer to
different meanings of the word “basis.”

## 9. Is one basis mathematically preferable?

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
   space. Its $-1$ eigenspace is the three-dimensional pair-difference space
   (9), while its $+1$ sector inside $\ker D$ has dimension two. The
   geometric basis respects this decomposition with
   $(P_g,L_g)$ even and $(\log R_\rho,\log R_\mu,\log R_G)$ odd. The
   total/contrast coordinates express the same decomposition without
   logarithms.

2. **Algebraic simplicity.** The dimensional dispersion relation contains
   $\rho_1+\rho_2$, and the canonical derivation naturally normalizes by
   $\rho_T,\mu_T,G_T$. Substitution therefore yields fractions
   $(1\pm A_x)/2$ and a single rate parameter $\Lambda$. Reference or
   geometric bases require repeated rational conversions to these sums. For
   this equation—not by dimensional analysis alone—the total coordinates are
   algebraically preferable.

3. **Boundedness.** A nonconstant positive monomial is unbounded on the full
   positive invariant space: along a suitable logarithmic kernel direction it
   is $\exp(ct)$, which approaches either zero or infinity. Hence bounded
   contrasts cannot be obtained by merely choosing another linear basis of
   $\ker D$; a nonlinear transformation such as (23) is necessary. The
   contrast variables achieve boundedness without losing invertibility.

4. **Physical interpretation.** The quotient decomposition gives exactly
   three material-partition coordinates plus two common-scale coordinates.
   Choosing the latter as $s$ and $\Lambda=(r_e/r_v)^2$ identifies the
   response and the elastic-to-viscous rate competition directly. This is more
   closely tied to the governing equation than the raw RREF groups.

Accordingly, the geometric-mean construction (18) is the mathematically
cleanest strict monomial basis when exchange symmetry is the sole priority.
The total/contrast variables (21) are the preferable coordinates when all four
criteria are considered together. This null-space analysis therefore supports
the existing recommendation; it does not reveal an inconsistency or a reason
to change the canonical nondimensional formulation.
