# Symbolic dispersion-relation verification

This directory contains the Mathematica/Wolfram Language symbolic audit
framework for the Kelvin--Voigt dispersion relation.

The current symbolic workflow is organized into three maintained phases:

1. **Phase I** verifies the dimensional-to-nondimensional transformation
   $D(\gamma)\mapsto D^*_\Lambda(s)$.
2. **Phase II** is the authoritative radical-elimination and polynomial
   derivation phase. It constructs $P_{14}(s)$ and verifies the forward
   no-root-loss chain.
3. **Phase III** is a thin polynomial-completeness audit built from the
   verified Phase II identities. It does not repeat the full derivation.

The Phase I target is

$$
D(\gamma)\longrightarrow D^*(s)
$$

The Phase II/III target is the one-way theorem

$$
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0.
$$

The reverse implication is not claimed.

## Directory layout

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
  derivation from $D^*_\Lambda(s)=0$ to the polynomial candidate
  $P_{14}(s)=0$.
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
- `../root_verification_workflow.md` defines the project-wide terminology and
  implementation contract for polynomial candidate roots, mathematically
  genuine roots, and physically admissible roots.
- `logs/` stores newly generated audit logs.
- `exports/canonical/` stores the latest frozen implementation inputs.
- older timestamped logs and exports remain in the upstream research archive;
  they are not required by this standalone solver release.

## Phase I scope

Implemented sections:

1. Section 0 — Definitions
2. Section 1 — Original dimensional dispersion relation
3. Section 2 — Nondimensionalization

Not yet implemented:

1. Phase IV — Production polynomial-root filtering
2. Phase V — Parameter-region and branch-structure analysis

## Verification target

The framework verifies the exact normalization used by the numerical module:

$$
D^*_\Lambda(s;A_\rho,A_\mu,A_G,\Lambda)
=
\frac{\rho_T}{2k^2}
D(\gamma)
$$

after applying the substitutions

$$
\gamma=\frac{s\mu_Tk^2}{\rho_T},
$$

$$
\rho_1=\frac{\rho_T(1+A_\rho)}{2},\qquad
\rho_2=\frac{\rho_T(1-A_\rho)}{2},
$$

$$
\mu_1=\frac{\mu_T(1+A_\mu)}{2},\qquad
\mu_2=\frac{\mu_T(1-A_\mu)}{2},
$$

$$
G_1=\frac{\Lambda\mu_T^2k^2(1+A_G)}{2\rho_T},\qquad
G_2=\frac{\Lambda\mu_T^2k^2(1-A_G)}{2\rho_T}.
$$

Equivalently,

$$
D^*_\Lambda-\frac{\rho_T}{2k^2}D=0
$$

under the Phase I assumptions.

The present symbolic verification is established under the explicit theorem
assumption $\Lambda>0$. The limiting case $\Lambda=0$ corresponds, within this
chosen parameterization, to both shear moduli vanishing simultaneously
($E_k=1$). That boundary case lies outside the scope of the current Phase I
theorem and may be analyzed separately if needed. Individual elastic moduli may
still vanish within the theorem through the boundary values $A_G=\pm1$ when
$\Lambda>0$.

The audit also verifies the nondimensionalization component by component:

1. Kelvin--Voigt effective-viscosity scaling
2. Radicand mapping
3. Principal square-root mapping
4. Reciprocal-denominator scaling
5. Constant-term scaling
6. Reciprocal-term scaling
7. Whole-expression equivalence

## Running the audit

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

## Phase II theorem scope

Phase II establishes only the forward implication

$$
D^*_\Lambda(s)=0\Longrightarrow P_{14}(s)=0
$$

on the recorded admissible domain. It does not claim the reverse implication.
Spurious roots may enter when reciprocal denominators are cleared and when the
two radical-elimination squarings are performed.

## Phase III theorem scope

Phase III packages the final polynomial-completeness theorem:

$$
\Omega_{\mathrm{definition}}
\land
\Omega_{\mathrm{no\ root\ loss}}
\land
D^*_\Lambda(s)=0
\Longrightarrow
P_{14}(s)=0.
$$

Phase II is authoritative for the symbolic identities and algebraic
transformations; Phase III reuses those results. It does not prove or claim
the reverse implication.

## Runtime candidate filtering

Polynomial roots are only candidates. The full mathematical contract for the
later runtime workflow is documented in
[`../root_verification_workflow.md`](../root_verification_workflow.md).

In summary, the workflow is:

1. solve $P_{14}(s)=0$;
2. reject candidates outside the mathematical definition domain;
3. substitute each remaining candidate into $D^*_\Lambda(s)$;
4. retain only true roots of the original dispersion relation;
5. compute $q_-$ and $q_+$ on the principal square-root branch; and
6. apply the spatial decay conditions

   $$
   \operatorname{Re}(q_-)>0,
   \qquad
   \operatorname{Re}(q_+)>0.
   $$

This runtime filtering workflow is separate from the symbolic no-root-loss
theorem.

Temporal behavior is classified by $\operatorname{Re}(s)$. Physical spatial
admissibility is classified by $\operatorname{Re}(q_\pm)$. These are different
concepts and must not be conflated.

## Mathematica GUI validation checklist

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

## Branch convention

The symbolic expressions use Mathematica's principal square-root branch,
matching the NumPy principal branch used by `src/dispersion.py`. Phase I checks
verify algebraic equivalence of the dimensional and nondimensional expressions
on the same branch convention. Later admissible-domain work will make branch
selection and possible spurious roots explicit.
