"""Full arbitrary-precision reference path used by the manuscript sweeps.

Acceptance arithmetic runs at the requested precision. Returned roots and q
values are converted to binary64 for compatibility with stored sweep outputs.
Use separate processes for parallel solves; mpmath's context is not thread-local.
"""
from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from .polynomial_backend import polynomial_roots_from_ek


@dataclass(frozen=True)
class ReferencePointResult:
    """Counts include multiplicity; accepted rows are deduplicated as in sweeps."""

    polynomial_candidate_count: int
    mathematically_genuine_count: int
    accepted: tuple[tuple[complex, complex, complex, float], ...]


def solve_reference_point(
    Arho: str | float, Amu: str | float, AG: str | float, Ek: str | float,
    *, decimal_precision: int = 80,
) -> ReferencePointResult:
    """Solve and verify canonical P14 in the interior Ek theorem scope.

    Pass decimal strings to avoid rounding inputs before coefficient evaluation.
    Each accepted row contains (s, q_plus, q_minus, relative_source_residual).
    A solver convergence failure raises; it is never reported as zero roots.
    """
    if (
        isinstance(decimal_precision, bool)
        or not isinstance(decimal_precision, int)
        or decimal_precision < 30
    ):
        raise ValueError("decimal_precision must be an integer >= 30")
    with mp.workdps(decimal_precision):
        values = tuple(mp.mpf(str(value)) for value in (Arho, Amu, AG, Ek))
        if not all(mp.isfinite(value) for value in values):
            raise ValueError("parameters must be finite real numbers")
        arho, amu, ag, ek = values
        if not (-1 < arho < 1 and -1 < amu < 1 and -1 <= ag <= 1 and -1 < ek < 1):
            raise ValueError("requires |Arho| < 1, |Amu| < 1, |AG| <= 1, |Ek| < 1")
        return _solve(str(Arho), str(Amu), str(AG), str(Ek), decimal_precision)


def _solve(arho_text, amu_text, ag_text, ek_text, precision):
    roots, lam = polynomial_roots_from_ek(
        arho_text, amu_text, ag_text, ek_text, backend="mp",
        decimal_precision=precision,
    )
    arho = mp.mpf(arho_text)
    amu = mp.mpf(amu_text)
    ag = mp.mpf(ag_text)
    accepted = []
    genuine = 0
    residual_cutoff = mp.power(10, -(precision // 2))
    domain_cutoff = mp.power(10, -(precision // 2))
    for s in roots:
        try:
            ep = s * (1 + amu) + lam * (1 + ag)
            em = s * (1 - amu) + lam * (1 - ag)
            if min(abs(s), abs(ep), abs(em)) <= domain_cutoff:
                continue
            qp = mp.sqrt(1 + (1 + arho) * s * s / ep)
            qm = mp.sqrt(1 + (1 - arho) * s * s / em)
            bm = ep + em * qm
            bp = em + ep * qp
            if min(abs(bp), abs(bm)) <= domain_cutoff:
                continue
            tm = s * s / bm
            tp = s * s / bp
            value = tm + tp + 2
            relative = abs(value) / max(mp.mpf(1), abs(tm) + abs(tp) + 2)
            if relative > residual_cutoff:
                continue
            genuine += 1
            # Match the production physical-admissibility criterion.  A fixed
            # MP-scale cutoff over-accepts roots that are numerically positive
            # but marginal relative to a large |q|.
            q_relative_tolerance = mp.mpf("1e-10")
            q_absolute_tolerance = mp.mpf("1e-12")
            qtol_plus = q_absolute_tolerance + q_relative_tolerance * max(
                mp.mpf(1), abs(qp)
            )
            qtol_minus = q_absolute_tolerance + q_relative_tolerance * max(
                mp.mpf(1), abs(qm)
            )
            if mp.re(qp) <= qtol_plus or mp.re(qm) <= qtol_minus:
                continue
            candidate = complex(float(mp.re(s)), float(mp.im(s)))
            if not any(abs(candidate - row[0]) <= 1e-12 * max(1.0, abs(candidate)) for row in accepted):
                accepted.append((candidate, complex(qp), complex(qm), float(relative)))
        except (ZeroDivisionError, ValueError):
            continue
    accepted.sort(key=lambda row: (row[0].real, row[0].imag))
    return ReferencePointResult(len(roots), genuine, tuple(accepted))
