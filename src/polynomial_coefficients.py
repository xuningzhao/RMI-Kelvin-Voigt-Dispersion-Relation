"""Numerical coefficients for the verified degree-14 polynomial P14.

This module is a mechanical Python translation of
``symbolic/exports/canonical/P14_coefficients.wl``.

Coefficient ordering is low-to-high:

    P14(s; Arho, Amu, AG, Lambda) = sum(c[j] * s**j for j in range(15))

Most numerical root solvers, including ``numpy.roots``, expect descending
order.  Conversion to descending order is handled in ``src.polynomial`` and is
not baked into this coefficient evaluator.

Do not edit the algebra here by hand.  Regenerate this file from the frozen
Wolfram export if the symbolic theorem is deliberately revised.
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray


CoefficientArray: TypeAlias = NDArray[np.complex128]


COEFFICIENT_COUNT = 15
GENERIC_DEGREE = 14
COEFFICIENT_ORDER = "ascending"
SOURCE_EXPORT = "symbolic/exports/canonical/P14_coefficients.wl"


def p14_coefficients(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    dtype: type[np.complexfloating] = np.complex128,
) -> CoefficientArray:
    """Return ``[c0, c1, ..., c14]`` for the verified polynomial ``P14``.

    Parameters are the nondimensional groups of the canonical viscous-clock
    formulation.  The return order matches the Mathematica export exactly:
    coefficient ``j`` multiplies ``s**j``.
    """

    Arho = np.asarray(Arho, dtype=dtype)
    Amu = np.asarray(Amu, dtype=dtype)
    AG = np.asarray(AG, dtype=dtype)
    Lambda = np.asarray(Lambda, dtype=dtype)

    coefficients = np.asarray(
        [1024*AG**4*Lambda**8, 4096*AG**4*Lambda**7 + 4096*AG**3*Amu*Lambda**7, 
         6144*AG**4*Lambda**6 + 16384*AG**3*Amu*Lambda**6 + 6144*AG**2*Amu**2*Lambda**6 + 
          2048*AG**4*Lambda**7 + 1024*AG**3*Arho*Lambda**7, 
         4096*AG**4*Lambda**5 + 24576*AG**3*Amu*Lambda**5 + 24576*AG**2*Amu**2*Lambda**5 + 
          4096*AG*Amu**3*Lambda**5 + 6144*AG**4*Lambda**6 + 8192*AG**3*Amu*Lambda**6 + 
          4096*AG**3*Arho*Lambda**6 + 3072*AG**2*Amu*Arho*Lambda**6, 
         1024*AG**4*Lambda**4 + 16384*AG**3*Amu*Lambda**4 + 36864*AG**2*Amu**2*Lambda**4 + 
          16384*AG*Amu**3*Lambda**4 + 1024*Amu**4*Lambda**4 + 6144*AG**4*Lambda**5 + 
          24576*AG**3*Amu*Lambda**5 + 12288*AG**2*Amu**2*Lambda**5 + 
          6144*AG**3*Arho*Lambda**5 + 12288*AG**2*Amu*Arho*Lambda**5 + 
          3072*AG*Amu**2*Arho*Lambda**5 + 64*AG**2*Lambda**6 + 1472*AG**4*Lambda**6 + 
          1536*AG**3*Arho*Lambda**6 + 384*AG**2*Arho**2*Lambda**6 - 
          128*AG**4*Arho**2*Lambda**6, 4096*AG**3*Amu*Lambda**3 + 
          24576*AG**2*Amu**2*Lambda**3 + 24576*AG*Amu**3*Lambda**3 + 4096*Amu**4*Lambda**3 + 
          2048*AG**4*Lambda**4 + 24576*AG**3*Amu*Lambda**4 + 36864*AG**2*Amu**2*Lambda**4 + 
          8192*AG*Amu**3*Lambda**4 + 4096*AG**3*Arho*Lambda**4 + 
          18432*AG**2*Amu*Arho*Lambda**4 + 12288*AG*Amu**2*Arho*Lambda**4 + 
          1024*Amu**3*Arho*Lambda**4 + 256*AG**2*Lambda**5 + 2944*AG**4*Lambda**5 + 
          128*AG*Amu*Lambda**5 + 5888*AG**3*Amu*Lambda**5 + 4608*AG**3*Arho*Lambda**5 + 
          4608*AG**2*Amu*Arho*Lambda**5 + 1536*AG**2*Arho**2*Lambda**5 - 
          256*AG**4*Arho**2*Lambda**5 + 768*AG*Amu*Arho**2*Lambda**5 - 
          512*AG**3*Amu*Arho**2*Lambda**5, 6144*AG**2*Amu**2*Lambda**2 + 
          16384*AG*Amu**3*Lambda**2 + 6144*Amu**4*Lambda**2 + 8192*AG**3*Amu*Lambda**3 + 
          36864*AG**2*Amu**2*Lambda**3 + 24576*AG*Amu**3*Lambda**3 + 2048*Amu**4*Lambda**3 + 
          1024*AG**3*Arho*Lambda**3 + 12288*AG**2*Amu*Arho*Lambda**3 + 
          18432*AG*Amu**2*Arho*Lambda**3 + 4096*Amu**3*Arho*Lambda**3 + 
          384*AG**2*Lambda**4 + 1472*AG**4*Lambda**4 + 512*AG*Amu*Lambda**4 + 
          11776*AG**3*Amu*Lambda**4 + 64*Amu**2*Lambda**4 + 8832*AG**2*Amu**2*Lambda**4 + 
          4608*AG**3*Arho*Lambda**4 + 13824*AG**2*Amu*Arho*Lambda**4 + 
          4608*AG*Amu**2*Arho*Lambda**4 + 2304*AG**2*Arho**2*Lambda**4 - 
          128*AG**4*Arho**2*Lambda**4 + 3072*AG*Amu*Arho**2*Lambda**4 - 
          1024*AG**3*Amu*Arho**2*Lambda**4 + 384*Amu**2*Arho**2*Lambda**4 - 
          768*AG**2*Amu**2*Arho**2*Lambda**4 + 128*AG**2*Lambda**5 + 448*AG**4*Lambda**5 + 
          32*AG*Arho*Lambda**5 + 800*AG**3*Arho*Lambda**5 + 384*AG**2*Arho**2*Lambda**5 - 
          128*AG**4*Arho**2*Lambda**5 + 64*AG*Arho**3*Lambda**5 - 64*AG**3*Arho**3*Lambda**5, 
         4096*AG*Amu**3*Lambda + 4096*Amu**4*Lambda + 12288*AG**2*Amu**2*Lambda**2 + 
          24576*AG*Amu**3*Lambda**2 + 6144*Amu**4*Lambda**2 + 
          3072*AG**2*Amu*Arho*Lambda**2 + 12288*AG*Amu**2*Arho*Lambda**2 + 
          6144*Amu**3*Arho*Lambda**2 + 256*AG**2*Lambda**3 + 768*AG*Amu*Lambda**3 + 
          5888*AG**3*Amu*Lambda**3 + 256*Amu**2*Lambda**3 + 17664*AG**2*Amu**2*Lambda**3 + 
          5888*AG*Amu**3*Lambda**3 + 1536*AG**3*Arho*Lambda**3 + 
          13824*AG**2*Amu*Arho*Lambda**3 + 13824*AG*Amu**2*Arho*Lambda**3 + 
          1536*Amu**3*Arho*Lambda**3 + 1536*AG**2*Arho**2*Lambda**3 + 
          4608*AG*Amu*Arho**2*Lambda**3 - 512*AG**3*Amu*Arho**2*Lambda**3 + 
          1536*Amu**2*Arho**2*Lambda**3 - 1536*AG**2*Amu**2*Arho**2*Lambda**3 - 
          512*AG*Amu**3*Arho**2*Lambda**3 + 384*AG**2*Lambda**4 + 448*AG**4*Lambda**4 + 
          256*AG*Amu*Lambda**4 + 1792*AG**3*Amu*Lambda**4 + 128*AG*Arho*Lambda**4 + 
          1600*AG**3*Arho*Lambda**4 + 32*Amu*Arho*Lambda**4 + 
          2400*AG**2*Amu*Arho*Lambda**4 + 1152*AG**2*Arho**2*Lambda**4 - 
          128*AG**4*Arho**2*Lambda**4 + 768*AG*Amu*Arho**2*Lambda**4 - 
          512*AG**3*Amu*Arho**2*Lambda**4 + 256*AG*Arho**3*Lambda**4 - 
          128*AG**3*Arho**3*Lambda**4 + 64*Amu*Arho**3*Lambda**4 - 
          192*AG**2*Amu*Arho**3*Lambda**4, 1024*Amu**4 + 8192*AG*Amu**3*Lambda + 
          6144*Amu**4*Lambda + 3072*AG*Amu**2*Arho*Lambda + 4096*Amu**3*Arho*Lambda + 
          64*AG**2*Lambda**2 + 512*AG*Amu*Lambda**2 + 384*Amu**2*Lambda**2 + 
          8832*AG**2*Amu**2*Lambda**2 + 11776*AG*Amu**3*Lambda**2 + 1472*Amu**4*Lambda**2 + 
          4608*AG**2*Amu*Arho*Lambda**2 + 13824*AG*Amu**2*Arho*Lambda**2 + 
          4608*Amu**3*Arho*Lambda**2 + 384*AG**2*Arho**2*Lambda**2 + 
          3072*AG*Amu*Arho**2*Lambda**2 + 2304*Amu**2*Arho**2*Lambda**2 - 
          768*AG**2*Amu**2*Arho**2*Lambda**2 - 1024*AG*Amu**3*Arho**2*Lambda**2 - 
          128*Amu**4*Arho**2*Lambda**2 + 384*AG**2*Lambda**3 + 768*AG*Amu*Lambda**3 + 
          1792*AG**3*Amu*Lambda**3 + 128*Amu**2*Lambda**3 + 2688*AG**2*Amu**2*Lambda**3 + 
          192*AG*Arho*Lambda**3 + 800*AG**3*Arho*Lambda**3 + 128*Amu*Arho*Lambda**3 + 
          4800*AG**2*Amu*Arho*Lambda**3 + 2400*AG*Amu**2*Arho*Lambda**3 + 
          1152*AG**2*Arho**2*Lambda**3 + 2304*AG*Amu*Arho**2*Lambda**3 - 
          512*AG**3*Amu*Arho**2*Lambda**3 + 384*Amu**2*Arho**2*Lambda**3 - 
          768*AG**2*Amu**2*Arho**2*Lambda**3 + 384*AG*Arho**3*Lambda**3 - 
          64*AG**3*Arho**3*Lambda**3 + 256*Amu*Arho**3*Lambda**3 - 
          384*AG**2*Amu*Arho**3*Lambda**3 - 192*AG*Amu**2*Arho**3*Lambda**3 + 
          80*AG**2*Lambda**4 + 48*AG**4*Lambda**4 + 48*AG*Arho*Lambda**4 + 
          176*AG**3*Arho*Lambda**4 + 4*Arho**2*Lambda**4 + 120*AG**2*Arho**2*Lambda**4 - 
          28*AG**4*Arho**2*Lambda**4 + 32*AG*Arho**3*Lambda**4 - 32*AG**3*Arho**3*Lambda**4 + 
          4*Arho**4*Lambda**4 - 8*AG**2*Arho**4*Lambda**4 + 4*AG**4*Arho**4*Lambda**4, 
         2048*Amu**4 + 1024*Amu**3*Arho + 128*AG*Amu*Lambda + 256*Amu**2*Lambda + 
          5888*AG*Amu**3*Lambda + 2944*Amu**4*Lambda + 4608*AG*Amu**2*Arho*Lambda + 
          4608*Amu**3*Arho*Lambda + 768*AG*Amu*Arho**2*Lambda + 
          1536*Amu**2*Arho**2*Lambda - 512*AG*Amu**3*Arho**2*Lambda - 
          256*Amu**4*Arho**2*Lambda + 128*AG**2*Lambda**2 + 768*AG*Amu*Lambda**2 + 
          384*Amu**2*Lambda**2 + 2688*AG**2*Amu**2*Lambda**2 + 1792*AG*Amu**3*Lambda**2 + 
          128*AG*Arho*Lambda**2 + 192*Amu*Arho*Lambda**2 + 
          2400*AG**2*Amu*Arho*Lambda**2 + 4800*AG*Amu**2*Arho*Lambda**2 + 
          800*Amu**3*Arho*Lambda**2 + 384*AG**2*Arho**2*Lambda**2 + 
          2304*AG*Amu*Arho**2*Lambda**2 + 1152*Amu**2*Arho**2*Lambda**2 - 
          768*AG**2*Amu**2*Arho**2*Lambda**2 - 512*AG*Amu**3*Arho**2*Lambda**2 + 
          256*AG*Arho**3*Lambda**2 + 384*Amu*Arho**3*Lambda**2 - 
          192*AG**2*Amu*Arho**3*Lambda**2 - 384*AG*Amu**2*Arho**3*Lambda**2 - 
          64*Amu**3*Arho**3*Lambda**2 + 160*AG**2*Lambda**3 + 160*AG*Amu*Lambda**3 + 
          192*AG**3*Amu*Lambda**3 + 144*AG*Arho*Lambda**3 + 176*AG**3*Arho*Lambda**3 + 
          48*Amu*Arho*Lambda**3 + 528*AG**2*Amu*Arho*Lambda**3 + 16*Arho**2*Lambda**3 + 
          240*AG**2*Arho**2*Lambda**3 + 240*AG*Amu*Arho**2*Lambda**3 - 
          112*AG**3*Amu*Arho**2*Lambda**3 + 96*AG*Arho**3*Lambda**3 - 
          32*AG**3*Arho**3*Lambda**3 + 32*Amu*Arho**3*Lambda**3 - 
          96*AG**2*Amu*Arho**3*Lambda**3 + 16*Arho**4*Lambda**3 - 
          16*AG**2*Arho**4*Lambda**3 - 16*AG*Amu*Arho**4*Lambda**3 + 
          16*AG**3*Amu*Arho**4*Lambda**3, 64*Amu**2 + 1472*Amu**4 + 1536*Amu**3*Arho + 
          384*Amu**2*Arho**2 - 128*Amu**4*Arho**2 + 256*AG*Amu*Lambda + 
          384*Amu**2*Lambda + 1792*AG*Amu**3*Lambda + 448*Amu**4*Lambda + 
          32*AG*Arho*Lambda + 128*Amu*Arho*Lambda + 2400*AG*Amu**2*Arho*Lambda + 
          1600*Amu**3*Arho*Lambda + 768*AG*Amu*Arho**2*Lambda + 
          1152*Amu**2*Arho**2*Lambda - 512*AG*Amu**3*Arho**2*Lambda - 
          128*Amu**4*Arho**2*Lambda + 64*AG*Arho**3*Lambda + 256*Amu*Arho**3*Lambda - 
          192*AG*Amu**2*Arho**3*Lambda - 128*Amu**3*Arho**3*Lambda + 80*AG**2*Lambda**2 + 
          320*AG*Amu*Lambda**2 + 80*Amu**2*Lambda**2 + 288*AG**2*Amu**2*Lambda**2 + 
          144*AG*Arho*Lambda**2 + 144*Amu*Arho*Lambda**2 + 528*AG**2*Amu*Arho*Lambda**2 + 
          528*AG*Amu**2*Arho*Lambda**2 + 24*Arho**2*Lambda**2 + 
          120*AG**2*Arho**2*Lambda**2 + 480*AG*Amu*Arho**2*Lambda**2 + 
          120*Amu**2*Arho**2*Lambda**2 - 168*AG**2*Amu**2*Arho**2*Lambda**2 + 
          96*AG*Arho**3*Lambda**2 + 96*Amu*Arho**3*Lambda**2 - 
          96*AG**2*Amu*Arho**3*Lambda**2 - 96*AG*Amu**2*Arho**3*Lambda**2 + 
          24*Arho**4*Lambda**2 - 8*AG**2*Arho**4*Lambda**2 - 32*AG*Amu*Arho**4*Lambda**2 - 
          8*Amu**2*Arho**4*Lambda**2 + 24*AG**2*Amu**2*Arho**4*Lambda**2 + 
          16*AG**2*Lambda**3 + 20*AG*Arho*Lambda**3 + 12*AG**3*Arho*Lambda**3 + 
          4*Arho**2*Lambda**3 + 12*AG**2*Arho**2*Lambda**3 + 4*AG*Arho**3*Lambda**3 - 
          4*AG**3*Arho**3*Lambda**3, 128*Amu**2 + 448*Amu**4 + 32*Amu*Arho + 
          800*Amu**3*Arho + 384*Amu**2*Arho**2 - 128*Amu**4*Arho**2 + 64*Amu*Arho**3 - 
          64*Amu**3*Arho**3 + 160*AG*Amu*Lambda + 160*Amu**2*Lambda + 
          192*AG*Amu**3*Lambda + 48*AG*Arho*Lambda + 144*Amu*Arho*Lambda + 
          528*AG*Amu**2*Arho*Lambda + 176*Amu**3*Arho*Lambda + 16*Arho**2*Lambda + 
          240*AG*Amu*Arho**2*Lambda + 240*Amu**2*Arho**2*Lambda - 
          112*AG*Amu**3*Arho**2*Lambda + 32*AG*Arho**3*Lambda + 96*Amu*Arho**3*Lambda - 
          96*AG*Amu**2*Arho**3*Lambda - 32*Amu**3*Arho**3*Lambda + 16*Arho**4*Lambda - 
          16*AG*Amu*Arho**4*Lambda - 16*Amu**2*Arho**4*Lambda + 
          16*AG*Amu**3*Arho**4*Lambda + 16*AG**2*Lambda**2 + 32*AG*Amu*Lambda**2 + 
          40*AG*Arho*Lambda**2 + 20*Amu*Arho*Lambda**2 + 36*AG**2*Amu*Arho*Lambda**2 + 
          12*Arho**2*Lambda**2 + 12*AG**2*Arho**2*Lambda**2 + 24*AG*Amu*Arho**2*Lambda**2 + 
          8*AG*Arho**3*Lambda**2 + 4*Amu*Arho**3*Lambda**2 - 12*AG**2*Amu*Arho**3*Lambda**2, 
         80*Amu**2 + 48*Amu**4 + 48*Amu*Arho + 176*Amu**3*Arho + 4*Arho**2 + 
          120*Amu**2*Arho**2 - 28*Amu**4*Arho**2 + 32*Amu*Arho**3 - 32*Amu**3*Arho**3 + 
          4*Arho**4 - 8*Amu**2*Arho**4 + 4*Amu**4*Arho**4 + 32*AG*Amu*Lambda + 
          16*Amu**2*Lambda + 20*AG*Arho*Lambda + 40*Amu*Arho*Lambda + 
          36*AG*Amu**2*Arho*Lambda + 12*Arho**2*Lambda + 24*AG*Amu*Arho**2*Lambda + 
          12*Amu**2*Arho**2*Lambda + 4*AG*Arho**3*Lambda + 8*Amu*Arho**3*Lambda - 
          12*AG*Amu**2*Arho**3*Lambda + AG**2*Lambda**2 + 2*AG*Arho*Lambda**2 + 
          Arho**2*Lambda**2, 16*Amu**2 + 20*Amu*Arho + 12*Amu**3*Arho + 4*Arho**2 + 
          12*Amu**2*Arho**2 + 4*Amu*Arho**3 - 4*Amu**3*Arho**3 + 2*AG*Amu*Lambda + 
          2*AG*Arho*Lambda + 2*Amu*Arho*Lambda + 2*Arho**2*Lambda, 
         Amu**2 + 2*Amu*Arho + Arho**2],
        dtype=dtype,
    )
    return coefficients


def p14_value(
    s: complex,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    dtype: type[np.complexfloating] = np.complex128,
) -> complex:
    """Evaluate ``P14(s)`` using low-to-high coefficients and Horner form."""

    coeffs = p14_coefficients(Arho, Amu, AG, Lambda, dtype=dtype)
    value = np.asarray(0.0, dtype=dtype)
    z = np.asarray(s, dtype=dtype)
    for coeff in coeffs[::-1]:
        value = value * z + coeff
    return complex(value.item() if value.ndim == 0 else value)
