"""Nondimensional Kelvin--Voigt dispersion-relation evaluators.

This module implements equations (6), (8), and (13) of
``THEORY.md#dispersion``.  It evaluates the
dispersion function only; it contains no root-finding logic.

NumPy's complex square root is used, so every radical is evaluated on the
principal branch.  Callers performing future continuation must track that
branch convention explicitly.

The dimensional and adopted nondimensional normalizations satisfy

    D_star_Lambda = (rho1 + rho2) / (2*k**2) * original_dimensional_D.

Multiplication of the dimensional equation by ``(rho1 + rho2) / k**2``
produces the pre-contrast equation with constant 4.  Contrast reconstruction
adds a factor 2 to its reciprocal terms, and the adopted final equation divides
the whole relation by 2, producing the prefactor above and constant 2.
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray


DispersionValue: TypeAlias = complex | NDArray[np.complex128]
GroupValue: TypeAlias = float | complex | NDArray[np.float64] | NDArray[np.complex128]


def _real_parameter(value: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return a finite real-valued parameter array."""

    raw = np.asarray(value)
    if np.iscomplexobj(raw) and np.any(np.imag(raw) != 0):
        raise ValueError(f"{name} must be real-valued")

    try:
        result = np.asarray(np.real(raw), dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a real numerical value") from exc

    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _contrast(value: ArrayLike, name: str) -> NDArray[np.float64]:
    """Validate a material contrast in its physical interval [-1, 1]."""

    result = _real_parameter(value, name)
    if np.any((result < -1.0) | (result > 1.0)):
        raise ValueError(f"{name} must lie in [-1, 1]")
    return result


def _complex_nonzero(value: ArrayLike, name: str) -> NDArray[np.complex128]:
    """Convert a value to a finite complex array and reject exact zero."""

    try:
        result = np.asarray(value, dtype=np.complex128)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a complex numerical value or array") from exc

    if not np.all(np.isfinite(result.real) & np.isfinite(result.imag)):
        raise ValueError(f"{name} must be finite")
    if np.any(result == 0):
        if name == "s":
            raise ValueError(
                "s = 0 is singular in the Kelvin--Voigt representation; "
                "evaluate an explicitly derived limit instead"
            )
        raise ValueError(
            f"{name} = 0 is singular in the dimensional Kelvin--Voigt "
            "representation"
        )
    return result


def _complex_response(s: ArrayLike) -> NDArray[np.complex128]:
    """Convert the nondimensional response to complex and reject s = 0."""

    return _complex_nonzero(s, "s")


def _positive_parameter(
    value: ArrayLike, name: str, *, allow_zero: bool = False
) -> NDArray[np.float64]:
    """Validate a positive, or optionally nonnegative, dimensional input."""

    result = _real_parameter(value, name)
    invalid = result < 0.0 if allow_zero else result <= 0.0
    if np.any(invalid):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _scalar_or_array(value: NDArray[np.complex128]) -> DispersionValue:
    """Return a Python complex for scalar inputs and an ndarray otherwise."""

    if value.ndim == 0:
        return complex(value.item())
    return value


def _native_scalar_or_array(value: NDArray) -> GroupValue:
    """Return a NumPy array, or its native scalar value for zero dimensions."""

    if value.ndim == 0:
        return value.item()
    return value


def _dimensional_inputs(
    gamma: ArrayLike,
    k: ArrayLike,
    rho1: ArrayLike,
    rho2: ArrayLike,
    mu1: ArrayLike,
    mu2: ArrayLike,
    G1: ArrayLike,
    G2: ArrayLike,
) -> tuple[NDArray, ...]:
    """Validate and broadcast the dimensional Kelvin--Voigt inputs."""

    values = (
        _complex_nonzero(gamma, "gamma"),
        _positive_parameter(k, "k"),
        _positive_parameter(rho1, "rho1"),
        _positive_parameter(rho2, "rho2"),
        _positive_parameter(mu1, "mu1"),
        _positive_parameter(mu2, "mu2"),
        _positive_parameter(G1, "G1", allow_zero=True),
        _positive_parameter(G2, "G2", allow_zero=True),
    )
    return tuple(np.broadcast_arrays(*values))


def original_dimensional_D(
    gamma: ArrayLike,
    k: ArrayLike,
    rho1: ArrayLike,
    rho2: ArrayLike,
    mu1: ArrayLike,
    mu2: ArrayLike,
    G1: ArrayLike,
    G2: ArrayLike,
) -> DispersionValue:
    """Evaluate the original dimensional dispersion relation from ``THEORY.md#specification``.

    The two radicals use NumPy's principal complex square-root branch, matching
    the branch convention used by the nondimensional evaluators.  Inputs obey
    NumPy broadcasting rules.  Densities, viscosities, and ``k`` must be
    positive; individual moduli may be zero.
    """

    gamma, k, rho1, rho2, mu1, mu2, G1, G2 = _dimensional_inputs(
        gamma, k, rho1, rho2, mu1, mu2, G1, G2
    )
    rho_total = rho1 + rho2

    try:
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            eta1 = mu1 + G1 / gamma
            eta2 = mu2 + G2 / gamma
            radical2 = np.sqrt(1.0 + rho2 * gamma / (eta2 * k**2))
            radical1 = np.sqrt(1.0 + rho1 * gamma / (eta1 * k**2))
            result = gamma * (
                1.0 / (eta1 + eta2 * radical2)
                + 1.0 / (eta2 + eta1 * radical1)
            ) + 4.0 * k**2 / rho_total
    except FloatingPointError as exc:
        raise FloatingPointError(
            "dimensional dispersion evaluation encountered a pole, overflow, "
            "or invalid intermediate value"
        ) from exc

    if not np.all(np.isfinite(result.real) & np.isfinite(result.imag)):
        raise FloatingPointError(
            "dimensional dispersion evaluation produced a non-finite value"
        )
    return _scalar_or_array(np.asarray(result, dtype=np.complex128))


def nondimensional_groups(
    gamma: ArrayLike,
    k: ArrayLike,
    rho1: ArrayLike,
    rho2: ArrayLike,
    mu1: ArrayLike,
    mu2: ArrayLike,
    G1: ArrayLike,
    G2: ArrayLike,
) -> tuple[GroupValue, GroupValue, GroupValue, GroupValue, GroupValue, GroupValue]:
    """Convert dimensional inputs to ``s, Arho, Amu, AG, Lambda, Ek``.

    ``AG`` is undefined when both moduli vanish, so this helper requires
    ``G1 + G2 > 0``.  An individual modulus may still be zero.
    """

    gamma, k, rho1, rho2, mu1, mu2, G1, G2 = _dimensional_inputs(
        gamma, k, rho1, rho2, mu1, mu2, G1, G2
    )
    rho_total = rho1 + rho2
    mu_total = mu1 + mu2
    G_total = G1 + G2
    if np.any(G_total == 0.0):
        raise ValueError("AG is undefined when G1 + G2 = 0")

    try:
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            s = rho_total * gamma / (mu_total * k**2)
            arho = (rho1 - rho2) / rho_total
            amu = (mu1 - mu2) / mu_total
            ag = (G1 - G2) / G_total
            lam = rho_total * G_total / (mu_total**2 * k**2)
            rate_ratio = np.sqrt(lam)
            ek = (1.0 - rate_ratio) / (1.0 + rate_ratio)
    except FloatingPointError as exc:
        raise FloatingPointError(
            "nondimensional conversion overflowed or produced an invalid value"
        ) from exc

    groups = (s, arho, amu, ag, lam, ek)
    if any(not np.all(np.isfinite(value)) for value in groups):
        raise FloatingPointError("nondimensional conversion produced non-finite data")
    return tuple(_native_scalar_or_array(np.asarray(value)) for value in groups)


def D_star_Lambda(
    s: ArrayLike,
    Arho: ArrayLike,
    Amu: ArrayLike,
    AG: ArrayLike,
    Lambda: ArrayLike,
) -> DispersionValue:
    """Evaluate the canonical Lambda-form dispersion function.

    Parameters are broadcast using NumPy rules.  ``s`` may be a scalar or an
    array and is always evaluated as complex, so negative radical arguments do
    not silently become real NaNs.  ``Arho``, ``Amu``, and ``AG`` must lie in
    [-1, 1], and ``Lambda`` must be nonnegative.

    The returned value is

        s * (1 / (eta1 + eta2*q2) + 1 / (eta2 + eta1*q1)) + 2,

    with the compact helper quantities defined exactly as in the final
    derivation.  No transformed equation is substituted for this definition.
    """

    s_array = _complex_response(s)
    arho = _contrast(Arho, "Arho")
    amu = _contrast(Amu, "Amu")
    ag = _contrast(AG, "AG")
    lam = _real_parameter(Lambda, "Lambda")
    if np.any(lam < 0.0):
        raise ValueError("Lambda must be nonnegative")

    s_array, arho, amu, ag, lam = np.broadcast_arrays(
        s_array, arho, amu, ag, lam
    )

    try:
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            eta1 = (1.0 + amu) + (lam / s_array) * (1.0 + ag)
            eta2 = (1.0 - amu) + (lam / s_array) * (1.0 - ag)

            q1 = np.sqrt(1.0 + ((1.0 + arho) * s_array) / eta1)
            q2 = np.sqrt(1.0 + ((1.0 - arho) * s_array) / eta2)

            denominator1 = eta1 + eta2 * q2
            denominator2 = eta2 + eta1 * q1
            result = s_array * (
                1.0 / denominator1 + 1.0 / denominator2
            ) + 2.0
    except FloatingPointError as exc:
        raise FloatingPointError(
            "dispersion evaluation encountered a pole, overflow, or invalid "
            "intermediate value"
        ) from exc

    if not np.all(np.isfinite(result.real) & np.isfinite(result.imag)):
        raise FloatingPointError("dispersion evaluation produced a non-finite value")

    return _scalar_or_array(np.asarray(result, dtype=np.complex128))


def D_star_Ek(
    s: ArrayLike,
    Arho: ArrayLike,
    Amu: ArrayLike,
    AG: ArrayLike,
    Ek: ArrayLike,
) -> DispersionValue:
    """Evaluate the bounded-Ek form of the dispersion function.

    The function uses the exact coordinate transformation

        Lambda = ((1 - Ek) / (1 + Ek))**2

    and then evaluates :func:`D_star_Lambda`.  Thus the two public evaluators
    share the same formula, normalization, and principal square-root branch.
    The physical finite-parameter interval is -1 < Ek <= 1.
    """

    ek = _real_parameter(Ek, "Ek")
    if np.any((ek <= -1.0) | (ek > 1.0)):
        raise ValueError("Ek must lie in the finite-parameter interval (-1, 1]")

    with np.errstate(divide="raise", invalid="raise", over="raise"):
        lam = ((1.0 - ek) / (1.0 + ek)) ** 2

    return D_star_Lambda(s, Arho, Amu, AG, lam)
