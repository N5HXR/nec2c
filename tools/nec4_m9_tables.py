#!/usr/bin/env python3
"""Documented NEC4 M9 table geometry and interpolation primitives.

This module intentionally contains no production dispatcher.  It freezes the
Section 5.1.3 coordinate conventions before direct field samples are stored.
"""

from dataclasses import dataclass
import cmath
import math


# Equation (5-10): l=1..3, m=0..3, n=0..2, followed by F1..F8.
MODEL_RAY_TERM_COUNT = 3 * 4 * 3
MODEL_SUPPLEMENTAL_TERM_COUNT = 8
MODEL_TOTAL_TERM_COUNT = MODEL_RAY_TERM_COUNT + MODEL_SUPPLEMENTAL_TERM_COUNT


def _beta(wave_number, lambda_value):
    """Outgoing beta=(k^2-lambda^2)^1/2 for the exp(-j phase) convention."""
    value = cmath.sqrt(wave_number * wave_number - lambda_value * lambda_value)
    if value.real < 0 or (abs(value.real) < 1e-14 * max(1.0, abs(value))
                          and value.imag > 0):
        value = -value
    return value


def saddle_lambda1(rho, z, source_depth, k_lower, k_upper,
                   tolerance=2e-13, iterations=40):
    """Solve F'(lambda1)=0 for the saddle connected to 0 <= lambda <= k2.

    This is the Section 5.1.3 / Appendix I air-dominant saddle.  Endpoint
    formulas avoid the indeterminate Newton equation at rho=0, z'=0, or z=0.
    General points use Newton's method on

        rho - lambda * (z/beta2 + source_depth/beta1) = 0.

    The caller supplies positive ``source_depth=|z'|``.
    """
    if min(rho, z, source_depth) < 0:
        raise ValueError("rho, z, and source depth must be nonnegative")
    if rho == 0:
        return 0j
    if source_depth == 0:
        return k_upper * rho / math.hypot(rho, z)
    if z == 0:
        return complex(k_upper)
    if (abs(k_lower.imag) <= 1e-14 * max(1.0, abs(k_lower)) and
            abs(k_upper.imag) <= 1e-14 * max(1.0, abs(k_upper))):
        # In lossless media the desired saddle is real and bracketed by the
        # upper-medium branch point.  Bisection avoids Newton landing exactly
        # on that branch point for grazing/deep sample combinations.
        low = 0.0
        high = k_upper.real * (1 - 2e-14)

        def residual_real(value):
            beta1 = math.sqrt(max(0.0, k_lower.real ** 2 - value ** 2))
            beta2 = math.sqrt(max(0.0, k_upper.real ** 2 - value ** 2))
            return rho - value * (z / beta2 + source_depth / beta1)

        for _ in range(80):
            midpoint = (low + high) / 2
            if residual_real(midpoint) > 0:
                low = midpoint
            else:
                high = midpoint
        return complex((low + high) / 2)
    scale = abs(k_upper / k_lower)
    estimate = k_upper * rho / math.hypot(rho, z + scale * source_depth)
    value = complex(estimate)
    for _ in range(iterations):
        beta1 = _beta(k_lower, value)
        beta2 = _beta(k_upper, value)
        if abs(beta1) == 0 or abs(beta2) == 0:
            value *= 1 - 1e-8
            continue
        residual = rho - value * (z / beta2 + source_depth / beta1)
        derivative = -(z * k_upper * k_upper / beta2 ** 3
                       + source_depth * k_lower * k_lower / beta1 ** 3)
        if derivative == 0:
            break
        update = residual / derivative
        value -= update
        if abs(update) <= tolerance * max(1.0, abs(value)):
            break
    beta1 = _beta(k_lower, value)
    beta2 = _beta(k_upper, value)
    residual = rho - value * (z / beta2 + source_depth / beta1)
    if abs(residual) > 2e-10 * max(1.0, rho):
        raise ArithmeticError(f"saddle did not converge: residual={residual!r}")
    if value.real > k_upper.real * (1 + 1e-10):
        return complex(k_upper)
    return value


def model_basis(rho, z, source_depth, k_lower, k_upper,
                include_supplemental=True, component=None):
    """Evaluate the documented equation (5-10) basis at one point.

    The scan places ``1/2`` after each supplemental exponential.  Reading it
    as a square root of either the phase or exponential conflicts with the
    preceding definition P2=exp(-F(lambda2)) and fits the guide's dielectric
    witness much worse.  This implementation therefore uses the ordinary
    second-ray phase exp[-j(beta3*z+k1*R2)]; the M9 model-selection artifact
    records the comparison rather than silently resolving the typography.

    ``component`` applies the two component-specific factors stated below
    equation (5-10): all terms receive rho/R1 for rho_v and z_h, while only
    l=1 ray terms receive it for z_v.
    """
    if min(rho, z, source_depth) < 0:
        raise ValueError("model coordinates must be nonnegative")
    z_prime = -source_depth
    ratio = abs(k_upper / k_lower)
    r1 = math.hypot(rho, z + ratio * source_depth)
    r2 = math.hypot(rho, source_depth)
    if r1 == 0:
        raise ValueError("equation (5-10) is undefined at R1=0")
    lambda1 = saddle_lambda1(rho, z, source_depth, k_lower, k_upper)
    beta1 = _beta(k_lower, lambda1)
    beta2 = _beta(k_upper, lambda1)
    phase1 = cmath.exp(-1j * (beta2 * z - beta1 * z_prime + lambda1 * rho))
    radial_factor = rho / r1
    values = []
    labels = []
    for ell in range(1, 4):
        for m in range(4):
            for n in range(3):
                value = (r1 ** -ell * (z / r1) ** m
                         * (z_prime / r1) ** n * phase1)
                if component in ("rho_v", "z_h") or (component == "z_v" and ell == 1):
                    value *= radial_factor
                values.append(value)
                labels.append(f"A_{ell}_{m}_{n}")
    if not include_supplemental:
        return tuple(values), tuple(labels)
    if r2 == 0:
        supplemental = [0j] * 8
    else:
        lambda2 = k_lower * rho / r2
        beta3 = _beta(k_upper, lambda2)
        phase2 = cmath.exp(-1j * (beta3 * z + k_lower * r2))
        f1 = r1 ** -1 * radial_factor * (z_prime / r1) * phase2
        f2 = r1 ** -1 * radial_factor * (z_prime / r1) ** 2 * phase2
        f3 = r1 ** -1 * r2 ** -1 * radial_factor * phase2
        f4 = f3 * (z_prime / r1)
        supplemental = [f1, f2, f3, f4,
                        f1 * z / r1, f2 * z / r1,
                        f3 * z / r1, f4 * z / r1]
        if component in ("rho_v", "z_h"):
            supplemental = [value * radial_factor for value in supplemental]
    values.extend(supplemental)
    labels.extend(f"B_{index}" for index in range(1, 9))
    return tuple(values), tuple(labels)


@dataclass(frozen=True)
class TableScales:
    wavelength0: float
    lg: float
    lmax: float


@dataclass(frozen=True)
class SampleAxis:
    start: float
    step: float
    count: int

    def values(self):
        return tuple(self.start + index * self.step
                     for index in range(self.count))


@dataclass(frozen=True)
class SampleLattice:
    rho: SampleAxis
    z: SampleAxis
    source_depth: SampleAxis

    def points(self):
        return tuple((rho, z, depth)
                     for rho in self.rho.values()
                     for z in self.z.values()
                     for depth in self.source_depth.values())


def table_scales(k_lower, k_upper, wavelength0):
    """Return Section 5.1.3 Lg and Lmax.

    The guide's medium 1 is lower and medium 2 is upper, hence the ratio
    |k2/k1| below.  ``wavelength0`` carries the units for all lengths.
    """
    if wavelength0 <= 0 or abs(k_lower) == 0:
        raise ValueError("positive wavelength and nonzero lower-medium k required")
    ratio = abs(k_upper / k_lower)
    return TableScales(float(wavelength0),
                       float(min(ratio, 0.3) * wavelength0),
                       float(min(15.0 * ratio, 1.0) * wavelength0))


def source_depth_region(source_z, scales):
    """Select the lowest-numbered table at a shared boundary.

    Section 5.1.3 specifies depths by |z'|.  A positive source height is not a
    lower-to-upper (12) table coordinate and is rejected here.
    """
    if source_z > 0:
        raise ValueError("the 12 tables require source_z <= 0")
    depth = abs(source_z)
    if depth <= 0.4 * scales.lg:
        return 1
    if depth <= scales.lg:
        return 2
    if depth <= scales.lmax:
        return 3
    return None


def near_table_contains(rho, observer_z, source_z, scales):
    return (rho >= 0 and observer_z >= 0 and source_z <= 0
            and rho <= scales.lg and observer_z <= scales.lg
            and source_depth_region(source_z, scales) is not None)


def separation(rho, observer_z, source_z):
    return math.hypot(rho, observer_z - source_z)


def storage_factor(region, rho, observer_z, source_z, k_lower, k_upper):
    """Section 5.1.3 F1, F2, and F3 storage transformations."""
    radius = separation(rho, observer_z, source_z)
    if radius == 0:
        if region == 1:
            return 0j
        raise ValueError("F2/F3 are undefined at zero separation")
    if region == 1:
        return radius * cmath.exp(1j * k_upper * radius)
    if region == 2:
        phase = (k_upper * math.hypot(rho, observer_z)
                 + k_lower * abs(source_z))
        return radius * radius * cmath.exp(1j * phase)
    if region == 3:
        phase = (k_upper * observer_z
                 + k_lower * math.hypot(rho, source_z))
        return radius * radius * cmath.exp(1j * phase)
    raise ValueError("region must be 1, 2, or 3")


def store_value(value, *factor_arguments):
    return value * storage_factor(*factor_arguments)


def restore_value(value, *factor_arguments):
    factor = storage_factor(*factor_arguments)
    if factor == 0:
        raise ValueError("the R=0 table-1 origin stores the analytic zero limit")
    return value / factor


def trilinear(bounds, corners, point):
    """Published Section 5.1.3 eight-corner linear interpolation.

    ``corners[(ir, iz, id)]`` uses 0 for the lower coordinate and 1 for the
    upper coordinate.  The algebraic barycentric form is equivalent to the
    printed eight-term equation while avoiding its sign-sensitive denominator.
    """
    weights = []
    for (low, high), coordinate in zip(bounds, point):
        if high <= low or coordinate < low or coordinate > high:
            raise ValueError("point must lie in a non-degenerate cell")
        fraction = (coordinate - low) / (high - low)
        weights.append((1.0 - fraction, fraction))
    result = 0j
    for ir in (0, 1):
        for iz in (0, 1):
            for idepth in (0, 1):
                result += (weights[0][ir] * weights[1][iz]
                           * weights[2][idepth]
                           * corners[(ir, iz, idepth)])
    return result


def _axis(start, step, count):
    return SampleAxis(float(start), float(step), count)


def model_sample_lattices(region, scales):
    """Return Figure 5-4/Table 1 general and surface sample lattices.

    Dimensionless constants in rho and z are multiplied by lambda0.  Constants
    containing Lg remain ground-scaled.  ``None`` denotes no surface lattice.
    """
    l0, lg = scales.wavelength0, scales.lg
    if region == 1:
        general = SampleLattice(_axis(0, 0.2 * l0, 4),
                                _axis(lg, 0.2 * l0, 4),
                                _axis(0, 0.4 * lg, 2))
        return general, None
    if region == 2:
        # Table 1 specifies five equally spaced rho samples from Lg to 0.6 l0.
        general = SampleLattice(_axis(lg, (0.6 * l0 - lg) / 4, 5),
                                _axis(0, 0.5 * lg, 4),
                                _axis(0, 0.4 * lg, 2))
        surface = SampleLattice(_axis(lg, 0.2 * lg, 6),
                                _axis(0, 0.4 * lg, 2),
                                _axis(0, 0.5 * lg, 3))
        return general, surface
    if region == 3:
        general = SampleLattice(_axis(0.6 * l0, 0.2 * l0, 5),
                                _axis(0, 0.2 * l0, 4),
                                _axis(0, 0.4 * lg, 2))
        surface = SampleLattice(_axis(0.6 * l0, 0.2 * lg, 8),
                                _axis(0, 0.4 * lg, 2),
                                _axis(0, 0.5 * lg, 3))
        return general, surface
    raise ValueError("region must be 1, 2, or 3")


def model_sample_count(region, scales):
    general, surface = model_sample_lattices(region, scales)
    return len(general.points()) + (0 if surface is None else len(surface.points()))
