#!/usr/bin/env python3
"""High-precision equation layer for the NEC4 M8 Sommerfeld oracle.

This module intentionally contains no quadrature and no production-solver
integration.  It locks the convention, spectral square-root sheet, equations
(5-2)/(5-7), and the cylindrical-to-Cartesian assembly in equation (5-1)
before M8 adds contour integration.

Convention: exp(+j omega t), with outgoing homogeneous waves exp(-j k R).
The medium numbering follows the theory guide: medium 1 is below z=0 and
medium 2 is above z=0.
"""

from dataclasses import dataclass

import mpmath as mp


C0 = mp.mpf("299792458")
EPS0 = mp.mpf("8.8541878128e-12")
MU0 = 1 / (EPS0 * C0 * C0)


@dataclass(frozen=True)
class Medium:
    """Homogeneous isotropic nonmagnetic medium at one frequency."""

    epsilon_relative: object
    wave_number: object


@dataclass(frozen=True)
class CrossGeometry:
    """Canonical source-below/observer-above geometry for equations (5-2)."""

    rho: object
    source_z: object
    observer_z: object

    def validate(self):
        if self.rho < 0:
            raise ValueError("rho must be nonnegative")
        if self.source_z > 0 or self.observer_z < 0:
            raise ValueError("canonical 12 geometry requires source_z <= 0 <= observer_z")
        if self.rho == 0 and self.source_z == 0 and self.observer_z == 0:
            raise ValueError("source and observer cannot coincide at the interface")


def material(frequency_hz, epsilon_real, conductivity=0, epsilon_loss=None):
    """Return a medium using the guide's exp(+j omega t) loss convention."""
    frequency_hz = mp.mpf(frequency_hz)
    epsilon_real = mp.mpf(epsilon_real)
    conductivity = mp.mpf(conductivity)
    if not frequency_hz > 0 or not epsilon_real > 0 or conductivity < 0:
        raise ValueError("frequency/epsilon must be positive and conductivity nonnegative")
    omega = 2 * mp.pi * frequency_hz
    if epsilon_loss is None:
        epsilon = mp.mpc(epsilon_real, -conductivity / (omega * EPS0))
    else:
        epsilon = mp.mpc(epsilon_real, -abs(mp.mpf(epsilon_loss)))
    # The principal root has Re(k) >= 0 and Im(k) <= 0 for passive media.
    wave_number = omega / C0 * mp.sqrt(epsilon)
    if mp.re(wave_number) < 0 or (mp.almosteq(mp.re(wave_number), 0) and
                                 mp.im(wave_number) > 0):
        wave_number = -wave_number
    return Medium(epsilon, wave_number)


def outgoing_gamma(lambda_value, wave_number):
    """Return gamma=(lambda^2-k^2)^1/2 on the outgoing Bessel-contour sheet.

    On the positive real lambda axis this is continuous from gamma(0)=+j*k
    to positive-real gamma above the branch point.  The first-quadrant Bessel
    contour used in Figure 5-2 remains on this principal continuation.  The
    multi-sheet Hankel contours are deliberately deferred to Gate B.
    """
    gamma = mp.sqrt(lambda_value * lambda_value - wave_number * wave_number)
    scale = max(mp.mpf(1), abs(gamma))
    tolerance = 32 * mp.eps * scale
    if mp.re(gamma) < -tolerance or (abs(mp.re(gamma)) <= tolerance and
                                    mp.im(gamma) < 0):
        gamma = -gamma
    return gamma


def _common(lambda_value, geometry, lower, upper):
    geometry.validate()
    gamma1 = outgoing_gamma(lambda_value, lower.wave_number)
    gamma2 = outgoing_gamma(lambda_value, upper.wave_number)
    exponential = mp.exp(gamma1 * geometry.source_z -
                         gamma2 * geometry.observer_z)
    denominator = (lower.wave_number ** 2 * gamma2 +
                   upper.wave_number ** 2 * gamma1)
    return gamma1, gamma2, exponential, denominator


def potential_densities(lambda_value, geometry, lower, upper):
    """Spectral densities of V12 and U12 in equations (5-3)/(5-4)."""
    gamma1, gamma2, exponential, denominator = _common(
        lambda_value, geometry, lower, upper)
    bessel0 = mp.besselj(0, lambda_value * geometry.rho)
    v12 = 2 * exponential * bessel0 * lambda_value / denominator
    u12 = 2 * exponential * bessel0 * lambda_value / (gamma1 + gamma2)
    return v12, u12


def appendix_h_difference_terms(gamma1, gamma2, source_z,
                                k1_squared, k2_squared):
    """Cancellation-resistant ``T1``, ``T2``, and ``T3`` from Appendix H.

    Equations (H-3)/(H-4) define the terms with
    ``exp(source_z*(gamma1-gamma2))``.  Direct subtraction loses digits when
    ``|lambda|`` is large.  ``gamma1-gamma2`` is therefore formed through the
    exact conjugate-product identity and the exponential through ``expm1``.
    This is algebraically equivalent to the guide's inverse-lambda series but
    remains useful away from its asymptotic region as well.
    """
    denominator = gamma1 + gamma2
    if denominator == 0:
        raise ZeroDivisionError("Appendix H gamma difference is singular")
    # gamma_i^2=lambda^2-k_i^2, hence gamma1-gamma2=(k2^2-k1^2)/(gamma1+gamma2).
    # Use the material wave numbers directly: reconstructing the numerator as
    # gamma1^2-gamma2^2 would reintroduce catastrophic cancellation.
    delta = (k2_squared - k1_squared) / denominator
    t1 = mp.expm1(source_z * delta)
    t2 = gamma2 * t1 - delta
    t3 = gamma1 * t1 + delta
    return t1, t2, t3


def appendix_h_vprime_densities(lambda_value, geometry, lower, upper):
    """Stable spectral densities for ``V'12`` and ``dV'12/dz'``.

    These are equations (H-3) and (H-4), including ``J0(lambda*rho)`` and
    the integration measure ``lambda``.  They are potential densities, not
    electric-field densities.
    """
    gamma1, gamma2, _, denominator = _common(
        lambda_value, geometry, lower, upper)
    k1_squared = lower.wave_number ** 2
    k2_squared = upper.wave_number ** 2
    sum_squared = k1_squared + k2_squared
    t1, t2, t3 = appendix_h_difference_terms(
        gamma1, gamma2, geometry.source_z, k1_squared, k2_squared)
    propagation = mp.exp(-gamma2 *
                         (geometry.observer_z - geometry.source_z))
    measure = mp.besselj(0, lambda_value * geometry.rho) * lambda_value
    vprime = (2 * (gamma2 * k1_squared * t1 + k2_squared * t2) /
              (gamma2 * sum_squared * denominator) * propagation * measure)
    source_derivative = (
        2 * (gamma1 * k2_squared * t1 + k1_squared * t3) /
        (sum_squared * denominator) * propagation * measure)
    return vprime, source_derivative


def appendix_h_remainder_field_densities(lambda_value, geometry, lower, upper,
                                         frequency_hz):
    """Five field densities after the first Appendix H subtraction.

    This is equations (H-1a-e) evaluated from the cancellation-resistant
    (H-3)/(H-4) coefficients. Integrating these densities gives ``E'``;
    adding :func:`appendix_h_v0_fields` recovers the total field.
    """
    gamma1, gamma2, _, denominator = _common(
        lambda_value, geometry, lower, upper)
    k1_squared = lower.wave_number ** 2
    k2_squared = upper.wave_number ** 2
    sum_squared = k1_squared + k2_squared
    t1, t2, t3 = appendix_h_difference_terms(
        gamma1, gamma2, geometry.source_z, k1_squared, k2_squared)
    propagation = mp.exp(-gamma2 *
                         (geometry.observer_z - geometry.source_z))
    v_coefficient = (2 * (gamma2 * k1_squared * t1 + k2_squared * t2) /
                     (gamma2 * sum_squared * denominator) * propagation)
    source_coefficient = (
        2 * (gamma1 * k2_squared * t1 + k1_squared * t3) /
        (sum_squared * denominator) * propagation)
    argument = lambda_value * geometry.rho
    j0 = mp.besselj(0, argument)
    j0_prime = -mp.besselj(1, argument)
    j0_second = (mp.besselj(2, argument) - j0) / 2
    u12 = (2 * mp.exp(gamma1 * geometry.source_z -
                      gamma2 * geometry.observer_z) * j0 * lambda_value /
           (gamma1 + gamma2))
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    q = 1j * omega * MU0 / (4 * mp.pi)
    return {
        "rho_v": q * gamma2 * v_coefficient * j0_prime * lambda_value ** 2,
        "z_v": -q * v_coefficient * j0 * lambda_value ** 3,
        "rho_h": -q * (v_coefficient * j0_second * lambda_value ** 3 + u12),
        "phi_h": q * (v_coefficient * j0_prime * lambda_value ** 2 /
                       geometry.rho + u12),
        "z_h": q * source_coefficient * j0_prime * lambda_value ** 2,
    }


def appendix_h_v0(frequency_hz, lower, upper, geometry):
    """Closed-form singular potential ``V0`` in Appendix H."""
    geometry.validate()
    dz = geometry.observer_z - geometry.source_z
    distance = mp.sqrt(geometry.rho ** 2 + dz ** 2)
    k2 = upper.wave_number
    return (2 / (lower.wave_number ** 2 + k2 ** 2) *
            mp.exp(-1j * k2 * distance) / distance)


def appendix_h_v0_fields(frequency_hz, lower, upper, geometry):
    """Analytic field terms (H-2a-e) generated by ``V0``."""
    geometry.validate()
    rho = geometry.rho
    dz = geometry.observer_z - geometry.source_z
    distance = mp.sqrt(rho ** 2 + dz ** 2)
    k2 = upper.wave_number
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    factor = (-1j * omega * MU0 / (4 * mp.pi) *
              2 / (lower.wave_number ** 2 + k2 ** 2) *
              mp.exp(-1j * k2 * distance) / distance ** 3)
    radial = 3 / distance ** 2 + 3j * k2 / distance - k2 ** 2
    return {
        "rho_v": factor * rho * dz * radial,                 # (H-2a)
        "z_v": factor * (dz ** 2 * radial - 1 -
                           1j * k2 * distance + k2 ** 2 * distance ** 2),
        "rho_h": factor * (rho ** 2 * radial - 1 -
                             1j * k2 * distance),             # (H-2c)
        "phi_h": factor * (1 + 1j * k2 * distance),          # (H-2d)
        "z_h": factor * rho * dz * radial,                   # (H-2e)
    }


def appendix_h_residual_singular_fields(frequency_hz, lower, upper, geometry):
    """The residual ``1/R`` field terms following equations (H-3)/(H-4).

    Subtracting these ``S'`` terms from the fields remaining after ``V0``
    subtraction produces Appendix H's finite ``E''`` interpolation values.
    The formula is evaluated only for ``rho>0``; the axial limiting forms are
    a separate M8 qualification case.
    """
    geometry.validate()
    if geometry.rho == 0:
        raise ValueError("Appendix H residual uses a directional rho->0 limit")
    rho = geometry.rho
    source_z = geometry.source_z
    dz = geometry.observer_z - source_z
    distance = mp.sqrt(rho ** 2 + dz ** 2)
    sine = dz / distance
    cosine = rho / distance
    k1_squared = lower.wave_number ** 2
    k2_squared = upper.wave_number ** 2
    c1 = (k1_squared - k2_squared) / (k1_squared + k2_squared)
    c2 = k2_squared * (k1_squared - k2_squared) / (k1_squared + k2_squared) ** 2
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    negative = (-1j * omega * MU0 / (4 * mp.pi) *
                mp.exp(-1j * upper.wave_number * distance) / distance)
    positive = -negative
    one_minus_sine = 1 - sine
    return {
        "rho_v": negative * (c2 * one_minus_sine / cosine -
                               c1 * source_z * cosine / distance),
        "z_v": negative * (c2 - c1 * source_z * sine / distance),
        "rho_h": negative * (
            c2 * (one_minus_sine / cosine ** 2 - 1) +
            c1 * source_z / distance *
            ((sine * (1 + cosine ** 2) - 1) / cosine ** 2) + 1),
        "phi_h": positive * (
            -c2 * one_minus_sine / cosine ** 2 +
            c1 * source_z * one_minus_sine /
            (distance * cosine ** 2) + 1),
        "z_h": positive * (
            k1_squared * c2 / k2_squared * one_minus_sine / cosine +
            c1 * source_z * cosine / distance),
    }


def expanded_field_densities(lambda_value, geometry, lower, upper,
                             frequency_hz):
    """Five canonical E-field spectral densities from equations (5-7a-e).

    Names use source orientation, not magnetic field: ``rho_h`` means the
    radial electric field due to a horizontal electric-current element.
    """
    gamma1, gamma2, exponential, denominator = _common(
        lambda_value, geometry, lower, upper)
    argument = lambda_value * geometry.rho
    j0 = mp.besselj(0, argument)
    j0_prime = -mp.besselj(1, argument)
    j0_second = (mp.besselj(2, argument) - j0) / 2
    _, u12 = potential_densities(lambda_value, geometry, lower, upper)
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    q = 1j * omega * MU0 / (4 * mp.pi)
    if geometry.rho == 0:
        phi_v_term = -exponential * lambda_value ** 3 / denominator
    else:
        phi_v_term = (2 * exponential * j0_prime * lambda_value ** 2 /
                      (geometry.rho * denominator))
    return {
        # Equations (5-7a) through (5-7e), in order.
        "rho_v": q * 2 * gamma2 * exponential * j0_prime *
                 lambda_value ** 2 / denominator,
        "z_v": -q * 2 * exponential * j0 * lambda_value ** 3 / denominator,
        "rho_h": -q * (2 * exponential * j0_second * lambda_value ** 3 /
                        denominator + u12),
        "phi_h": q * (phi_v_term + u12),
        "z_h": q * 2 * gamma1 * exponential * j0_prime *
               lambda_value ** 2 / denominator,
    }


def appendix_f_21_field_densities(lambda_value, rho, source_z, observer_z,
                                  lower, upper, frequency_hz):
    """Direct Appendix F density for source above and observer below.

    This deliberately does not use Appendix G reciprocity.  ``source_z`` is
    positive, ``observer_z`` negative, medium 1 is ``lower``, and medium 2 is
    ``upper``.  The returned component convention omits the printed
    ``cos(phi)``/``sin(phi)`` orientation factors, matching the canonical
    component dictionaries used throughout the M8 oracle.
    """
    rho = mp.mpf(rho)
    source_z = mp.mpf(source_z)
    observer_z = mp.mpf(observer_z)
    if not rho > 0 or not source_z > 0 or not observer_z < 0:
        raise ValueError("Appendix F 21 requires rho>0, source_z>0, observer_z<0")
    gamma1 = outgoing_gamma(lambda_value, lower.wave_number)
    gamma2 = outgoing_gamma(lambda_value, upper.wave_number)
    exponential = mp.exp(gamma1 * observer_z - gamma2 * source_z)
    denominator = (lower.wave_number ** 2 * gamma2 +
                   upper.wave_number ** 2 * gamma1)
    argument = lambda_value * rho
    j0 = mp.besselj(0, argument)
    j0_prime = -mp.besselj(1, argument)
    j0_second = (mp.besselj(2, argument) - j0) / 2
    u21 = 2 * exponential * j0 * lambda_value / (gamma1 + gamma2)
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    q = 1j * omega * MU0 / (4 * mp.pi)
    return {
        "rho_v": -q * 2 * gamma1 * exponential * j0_prime *
                 lambda_value ** 2 / denominator,
        "z_v": -q * 2 * exponential * j0 * lambda_value ** 3 / denominator,
        "rho_h": -q * (2 * exponential * j0_second * lambda_value ** 3 /
                        denominator + u21),
        "phi_h": q * (2 * exponential * j0_prime * lambda_value ** 2 /
                       (rho * denominator) + u21),
        "z_h": -q * 2 * gamma2 * exponential * j0_prime *
               lambda_value ** 2 / denominator,
    }


def hankel_field_densities_with_roots(lambda_value, geometry, lower, upper,
                                      frequency_hz, gamma1, gamma2):
    """Full-line Hankel densities with caller-supplied continuous roots."""
    geometry.validate()
    exponential = mp.exp(gamma1 * geometry.source_z -
                         gamma2 * geometry.observer_z)
    denominator = (lower.wave_number ** 2 * gamma2 +
                   upper.wave_number ** 2 * gamma1)
    argument = lambda_value * geometry.rho
    h0 = mp.hankel2(0, argument)
    h0_prime = -mp.hankel2(1, argument)
    # Bessel's equation avoids a third expensive Hankel evaluation.
    h0_second = -h0 - h0_prime / argument
    u12 = exponential * h0 * lambda_value / (gamma1 + gamma2)
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    q = 1j * omega * MU0 / (4 * mp.pi)
    return {
        # The factor 2 in (5-7) is absent because (5-5)/(5-6) span -inf..inf.
        "rho_v": q * gamma2 * exponential * h0_prime *
                 lambda_value ** 2 / denominator,
        "z_v": -q * exponential * h0 * lambda_value ** 3 / denominator,
        "rho_h": -q * (exponential * h0_second * lambda_value ** 3 /
                        denominator + u12),
        "phi_h": q * (exponential * h0_prime * lambda_value ** 2 /
                       (geometry.rho * denominator) + u12),
        "z_h": q * gamma1 * exponential * h0_prime *
               lambda_value ** 2 / denominator,
    }


def hankel_field_densities(lambda_value, geometry, lower, upper,
                           frequency_hz):
    """Pointwise principal-sheet form of equations (5-5)/(5-7).

    This form is suitable for identities on the real contour.  A deformed
    contour must track both square roots continuously and call
    :func:`hankel_field_densities_with_roots` instead.
    """
    gamma1, gamma2, _, _ = _common(lambda_value, geometry, lower, upper)
    return hankel_field_densities_with_roots(
        lambda_value, geometry, lower, upper, frequency_hz, gamma1, gamma2)


def differentiated_field_densities(lambda_value, geometry, lower, upper,
                                   frequency_hz):
    """Independent equations (5-2a-e) evaluated by numerical differentiation."""
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    q = 1j * omega * MU0 / (4 * mp.pi)

    def v(rho, source_z, observer_z):
        point = CrossGeometry(rho, source_z, observer_z)
        return potential_densities(lambda_value, point, lower, upper)[0]

    u12 = potential_densities(lambda_value, geometry, lower, upper)[1]
    rho = geometry.rho
    source_z = geometry.source_z
    observer_z = geometry.observer_z
    d_rho = mp.diff(lambda r: v(r, source_z, observer_z), rho)
    d_rho_rho = mp.diff(lambda r: v(r, source_z, observer_z), rho, 2)
    d_z_z = mp.diff(lambda z: v(rho, source_z, z), observer_z, 2)
    d_rho_z = mp.diff(
        lambda r: mp.diff(lambda z: v(r, source_z, z), observer_z), rho)
    d_rho_source_z = mp.diff(
        lambda r: mp.diff(lambda z: v(r, z, observer_z), source_z), rho)
    return {
        "rho_v": -q * d_rho_z,
        "z_v": -q * (d_z_z + upper.wave_number ** 2 *
                       v(rho, source_z, observer_z)),
        "rho_h": -q * (d_rho_rho + u12),
        "phi_h": q * (d_rho / rho + u12),
        "z_h": q * d_rho_source_z,
    }


def assemble_cartesian(components, phi):
    """Assemble the Cartesian electric dyadic using equation (5-1)."""
    cosine = mp.cos(phi)
    sine = mp.sin(phi)
    rho_h = components["rho_h"]
    phi_h = components["phi_h"]
    rho_v = components["rho_v"]
    z_h = components["z_h"]
    z_v = components["z_v"]
    return mp.matrix([
        [cosine ** 2 * rho_h - sine ** 2 * phi_h,
         sine * cosine * (rho_h + phi_h), cosine * rho_v],
        [sine * cosine * (rho_h + phi_h),
         sine ** 2 * rho_h - cosine ** 2 * phi_h, sine * rho_v],
        [cosine * z_h, sine * z_h, z_v],
    ])


def homogeneous_components(frequency_hz, medium, rho, source_z, observer_z):
    """Exact finite-medium electric field in Appendix G, equations (G-1)-(G-5)."""
    rho = mp.mpf(rho)
    source_z = mp.mpf(source_z)
    observer_z = mp.mpf(observer_z)
    dz = observer_z - source_z
    distance = mp.sqrt(rho * rho + dz * dz)
    if not distance > 0:
        raise ValueError("source and observer must not coincide")
    wave_number = medium.wave_number
    omega = 2 * mp.pi * mp.mpf(frequency_hz)
    common = (-1j * omega * MU0 / (4 * mp.pi * wave_number ** 2) *
              mp.exp(-1j * wave_number * distance) / distance ** 3)
    radial = 3 + 3j * wave_number * distance - wave_number ** 2 * distance ** 2
    transverse = 1 + 1j * wave_number * distance - wave_number ** 2 * distance ** 2
    mixed = rho * dz / distance ** 2 * radial
    return {
        "rho_v": common * mixed,                                      # (G-1)
        "z_v": common * (dz ** 2 / distance ** 2 * radial - transverse),  # (G-2)
        "rho_h": common * (rho ** 2 / distance ** 2 * radial - transverse),  # (G-3)
        "phi_h": common * transverse,                                # (G-4)
        "z_h": common * mixed,                                       # (G-5)
    }


def relative_error(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))
