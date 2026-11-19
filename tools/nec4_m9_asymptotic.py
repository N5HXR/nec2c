#!/usr/bin/env python3
"""Documented NEC4 first-order cross-interface asymptotic fields.

The implementation is an equation-level diagnostic for M9 Gate D.  It covers
the generic nonzero-rho saddle contribution in equations (5-21)--(5-23).  Pole
and second-order corrections are intentionally separate so their contribution
can be measured rather than folded into an opaque expression.
"""

import cmath
import math

import mpmath as mp
from scipy import special

from nec4_m8_equations import MU0
from nec4_m9_tables import _beta, saddle_lambda1


COMPONENTS = ("rho_v", "z_v", "rho_h", "phi_h", "z_h")


def outgoing_gamma(value, wave_number):
    """Return Gamma=(lambda^2-k^2)^1/2 on the continued outgoing sheet.

    Selecting each complex square root by its instantaneous real part causes
    a sheet jump when a lossy saddle moves just below the real lambda axis.
    ``beta`` is continued from sqrt(k^2-lambda^2), so Gamma=j*beta preserves
    the outgoing real-axis branch through that displacement.
    """
    return 1j * _beta(wave_number, value)


def first_order_cross(frequency_hz, k_lower, k_upper,
                      rho, observer_z, source_depth):
    """Evaluate equations (5-21)--(5-23) away from limiting coordinates.

    Coordinates follow the guide: lower-medium source ``z'=-source_depth``
    and upper-medium observer ``z=observer_z``.  The limiting forms printed
    below equation (5-22) are deferred until the generic branch is validated.
    """
    if frequency_hz <= 0 or min(rho, observer_z, source_depth) < 0:
        raise ValueError("positive frequency and nonnegative coordinates required")
    if rho == 0:
        raise ValueError("generic first-order formula requires rho > 0")
    source_z = -source_depth
    saddle = saddle_lambda1(rho, observer_z, source_depth,
                            k_lower, k_upper)
    gamma1 = outgoing_gamma(saddle, k_lower)
    gamma2 = outgoing_gamma(saddle, k_upper)
    if saddle == 0 or gamma1 == 0 or gamma2 == 0:
        raise ValueError("generic first-order formula reached a limiting saddle")
    denominator = k_lower * k_lower * gamma2 + k_upper * k_upper * gamma1
    curvature = (k_lower * k_lower * gamma2 * gamma2 * source_z / gamma1
                 - k_upper * k_upper * gamma1 * gamma1 * observer_z / gamma2)
    # Equation (5-22) writes sqrt(rho/lambda_s) in the denominator;
    # equivalently the original (5-20) form carries sqrt(lambda_s/rho)
    # in the numerator.  Keeping the printed (5-22) arrangement avoids a
    # three-order scale error at wavelength-sized radial separations.
    # The sign of F''**(-1/2) follows the tangent of the deformed SDP, not the
    # principal square root.  For the outgoing roots above that continuation
    # is the negative principal root; the lossless far-ray limit fixes it.
    spreading = -cmath.sqrt(rho / saddle) * cmath.sqrt(curvature)
    psi = (gamma1 * gamma2
           * cmath.exp(gamma1 * source_z - gamma2 * observer_z
                       - 1j * saddle * rho) / spreading)
    omega = 2 * math.pi * frequency_hz
    scale = omega * float(MU0) / (2 * math.pi)
    plus = cmath.exp(1j * math.pi / 4) * scale
    minus = cmath.exp(-1j * math.pi / 4) * scale
    return {
        "rho_v": plus * saddle * gamma2 / denominator * psi,
        "z_v": minus * saddle ** 3 / denominator * psi,
        "rho_h": minus * (-saddle ** 2 / denominator
                           + 1 / (gamma1 + gamma2)) * psi,
        "phi_h": minus * (1j * saddle / rho / denominator
                           - 1 / (gamma1 + gamma2)) * psi,
        "z_h": plus * saddle * gamma1 / denominator * psi,
    }


def q_function(value):
    """Modified saddle-point Q(z) from equation (5-25)."""
    return -math.sqrt(math.pi) * value * cmath.exp(value * value) * special.erfc(value)


def _pole_geometry(k_lower, k_upper, rho, observer_z, source_depth):
    if rho <= 0:
        raise ValueError("pole formula requires rho > 0")
    source_z = -source_depth
    saddle = saddle_lambda1(rho, observer_z, source_depth,
                            k_lower, k_upper)
    sum_squared = k_lower * k_lower + k_upper * k_upper
    pole = k_lower * k_upper / cmath.sqrt(sum_squared)

    def phase(value):
        return (outgoing_gamma(value, k_upper) * observer_z
                - outgoing_gamma(value, k_lower) * source_z
                + 1j * value * rho)

    saddle_phase = phase(saddle)
    difference = phase(pole) - saddle_phase
    b = cmath.sqrt(difference)
    # The modified-SDP mapping chooses b continuously from the saddle toward
    # the pole; the guide's Q(j*b) form uses the branch with Im(b) <= 0.
    if b.imag > 0 or (abs(b.imag) < 1e-14 * max(1.0, abs(b)) and b.real < 0):
        b = -b
    if abs(b) == 0:
        raise ValueError("pole and saddle coalescence requires a limiting form")
    pole_gamma1 = outgoing_gamma(pole, k_lower)
    pole_gamma2 = outgoing_gamma(pole, k_upper)
    denominator_derivative = (k_lower * k_lower * pole / pole_gamma2
                              + k_upper * k_upper * pole / pole_gamma1)
    # Evaluate the defining residue instead of separately taking the many
    # fractional powers in the printed closed form; this keeps all roots on
    # the same outgoing sheet as F(lambda).
    a = (cmath.sqrt(2j / (math.pi * rho)) * cmath.sqrt(pole)
         / denominator_derivative)
    return pole, pole_gamma1, pole_gamma2, a, b, saddle_phase


def first_order_pole(frequency_hz, k_lower, k_upper,
                     rho, observer_z, source_depth):
    """Surface-wave/pole correction in equations (5-25)--(5-27)."""
    pole, _, _, a, b, saddle_phase = _pole_geometry(
        k_lower, k_upper, rho, observer_z, source_depth)
    sum_squared = k_lower * k_lower + k_upper * k_upper
    potential = (math.sqrt(math.pi) * a / b * cmath.exp(-saddle_phase)
                 * (q_function(1j * b) + 1))
    omega = 2 * math.pi * frequency_hz
    common = omega * float(MU0) / (4 * math.pi)
    return {
        "rho_v": -1j * common * k_lower * k_upper ** 3 / sum_squared * potential,
        "z_v": 1j * common * k_lower ** 2 * k_upper ** 2 / sum_squared * potential,
        "rho_h": 1j * common * k_lower ** 2 * k_upper ** 2 / sum_squared * potential,
        "phi_h": common / rho * k_lower * k_upper / cmath.sqrt(sum_squared) * potential,
        "z_h": -1j * common * k_lower ** 3 * k_upper / sum_squared * potential,
    }


def second_order_pole(frequency_hz, k_lower, k_upper,
                      rho, observer_z, source_depth):
    """Second-order uniform correction in equations (5-26)/(5-28)."""
    pole, gamma1, gamma2, a, b, saddle_phase = _pole_geometry(
        k_lower, k_upper, rho, observer_z, source_depth)
    potential = (math.sqrt(math.pi) * a / b * cmath.exp(-saddle_phase)
                 * (q_function(1j * b) + 1 + 1 / (2 * b * b)))
    omega = 2 * math.pi * frequency_hz
    common = omega * float(MU0) / (4 * math.pi)
    radial_first = pole - 3j / (8 * rho)
    return {
        "rho_v": common * radial_first * gamma2 * potential,
        "z_v": -1j * common * (pole * pole + 1j * pole / (8 * rho)) * potential,
        "rho_h": 1j * common * (pole * pole - 7j * pole / (8 * rho)) * potential,
        "phi_h": common / rho * radial_first * potential,
        "z_h": common * radial_first * gamma1 * potential,
    }


def first_order_uniform_cross(frequency_hz, k_lower, k_upper,
                              rho, observer_z, source_depth):
    saddle = first_order_cross(frequency_hz, k_lower, k_upper,
                               rho, observer_z, source_depth)
    pole = first_order_pole(frequency_hz, k_lower, k_upper,
                            rho, observer_z, source_depth)
    return {name: saddle[name] + pole[name] for name in COMPONENTS}


def second_order_cross_reference(frequency_hz, k_lower, k_upper,
                                 rho, observer_z, source_depth,
                                 precision_digits=60):
    """High-precision Appendix-J saddle contribution before the pole term.

    This evaluates the universal equations (J-1)/(J-2) by differentiating the
    component amplitudes formed with the printed two-term Hankel expansions.
    It is intentionally an oracle; optimized production code will use the
    algebraic J-9--J-13 forms only after matching this implementation.
    """
    if rho <= 0:
        raise ValueError("second-order saddle formula requires rho > 0")
    with mp.workdps(precision_digits):
        k1 = mp.mpc(k_lower)
        k2 = mp.mpc(k_upper)
        radial = mp.mpf(rho)
        z = mp.mpf(observer_z)
        source_z = -mp.mpf(source_depth)
        saddle = mp.mpc(saddle_lambda1(
            rho, observer_z, source_depth, k_lower, k_upper))

        def gamma(value, wave_number):
            beta = mp.sqrt(wave_number * wave_number - value * value)
            if (mp.re(beta) < 0 or
                    (abs(mp.re(beta)) < mp.mpf("1e-40") and mp.im(beta) > 0)):
                beta = -beta
            return 1j * beta

        def phase(value):
            return (gamma(value, k2) * z - gamma(value, k1) * source_z
                    + 1j * value * radial)

        def amplitude(value, component, two_term=True):
            gamma1 = gamma(value, k1)
            gamma2 = gamma(value, k2)
            denominator = k1 * k1 * gamma2 + k2 * k2 * gamma1
            root = mp.sqrt(2 / (mp.pi * value * radial))
            correction = 1 / (8 * value * radial) if two_term else 0
            h0 = root * (1 + 1j * correction) * mp.exp(1j * mp.pi / 4)
            h0_prime = root * (1 - 3j * correction) * mp.exp(-1j * mp.pi / 4)
            h0_second = -root * (1 - 7j * correction) * mp.exp(1j * mp.pi / 4)
            omega = 2 * mp.pi * mp.mpf(frequency_hz)
            q = 1j * omega * MU0 / (4 * mp.pi)
            u = h0 * value / (gamma1 + gamma2)
            if component == "rho_v":
                return q * gamma2 * h0_prime * value ** 2 / denominator
            if component == "z_v":
                return -q * h0 * value ** 3 / denominator
            if component == "rho_h":
                return -q * (h0_second * value ** 3 / denominator + u)
            if component == "phi_h":
                return q * (h0_prime * value ** 2 /
                            (radial * denominator) + u)
            if component == "z_h":
                return q * gamma1 * h0_prime * value ** 2 / denominator
            raise ValueError(component)

        f0 = phase(saddle)
        f2 = mp.diff(phase, saddle, 2)
        f3 = mp.diff(phase, saddle, 3)
        f4 = mp.diff(phase, saddle, 4)
        prefactor = mp.sqrt(2 * mp.pi / f2) * mp.exp(-f0)
        # Select the SDP square-root continuation by matching the leading
        # term to the independently implemented equations (5-21)--(5-23).
        first = first_order_cross(frequency_hz, k_lower, k_upper,
                                  rho, observer_z, source_depth)
        probe = prefactor * amplitude(saddle, "rho_v", two_term=False)
        if abs(complex(probe) - first["rho_v"]) > abs(complex(-probe) - first["rho_v"]):
            prefactor = -prefactor
        result = {}
        for component in COMPONENTS:
            function = lambda value: amplitude(value, component, two_term=True)
            g0 = function(saddle)
            g1 = mp.diff(function, saddle, 1)
            g2 = mp.diff(function, saddle, 2)
            q2 = (g0 * (5 * f3 * f3 - 3 * f2 * f4)
                  - 12 * g1 * f2 * f3 + 12 * g2 * f2 * f2) / (24 * f2 ** 3)
            result[component] = complex(prefactor * (g0 + q2))
        return result


def second_order_uniform_cross_reference(frequency_hz, k_lower, k_upper,
                                         rho, observer_z, source_depth,
                                         precision_digits=60):
    saddle = second_order_cross_reference(
        frequency_hz, k_lower, k_upper, rho, observer_z, source_depth,
        precision_digits=precision_digits)
    pole = second_order_pole(frequency_hz, k_lower, k_upper,
                             rho, observer_z, source_depth)
    return {name: saddle[name] + pole[name] for name in COMPONENTS}
