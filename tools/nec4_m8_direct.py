#!/usr/bin/env python3
"""Direct high-precision Sommerfeld integration for NEC4 M8.

Only the Figure 5-2 Bessel contour is enabled in this first Gate B slice.  The
documented selection rule is enforced: |z|+|z'| must exceed 0.5*rho.  Callers
cannot accidentally apply it in the Hankel-contour region.
"""

from dataclasses import dataclass

import mpmath as mp

from nec4_m8_equations import (
    CrossGeometry,
    appendix_h_remainder_field_densities,
    appendix_h_residual_singular_fields,
    appendix_h_v0_fields,
    expanded_field_densities,
    hankel_field_densities_with_roots,
    homogeneous_components,
    outgoing_gamma,
)


COMPONENTS = ("rho_v", "z_v", "rho_h", "phi_h", "z_h")


@dataclass(frozen=True)
class BesselResult:
    components: dict
    contour_break: object
    cutoff: object
    panel_width: object
    panel_count: int
    quadrature_order: int
    precision_digits: int


@dataclass(frozen=True)
class HankelResult:
    components: dict
    contour_points: tuple
    ray_parameter: object
    panel_count: int
    quadrature_order: int
    precision_digits: int
    principal_root_sign_changes: tuple


@dataclass(frozen=True)
class AppendixHResult:
    """Two-stage Appendix H evaluation and its finite interpolation value."""

    components: dict
    first_remainder: dict
    finite_remainder: dict
    panel_count: int
    quadrature_order: int
    precision_digits: int


def bessel_contour(frequency_hz, lower, upper, geometry, precision_digits=60,
                   cutoff_decay=60, panel_subdivisions=1,
                   quadrature_order=32, permit_outside_region=False):
    """Evaluate equations (5-7) on the Figure 5-2 deformed contour.

    The path is 0 -> p0 -> infinity+j*Im(p0), where
    p0=(1+j)*min(1/rho,1/(|z|+|z'|)).  The infinite leg is truncated only
    after ``cutoff_decay`` exponential decay lengths.  Independent calls with
    different precision, cutoff, panel subdivision, and quadrature order
    provide the external error estimate required by the M8 contract.
    """
    geometry.validate()
    height_sum = abs(geometry.source_z) + abs(geometry.observer_z)
    if (not permit_outside_region and
            not height_sum > mp.mpf("0.5") * geometry.rho):
        raise ValueError("Figure 5-2 Bessel contour requires |z|+|z'| > 0.5*rho")
    if precision_digits < 30:
        raise ValueError("M8 reference evaluation requires at least 30 digits")
    if cutoff_decay < 30 or panel_subdivisions < 1 or quadrature_order < 12:
        raise ValueError("insufficient tail cutoff, panel subdivision, or order")

    with mp.workdps(precision_digits):
        rho = +geometry.rho
        height_sum = +(abs(geometry.source_z) + abs(geometry.observer_z))
        scale = 1 / height_sum if rho == 0 else min(1 / rho, 1 / height_sum)
        p0 = mp.mpc(scale, scale)
        panel_width = mp.pi * scale / panel_subdivisions
        cutoff = mp.re(p0) + mp.mpf(cutoff_decay) / height_sum
        panel_count = int(mp.ceil((cutoff - mp.re(p0)) / panel_width))

        nodes, weights = mp.gauss_quadrature(quadrature_order, "legendre")
        values = {name: mp.mpc(0) for name in COMPONENTS}

        def add_segment(start, end):
            midpoint = (start + end) / 2
            half_width = (end - start) / 2
            for node, weight in zip(nodes, weights):
                lambda_value = midpoint + half_width * node
                densities = expanded_field_densities(
                    lambda_value, geometry, lower, upper, frequency_hz)
                scaled_weight = half_width * weight
                for name in COMPONENTS:
                    values[name] += scaled_weight * densities[name]

        add_segment(mp.mpc(0), p0)
        left = mp.re(p0)
        for panel in range(panel_count):
            right = min(left + panel_width, cutoff)
            add_segment(mp.mpc(left, mp.im(p0)),
                        mp.mpc(right, mp.im(p0)))
            left = right
        values = {name: +value for name, value in values.items()}

        return BesselResult(
            components=values,
            contour_break=+p0,
            cutoff=+cutoff,
            panel_width=+panel_width,
            panel_count=panel_count,
            quadrature_order=quadrature_order,
            precision_digits=precision_digits,
        )


def appendix_h_bessel_contour(frequency_hz, lower, upper, geometry,
                              precision_digits=70, cutoff_decay=80,
                              panel_subdivisions=2, quadrature_order=52):
    """Evaluate a near-singular point with both Appendix H subtractions.

    The cancellation-resistant H-1 remainder is integrated, H-2 is added back
    analytically to recover the total field, and the residual 1/R term is
    removed to expose the finite interpolation quantity described on page 129.
    """
    geometry.validate()
    if geometry.rho == 0:
        raise ValueError("Appendix H residual uses a directional rho->0 limit")
    height_sum = abs(geometry.source_z) + abs(geometry.observer_z)
    if not height_sum > mp.mpf("0.5") * geometry.rho:
        raise ValueError("Appendix H Bessel path requires |z|+|z'| > 0.5*rho")
    with mp.workdps(precision_digits):
        scale = min(1 / geometry.rho, 1 / height_sum)
        p0 = mp.mpc(scale, scale)
        panel_width = mp.pi * scale / panel_subdivisions
        cutoff = mp.re(p0) + mp.mpf(cutoff_decay) / height_sum
        panel_count = int(mp.ceil((cutoff - mp.re(p0)) / panel_width))
        nodes, weights = mp.gauss_quadrature(quadrature_order, "legendre")
        values = {name: mp.mpc(0) for name in COMPONENTS}

        def add_segment(start, end):
            midpoint = (start + end) / 2
            half_width = (end - start) / 2
            for node, weight in zip(nodes, weights):
                lambda_value = midpoint + half_width * node
                densities = appendix_h_remainder_field_densities(
                    lambda_value, geometry, lower, upper, frequency_hz)
                for name in COMPONENTS:
                    values[name] += half_width * weight * densities[name]

        add_segment(mp.mpc(0), p0)
        left = mp.re(p0)
        for _ in range(panel_count):
            right = min(left + panel_width, cutoff)
            add_segment(mp.mpc(left, mp.im(p0)),
                        mp.mpc(right, mp.im(p0)))
            left = right
        values = {name: +value for name, value in values.items()}
        singular = appendix_h_v0_fields(
            frequency_hz, lower, upper, geometry)
        residual = appendix_h_residual_singular_fields(
            frequency_hz, lower, upper, geometry)
        total = {name: +(values[name] + singular[name]) for name in COMPONENTS}
        finite = {name: +(values[name] - residual[name]) for name in COMPONENTS}
        return AppendixHResult(
            components=total,
            first_remainder=values,
            finite_remainder=finite,
            panel_count=panel_count,
            quadrature_order=quadrature_order,
            precision_digits=precision_digits,
        )


def hankel_contour(frequency_hz, lower, upper, geometry, precision_digits=70,
                   cutoff_decay=80, panel_subdivisions=3,
                   quadrature_order=48):
    """Evaluate equations (5-5)/(5-7) on the Figure 5-3a contour.

    Unlike a pointwise complex square root, ``gamma1`` and ``gamma2`` are
    continued from the left asymptotic ray through every ordered quadrature
    node.  This is essential: independently selecting the principal root
    crosses a Riemann sheet near the upper-medium branch point.
    """
    geometry.validate()
    if geometry.rho == 0:
        raise ValueError("Figure 5-3 Hankel form is singular on axis; use Bessel form")
    height_sum = abs(geometry.source_z) + abs(geometry.observer_z)
    if not height_sum <= mp.mpf("0.6") * geometry.rho:
        raise ValueError("Figure 5-3 Hankel contour is restricted to the overlap/Hankel region")
    if precision_digits < 30 or cutoff_decay < 30:
        raise ValueError("insufficient reference precision or tail cutoff")
    if panel_subdivisions < 1 or quadrature_order < 12:
        raise ValueError("insufficient panel subdivision or quadrature order")

    with mp.workdps(precision_digits):
        rho = +geometry.rho
        height_sum = +(abs(geometry.source_z) + abs(geometry.observer_z))
        radial_distance = mp.sqrt(rho ** 2 + height_sum ** 2)
        k1 = lower.wave_number
        k2 = upper.wave_number
        if abs(mp.im(k2)) > 128 * mp.eps * max(1, abs(k2)):
            raise ValueError("M8 Figure 5-3 implementation currently requires a lossless upper medium")
        k2_real = mp.re(k2)
        d_min = k2_real * min(
            mp.mpf("0.06") / (k2_real * radial_distance) + mp.mpf("0.04"),
            mp.mpf("0.8") / (1 + (geometry.observer_z - geometry.source_z) / rho))
        # Figure 5-3 limits Re(p2) to Re(k1+k2)/2.  Exact identical media
        # collapse this detour and are evaluated analytically/Bessel-wise.
        d_cap = mp.re(k1 + k2) / 2 - k2_real
        d_min = min(d_min, d_cap)
        if not d_min > 0:
            raise ValueError("Hankel branch detour collapses for identical/reversed-index media")
        p1 = k2 + (-1 + 1j) * d_min
        p2 = k2 + (1 + 1j) * d_min
        p3 = (max(mp.mpf("1.01") * mp.re(k1),
                  mp.mpf("1.2") * k2_real) +
              1j * mp.mpf("0.95") * mp.im(k1))

        slope = rho / height_sum
        left_direction = 1 + 1j * slope
        right_direction = 1 - 1j * slope
        asymptotic_decay = height_sum + rho ** 2 / height_sum
        ray_parameter = mp.mpf(cutoff_decay) / asymptotic_decay
        left_end = p1 - ray_parameter * left_direction
        right_end = p3 + ray_parameter * right_direction

        nodes, weights = mp.gauss_quadrature(quadrature_order, "legendre")
        node_weights = sorted(zip(nodes, weights), key=lambda pair: pair[0])
        values = {name: mp.mpc(0) for name in COMPONENTS}
        previous_roots = [None, None]
        sign_changes = [0, 0]
        total_panels = 0

        def continuous_root(lambda_value, wave_number, index):
            principal = mp.sqrt(lambda_value ** 2 - wave_number ** 2)
            if previous_roots[index] is None:
                selected = outgoing_gamma(lambda_value, wave_number)
            elif abs(-principal - previous_roots[index]) < abs(
                    principal - previous_roots[index]):
                selected = -principal
            else:
                selected = principal
            if selected == -principal:
                sign_changes[index] += 1
            previous_roots[index] = selected
            return selected

        def add_segment(start, end):
            nonlocal total_panels
            midpoint = (start + end) / 2
            half_width = (end - start) / 2
            for node, weight in node_weights:
                lambda_value = midpoint + half_width * node
                gamma1 = continuous_root(lambda_value, k1, 0)
                gamma2 = continuous_root(lambda_value, k2, 1)
                densities = hankel_field_densities_with_roots(
                    lambda_value, geometry, lower, upper, frequency_hz,
                    gamma1, gamma2)
                scaled_weight = half_width * weight
                for name in COMPONENTS:
                    values[name] += scaled_weight * densities[name]
            total_panels += 1

        panel_length = mp.pi / (max(rho, height_sum) * panel_subdivisions)

        def add_line(start, end, minimum_panels=1, fixed_panels=None):
            count = (fixed_panels if fixed_panels is not None else
                     max(minimum_panels,
                         int(mp.ceil(abs(end - start) / panel_length))))
            for panel in range(count):
                add_segment(start + (end - start) * panel / count,
                            start + (end - start) * (panel + 1) / count)

        # Along the asymptotic rays the phase is stationary and ``t`` is the
        # meaningful smoothness variable.  Counting panels by complex path
        # length would spuriously multiply work by sqrt(1+slope^2) for shallow
        # geometries even though the integrand is exponentially damped.
        ray_scale = max(abs(k1), abs(k2), 1 / radial_distance)
        ray_panels = max(4 * panel_subdivisions,
                         int(mp.ceil(ray_parameter * ray_scale *
                                     panel_subdivisions)))
        # The principal sqrt representation has artificial cuts on the
        # imaginary axis and on the real interval around the lossless k2
        # branch point.  The physical root is continuous, but explicitly
        # splitting at both axis crossings prevents a coarse quadrature gap
        # from making the nearest-root continuation ambiguous.
        left_breaks = [mp.mpf(0), mp.mpf(1)]
        delta = p1 - left_end
        if mp.re(delta):
            fraction = -mp.re(left_end) / mp.re(delta)
            if 0 < fraction < 1:
                left_breaks.append(fraction)
        if mp.im(delta):
            fraction = -mp.im(left_end) / mp.im(delta)
            if 0 < fraction < 1:
                left_breaks.append(fraction)
        left_breaks = sorted(set(left_breaks))
        for start_fraction, end_fraction in zip(left_breaks, left_breaks[1:]):
            fraction_span = end_fraction - start_fraction
            add_line(left_end + delta * start_fraction,
                     left_end + delta * end_fraction,
                     fixed_panels=max(2 * panel_subdivisions,
                                      int(mp.ceil(ray_panels * fraction_span))))
        add_line(p1, p2, 2 * panel_subdivisions)
        add_line(p2, p3, 2 * panel_subdivisions)
        add_line(p3, right_end, fixed_panels=ray_panels)

        return HankelResult(
            components={name: +value for name, value in values.items()},
            contour_points=(+left_end, +p1, +p2, +p3, +right_end),
            ray_parameter=+ray_parameter,
            panel_count=total_panels,
            quadrature_order=quadrature_order,
            precision_digits=precision_digits,
            principal_root_sign_changes=tuple(sign_changes),
        )


def canonical_cross(frequency_hz, lower, upper, geometry, method="auto",
                    precision_digits=55, cutoff_decay=65,
                    panel_subdivisions=2, quadrature_order=40):
    """Evaluate canonical lower-source/upper-observer components."""
    height_sum = abs(geometry.source_z) + abs(geometry.observer_z)
    if method == "auto":
        method = "bessel" if height_sum > mp.mpf("0.5") * geometry.rho else "hankel"
    if method == "bessel":
        return bessel_contour(
            frequency_hz, lower, upper, geometry,
            precision_digits=precision_digits, cutoff_decay=cutoff_decay,
            panel_subdivisions=panel_subdivisions,
            quadrature_order=quadrature_order,
            permit_outside_region=True).components
    if method == "hankel":
        return hankel_contour(
            frequency_hz, lower, upper, geometry,
            precision_digits=max(precision_digits, 55),
            cutoff_decay=max(cutoff_decay, 80),
            panel_subdivisions=max(panel_subdivisions, 3),
            quadrature_order=max(quadrature_order, 48)).components
    raise ValueError("method must be auto, bessel, or hankel")


def transformed_pairing(frequency_hz, lower, upper, rho, source_z, observer_z,
                        method="auto", **integration_options):
    """Evaluate all source/observer pairings using Appendix G (G-6)-(G-20)."""
    rho = mp.mpf(rho)
    source_z = mp.mpf(source_z)
    observer_z = mp.mpf(observer_z)
    if rho < 0 or source_z == 0 or observer_z == 0:
        raise ValueError("pairing evaluation requires rho>=0 and off-interface endpoints")

    def cross(z_source, z_observer):
        return canonical_cross(
            frequency_hz, lower, upper,
            CrossGeometry(rho, z_source, z_observer), method=method,
            **integration_options)

    if source_z < 0 < observer_z:                              # 12
        return cross(source_z, observer_z)

    if source_z > 0 > observer_z:                              # G-6..G-10
        base = cross(observer_z, source_z)
        return {
            "rho_v": -base["z_h"],
            "z_v": base["z_v"],
            "rho_h": base["rho_h"],
            "phi_h": base["phi_h"],
            "z_h": -base["rho_v"],
        }

    if source_z < 0 and observer_z < 0:                        # G-11..G-15
        a = -source_z
        b = -observer_z
        base = cross(-(a + b), mp.mpf(0))
        direct = homogeneous_components(
            frequency_hz, lower, rho, -a, -b)
        image = homogeneous_components(
            frequency_hz, lower, rho, -a, b)
        # Appendix G prints epsilon_1 factors because NEC's upper medium is
        # free space (epsilon_2=1).  Retain the dimensionally general ratio so
        # the oracle also has a valid identical-nonvacuum null test.
        epsilon_ratio = lower.epsilon_relative / upper.epsilon_relative
        return {
            "rho_v": base["rho_v"] + direct["rho_v"] - image["rho_v"],
            "z_v": base["z_v"] / epsilon_ratio + direct["z_v"] - image["z_v"],
            "rho_h": base["rho_h"] + direct["rho_h"] - image["rho_h"],
            "phi_h": base["phi_h"] + direct["phi_h"] - image["phi_h"],
            "z_h": base["z_h"] / epsilon_ratio + direct["z_h"] - image["z_h"],
        }

    if source_z > 0 and observer_z > 0:                        # G-16..G-20
        a = source_z
        b = observer_z
        base = cross(mp.mpf(0), a + b)
        direct = homogeneous_components(
            frequency_hz, upper, rho, a, b)
        image = homogeneous_components(
            frequency_hz, upper, rho, a, -b)
        epsilon_ratio = lower.epsilon_relative / upper.epsilon_relative
        return {
            "rho_v": -base["z_h"] + direct["rho_v"] - image["rho_v"],
            "z_v": epsilon_ratio * base["z_v"] + direct["z_v"] - image["z_v"],
            "rho_h": base["rho_h"] + direct["rho_h"] - image["rho_h"],
            "phi_h": base["phi_h"] + direct["phi_h"] - image["phi_h"],
            "z_h": -epsilon_ratio * base["rho_v"] + direct["z_h"] - image["z_h"],
        }

    raise ValueError("unclassified source/observer pairing")
