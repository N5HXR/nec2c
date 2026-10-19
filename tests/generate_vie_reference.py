#!/usr/bin/env python3
"""Generate M4 reference tensors with an independent high-order oracle.

This deliberately does not call the C kernel.  It evaluates the published
PWC face-pair reduction with vectorized high-order Gauss quadrature and an
independent Duffy implementation. Output is written to stdout for review;
redirect it to tests/vie-reference.csv only after the order-convergence check
passes. Singular face-pair results were also cross-checked with DIRECTFN.
"""

import itertools
import sys

import numpy as np


CASES = (
    ("self_dynamic", 1.3, (0, 0, 0), (0.2, 0.3, 0.4),
     (0, 0, 0), (0.2, 0.3, 0.4)),
    ("face", 2.1, (0, 0, 0), (0.2, 0.3, 0.4),
     (0.2, 0, 0), (0.2, 0.3, 0.4)),
    ("edge", 1.7, (0, 0, 0), (0.2, 0.3, 0.4),
     (0.2, 0.3, 0), (0.2, 0.3, 0.4)),
    ("corner", 0.9, (0, 0, 0), (0.2, 0.3, 0.4),
     (0.2, 0.3, 0.4), (0.2, 0.3, 0.4)),
    ("near_unequal", 3.2, (0.1, -0.2, 0.3), (0.3, 0.2, 0.1),
     (0.48, 0.03, 0.12), (0.2, 0.25, 0.15)),
    ("far_unequal", 3.0, (0.1, -0.2, 0.3), (0.3, 0.2, 0.1),
     (1.7, 1.1, -0.8), (0.2, 0.25, 0.15)),
)


def segments(a0, a1, av, b0, b1, bv):
    if not av and not bv:
        return (), a0 - b0
    if av and bv:
        points = {a0 - b1, a0 - b0, a1 - b1, a1 - b0}
    elif av:
        points = {a0 - b0, a1 - b0}
    else:
        points = {a0 - b1, a0 - b0}
    if min(points) < 0 < max(points):
        points.add(0.0)
    points = sorted(points)
    return tuple((points[i], points[i + 1], av and bv, a0, a1, b0, b1)
                 for i in range(len(points) - 1)), None


def segment_weight(segment, x):
    if not segment[2]:
        return np.ones_like(x)
    return np.maximum(0.0, np.minimum(segment[4], segment[6] + x) -
                      np.maximum(segment[3], segment[5] + x))


def reduced_box(selected, fixed_squared, k0, nodes, weights):
    singular = fixed_squared == 0 and all(lo == 0 or hi == 0
                                          for lo, hi, *_ in selected)
    dimension = len(selected)
    grids = np.meshgrid(*([nodes] * dimension), indexing="ij")
    weight_grids = np.meshgrid(*([weights] * dimension), indexing="ij")
    quadrature_weight = np.prod(weight_grids, axis=0)
    if singular:
        total = 0j
        for sector in range(dimension):
            extents = np.array([abs(item[1] - item[0]) for item in selected])
            signs = np.array([-1 if item[1] == 0 else 1 for item in selected])
            t = grids[sector]
            mapped = [t if axis == sector else t * grids[axis]
                      for axis in range(dimension)]
            coordinates = [signs[axis] * extents[axis] * mapped[axis]
                           for axis in range(dimension)]
            physical_weight = np.prod(
                [segment_weight(selected[axis], coordinates[axis])
                 for axis in range(dimension)], axis=0)
            radius = np.sqrt(np.sum([value * value for value in coordinates], axis=0))
            jacobian = np.prod(extents) * t ** (dimension - 1)
            integrand = (jacobian * physical_weight * np.exp(-1j * k0 * radius) /
                         (4 * np.pi * radius))
            total += np.sum(quadrature_weight * integrand)
        return total
    coordinates = [item[0] + (item[1] - item[0]) * grids[axis]
                   for axis, item in enumerate(selected)]
    jacobian = np.prod([item[1] - item[0] for item in selected])
    physical_weight = np.prod(
        [segment_weight(selected[axis], coordinates[axis])
         for axis in range(dimension)], axis=0)
    radius = np.sqrt(fixed_squared + np.sum([value * value for value in coordinates], axis=0))
    integrand = jacobian * physical_weight * np.exp(-1j * k0 * radius) / (4 * np.pi * radius)
    return np.sum(quadrature_weight * integrand)


def face_pair(tc, ts, sc, ss, ta, tsgn, sa, ssgn, k0, nodes, weights):
    lists = []
    fixed_squared = 0.0
    for axis in range(3):
        tv = axis != ta
        sv = axis != sa
        t0 = -ts[axis] / 2 if tv else tsgn * ts[axis] / 2
        t1 = ts[axis] / 2 if tv else t0
        source_center = sc[axis] - tc[axis]
        s0 = source_center - ss[axis] / 2 if sv else source_center + ssgn * ss[axis] / 2
        s1 = source_center + ss[axis] / 2 if sv else s0
        pieces, fixed = segments(t0, t1, tv, s0, s1, sv)
        if pieces:
            lists.append(pieces)
        else:
            fixed_squared += fixed * fixed
    return sum(reduced_box(selected, fixed_squared, k0, nodes, weights)
               for selected in itertools.product(*lists))


def tensor(case, order):
    _, k0, tc, ts, sc, ss = case
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights
    result = np.zeros((3, 3), dtype=np.complex128)
    volume = np.prod(ts)
    for ta, sa in itertools.product(range(3), repeat=2):
        for tsgn, ssgn in itertools.product((-1, 1), repeat=2):
            integral = face_pair(tc, ts, sc, ss, ta, tsgn, sa, ssgn,
                                 k0, nodes, weights) / volume
            if ta == sa:
                for component in range(3):
                    if component != ta:
                        result[component][component] += tsgn * ssgn * integral
            else:
                result[sa][ta] -= tsgn * ssgn * integral
    return result


def emit(case, values):
    fields = [case[0], f"{case[1]:.17g}"]
    for group in case[2:]:
        fields.extend(f"{value:.17g}" for value in group)
    for row in values:
        for value in row:
            fields.extend((f"{value.real:.17g}", f"{value.imag:.17g}"))
    print(",".join(fields))


def main():
    print("# name,k0,test-center(3),test-size(3),source-center(3),source-size(3),N00.re,N00.im,...,N22.re,N22.im")
    for case in CASES:
        lower = tensor(case, 36)
        values = tensor(case, 52)
        error = np.max(np.abs(lower - values))
        if error > 2e-11:
            raise RuntimeError(f"{case[0]} oracle did not converge: {error:.3e}")
        emit(case, values.tolist())
        print(f"generated {case[0]} (36/52 max delta {error:.3e})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
