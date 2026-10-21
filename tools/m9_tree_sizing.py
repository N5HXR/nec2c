#!/usr/bin/env python3
"""Geometric/cost inventory for M9; estimates only, never allocates a matrix.

Compare occupied uniform-grid cells, padded FFT embedding, and occupied octree
nodes. The latter counts geometry, not low-rank storage or achieved accuracy.
All lengths are metres. The idealized tree is explicitly suspended in free space.
"""
import argparse
import csv
import math
from pathlib import Path
import sys

# Disjoint boxes: trunk, two branches, and canopy. All coordinates align to 0.25 m.
TREE = (
    (-0.25, -0.25, 0.0, 0.25, 0.25, 8.0),
    (0.25, -0.25, 5.5, 3.0, 0.25, 6.0),
    (-0.25, 0.25, 6.5, 0.25, 3.0, 7.0),
    (-3.0, -3.0, 8.0, 3.0, 3.0, 12.0),
)


def product(values):
    result = 1
    for value in values:
        result *= value
    return result


def power_two(n):
    return 1 << (n-1).bit_length()


def inventory(boxes, pitch, leaf_capacity=32):
    if not math.isfinite(pitch) or pitch <= 0 or leaf_capacity <= 0:
        raise ValueError('pitch and leaf capacity must be positive')
    integer_boxes = []
    for box in boxes:
        scaled = [coordinate/pitch for coordinate in box]
        if any(not math.isfinite(v) or abs(v-round(v)) > 1e-8 for v in scaled):
            raise ValueError('all box bounds must align to the requested uniform pitch')
        bounds = tuple(round(v) for v in scaled)
        if any(bounds[d] >= bounds[d+3] for d in range(3)):
            raise ValueError('every box must have positive dimensions')
        for other in integer_boxes:
            if all(min(bounds[d+3], other[d+3]) > max(bounds[d], other[d]) for d in range(3)):
                raise ValueError('box interiors must not overlap')
        integer_boxes.append(bounds)
    if not integer_boxes:
        raise ValueError('at least one physical box is required')
    lo = [min(box[d] for box in integer_boxes) for d in range(3)]
    hi = [max(box[d+3] for box in integer_boxes) for d in range(3)]
    shape = [hi[d]-lo[d] for d in range(3)]
    cells = sum(product(box[d+3]-box[d] for d in range(3)) for box in integer_boxes)
    nodes = leaves = depth = 0

    def visit(origin, side, level):
        nonlocal nodes, leaves, depth
        occupied = sum(product(max(0, min(origin[d]+side, box[d+3])-
                                           max(origin[d], box[d])) for d in range(3))
                       for box in integer_boxes)
        if occupied == 0:
            return
        nodes += 1
        depth = max(depth, level)
        if occupied <= leaf_capacity or side == 1:
            leaves += 1
            return
        half = side//2
        for mask in range(8):
            visit([origin[d]+(half if mask & (1 << d) else 0) for d in range(3)], half, level+1)

    visit(lo, power_two(max(shape)), 0)
    padded = [power_two(2*n-1) for n in shape]
    padded_cells = product(padded)
    order = 3*cells
    restart = 30
    # 6 symmetric tensor spectra + 3 input and 3 output complex grids.
    # FFT-library workspace, conductor coupling, maps, and allocator overhead excluded.
    fft_array_bytes = 12*16*padded_cells
    krylov_bytes = 16*((2*restart+3)*order+(restart+1)*restart+3*restart+1)+8*restart
    return dict(
        pitch_m=pitch, occupied_cells=cells, volume_unknowns=order,
        grid_x=shape[0], grid_y=shape[1], grid_z=shape[2],
        occupied_fraction=cells/product(shape),
        fft_x=padded[0], fft_y=padded[1], fft_z=padded[2],
        padded_fft_cells=padded_cells, fft_array_bytes=fft_array_bytes,
        gmres30_workspace_bytes=krylov_bytes,
        dense_one_matrix_bytes=16*order*order,
        dense_two_matrix_bytes=32*order*order,
        occupied_directed_cell_pairs=cells*cells,
        translation_table_upper_bound=product(2*n-1 for n in shape),
        octree_nodes=nodes, octree_leaves=leaves, octree_depth=depth,
        octree_leaf_capacity=leaf_capacity,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pitch', type=float, nargs='+', default=[0.25, 0.125, 0.0625])
    parser.add_argument('--frequency-mhz', type=float, default=28.0)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--geometry', choices=('full', 'skeleton', 'both'), default='both')
    args = parser.parse_args()
    if not math.isfinite(args.frequency_mhz) or args.frequency_mhz <= 0:
        parser.error('frequency must be finite and positive')
    rows = []
    try:
        for geometry in ('full', 'skeleton') if args.geometry == 'both' else (args.geometry,):
            for pitch in args.pitch:
                row = dict(geometry=geometry)
                row.update(inventory(TREE if geometry == 'full' else TREE[:3], pitch))
                row.update(frequency_mhz=args.frequency_mhz, epsilon_r=3.0,
                           max_kh=2*math.pi*args.frequency_mhz*1e6/299792458*math.sqrt(3)*pitch)
                rows.append(row)
    except ValueError as error:
        parser.error(str(error))
    stream = args.output.open('w', newline='') if args.output else sys.stdout
    try:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.output:
            stream.close()


if __name__ == '__main__':
    main()
