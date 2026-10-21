#!/usr/bin/env python3
"""GN 2 wires/volume integration: compatibility, fields, power and sweep state."""
import math
from pathlib import Path
import subprocess
import tempfile

import dielectric_ground as ground
import dielectric_wire as wire

GN = 'GN 2 0 0 0 13 .02'
EX = 'EX 0 1 6 0 1 0'


def solve(d, **kwargs):
    options = dict(ground=GN, geometry=ground.WIRE, excitation=EX)
    options.update(kwargs)
    return ground.solve(d, **options)


def test_vacuum_and_solvers(d):
    vacuum = solve(d, er=1)
    deck = (d/'model.nec').read_text().replace('DF "model.diel"\n', '')
    (d/'model.nec').write_text(deck)
    result = subprocess.run([str(ground.BINARY), '-i', str(d/'model.nec'),
                             '-o', str(d/'legacy.out')], capture_output=True)
    legacy = (d/'legacy.out').read_text()
    assert result.returncode == 0, legacy
    ground.close([wire.coupled_impedance(vacuum)], [wire.legacy_impedance(legacy)], rtol=6e-5)
    for j in range(1, 12):
        ground.close([wire.coupled_current(vacuum, j)], [wire.legacy_current(legacy, j)], rtol=6e-5)
    dense = solve(d)
    iterative = solve(d, mode='gmres')
    ground.close(ground.currents(dense, 1, True)[0], ground.currents(iterative, 1, True)[0], rtol=2e-8)
    ground.close([wire.coupled_impedance(dense)], [wire.coupled_impedance(iterative)], rtol=2e-10)
    for a, b in zip(ground.rows(dense, 'COUPLED SCATTERED FAR FIELD', 6),
                    ground.rows(iterative, 'COUPLED SCATTERED FAR FIELD', 6)):
        ground.close(a[2:], b[2:], rtol=2e-8)
    defect = ground.scalar(dense, 'RELATIVE POWER DEFECT:')
    assert defect < 2e-3, defect
    assert ground.scalar(dense, 'LOWER-HALF-SPACE POWER:') > 0
    print('GN 2: vacuum recovers legacy; loaded dense/GMRES agree; driven power defect', defect, flush=True)


def test_plane_waves(d):
    for kind in (1, 2, 3):
        ex = 'EX %d 1 1 0 37 23 19 0 0 .6' % kind
        volume = solve(d, geometry='', excitation=ex)
        defect = ground.scalar(volume, 'RELATIVE IMBALANCE:')
        assert defect < 2e-4, defect
        coupled = solve(d, excitation=ex)
        defect = ground.scalar(coupled, 'RELATIVE POWER DEFECT:')
        assert defect < 5e-3, defect
        print('GN 2 EX', kind, 'volume and coupled extinction/air/ground power pass; coupled defect', defect, flush=True)
    lossless = solve(d, geometry='', ground='GN 2 0 0 0 4 0', excitation=ex)
    assert ground.scalar(lossless, 'LOWER-HALF-SPACE SCATTERED POWER:') > 0
    assert ground.scalar(lossless, 'RELATIVE IMBALANCE:') < 2e-4


def test_matched_and_sweeps(d):
    for geometry in ('', ground.WIRE):
        ex = EX if geometry else 'EX 1 1 1 0 37 23 19 0 0 .6'
        matched = solve(d, geometry=geometry, excitation=ex, ground='GN 2 0 0 0 1 0')
        free = solve(d, geometry=geometry, excitation=ex, ground='')
        ground.close(ground.currents(matched, 1, bool(geometry))[0],
                     ground.currents(free, 1, bool(geometry))[0], rtol=2e-10)
        if geometry:
            ground.close([wire.coupled_impedance(matched)], [wire.coupled_impedance(free)], rtol=2e-10)
    # Negative GN sigma specifies imaginary permittivity, not conductivity.
    for gn in (GN, 'GN 2 0 0 0 13 -.4'):
        sweep = solve(d, ground=gn, frequencies=2, frequency=280, step=40)
        parts = sweep.split('---------- COUPLED CONDUCTOR-DIELECTRIC SOLUTION (M7) ----------')[1:]
        assert len(parts) == 2, sweep
        for frequency, part in zip((280, 320), parts):
            separate = solve(d, ground=gn, frequency=frequency)
            ground.close(ground.currents(part, 1, True)[0], ground.currents(separate, 1, True)[0], rtol=2e-10)
    print('GN 2: matched half-space recovers free space; conductivity and imaginary-epsilon sweeps restore state', flush=True)


def test_maxwell_near_fields(d):
    point = [.4, .2, .5]
    step = 1e-4
    derivative = []
    for axis in range(3):
        start = point.copy()
        start[axis] -= step
        counts = [1, 1, 1]
        counts[axis] = 3
        spacing = [0, 0, 0]
        spacing[axis] = step
        card = 'NE 0 ' + ' '.join(map(str, counts + start + spacing))
        text = solve(d, output=card)
        values = ground.rows(text, 'DIELECTRIC SCATTERED NEAR ELECTRIC', 3)
        fields = [[complex(*row[j:j+2]) for j in (3, 5, 7)] for row in values]
        derivative.append([(b-a)/(2*step) for a, b in zip(fields[0], fields[2])])
    curl = [derivative[1][2]-derivative[2][1], derivative[2][0]-derivative[0][2],
            derivative[0][1]-derivative[1][0]]
    magnetic = ground.field(solve(d, output='NH 0 1 1 1 .4 .2 .5 0 0 0'), True)
    from_curl = [-v/(1j*2*math.pi*300e6*1.25663706212e-6) for v in curl]
    ground.close(magnetic, from_curl, rtol=3e-4, atol=1e-10)
    print('GN 2: total scattered E/H satisfy curl E = -i omega mu H', flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix='nec2c-m10-real-') as tmp:
        d = Path(tmp)
        test_vacuum_and_solvers(d)
        test_plane_waves(d)
        test_matched_and_sweeps(d)
        test_maxwell_near_fields(d)


if __name__ == '__main__':
    main()
