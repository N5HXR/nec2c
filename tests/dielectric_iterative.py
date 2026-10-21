#!/usr/bin/env python3
"""M9: compare iterative currents and fields with the dense solve on one mesh."""
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(os.environ.get('TEST_SRCDIR', '.')).resolve()
BINARY = Path(os.environ.get('NEC2C', './nec2c')).resolve()


def solve(directory, deck, mode):
    # Strict per-component overlap checks need a tighter target than routine runs.
    result = subprocess.run([str(BINARY), '-s', mode, '-r', '1e-13', '-i', str(deck),
                             '-o', str(directory/'result.out')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    text = (directory/'result.out').read_text()
    assert result.returncode == 0, text+result.stderr.decode()
    if mode in ('gmres', 'fft'):
        assert ('FFT-OPERATOR GMRES' if mode == 'fft' else 'DENSE-OPERATOR GMRES (EXPERIMENTAL; NO ACCELERATION)') in text
        residual = float(re.search(r'GMRES TRUE RELATIVE RESIDUAL:\s*(\S+)', text)[1])
        assert residual <= 1e-13, residual
    return text


def numeric_rows(text):
    section = text.split('SOLUTION', 1)[1]
    rows = []
    for line in section.splitlines():
        if not re.match(r'^\s*[-+0-9]', line):
            continue
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            continue
    return rows


def main():
    # Validate the sizing inventory against a fully occupied 2x2x2 grid and
    # disjoint cells separated by holes; invalid geometry must fail explicitly.
    sys.path.insert(0, str(ROOT/'tools'))
    from m9_tree_sizing import inventory, TREE
    cube = inventory([(0, 0, 0, 2, 2, 2)], 1.0, leaf_capacity=1)
    assert cube['occupied_cells'] == 8
    assert cube['octree_nodes'] == 9 and cube['octree_leaves'] == 8
    assert cube['padded_fft_cells'] == 64
    sparse = inventory([(0, 0, 0, 1, 1, 1), (3, 0, 0, 4, 1, 1)], 1.0)
    assert sparse['occupied_cells'] == 2 and sparse['occupied_fraction'] == 0.5
    assert inventory(TREE, 0.25)['occupied_cells'] == 9432
    for boxes, pitch in ((TREE, 0.3), (TREE, float('nan')),
                         ([(0, 0, 0, 2, 2, 2), (1, 1, 1, 3, 3, 3)], 1.0)):
        try:
            inventory(boxes, pitch)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid sizing geometry accepted')
    cases = [
        ('m5-cube.nec', 'm5-cube.diel', None),
        ('m6-cube.nec', 'm6-cube.diel', None),
        ('m7-patch-cube-coarse.nec', 'm7-cube-coarse.diel',
         'BOX 1 0.075 -0.025 0.345 0.125 0.025 0.395'),
        ('m7-mixed-cube.nec', 'm7-mixed-cube.diel',
         'BOX 1 -0.125 0.325 -0.025 -0.075 0.375 0.025'),
    ]
    with tempfile.TemporaryDirectory(prefix='nec2c-m9-comparison-') as temporary:
        directory = Path(temporary)
        for name, resource_name, box in cases:
            resource = (ROOT/'tests/dielectric'/resource_name).read_text()
            if box:
                resource = re.sub(r'^BOX .*$', box, resource, flags=re.M)
            # Eight equal cubes exercise actual cross-cell convolution and coupling.
            resource = re.sub(r'MAX_KH \S+', 'MAX_KH '+('0.05' if name.startswith('m5') else '0.30'), resource)
            (directory/resource_name).write_text(resource)
            source = (ROOT/'tests/dielectric'/name).read_text()
            # Principal-cut complex fields exercise the complete solved vector.
            source = source.split('RP ')[0] if 'RP ' in source else source.split('XQ')[0]
            source += 'RP 0 5 2 0000 0 0 45 90 0 0\nEN\n'
            deck = directory/'model.nec'
            deck.write_text(source)
            dense = solve(directory, deck, 'dense')
            for mode in ('gmres', 'fft'):
                iterative = solve(directory, deck, mode)
                compare(name+'/'+mode, dense, iterative)
            print(name+': dense/GMRES/FFT currents and complex pattern agree')
        # Solver selection must leave a no-DF legacy deck on its existing path.
        legacy_results = []
        from dielectric_wire import legacy_current, legacy_impedance
        for mode in ('dense', 'gmres', 'fft'):
            subprocess.run([str(BINARY), '-s', mode, '-i',
                            str(ROOT/'tests/dielectric/m6-legacy.nec'),
                            '-o', str(directory/'legacy.out')], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            text = (directory/'legacy.out').read_text()
            assert 'OPERATOR GMRES' not in text
            legacy_results.append((legacy_current(text, 6), legacy_impedance(text)))
        assert all(item == legacy_results[0] for item in legacy_results)
        for arguments in (['-s', 'unknown'], ['-r', 'nan'], ['-r', '0'],
                          ['-n', '-1'], ['-n', '0']):
            bad = subprocess.run([str(BINARY)]+arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert bad.returncode != 0
        # An iteration limit must fail before solved currents or patterns appear.
        failed = subprocess.run([str(BINARY), '-s', 'fft', '-n', '1',
                                 '-r', '1e-15', '-i', str(deck),
                                 '-o', str(directory/'failed.out')],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        failure = (directory/'failed.out').read_text()
        assert failed.returncode != 0 and 'iteration limit' in failure
        assert 'COUPLED CONDUCTOR-DIELECTRIC SOLUTION' not in failure

    check_vacuum()
    print('M9 dense/iterative/FFT comparison passed')


def check_vacuum():
    with tempfile.TemporaryDirectory(prefix='nec2c-m9-vacuum-') as temporary:
        directory = Path(temporary)
        for name in ('m6-vacuum.nec', 'm7-patch-vacuum.nec', 'm7-mixed-vacuum.nec'):
            deck = ROOT/'tests/dielectric'/name
            compare(name, solve(directory, deck, 'dense'), solve(directory, deck, 'fft'))
    print('FFT zero-volume wire/patch/mixed identities passed', flush=True)


def compare(name,dense,iterative):
    a, b = numeric_rows(dense), numeric_rows(iterative)
    assert a and len(a) == len(b), (name, len(a), len(b))
    for left, right in zip(a, b):
        assert len(left) == len(right)
        for expected, actual in zip(left, right):
            assert math.isclose(actual, expected, rel_tol=1e-8, abs_tol=1e-12), (
                name, expected, actual)
    for label in ('RELATIVE POWER DEFECT:', 'RADIATED POWER:'):
        if label in dense:
            expected = float(re.search(re.escape(label)+r'\s*(\S+)', dense)[1])
            actual = float(re.search(re.escape(label)+r'\s*(\S+)', iterative)[1])
            assert math.isclose(actual, expected, rel_tol=1e-8, abs_tol=1e-12)

if __name__ == '__main__':
    main()
