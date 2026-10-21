#!/usr/bin/env python3
"""M8 analytical field validation; no external EM solver or NumPy required."""
import cmath
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(os.environ.get('TEST_SRCDIR', '.')).resolve()
BINARY = Path(os.environ.get('NEC2C', './nec2c')).resolve()
EPS = 8.8541878128e-12
C0 = 299792458.0
K = 2*math.pi*100e6/C0


def run(directory, resource, output_card, conductor=False, success=True):
    (directory/'model.diel').write_text(resource)
    geometry = 'GW 1 11 0 0 -0.24 0 0 0.24 0.001\n' if conductor else ''
    excitation = 'EX 0 1 6 0 1 0' if conductor else 'EX 1 1 1 0 90 0 90 0 0 0'
    (directory/'model.nec').write_text(
        'CE\n'+geometry+'GE\nFR 0 1 0 0 100\n'+excitation+
        '\nDF "model.diel"\n'+output_card+'\nEN\n')
    result = subprocess.run([str(BINARY), '-i', str(directory/'model.nec'),
                             '-o', str(directory/'model.out')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    text = (directory/'model.out').read_text()
    assert (result.returncode == 0) == success, text
    return text


def row(text, heading):
    section = text.split(heading, 1)[1]
    for line in section.splitlines()[2:]:
        values = line.split()
        try:
            return [float(x) for x in values]
        except ValueError:
            continue
    raise AssertionError('missing data: '+heading)


def near(text, magnetic=False):
    values = row(text, 'DIELECTRIC SCATTERED NEAR '+('MAGNETIC' if magnetic else 'ELECTRIC'))
    return [complex(values[i], values[i+1]) for i in (3, 5, 7)]


def close(actual, expected, tolerance):
    error = math.sqrt(sum(abs(a-b)**2 for a, b in zip(actual, expected)))
    scale = math.sqrt(sum(abs(b)**2 for b in expected))
    assert error <= tolerance*max(scale, 1e-30), (actual, expected, error/scale)


def main():
    resource = (ROOT/'tests/dielectric/m5-cube.diel').read_text()
    with tempfile.TemporaryDirectory(prefix='nec2c-m8-fields-') as temporary:
        directory = Path(temporary)
        output = run(directory, resource, 'RP 0 1 1 0000 90 0 0 0 0 0')
        values = row(output, 'TOTAL ELECTRIC DIPOLE MOMENT (C M)')
        dipole = [complex(values[i], values[i+1]) for i in (0, 2, 4)]
        values = row(output, 'DIELECTRIC SCATTERED FAR FIELD (R TIMES E, V)')
        far = [complex(values[i], values[i+1]) for i in (2, 4)]
        # Observation +x: theta basis is -z and phi basis is +y.
        exact = [-K*K*dipole[2]/(4*math.pi*EPS), K*K*dipole[1]/(4*math.pi*EPS)]
        # The finite constant-current cube has an exact sinc form factor.
        sinc = math.sin(K*0.01)/(K*0.01)
        close(far, [sinc*v for v in exact], 1e-4)
        for distance in (1.0, 10.0, 1000.0):
            electric = near(run(directory, resource, f'NE 0 1 1 1 {distance} 0 0 0 0 0'))
            magnetic = near(run(directory, resource, f'NH 0 1 1 1 {distance} 0 0 0 0 0 0'), True)
            phase = cmath.exp(-1j*K*distance)/(4*math.pi*EPS)
            expected = [phase*2*dipole[0]*(1/distance**3+1j*K/distance**2)]
            expected += [phase*p*(K*K/distance-1/distance**3-1j*K/distance**2) for p in dipole[1:]]
            close(electric, expected, 4e-4)
            factor = (2*math.pi*100e6)*(K-1j/distance)*cmath.exp(-1j*K*distance)/(4*math.pi*distance)
            close(magnetic, [0, -factor*dipole[2], factor*dipole[1]], 4e-4)
            if distance == 1000:
                close([-electric[2]*distance*cmath.exp(1j*K*distance),
                       electric[1]*distance*cmath.exp(1j*K*distance)], far, 6e-4)
        for card, magnetic in (('NE', False), ('NH', True)):
            cartesian = near(run(directory, resource, f'{card} 0 1 1 1 1 0 0 0 0 0'), magnetic)
            spherical = near(run(directory, resource, f'{card} 1 1 1 1 1 0 90 0 0 0'), magnetic)
            close(spherical, cartesian, 1e-10)
        vacuum = resource.replace('ER 3.0', 'ER 1.0')
        assert near(run(directory, vacuum, 'NE 0 1 1 1 1 0 0 0 0 0')) == [0j]*3
        rejected = run(directory, resource, 'NE 0 1 1 1 0 0 0 0 0 0', success=False)
        assert 'one cell longest-edge clearance' in rejected
        # The audited network paths remain explicitly gated even in vacuum.
        for card in ('NT 1 3 1 9 0.02 0 0.01 0 0.02 0',
                     'TL 1 3 1 9 50 0.2 0 0 0 0'):
            rejected = run(directory, vacuum, card+'\nXQ', conductor=True, success=False)
            assert 'NETWORK AND TRANSMISSION-LINE CARDS ARE NOT YET SUPPORTED' in rejected
        rejected = run(directory, vacuum, 'RP 1 1 1 0000 90 0 0 0 0 0',
                       conductor=True, success=False)
        assert 'ONLY RP MODE 0 IS SUPPORTED' in rejected
        # Deferred multi-frequency NE requests must survive until EN.
        run(directory, resource, 'NE 0 1 1 1 1 0 0 0 0 0')
        deck = directory/'model.nec'
        deck.write_text(deck.read_text().replace('FR 0 1 0 0 100', 'FR 0 2 0 0 100 10'))
        subprocess.run([str(BINARY), '-i', str(deck), '-o', str(directory/'model.out')],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (directory/'model.out').read_text().count('SCATTERED NEAR ELECTRIC FIELD') == 2

    # Vacuum identity against the independent legacy output dispatch. Both
    # wire and patch currents must reach ffld/nefld/nhfld with their NEC units.
    with tempfile.TemporaryDirectory(prefix='nec2c-m8-vacuum-') as temporary:
        directory = Path(temporary)
        for fixture in ('m6-legacy.nec', 'm7-patch-legacy.nec', 'm7-mixed-legacy.nec'):
            source = (ROOT/'tests/dielectric'/fixture).read_text()
            prefix = source.split('XQ')[0]
            vacuum_path = directory/'vacuum.diel'
            vacuum_path.write_text('NEC2C-DIELECTRIC 1\nMATERIAL 1 ER 1\nREGION 1 MATERIAL 1\nBOX 1 2 2 2 2.1 2.1 2.1\nEND\n')
            for card, heading, magnetic in (
                    ('RP 0 1 1 0000 60 30 0 0 0 0', 'RADIATION PATTERNS', False),
                    ('NE 0 1 1 1 1 1 1 0 0 0', 'NEAR ELECTRIC FIELDS', False),
                    ('NH 0 1 1 1 1 1 1 0 0 0', 'NEAR MAGNETIC FIELDS', True)):
                results = []
                for df in ('', 'DF "vacuum.diel"\n'):
                    (directory/'case.nec').write_text(prefix+df+card+'\nEN\n')
                    subprocess.run([str(BINARY), '-i', str(directory/'case.nec'),
                                    '-o', str(directory/'case.out')], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    results.append((directory/'case.out').read_text())
                if card.startswith('RP'):
                    # Legacy prints magnitudes/phases at SI r=1 m when RP distance is zero.
                    section = results[0].split('RADIATION PATTERNS', 1)[1]
                    line = next(line for line in section.splitlines() if 'LINEAR' in line or 'RIGHT' in line or 'LEFT' in line)
                    values = [float(v) for v in line.split()[-4:]]
                    expected = [cmath.rect(values[i], math.radians(values[i+1])) for i in (0, 2)]
                    values = row(results[1], 'COUPLED SCATTERED FAR FIELD')
                    actual = [complex(values[i], values[i+1]) for i in (2, 4)]
                else:
                    values = row(results[0], heading)
                    expected = [cmath.rect(values[i], math.radians(values[i+1])) for i in (3, 5, 7)]
                    actual = near(results[1], magnetic)
                close(actual, expected, 3e-4)

    # Exercise combined fields and accepted/radiated power on retained M6 decks.
    with tempfile.TemporaryDirectory(prefix='nec2c-m8-power-') as temporary:
        output_path = Path(temporary)/'power.out'
        for name in ('m6-vacuum.nec', 'm6-cube.nec'):
            subprocess.run([str(BINARY), '-i', str(ROOT/'tests/dielectric'/name),
                            '-o', str(output_path)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            text = output_path.read_text()
            defect = float(re.search(r'RELATIVE POWER DEFECT:\s*(\S+)', text)[1])
            angular = float(re.search(r'ANGULAR QUADRATURE CHANGE:\s*(\S+)', text)[1])
            assert defect < 0.005, defect
            assert angular < 1e-5, angular
    print('M8 analytical fields and driven lossless power validation passed')


if __name__ == '__main__':
    main()
