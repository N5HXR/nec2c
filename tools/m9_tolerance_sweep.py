#!/usr/bin/env python3
"""Compare feed impedance and complex patterns to dense solutions on one mesh."""
import argparse
import csv
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tests'))
from dielectric_wire import coupled_impedance


def fields(text):
    section = text.split('COUPLED SCATTERED FAR FIELD', 1)[1]
    result = []
    for line in section.splitlines():
        try:
            row = [float(v) for v in line.split()]
        except ValueError:
            continue
        if len(row) == 6:
            result.extend((complex(row[2], row[3]), complex(row[4], row[5])))
    assert result
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT/'nec2c')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true', help='append missing points to an existing CSV')
    args = parser.parse_args()
    completed = set()
    if args.resume and args.output.exists():
        with args.output.open() as prior:
            completed = {(row['fixture'], float(row['rtol'])) for row in csv.DictReader(prior)}
    columns = ['fixture', 'rtol', 'true_residual', 'iterations', 'impedance_error_ohm',
               'pattern_scaled_error']
    with args.output.open('a' if completed else 'w', newline='') as output, tempfile.TemporaryDirectory(prefix='m9-tol-') as temp:
        writer = csv.DictWriter(output, fieldnames=columns)
        if not completed:
            writer.writeheader()
        directory = Path(temp)
        for name, resource_name, box in [
                ('m6-cube.nec', 'm6-cube.diel', None),
                ('m7-mixed-cube.nec', 'm7-mixed-cube.diel',
                 'BOX 1 -0.125 0.325 -0.025 -0.075 0.375 0.025')]:
            if all((name, tol) in completed for tol in (1e-4, 1e-6, 1e-8, 1e-10)):
                continue
            material = (ROOT/'tests/dielectric'/resource_name).read_text()
            if box:
                material = re.sub(r'^BOX .*$', box, material, flags=re.M)
            material = re.sub(r'MAX_KH \S+', 'MAX_KH 0.30', material)
            (directory/resource_name).write_text(material)
            deck = (ROOT/'tests/dielectric'/name).read_text()
            deck = deck.split('RP ')[0] if 'RP ' in deck else deck.split('XQ')[0]
            deck += 'RP 0 5 2 0000 0 0 45 90 0 0\nEN\n'
            (directory/'model.nec').write_text(deck)

            def solve(mode, tol):
                result = subprocess.run([str(args.binary.resolve()), '-s', mode, '-r', str(tol),
                                         '-i', str(directory/'model.nec'), '-o', str(directory/'model.out')],
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                text = (directory/'model.out').read_text()
                assert result.returncode == 0, text+result.stderr.decode()
                return text

            dense = solve('dense', 1e-12)
            reference_z, reference_fields = coupled_impedance(dense), fields(dense)
            for tol in (1e-4, 1e-6, 1e-8, 1e-10):
                if (name, tol) in completed:
                    continue
                text = solve('fft', tol)
                residual = float(re.search(r'GMRES TRUE RELATIVE RESIDUAL:\s*(\S+)', text)[1])
                current = fields(text)
                assert len(current) == len(reference_fields)
                z_error = abs(coupled_impedance(text)-reference_z)
                pattern_error = max(abs(a-b) for a,b in zip(current,reference_fields))/max(map(abs,reference_fields))
                row = dict(fixture=name, rtol=tol, true_residual=residual,
                           iterations=int(re.search(r'GMRES ITERATIONS:\s*(\d+)', text)[1]),
                           impedance_error_ohm=z_error, pattern_scaled_error=pattern_error)
                writer.writerow(row)
                output.flush()
                print(row, flush=True)
                assert residual <= tol
                if tol <= 1e-10:
                    assert z_error < 1e-4 and pattern_error < 1e-6, row
    print('Dense-reference impedance/pattern tolerance sweeps passed', flush=True)


if __name__ == '__main__':
    main()
