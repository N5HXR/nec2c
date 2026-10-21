#!/usr/bin/env python3
"""Run the production FFT antenna/tree path and retain auditable inputs/results."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import subprocess
import time
from m9_tree_sizing import TREE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['tree', 'skeleton'], default='tree')
    parser.add_argument('--ground',choices=('free','pec'),default='free')
    parser.add_argument('--ground-gap',type=float,default=0,help='translate dielectric upward in metres; PEC only')
    parser.add_argument('--binary', type=Path, default=Path('./nec2c'))
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--max-iterations', type=int, default=300)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if not math.isfinite(args.ground_gap) or args.ground_gap<0 or (args.ground=='free' and args.ground_gap):
        parser.error('--ground-gap requires PEC and a finite nonnegative gap')
    boxes = TREE if args.model == 'tree' else TREE[:3]
    boxes = [tuple(v+(args.ground_gap if i in (2,5) else 0) for i,v in enumerate(box)) for box in boxes]
    q = 2*math.pi*28e6/299792458*math.sqrt(3)*.25
    material = 'nec2c-dielectric 1\nMATERIAL 1 ER 3\nREGION 1 MATERIAL 1\n'
    material += ''.join('BOX 1 '+' '.join(str(v) for v in box)+'\n' for box in boxes)
    material += 'MESH MAX_KH %.17g MAX_CELLS 10000 MAX_ASPECT 1\nEND\n' % q
    (root/'model.diel').write_text(material)
    ground_card='GN 1\n' if args.ground=='pec' else ''
    pattern='RP 0 3 2 0000 0 0 45 90 0 0' if args.ground=='pec' else 'RP 0 5 2 0000 0 0 45 90 0 0'
    deck = ('CM M9/M11 N3 ANTENNA TREE\nCE\n'
            'GW 1 11 -5 0 3.5 -5 0 8.5 .002\nGE\n'
            'FR 0 1 0 0 28\n'+ground_card+'EX 0 1 6 0 1 0\nDF "model.diel"\n'+pattern+'\nEN\n')
    (root/'model.nec').write_text(deck)
    command = [str(args.binary.resolve()), '-s', 'fft', '-r', str(args.rtol),
               '-n', str(args.max_iterations), '-i', str(root/'model.nec'),
               '-o', str(root/'model.out')]
    start = time.monotonic()
    with (root/'stdout.txt').open('w') as stdout, (root/'stderr.txt').open('w') as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr)
    text = (root/'model.out').read_text()
    record = dict(model=args.model, ground=args.ground, dielectric_bottom_gap_m=args.ground_gap,
                  binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                  target_relative_residual=args.rtol,
                  wall_seconds=time.monotonic()-start,
                  peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                  returncode=result.returncode, command=command)
    for label, key in [('GMRES ITERATIONS:', 'iterations'),
                       ('GMRES TRUE RELATIVE RESIDUAL:', 'true_relative_residual'),
                       ('FFT ACCOUNTED SOLVER BYTES:', 'accounted_solver_bytes'),
                       ('RELATIVE POWER DEFECT:', 'power_defect'),
                       ('TOTAL RUN TIME:', 'solver_cpu_milliseconds')]:
        match = re.search(re.escape(label)+r'\s*(\S+)', text)
        if match:
            record[key] = float(match[1])
    match = re.search(r'COUPLED ANTENNA INPUT PARAMETERS.*?TAG SEG[^\n]*\n\s*([^\n]+)', text, re.S)
    if match:
        row = match[1].split()
        record['feed_impedance_real_ohm'] = float(row[6])
        record['feed_impedance_imag_ohm'] = float(row[7])
    (root/'benchmark.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2), flush=True)
    if result.returncode:
        print(text[-2000:])
        raise SystemExit(1)


if __name__ == '__main__':
    main()
