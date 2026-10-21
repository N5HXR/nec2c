#!/usr/bin/env python3
"""Reproducible VM-only PEC impedance-refinement checkpoints.

Use a fresh workdir for each binary build. Completed outputs are reused only
when the binary manifest, generated inputs and completion marker agree.
Groups can be run separately; the three groups cover volume, wire and
frequency refinement without requiring OpenEMS.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time


def impedance(text):
    section = text.split('COUPLED ANTENNA INPUT PARAMETERS', 1)[1]
    for line in section.splitlines():
        try:
            row = [float(x) for x in line.split()]
        except ValueError:
            continue
        if len(row) == 8:
            return complex(row[-2], row[-1])
    raise ValueError('missing feed impedance')


def scalar(text, label):
    return float(re.search(re.escape(label)+r'\s*(\S+)', text)[1])


def run_case(binary, root, frequency, kh, segments):
    case = root / ('f%g-kh%g-n%d' % (frequency, kh, segments))
    case.mkdir(parents=True, exist_ok=True)
    records = []
    for er in (1, 3):
        stem = 'vacuum' if er == 1 else 'loaded'
        deck = ('CM M10 PEC REFINEMENT; 5MM RADIUS; EK; POSITIVE GAP\nCE\n'
                f'GW 1 {segments} 0 0 .3625 0 0 .8375 .005\nGE\n'
                f'FR 0 1 0 0 {frequency} 0\nGN 1\nEK\n'
                f'EX 0 1 {(segments+1)//2} 0 1 0\nDF "{stem}.diel"\nXQ\nEN\n')
        resource = (f'nec2c-dielectric 1\nMATERIAL 1 ER {er}\nREGION 1 MATERIAL 1\n'
                    f'BOX 1 .05 -.05 .55 .15 .05 .65\nMESH MAX_KH {kh} MAX_ASPECT 2 MAX_CELLS 1000\nEND\n')
        same_inputs = ((case/(stem+'.nec')).exists() and
                       (case/(stem+'.nec')).read_text()==deck and
                       (case/(stem+'.diel')).exists() and
                       (case/(stem+'.diel')).read_text()==resource)
        (case/(stem+'.nec')).write_text(deck)
        (case/(stem+'.diel')).write_text(resource)
        output = case/(stem+'.out')
        log = case/(stem+'.log')
        # Cache is local to a run directory and its pinned binary, not a global
        # cross-version numerical cache. Reuse only normally completed results.
        valid = same_inputs and output.exists() and 'TOTAL RUN TIME:' in output.read_text()
        start = time.monotonic()
        if not valid:
            with log.open('w') as stream:
                subprocess.run([str(binary), '-s', 'gmres', '-r', '1e-10',
                                '-i', str(case/(stem+'.nec')), '-o', str(output)],
                               stdout=stream, stderr=subprocess.STDOUT, check=True)
        text = output.read_text()
        z = impedance(text)
        records.append(dict(er=er, z_ohm=[z.real, z.imag],
                            active_cells=int(scalar(text, 'ACTIVE CELLS:')),
                            true_residual=scalar(text, 'GMRES TRUE RELATIVE RESIDUAL:'),
                            power_defect=scalar(text, 'RELATIVE POWER DEFECT:'),
                            cpu_seconds=scalar(text, 'TOTAL RUN TIME:')/1000,
                            wall_seconds=None if valid else time.monotonic()-start))
    delta = complex(*records[1]['z_ohm'])-complex(*records[0]['z_ohm'])
    result = dict(frequency_mhz=frequency, mesh_max_kh=kh, wire_segments=segments,
                  delta_z_ohm=[delta.real, delta.imag], runs=records)
    (case/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=Path('./nec2c'))
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--group', choices=('volume', 'frequency', 'wire'), required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    provenance = {'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
                  'solver':'gmres','requested_relative_tolerance':1e-10}
    manifest=root/'run-manifest.json'
    if not manifest.exists() and any(root.rglob('*.out')):
        parser.error('unprovenanced existing outputs; choose a fresh workdir')
    try:
        with manifest.open('x') as stream: stream.write(json.dumps(provenance,indent=2)+'\n')
    except FileExistsError:
        if json.loads(manifest.read_text())!=provenance:
            parser.error('workdir belongs to a different binary/configuration; choose a fresh directory')
    cases = {'volume': [(300, .4, 11), (300, .24, 11)],
             'frequency': [(280, .3, 11), (320, .3, 11)],
             'wire': [(300, .3, 21)]}[args.group]
    results = [run_case(binary, root, *case) for case in cases]
    (root/(args.group+'-results.json')).write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()
