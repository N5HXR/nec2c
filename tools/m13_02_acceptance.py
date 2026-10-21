#!/usr/bin/env python3
"""Cached VM counterparts to the M13-02 lossy laptop grounded-pillar benchmark."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import time
from m10_pec_acceptance import impedance,scalar


def run_case(binary,root,frequency,pitch,segments):
    er_loaded=10.5-frequency/40
    sigma=.01+(frequency-260)*.00025
    if not 260<=frequency<=340: raise ValueError('table range is 260..340 MHz')
    index=math.sqrt(math.hypot(4,.03/(2*math.pi*340e6*8.8541878128e-12)))
    kh=2*math.pi*frequency*1e6/299792458*index*pitch
    name='f%g-h%g-n%d-s%g'%(frequency,pitch,segments,sigma)
    case=root/name;case.mkdir(parents=True,exist_ok=True)
    runs=[]
    for er in (1,er_loaded):
        stem='vacuum' if er==1 else 'loaded'
        deck=('CM M13-02 LOSSY PEC CONTACT PILLAR; POSITIVE WIRE GAP; EK\nCE\n'
              f'GW 1 {segments} 0 0 .3625 0 0 .8375 .005\nGE\n'
              f'FR 0 1 0 0 {frequency} 0\nGN 1\nEK\nEX 0 1 {(segments+1)//2} 0 1 0\n'
              f'DF "{stem}.diel"\nXQ\nEN\n')
        material=(f'nec2c-dielectric 2\nMATERIAL 1 ER {er} SIGMA {sigma if er!=1 else 0}\nREGION 1 MATERIAL 1\n'
                  f'BOX 1 .05 -.05 0 .15 .05 .65\nMESH MAX_KH {kh:.17g} MAX_ASPECT 1 MAX_CELLS 10000\nEND\n')
        if er!=1:
            material=('NEC2C-DIELECTRIC 3\nMATERIAL 1 TABLE\n'
                'SAMPLE 1 FREQ_MHZ 260 ER 4 SIGMA .01\n'
                'SAMPLE 1 FREQ_MHZ 300 ER 3 SIGMA .02\n'
                'SAMPLE 1 FREQ_MHZ 340 ER 2 SIGMA .03\n'
                'REGION 1 MATERIAL 1\nBOX 1 .05 -.05 0 .15 .05 .65\n'
                f'MESH MAX_KH {kh:.17g} MAX_ASPECT 1 MAX_CELLS 10000\nEND\n')
        spec=dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),deck=deck,
                  material=material,solver='fft',relative_tolerance=1e-10)
        specfile=case/(stem+'-spec.json');out=case/(stem+'.out')
        cached=(specfile.exists() and json.loads(specfile.read_text())==spec and
                out.exists() and 'TOTAL RUN TIME:' in out.read_text())
        (case/(stem+'.nec')).write_text(deck);(case/(stem+'.diel')).write_text(material)
        specfile.write_text(json.dumps(spec,indent=2)+'\n')
        started=time.monotonic()
        if not cached:
            with (case/(stem+'.log')).open('w') as log:
                subprocess.run([str(binary),'-s','fft','-r','1e-10','-i',str(case/(stem+'.nec')),
                                '-o',str(out)],stdout=log,stderr=subprocess.STDOUT,check=True)
        text=out.read_text();z=impedance(text)
        runs.append(dict(epsilon_r=er,z_ohm=[z.real,z.imag],
                         active_cells=int(scalar(text,'ACTIVE CELLS:')),
                         true_relative_residual=scalar(text,'GMRES TRUE RELATIVE RESIDUAL:'),
                         power_defect=scalar(text,'RELATIVE POWER DEFECT:'),
                         absorbed_power_w=scalar(text,'DIELECTRIC ABSORBED POWER (W):'),
                         accounted_solver_bytes=int(scalar(text,'FFT ACCOUNTED SOLVER BYTES:')),
                         wall_seconds=None if cached else time.monotonic()-started,
                         solver_cpu_seconds=scalar(text,'TOTAL RUN TIME:')/1000))
    delta=complex(*runs[1]['z_ohm'])-complex(*runs[0]['z_ohm'])
    result=dict(experiment="M13-02",conductivity_s_m=sigma,evaluated_epsilon_r=er_loaded,frequency_mhz=frequency,pitch_m=pitch,wire_segments=segments,
                delta_z_ohm=[delta.real,delta.imag],runs=runs,
                peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (case/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('./nec2c'))
    p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--frequencies',type=float,nargs='+',default=[300])
    p.add_argument('--pitches',type=float,nargs='+',default=[.025,.0125])
    p.add_argument('--segments',type=int,default=11)
    args=p.parse_args()
    if args.segments<3 or args.segments%2!=1: p.error('use an odd segment count >=3')
    if any(not math.isfinite(f) or f<=0 for f in args.frequencies): p.error('positive finite frequencies required')
    if any(h not in (.025,.0125) for h in args.pitches): p.error('supported uniform pitches: .025, .0125 m')
    root=args.workdir.resolve();root.mkdir(parents=True,exist_ok=True)
    results=[run_case(args.binary.resolve(),root,f,h,args.segments)
             for f in args.frequencies for h in args.pitches]
    (root/'m13-02-results.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__': main()
