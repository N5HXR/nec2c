#!/usr/bin/env python3
"""Cached native M13-03 split-pillar OpenEMS counterparts and refinements."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import time
from m10_pec_acceptance import impedance,scalar
from openems_m10_ground import material_boxes


def material_resource(split,frequency,pitch,loaded):
    strongest=math.sqrt(abs(complex(4,-.03/(2*math.pi*frequency*1e6*8.8541878128e-12))))
    kh=2*math.pi*frequency*1e6/299792458*strongest*pitch
    lines=['NEC2C-DIELECTRIC 2']
    boxes=material_boxes(dict(epsilon_r=2,dielectric_sigma_s_m=.01,cube_bottom_mm=0,material_split=split))
    if loaded:
        for region,(_,er,sigma,start,stop) in enumerate(boxes,1):
            lines+=['MATERIAL %d ER %g SIGMA %g'%(region,er,sigma),'REGION %d MATERIAL %d'%(region,region)]
            counts=[round((b-a)*.001/pitch) for a,b in zip(start,stop)]
            # Form both bounds from integer global indices, avoiding roundoff gaps.
            origins=[round(a*.001/pitch) for a in start]
            for x in range(counts[0]):
                for y in range(counts[1]):
                    for z in range(counts[2]):
                        lo=[(a+i)*pitch for a,i in zip(origins,(x,y,z))]
                        hi=[(a+i+1)*pitch for a,i in zip(origins,(x,y,z))]
                        lines+=['BOX %d '%region+' '.join(format(v,'.17g') for v in lo+hi)]
    else:lines+=['MATERIAL 1 ER 1','REGION 1 MATERIAL 1','BOX 1 .05 -.05 0 .15 .05 .65']
    lines+=['MESH MAX_KH %.17g MAX_ASPECT 1 MAX_CELLS 10000'%kh,'END']
    return '\n'.join(lines)+'\n'


def run_case(binary,root,split,frequency,pitch,segments):
    case=root/('split-%s-f%g-h%g-n%d'%(split,frequency,pitch,segments));case.mkdir(parents=True,exist_ok=True)
    runs=[];binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
    for loaded in (False,True):
        stem='loaded' if loaded else 'vacuum'
        deck=('CM M13-03 TWO TOUCHING LOSSY MATERIALS OVER PEC; POSITIVE WIRE GAP; EK\nCE\n'
              f'GW 1 {segments} 0 0 .3625 0 0 .8375 .005\nGE\n'
              f'FR 0 1 0 0 {frequency} 0\nGN 1\nEK\nEX 0 1 {(segments+1)//2} 0 1 0\n'
              f'DF "{stem}.diel"\nXQ\nEN\n')
        material=material_resource(split,frequency,pitch,loaded)
        spec=dict(binary_sha256=binary_hash,deck=deck,material=material,solver='fft',relative_tolerance=1e-10)
        specfile=case/(stem+'-spec.json');out=case/(stem+'.out');resultfile=case/(stem+'-result.json')
        cached=(specfile.exists() and json.loads(specfile.read_text())==spec and resultfile.exists() and out.exists() and 'TOTAL RUN TIME:' in out.read_text())
        (case/(stem+'.nec')).write_text(deck);(case/(stem+'.diel')).write_text(material)
        specfile.write_text(json.dumps(spec,indent=2)+'\n')
        started=time.monotonic()
        if not cached:
            with (case/(stem+'.log')).open('w') as log:
                subprocess.run([str(binary),'-s','fft','-r','1e-10','-i',str(case/(stem+'.nec')),'-o',str(out)],stdout=log,stderr=subprocess.STDOUT,check=True)
        text=out.read_text();z=impedance(text)
        old=json.loads(resultfile.read_text()) if cached else {}
        run=dict(loaded=loaded,z_ohm=[z.real,z.imag],active_cells=int(scalar(text,'ACTIVE CELLS:')),
                 true_relative_residual=scalar(text,'GMRES TRUE RELATIVE RESIDUAL:'),
                 power_defect=scalar(text,'RELATIVE POWER DEFECT:'),absorbed_power_w=scalar(text,'DIELECTRIC ABSORBED POWER (W):'),
                 accounted_solver_bytes=int(scalar(text,'FFT ACCOUNTED SOLVER BYTES:')),
                 wall_seconds=old.get('wall_seconds') if cached else time.monotonic()-started,
                 solver_cpu_seconds=scalar(text,'TOTAL RUN TIME:')/1000)
        assert run['true_relative_residual']<=1e-10 and run['absorbed_power_w']>=0 and run['power_defect']<.003,run
        resultfile.write_text(json.dumps(run,indent=2)+'\n');runs.append(run)
    delta=complex(*runs[1]['z_ohm'])-complex(*runs[0]['z_ohm'])
    result=dict(experiment='M13-03',material_split=split,frequency_mhz=frequency,pitch_m=pitch,wire_segments=segments,
                delta_z_ohm=[delta.real,delta.imag],runs=runs,binary_sha256=binary_hash,
                peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (case/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('./nec2c'))
    p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--splits',choices=('x','z'),nargs='+',default=['x','z'])
    p.add_argument('--frequencies',type=float,nargs='+',default=[280,300,320])
    p.add_argument('--pitches',type=float,nargs='+',default=[.025])
    p.add_argument('--segments',type=int,default=11)
    args=p.parse_args()
    if args.segments<3 or args.segments%2!=1:p.error('odd segment count >=3 required')
    if any(not math.isfinite(f) or f<=0 for f in args.frequencies):p.error('positive finite frequencies required')
    if any(h not in (.025,.0125) for h in args.pitches):p.error('supported pitch is .025 or .0125 m')
    root=args.workdir.resolve();root.mkdir(parents=True,exist_ok=True)
    results=[run_case(args.binary.resolve(),root,s,f,h,args.segments) for s in args.splits for f in args.frequencies for h in args.pitches]
    (root/'m13-03-results.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':main()
