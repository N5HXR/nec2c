#!/usr/bin/env python3
"""Cached M13-04 GN 2 wire/cube benchmarks, including independent power channels."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import time
from m10_pec_acceptance import impedance,scalar


def run_case(binary,root,frequency,divisions,segments,sigma):
    case=root/('f%g-n%d-w%d-s%g'%(frequency,divisions,segments,sigma));case.mkdir(parents=True,exist_ok=True)
    k=2*math.pi*frequency*1e6/299792458
    index=math.sqrt(abs(complex(3,-sigma/(2*math.pi*frequency*1e6*8.8541878128e-12))))
    kh=k*index*.1/divisions*(1+1e-12)
    runs=[];binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
    for loaded in (False,True):
        stem='loaded' if loaded else 'vacuum'
        deck=('CM M13-04 GN 2 WIRES; POSITIVE CUBE/GROUND/CONDUCTOR GAPS; EK\nCE\n'
              f'GW 1 {segments} 0 0 .3625 0 0 .8375 .005\nGE\n'
              f'FR 0 1 0 0 {frequency} 0\nGN 2 0 0 0 13 .02\nEK\nEX 0 1 {(segments+1)//2} 0 1 0\nDF "{stem}.diel"\nXQ\nEN\n')
        er=3 if loaded else 1;loss=sigma if loaded else 0
        material=(f'NEC2C-DIELECTRIC 2\nMATERIAL 1 ER {er} SIGMA {loss:.17g}\nREGION 1 MATERIAL 1\n'
                  'BOX 1 .05 -.05 .55 .15 .05 .65\n'
                  f'MESH MAX_KH {kh:.17g} MAX_ASPECT 1 MAX_CELLS 1000\nEND\n')
        spec=dict(binary_sha256=binary_hash,deck=deck,material=material,solver='gmres',relative_tolerance=1e-10)
        specfile=case/(stem+'-spec.json');out=case/(stem+'.out');resultfile=case/(stem+'-result.json')
        cached=(specfile.exists() and json.loads(specfile.read_text())==spec and resultfile.exists() and out.exists() and 'TOTAL RUN TIME:' in out.read_text())
        (case/(stem+'.nec')).write_text(deck);(case/(stem+'.diel')).write_text(material);specfile.write_text(json.dumps(spec,indent=2)+'\n')
        started=time.monotonic()
        if not cached:
            with (case/(stem+'.log')).open('w') as log:
                subprocess.run([str(binary),'-s','gmres','-r','1e-10','-i',str(case/(stem+'.nec')),'-o',str(out)],stdout=log,stderr=subprocess.STDOUT,check=True)
        text=out.read_text();z=impedance(text);old=json.loads(resultfile.read_text()) if cached else {}
        run=dict(loaded=loaded,z_ohm=[z.real,z.imag],active_cells=int(scalar(text,'ACTIVE CELLS:')),
                 true_relative_residual=scalar(text,'GMRES TRUE RELATIVE RESIDUAL:'),power_defect=scalar(text,'RELATIVE POWER DEFECT:'),
                 accepted_power_w=scalar(text,'ACCEPTED POWER:'),upper_radiated_power_w=scalar(text,'UPPER-HEMISPHERE RADIATED POWER:'),
                 lower_halfspace_power_w=scalar(text,'LOWER-HALF-SPACE POWER:'),absorbed_power_w=scalar(text,'DIELECTRIC ABSORBED POWER (W):'),
                 ground_flux_quadrature_change_w=scalar(text,'GROUND FLUX QUADRATURE CHANGE (W):'),
                 wall_seconds=old.get('wall_seconds') if cached else time.monotonic()-started,solver_cpu_seconds=scalar(text,'TOTAL RUN TIME:')/1000)
        assert run['active_cells']==(divisions**3 if loaded else 0),run
        assert run['true_relative_residual']<=1e-10 and run['power_defect']<.003,run
        assert run['absorbed_power_w']>=0 and run['lower_halfspace_power_w']>0 and run['upper_radiated_power_w']>=0,run
        independent_defect=abs(run['accepted_power_w']-run['upper_radiated_power_w']-run['lower_halfspace_power_w']-run['absorbed_power_w'])/run['accepted_power_w']
        assert abs(independent_defect-run['power_defect'])<1e-10
        resultfile.write_text(json.dumps(run,indent=2)+'\n');runs.append(run)
    delta=complex(*runs[1]['z_ohm'])-complex(*runs[0]['z_ohm'])
    result=dict(experiment='M13-04',frequency_mhz=frequency,divisions_per_axis=divisions,pitch_m=.1/divisions,wire_segments=segments,
                dielectric_sigma_s_m=sigma,ground_epsilon_r=13,ground_sigma_s_m=.02,delta_z_ohm=[delta.real,delta.imag],runs=runs,
                binary_sha256=binary_hash,peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (case/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('./nec2c'));p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--frequencies',type=float,nargs='+',default=[280,300,320]);p.add_argument('--divisions',type=int,choices=(2,3,4),default=3)
    p.add_argument('--segments',type=int,default=11);p.add_argument('--sigma',type=float,default=.02)
    args=p.parse_args()
    if args.segments<3 or args.segments%2!=1:p.error('odd segment count >=3 required')
    if any(not math.isfinite(f) or f<=0 for f in args.frequencies):p.error('positive finite frequencies required')
    if not math.isfinite(args.sigma) or args.sigma<0:p.error('passive finite conductivity required')
    root=args.workdir.resolve();root.mkdir(parents=True,exist_ok=True)
    results=[run_case(args.binary.resolve(),root,f,args.divisions,args.segments,args.sigma) for f in args.frequencies]
    (root/'m13-04-results.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':main()
