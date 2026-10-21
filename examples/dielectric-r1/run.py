#!/usr/bin/env python3
"""Run the small revision-1 reference examples without OpenEMS or extra packages."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess


def scalars(text,label):
    return [float(x) for x in re.findall(re.escape(label)+r'\s*(\S+)',text)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('./nec2c'))
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();binary=args.binary.resolve();args.output_dir.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).resolve().parent;results=[]
    for name,solver,count,solves in [('free-space','dense',1,3),('pec-materials','fft',3,3),('real-ground','gmres',1,1),('plane-wave','dense',1,1)]:
        out=(args.output_dir/(name+'.out')).resolve()
        with (args.output_dir/(name+'.log')).open('w') as log:
            subprocess.run([str(binary),'-s',solver,'-r','1e-10','-i',str(source/(name+'.nec')),'-o',str(out)],stdout=log,stderr=subprocess.STDOUT,check=True)
        text=out.read_text()
        if 'TOTAL RUN TIME:' not in text:raise ValueError('incomplete output: '+name)
        residuals=scalars(text,'GMRES TRUE RELATIVE RESIDUAL:' if solver!='dense' else 'RELATIVE RESIDUAL:')
        defects=scalars(text,'RELATIVE IMBALANCE:' if name=='plane-wave' else 'RELATIVE POWER DEFECT:')
        absorbed=scalars(text,'DIELECTRIC ABSORBED POWER (W):')
        if len(defects)!=solves or len(absorbed)!=solves or len(residuals)!=solves or any(not math.isfinite(x) for x in defects+absorbed+residuals):raise ValueError('missing/nonfinite diagnostics: '+name)
        if any(x<0 or x>1e-9 for x in residuals):raise ValueError('residual check failed: '+name)
        if any(x<0 or x>=.003 for x in defects) or any(x<0 for x in absorbed):raise ValueError('power check failed: '+name)
        rows=[]
        for section in text.split('COUPLED ANTENNA INPUT PARAMETERS')[1:]:
            for line in section.splitlines():
                try:values=[float(x) for x in line.split()]
                except ValueError:continue
                if len(values)==8:rows.append(values[-2:]);break
        results.append(dict(example=name,solver=solver,frequency_count=count,z_ohm=rows,absorbed_power_w=absorbed,power_defects=defects,residuals=residuals,output_file=out.name))
        print(name+': completed; '+str(out),flush=True)
    (args.output_dir/'summary.json').write_text(json.dumps(dict(feature_revision=1,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),examples=results),indent=2)+'\n')


if __name__=='__main__':main()
