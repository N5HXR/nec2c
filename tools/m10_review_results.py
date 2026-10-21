#!/usr/bin/env python3
"""Offline OpenEMS port/log review, including explicitly provisional traces.

Never modifies imported data or promotes an incomplete run to converged.
DFTs use each probe's native timestamps, preserving the H half-step offset.
Tail-window sensitivity is a diagnostic, not a bound on truncation error.
"""
import argparse
import cmath
import csv
import json
import math
from pathlib import Path
import re

FREQUENCIES=(280.0,300.0,320.0)


def trace(path):
    data=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('%'): continue
        fields=line.split(); data.append((float(fields[0]),float(fields[1])))
    if len(data)<3: raise ValueError('insufficient probe samples: '+str(path))
    if any(not math.isfinite(t) or not math.isfinite(v) for t,v in data):
        raise ValueError('nonfinite probe data: '+str(path))
    if any(b[0]<=a[0] for a,b in zip(data,data[1:])): raise ValueError('nonmonotonic probe time')
    return data


def impedance(voltage,current,frequency_mhz,fraction=1.0):
    values=[]
    for data in (voltage,current):
        samples=data[:max(3,int(len(data)*fraction))]
        dt=samples[1][0]-samples[0][0]
        values.append(dt*sum(v*cmath.exp(-2j*math.pi*frequency_mhz*1e6*t) for t,v in samples))
    if abs(values[1])==0: raise ValueError('zero spectral port current')
    return values[0]/values[1]


def review(root):
    runs=[]; grouped={}
    for file in sorted(Path(root).rglob('spec.json')):
        spec=json.loads(file.read_text())
        if spec.get('action')!='verify': continue
        directory=file.parent
        log=(directory/'worker.log').read_text(errors='replace')
        complete=(directory/'result.json').exists()
        result=json.loads((directory/'result.json').read_text()) if complete else None
        if complete and result['spec']!=spec: raise ValueError('result/spec mismatch: '+directory.name)
        final=re.search(r'Time for\s+(\d+) iterations[^\n]*?:\s*([0-9.eE+-]+) sec',log)
        energy=re.findall(r'Energy:[^\n]*?\(\s*([-+0-9.]+)dB\)',log)
        step=re.search(r'FDTD timestep is:\s*(\S+)',log)
        cap=('before the end-criteria' in log or 'Cutting to max number of timesteps' in log)
        early=bool(complete and final and int(final[1])<spec['steps'] and not cap)
        run={'job_id':directory.name,'material_split':spec.get('material_split',''),'witness_frequency_mhz':spec.get('witness_frequency_mhz',0),'ground':spec['ground'],'epsilon_r':spec['epsilon_r'],'dielectric_sigma_s_m':spec.get('dielectric_sigma_s_m',0),
             'mesh_mm':spec['mesh_mm'],'wire_cell_mm':spec['wire_cell_mm'],
             'wire_radius_mm':spec['wire_radius_mm'],'domain_margin_mm':spec.get('domain_margin_mm',0),
             'feed_gap_mm':spec.get('feed_gap_mm',10),'cube_bottom_mm':spec.get('cube_bottom_mm',550),'result_present':complete,
             'normal_early_stop':early,'last_logged_energy_db':float(energy[-1]) if energy else None,
             'engine_seconds':float(final[2]) if final else None,
             'engine_timestep_seconds':float(step[1]) if step else None,
             'source_clipped_or_step_cap':cap,
             'status':'completed_early' if early else 'provisional_incomplete_or_unconfirmed',
             'impedance':[]}
        voltage=trace(directory/'solver/port_ut_1'); current=trace(directory/'solver/port_it_1')
        for frequency in FREQUENCIES:
            z=impedance(voltage,current,frequency)
            if result:
                saved=next(a for a in result['impedance'] if a['frequency_mhz']==frequency)
                error=abs(z-complex(saved['r_ohm'],saved['x_ohm']))
                if error>1e-8: raise ValueError('native DFT reproduction mismatch: '+directory.name)
            run['impedance'].append({'frequency_mhz':frequency,'r_ohm':z.real,'x_ohm':z.imag,
                                      'last_10_percent_sensitivity_ohm':abs(z-impedance(voltage,current,frequency,.9))})
        runs.append(run)
        key=tuple(spec.get(k) for k in ('ground','mesh_mm','wire_cell_mm','wire_radius_mm',
                                       'bandwidth_mhz','steps','end_criterion','mesh_revision','feed_model'))
        key+=(spec.get('domain_margin_mm',0),spec.get('feed_gap_mm',10),spec.get('cube_bottom_mm',550),spec.get('witness_frequency_mhz',0))
        key+=(spec.get('material_split',''),)
        group=grouped.setdefault(key,{})
        if spec['epsilon_r'] in group: raise ValueError('duplicate matching physics: '+str(key))
        group[spec['epsilon_r']]=run
    pairs=[]
    for key,group in sorted(grouped.items()):
        if 1 not in group or len(group)!=2: continue
        a=group[1];b=next(v for k,v in group.items() if k!=1)
        for vacuum,loaded in zip(a['impedance'],b['impedance']):
            delta=complex(loaded['r_ohm']-vacuum['r_ohm'],loaded['x_ohm']-vacuum['x_ohm'])
            pairs.append({'material_split':b['material_split'],'loaded_epsilon_r':b['epsilon_r'],'witness_frequency_mhz':b['witness_frequency_mhz'],'dielectric_sigma_s_m':b['dielectric_sigma_s_m'],'ground':a['ground'],'mesh_mm':a['mesh_mm'],'wire_cell_mm':a['wire_cell_mm'],
                          'domain_margin_mm':a['domain_margin_mm'],'feed_gap_mm':a['feed_gap_mm'],'cube_bottom_mm':a['cube_bottom_mm'],
                          'frequency_mhz':vacuum['frequency_mhz'],'delta_r_ohm':delta.real,'delta_x_ohm':delta.imag,
                          'status':'completed_pair' if a['normal_early_stop'] and b['normal_early_stop'] else 'provisional_incomplete_pair',
                          'vacuum_job':a['job_id'],'loaded_job':b['job_id']})
    if not runs: raise ValueError('no verification jobs found')
    return {'runs':runs,'pairs':pairs,'acceptance':'unreviewed physical/source/domain accuracy; provisional entries are not converged evidence'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('results',type=Path); p.add_argument('--output-dir',type=Path)
    args=p.parse_args(); result=review(args.results)
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        (args.output_dir/'review.json').write_text(json.dumps(result,indent=2)+'\n')
        with (args.output_dir/'pairs.csv').open('w',newline='') as stream:
            fields=('material_split','loaded_epsilon_r','witness_frequency_mhz','dielectric_sigma_s_m','ground','mesh_mm','wire_cell_mm','domain_margin_mm','feed_gap_mm','cube_bottom_mm','frequency_mhz','delta_r_ohm','delta_x_ohm','status','vacuum_job','loaded_job')
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(result['pairs'])
    else: print(json.dumps(result,indent=2))

if __name__=='__main__': main()
