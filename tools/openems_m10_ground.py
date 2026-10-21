#!/usr/bin/env python3
"""Portable M10 ground-gap jobs with a shared wall-time budget (POSIX).

Default --mode plan runs no simulation and needs only the Python standard
library. The OpenEMS/CSXCAD/NumPy imports occur only inside laptop workers.
Pilot results measure cost, not accuracy. Verification output is unreviewed
until decay, mesh and boundary-domain convergence have been checked.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import signal
import subprocess
import sys
import time


def run_worker(command, seconds, log):
    """Timeout/interrupt terminates the whole POSIX worker process group."""
    started=time.monotonic()
    with Path(log).open('w') as stream:
        process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            code=process.wait(timeout=seconds)
            return {'status':'returned' if code==0 else 'failed','returncode':code,
                    'wall_seconds':time.monotonic()-started}
        except (subprocess.TimeoutExpired,KeyboardInterrupt) as error:
            # Kill the group even if the worker itself exits before its solver.
            try: os.killpg(process.pid,signal.SIGTERM)
            except ProcessLookupError: pass
            try: process.wait(timeout=2)
            except subprocess.TimeoutExpired: pass
            try: os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError: pass
            process.wait()
            if isinstance(error,KeyboardInterrupt): raise
            return {'status':'timed_out','wall_seconds':time.monotonic()-started}


MESH_REVISION = 2


def graded_axis(anchors, maximum_step, ratio=2.0):
    """Preserve anchors, balance each interval, then bound adjacent growth.

    Bisecting only intervals larger than twice a neighbor cannot create a
    cell smaller than the input minimum. This avoids residual end slivers
    and guarantees a 2:1 bound instead of relying on a smoother's hints.
    """
    anchors=sorted(set(float(v) for v in anchors))
    points=[anchors[0]]
    for a,b in zip(anchors,anchors[1:]):
        count=max(1,int(math.ceil((b-a)/maximum_step)))
        points.extend(a+(b-a)*i/count for i in range(1,count+1))
    for _ in range(64):
        widths=[b-a for a,b in zip(points,points[1:])]
        split=[False]*len(widths)
        for i,width in enumerate(widths):
            near=widths[max(0,i-1):min(len(widths),i+2)]
            split[i]=width>ratio*min(near)*(1+1e-12)
        if not any(split): return points
        refined=[points[0]]
        for i,do_split in enumerate(split):
            if do_split: refined.append((points[i]+points[i+1])/2)
            refined.append(points[i+1])
        points=refined
    raise RuntimeError('mesh grading did not settle')


def build_axes(mesh_mm,wire_cell_mm,ground,wire_radius_mm=5,domain_margin_mm=0,feed_gap_mm=10,cube_bottom_mm=550,material_split=""):
    if not all(math.isfinite(v) and v>0 for v in (mesh_mm,wire_cell_mm,wire_radius_mm)):
        raise ValueError('mesh sizes must be positive and finite')
    # Wire radius and grid spacing are separate refinement parameters.
    # Include the finite feed boundaries, axis and requested near-wire spacing.
    zmin=0 if ground=='pec' else -500-domain_margin_mm
    extent=650+domain_margin_mm
    divisions=max(1,int(math.ceil(wire_radius_mm/wire_cell_mm-1e-12)))
    wire_anchors=[wire_radius_mm*i/divisions for i in range(-divisions,divisions+1)]
    return [graded_axis([-extent,50,150,extent]+wire_anchors+([100] if material_split=="x" else []),mesh_mm),
            graded_axis([-extent,-50,50,extent]+wire_anchors,mesh_mm),
            graded_axis([zmin,0,362.5,cube_bottom_mm,550,600-feed_gap_mm/2,600+feed_gap_mm/2,650,837.5,1400+domain_margin_mm]+([325] if material_split=="z" else []),mesh_mm)]


def mesh_metrics(lines):
    widths=[[b-a for a,b in zip(axis,axis[1:])] for axis in lines]
    minimum=[min(w) for w in widths]
    maximum_ratio=[max([max(a/b,b/a) for a,b in zip(w,w[1:])] or [1]) for w in widths]
    counts=[len(w) for w in widths]
    cells=counts[0]*counts[1]*counts[2]
    dt=1/(299792458*math.sqrt(sum(1/(v*.001)**2 for v in minimum)))
    return {'cells':cells,'grid_cells':counts,'minimum_step_mm':minimum,
            'maximum_adjacent_ratio':maximum_ratio,'vacuum_cfl_seconds':dt,
            'epsilon13_cfl_lower_estimate_seconds':dt/math.sqrt(13)}


def inspect_results(directory):
    """Audit returned XML and engine logs without importing OpenEMS."""
    reports=[]
    for xml in sorted(Path(directory).rglob('geometry.xml')):
        root=ET.parse(xml).getroot()
        axes=[[float(v) for v in root.find('.//'+axis+'Lines').text.split(',')] for axis in 'XYZ']
        unit=float(root.find('.//RectilinearGrid').get('DeltaUnit'))
        axes=[[v*unit*1000 for v in line] for line in axes]
        report={'directory':str(xml.parent),**mesh_metrics(axes)}
        log=xml.with_name('worker.log')
        text=log.read_text(errors='replace') if log.exists() else ''
        timing=re.search(r'Time for\s+(\d+) iterations[^\n]*?:\s*([0-9.eE+-]+) sec',text)
        step=re.search(r'Timestep \(s\)\s*:\s*(\S+)',text)
        pulse=re.search(r'Requested excitation (?:pusle|pulse) would be (\d+) timesteps or ([0-9.eE+-]+) s',text)
        if step: report['engine_timestep_seconds']=float(step[1])
        if timing:
            per_step=float(timing[2])/int(timing[1])
            report.update({'measured_seconds_per_step':per_step,
                           'projected_120000_step_seconds':120000*per_step,
                           'projected_pair_120000_step_seconds':240000*per_step})
            if step: report['pilot_simulated_seconds']=int(timing[1])*float(step[1])
        if pulse:
            report['requested_pulse_steps']=int(pulse[1]); report['requested_pulse_seconds']=float(pulse[2])
        report['source_was_clipped']='Cutting to max number of timesteps' in text
        report['step_cap_before_decay']='before the end-criteria' in text
        load=re.search(r'--- Series/Parallel Lumped RLC load ---\s*Active cells\s*:\s*(\d+)',text)
        if load: report['reported_lumped_load_active_cells']=int(load[1])
        reports.append(report)
    if not reports: raise ValueError('no geometry.xml found in '+str(directory))
    return reports


def material_boxes(spec):
    """M13-03 split pillar: near/lower ER=2 sigma=.01, far/upper ER=4 sigma=.03.

    The first material still uses the ordinary runner controls. Vacuum keeps
    identical interface mesh anchors, but contains neither material.
    """
    if spec['epsilon_r']==1: return []
    start=[50,-50,spec.get('cube_bottom_mm',550)];stop=[150,50,650]
    first=('cube',spec['epsilon_r'],spec.get('dielectric_sigma_s_m',0),start,stop)
    split=spec.get('material_split','')
    if not split: return [first]
    axis=0 if split=='x' else 2
    cut=100 if split=='x' else 325
    stop1=stop.copy();stop1[axis]=cut
    start2=start.copy();start2[axis]=cut
    return [('cube',spec['epsilon_r'],spec.get('dielectric_sigma_s_m',0),start,stop1),
            ('cube2',4,.03,start2,stop)]


def write_nec_reference(root,spec):
    """Keep the benchmark wire radius/materials in the NEC counterparts too."""
    ground='GN 1' if spec['ground']=='pec' else 'GN 2 0 0 0 13 .02'
    bottom=spec.get('cube_bottom_mm',550)*.001
    contact=bottom==0
    reference_frequency=spec.get('witness_frequency_mhz',300)
    frequency=('FR 0 1 0 0 %g 0'%reference_frequency) if contact else 'FR 0 3 0 0 280 20'
    kh=2*math.pi*reference_frequency*1e6/299792458*math.sqrt(abs(complex(spec['epsilon_r'],-spec.get('dielectric_sigma_s_m',0)/(2*math.pi*reference_frequency*1e6*8.8541878128e-12))))*.025 if contact else .3
    feed_label=str(spec.get('feed_gap_mm',10))
    before=('CM M10 LAPTOP BENCHMARK; NEC DELTA-GAP VERSUS OPENEMS '+feed_label+'MM FEED\nCE\n'
            'GW 1 11 0 0 .3625 0 0 .8375 '+str(spec['wire_radius_mm']*.001)+'\n'
            'GE\n'+frequency+'\n'+ground+'\nEK\nEX 0 1 6 0 1 0\n')
    (root/'reference.nec').write_text(before+'DF "reference.diel"\nXQ\nEN\n')
    (root/'reference-vacuum-legacy.nec').write_text(before+'XQ\nEN\n')
    (root/'reference.diel').write_text(('nec2c-dielectric 2' if spec.get('dielectric_sigma_s_m',0) else 'nec2c-dielectric 1')+'\nMATERIAL 1 ER '+str(spec['epsilon_r'])+
        (' SIGMA '+str(spec['dielectric_sigma_s_m']) if spec.get('dielectric_sigma_s_m',0) else '')+
        '\nREGION 1 MATERIAL 1\nBOX 1 .05 -.05 '+str(bottom)+' .15 .05 .65\n'
        'MESH MAX_KH '+format(kh,'.17g')+' MAX_ASPECT 2 MAX_CELLS 1000\nEND\n')
    if spec.get('material_split'):
        # Explicit cubic physical boxes preserve one common pitch across materials.
        # The legacy per-box mesher may otherwise choose a different pitch for ER=2/4.
        frequency=spec.get('witness_frequency_mhz',300)
        pitch=.025
        strongest=math.sqrt(abs(complex(4,-.03/(2*math.pi*frequency*1e6*8.8541878128e-12))))
        lines=['NEC2C-DIELECTRIC 2']
        for index,(_,er,sigma,start,stop) in enumerate(material_boxes(spec),1):
            lines += ['MATERIAL %d ER %g SIGMA %g'%(index,er,sigma),'REGION %d MATERIAL %d'%(index,index)]
            counts=[round((b-a)/25) for a,b in zip(start,stop)]
            for x in range(counts[0]):
                for y in range(counts[1]):
                    for z in range(counts[2]):
                        lo=[(a+25*i)*.001 for a,i in zip(start,(x,y,z))]
                        hi=[(a+25*(i+1))*.001 for a,i in zip(start,(x,y,z))]
                        lines.append('BOX %d '%index+' '.join(format(v,'.17g') for v in lo+hi))
        if not material_boxes(spec):
            lines+=['MATERIAL 1 ER 1','REGION 1 MATERIAL 1','BOX 1 .05 -.05 0 .15 .05 .65']
        lines+=['MESH MAX_KH %.17g MAX_ASPECT 1 MAX_CELLS 10000'%(2*math.pi*frequency*1e6/299792458*strongest*pitch),'END']
        (root/'reference.diel').write_text('\n'.join(lines)+'\n')


def worker(spec):
    import numpy as np
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    root=Path(spec['directory']); root.mkdir(parents=True,exist_ok=True)
    write_nec_reference(root,spec)
    fdtd=openEMS(NrTS=spec['steps'],EndCriteria=spec['end_criterion'],CellConstantMaterial=True)
    fdtd.SetGaussExcite(300e6,spec['bandwidth_mhz']*1e6)
    boundaries=['PML_8']*6
    if spec['ground']=='pec': boundaries[4]='PEC'
    fdtd.SetBoundaryCond(boundaries)
    csx=ContinuousStructure(); fdtd.SetCSX(csx)
    grid=csx.GetGrid(); grid.SetDeltaUnit(.001)
    margin=spec.get('domain_margin_mm',0)
    gap=spec.get('feed_gap_mm',10)
    zmin=0 if spec['ground']=='pec' else -500-margin
    wire=csx.AddMetal('dipole')
    # Preserve the 475 mm outer span and 600 mm center height. The optional
    # feed variant changes only the source gap and adjoining conductor ends.
    wire.AddCylinder([0,0,362.5],[0,0,600-gap/2],radius=spec['wire_radius_mm'],priority=10)
    wire.AddCylinder([0,0,600+gap/2],[0,0,837.5],radius=spec['wire_radius_mm'],priority=10)
    radius=spec['wire_radius_mm']
    port=fdtd.AddLumpedPort(port_nr=1,R=50,start=[-radius,-radius,600-gap/2],stop=[radius,radius,600+gap/2],
                          p_dir='z',excite=True,priority=100)
    if spec['ground']=='real':
        soil=csx.AddMaterial('soil',epsilon=13,kappa=.02)
        soil.AddBox([-650-margin,-650-margin,zmin],[650+margin,650+margin,0],priority=1)
    for name,er,sigma,start,stop in material_boxes(spec):
        cube=csx.AddMaterial(name,epsilon=er,kappa=sigma)
        cube.AddBox(start,stop,priority=20)
    axes=build_axes(spec['mesh_mm'],spec['wire_cell_mm'],spec['ground'],spec['wire_radius_mm'],margin,gap,spec.get('cube_bottom_mm',550),spec.get('material_split',''))
    for direction,axis in zip('xyz',axes): grid.AddLine(direction,axis)
    lines=[list(map(float,grid.GetLines(axis))) for axis in 'xyz']
    metrics=mesh_metrics(lines)
    cells=metrics['cells']
    report={'spec':spec,**metrics,
            'requested_max_steps':spec['steps'],
            'estimated_field_storage_mib':cells*64*8/1024**2,
            'storage_estimate_note':'64 doubles/cell heuristic, not a total-memory guarantee',
            'validation_status':'unreviewed; check solver decay and mesh/domain convergence'}
    csx.Write2XML(str(root/'geometry.xml'))
    (root/'setup.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if cells>spec['max_cells']:
        raise RuntimeError('mesh exceeds --max-cells; inspect setup.json before choosing another limit')
    if spec['action']=='prepare': return
    (root/'solver').mkdir(exist_ok=True)
    code=fdtd.Run(str(root/'solver'),cleanup=True,verbose=1,numThreads=spec['threads'])
    if code not in (None,0): raise RuntimeError('OpenEMS returned failure: '+str(code))
    if spec['action']=='pilot':
        report['validation_status']='pilot_only; timing/setup diagnostic, no impedance acceptance'
    else:
        frequencies=np.array([280e6,300e6,320e6])
        port.CalcPort(str(root/'solver'),frequencies)
        impedance=port.uf_tot/port.if_tot
        if not np.all(np.isfinite(impedance)): raise RuntimeError('nonfinite port impedance')
        report['impedance']=[{'frequency_mhz':float(f/1e6),'r_ohm':float(z.real),'x_ohm':float(z.imag)}
                             for f,z in zip(frequencies,impedance)]
    (root/'result.json').write_text(json.dumps(report,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('plan','prepare','pilot','verify'),default='plan')
    p.add_argument('--ground',choices=('pec','real'),default='real')
    p.add_argument('--mesh-mm',type=float)
    p.add_argument('--wire-radius-mm',type=float,default=5,help='physical radius, default 5 mm (10 mm diameter)')
    p.add_argument('--wire-cell-mm',type=float,help='near-wire transverse spacing; defaults to wire radius')
    p.add_argument('--cube-bottom-mm',type=float,choices=(0,550),default=550,help='550: M10 cube; 0: M11 grounded pillar (PEC only)')
    p.add_argument('--material-split',choices=('', 'x','z'),default='',help='M13-03: split contact pillar at x=100 or z=325 mm; second material ER=4, sigma=.03 S/m')
    p.add_argument('--loaded-epsilon-r',type=float,default=3,help='loaded scalar permittivity; default 3 preserves previous identities')
    p.add_argument('--witness-frequency-mhz',type=float,default=0,help='M13-02 witness identity; only this frequency is a material-table comparison')
    p.add_argument('--dielectric-sigma',type=float,default=0,help='loaded dielectric conductivity in S/m (M13-01); vacuum remains zero')
    p.add_argument('--pilot-epsilon-r',type=float,default=1,help='pilot/prepare material, default vacuum')
    p.add_argument('--domain-margin-mm',type=float,default=0,help='additional side/top clearance, default 0 mm')
    p.add_argument('--feed-gap-mm',type=float,default=10,help='finite feed gap, default 10 mm')
    p.add_argument('--bandwidth-mhz',type=float,default=240,help='Gaussian bandwidth, default 240 MHz')
    p.add_argument('--seconds',type=float,default=60,help='total wall budget for this invocation, including setup')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--max-cells',type=int,default=2000000)
    p.add_argument('--workdir',type=Path,default=Path('m10-openems-results'))
    p.add_argument('--analyze',type=Path,help='inspect returned XML and logs without OpenEMS')
    p.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=p.parse_args()
    if args.analyze:
        print(json.dumps(inspect_results(args.analyze),indent=2)); return 0
    if args.worker:
        worker(json.loads(args.worker.read_text())); return 0
    if not math.isfinite(args.seconds) or args.seconds<=0: p.error('--seconds must be positive and finite')
    if args.threads<1 or args.max_cells<1: p.error('--threads and --max-cells must be positive')
    mesh=args.mesh_mm if args.mesh_mm is not None else (30 if args.mode in ('plan','pilot','prepare') else 20)
    if not math.isfinite(mesh) or mesh<=0: p.error('--mesh-mm must be positive and finite')
    wire_cell=args.wire_cell_mm if args.wire_cell_mm is not None else args.wire_radius_mm
    if not math.isfinite(args.wire_radius_mm) or not 0<args.wire_radius_mm<=20: p.error('--wire-radius-mm must be in (0,20]')
    if not math.isfinite(wire_cell) or not 0<wire_cell<=50: p.error('--wire-cell-mm must be in (0,50]')
    if not math.isfinite(args.bandwidth_mhz) or not 20<=args.bandwidth_mhz<=300: p.error('--bandwidth-mhz must be between 20 and 300')
    if not math.isfinite(args.domain_margin_mm) or not 0<=args.domain_margin_mm<=1000: p.error('--domain-margin-mm must be in [0,1000]')
    if not math.isfinite(args.feed_gap_mm) or not 2.5<=args.feed_gap_mm<=20: p.error('--feed-gap-mm must be in [2.5,20]')
    if args.cube_bottom_mm==0 and args.ground!='pec': p.error('contact fixture requires --ground pec')
    if args.material_split and (args.ground!='pec' or args.cube_bottom_mm!=0): p.error('split fixture requires PEC contact pillar')
    metrics=mesh_metrics(build_axes(mesh,wire_cell,args.ground,args.wire_radius_mm,args.domain_margin_mm,args.feed_gap_mm,args.cube_bottom_mm,args.material_split))
    plan={'mesh_revision':MESH_REVISION,'wire_radius_mm':args.wire_radius_mm,'wire_cell_mm':wire_cell,'bandwidth_mhz':args.bandwidth_mhz,
          'cube_bottom_mm':args.cube_bottom_mm,'pilot_epsilon_r':args.pilot_epsilon_r,'domain_margin_mm':args.domain_margin_mm,'feed_gap_mm':args.feed_gap_mm,'predicted_mesh':metrics,'ground':args.ground,'mesh_mm':mesh,'total_wall_budget_seconds':args.seconds,
          'threads':args.threads,'max_cells':args.max_cells,'dielectric_sigma_s_m':args.dielectric_sigma,'loaded_epsilon_r':args.loaded_epsilon_r,'witness_frequency_mhz':args.witness_frequency_mhz,
          'pilot':'one epsilon_r='+str(args.pilot_epsilon_r)+' run, full source with 120000-step ceiling and external wall budget; not accuracy evidence',
          'verify':'vacuum and er='+str(args.loaded_epsilon_r)+' runs, at most 120000 steps each; results require convergence review',
          'remaining_gate':('M13-02: compare only designated witness frequency; audit decay and mesh/domain convergence' if args.witness_frequency_mhz else ('M13-01: review loss, decay, mesh/domain refinement and native comparisons' if args.dielectric_sigma else ('M11: review contact decay, mesh refinement and native comparisons' if args.cube_bottom_mm==0 else 'M10 fixture: refine PEC meshes and assess feed/domain error; real ground is preliminary M13'))),
          'runtime_estimate':'returned 5mm-radius laptop results: pilot 7.63 s, coarse 20mm about 21 s/run, fine 15mm/2.5mm about 90 s/run; other setups unmeasured' }
    if args.material_split:
        plan.update(material_split=args.material_split,secondary_material={'epsilon_r':4,'conductivity_s_m':.03},remaining_gate='M13-03: audit both materials, decay, interface mesh, domain and native impedance convergence')
    print(json.dumps(plan,indent=2),flush=True)
    if not math.isfinite(args.dielectric_sigma) or args.dielectric_sigma<0: p.error('--dielectric-sigma must be finite and nonnegative')
    if not math.isfinite(args.loaded_epsilon_r) or args.loaded_epsilon_r<=1: p.error('--loaded-epsilon-r must be finite and >1')
    if not math.isfinite(args.pilot_epsilon_r) or args.pilot_epsilon_r<1: p.error('--pilot-epsilon-r must be finite and >=1')
    if not math.isfinite(args.witness_frequency_mhz) or args.witness_frequency_mhz<0: p.error('invalid witness frequency')
    if args.mode=='plan': return 0
    if os.name!='posix': p.error('bounded worker termination currently requires POSIX (Linux/macOS)')
    root=args.workdir.resolve(); root.mkdir(parents=True,exist_ok=True)
    started=time.monotonic(); summary=[]
    for epsilon in ((args.pilot_epsilon_r,) if args.mode in ('pilot','prepare') else (1,args.loaded_epsilon_r)):
        if float(epsilon).is_integer(): epsilon=int(epsilon)
        spec={'mesh_revision':MESH_REVISION,'feed_model':'finite_box',
              'wire_radius_mm':args.wire_radius_mm,
              'wire_cell_mm':wire_cell,'bandwidth_mhz':args.bandwidth_mhz,'ground':args.ground,'epsilon_r':epsilon,'mesh_mm':mesh,'threads':args.threads,
              'max_cells':args.max_cells,'steps':120000,
              'end_criterion':1e-4 if args.mode=='pilot' else 1e-6,'action':args.mode}
        # Preserve existing default configuration identities and completed caches.
        if args.material_split: spec['material_split']=args.material_split
        if args.witness_frequency_mhz: spec['witness_frequency_mhz']=args.witness_frequency_mhz
        if args.dielectric_sigma and epsilon!=1: spec['dielectric_sigma_s_m']=args.dielectric_sigma
        if args.cube_bottom_mm!=550: spec['cube_bottom_mm']=args.cube_bottom_mm
        if args.domain_margin_mm: spec['domain_margin_mm']=args.domain_margin_mm
        if args.feed_gap_mm!=10: spec['feed_gap_mm']=args.feed_gap_mm
        identity=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
        directory=root/(args.ground+'-'+args.mode+'-'+identity)
        directory.mkdir(parents=True,exist_ok=True)
        spec['directory']=str(directory)
        specfile=directory/'spec.json'; specfile.write_text(json.dumps(spec,indent=2)+'\n')
        resultfile=directory/'result.json'
        # Cache only an exact configuration, never a previous timed-out attempt.
        if resultfile.exists() and json.loads(resultfile.read_text())['spec']==spec:
            status={'status':'cached','directory':str(directory)}
        else:
            remaining=args.seconds-(time.monotonic()-started)
            if remaining<=0:
                summary.append({'status':'budget_exhausted','directory':str(directory)}); break
            status=run_worker([sys.executable,str(Path(__file__).resolve()),'--worker',str(specfile)],
                              remaining,directory/'worker.log')
            status['directory']=str(directory)
        summary.append(status)
        if status['status'] not in ('returned','cached'): break
    result={'plan':plan,'jobs':summary,'wall_seconds':time.monotonic()-started}
    (root/'last-invocation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    # Export paired impedance differences only when both exact configurations exist.
    reports=[]
    for job in summary:
        file=Path(job['directory'])/'result.json'
        if file.exists() and job['status'] in ('returned','cached'): reports.append(json.loads(file.read_text()))
    if args.mode=='verify' and len(reports)==2:
        suffix='' if args.domain_margin_mm==0 and args.feed_gap_mm==10 and args.cube_bottom_mm==550 else '-bottom'+str(args.cube_bottom_mm)+'-margin'+str(args.domain_margin_mm)+'-gap'+str(args.feed_gap_mm)
        out=root/('impedance-'+args.ground+'-'+str(mesh)+'mm'+suffix+'.csv')
        rows=[]
        for a,b in zip(reports[0]['impedance'],reports[1]['impedance']):
            rows.append({'frequency_mhz':a['frequency_mhz'],'vacuum_r_ohm':a['r_ohm'],
                         'vacuum_x_ohm':a['x_ohm'],'loaded_r_ohm':b['r_ohm'],'loaded_x_ohm':b['x_ohm'],
                         'delta_r_ohm':b['r_ohm']-a['r_ohm'],'delta_x_ohm':b['x_ohm']-a['x_ohm'],
                         'mesh_mm':mesh,'ground':args.ground,'validation_status':'unreviewed'})
        with out.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    return 0 if summary and all(s['status'] in ('returned','cached') for s in summary) else 2

if __name__=='__main__':
    sys.exit(main())
