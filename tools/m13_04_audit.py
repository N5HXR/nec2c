#!/usr/bin/env python3
"""Audit M13-04 returned physics and compare completed native benchmark sets."""
import argparse,hashlib,json,re,sys,xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from m13_04_review import review_controls
from openems_m10_ground import inspect_results
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('results',type=Path)
parser.add_argument('--archive',type=Path,required=True)
parser.add_argument('--native',type=Path,nargs='+',required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
root=args.results
r=review_controls(root);meshes=inspect_results(root);audits=[]
for p in sorted(root.rglob('spec.json')):
 s=json.loads(p.read_text())
 if s['action']!='verify':continue
 log=p.with_name('worker.log').read_text();match=re.search(r'end-criteria of\s+(-?[\d.]+)dB reached after\s+(\d+) timesteps\s+\((-?[\d.]+)dB\)',log)
 assert match and float(match[1])<=-60 and float(match[3])<=-60 and int(match[2])<s['steps'],str(p)
 assert 'Cutting to max number of timesteps' not in log and 'before the end-criteria' not in log
 x=ET.parse(p.with_name('geometry.xml')).getroot();materials={m.get('Name'):m for m in x.findall('.//Material')}
 def check_material(name,er,sigma,lo,hi):
  m=materials[name];assert m.get('Isotropy')=='1';prop=m.find('Property')
  assert float(prop.get('Epsilon').split(',')[0])==er and float(prop.get('Kappa').split(',')[0])==sigma
  assert float(prop.get('Mue').split(',')[0])==1
  b=m.find('Primitives/Box');actual=[[float(b.find(q).get(a)) for a in 'XYZ'] for q in ('P1','P2')]
  assert actual==[lo,hi],(actual,lo,hi)
 margin=s.get('domain_margin_mm',0);side=650+margin;depth=500+margin
 check_material('soil',13,.02,[-side,-side,-depth],[side,side,0])
 if s['epsilon_r']==3:check_material('cube',3,s.get('dielectric_sigma_s_m',0),[50,-50,550],[150,50,650])
 else:assert 'cube' not in materials
 cylinders=x.findall('.//Metal[@Name="dipole"]/Primitives/Cylinder');assert len(cylinders)==2
 assert all(float(c.get('Radius'))==5 for c in cylinders)
 endpoints=[[float(c.find(q).get(a)) for a in 'XYZ'] for c in cylinders for q in ('P1','P2')]
 assert endpoints==[[0,0,362.5],[0,0,595],[0,0,605],[0,0,837.5]]
 audits.append(dict(directory=str(p.parent.relative_to(root)),target_db=float(match[1]),final_db=float(match[3]),steps=int(match[2]),xml_material_geometry_pass=True))
assert len(audits)==12 and len(meshes)==14
for m in meshes:
 assert m['cells']<2000000 and min(m['minimum_step_mm'])>=2.5-1e-9 and max(m['maximum_adjacent_ratio'])<=2+1e-9
 assert not m['source_was_clipped'] and not m['step_cap_before_decay']
assert all(p['status']=='completed_pair' for p in r['pairs'])
def pair(control,f,mesh,margin):return next(p for p in r['pairs'] if (p['experiment_case'],p['frequency_mhz'],p['mesh_mm'],p['domain_margin_mm'])==(control,f,mesh,margin))
convergence=[]
for c in ('lossless','lossy'):
 for f in (280,300,320):
  coarse=pair(c,f,20,0);fine=pair(c,f,15,0);large=pair(c,f,15,150)
  mesh=[abs(coarse[k]-fine[k]) for k in ('delta_r_ohm','delta_x_ohm')];domain=[abs(large[k]-fine[k]) for k in ('delta_r_ohm','delta_x_ohm')]
  assert max(mesh+domain)<2
  convergence.append(dict(control=c,frequency_mhz=f,mesh_component_change_ohm=mesh,domain_component_change_ohm=domain))
native=[json.loads(p.read_text()) for d in args.native if d.is_dir() for p in sorted(d.glob('*/result.json'))]
assert len(native)==6
comparisons=[]
for n in native:
 reference=pair('lossy' if n['dielectric_sigma_s_m'] else 'lossless',n['frequency_mhz'],15,150)
 ref=[reference['delta_r_ohm'],reference['delta_x_ohm']];errors=[abs(a-b) for a,b in zip(n['delta_z_ohm'],ref)];limits=[max(2,.02*abs(a)) for a in ref]
 assert all(a<=b for a,b in zip(errors,limits))
 assert all(a['true_relative_residual']<=1e-10 and a['power_defect']<.003 for a in n['runs'])
 comparisons.append(dict(frequency_mhz=n['frequency_mhz'],cells=n['divisions_per_axis']**3,wire_segments=n['wire_segments'],sigma_s_m=n['dielectric_sigma_s_m'],native_delta_z_ohm=n['delta_z_ohm'],openems_delta_z_ohm=ref,component_errors_ohm=errors,component_limits_ohm=limits,passed=True))
r.update(acceptance='PASS: all 12 verification decay/material/geometry audits, mesh/domain checks, and six native component comparisons pass',archive_sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(),physical_audits=audits,mesh_audits=meshes,convergence=convergence,native_comparisons=comparisons)
args.output.write_text(json.dumps(r,indent=2)+'\n')
print('PASS worst component error',max(max(c['component_errors_ohm']) for c in comparisons))
print('max mesh/domain change',max(max(c['mesh_component_change_ohm']) for c in convergence),max(max(c['domain_component_change_ohm']) for c in convergence))
