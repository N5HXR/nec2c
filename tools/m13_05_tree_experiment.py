#!/usr/bin/env python3
"""Reproduce the small synthetic PEC-tree material-sensitivity witness."""
import argparse,hashlib,json,math,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from m10_pec_acceptance import scalar,impedance
from m13_05_tree_report import report
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--binary',type=Path,default=Path('./nec2c'));parser.add_argument('--workdir',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True);args=parser.parse_args()
binary=args.binary.resolve();root=args.workdir.resolve();root.mkdir(parents=True,exist_ok=True);args.output_dir.mkdir(parents=True,exist_ok=True)
rows=[]
for name,er,sigma in [('nominal',3,.01),('er-low',2,.01),('er-high',4,.01),('sigma-low',3,.005),('sigma-high',3,.02),('refined',3,.01)]:
 d=root/name;d.mkdir(exist_ok=True);pitch=.05 if name=='refined' else .1
 kh=2*math.pi*100e6/299792458*math.sqrt(abs(complex(er,-sigma/(2*math.pi*100e6*8.8541878128e-12))))*pitch*(1+1e-12)
 (d/'tree.diel').write_text(f'NEC2C-DIELECTRIC 2\nMATERIAL 1 ER {er} SIGMA {sigma}\nREGION 1 MATERIAL 1\nBOX 1 -.05 -.05 0 .05 .05 .8\nBOX 1 .05 -.05 .5 .25 .05 .6\nBOX 1 -.05 .05 .6 .05 .25 .7\nBOX 1 -.25 -.25 .8 .35 .35 1.2\nMESH MAX_KH {kh:.17g} MAX_ASPECT 1 MAX_CELLS 2000\nEND\n')
 (d/'tree.nec').write_text('CE\nGW 1 11 -.5 0 .35 -.5 0 .85 .002\nGE\nFR 0 1 0 0 100 0\nGN 1\nEX 0 1 6 0 1 0\nDF "tree.diel"\nXQ\nEN\n')
 with (d/'run.log').open('w') as f:subprocess.run([str(binary),'-s','fft','-r','1e-10','-i',str(d/'tree.nec'),'-o',str(d/'tree.out')],stdout=f,stderr=subprocess.STDOUT,check=True)
 t=(d/'tree.out').read_text();z=impedance(t)
 row=dict(name=name,epsilon_r=er,sigma_s_m=sigma,frequency_mhz=100,pitch_m=pitch,z_ohm=[z.real,z.imag],true_relative_residual=scalar(t,'GMRES TRUE RELATIVE RESIDUAL:'),power_defect=scalar(t,'RELATIVE POWER DEFECT:'),absorbed_power_w=scalar(t,'DIELECTRIC ABSORBED POWER (W):'),active_cells=int(scalar(t,'ACTIVE CELLS:')))
 rows.append(row);print(row,flush=True)
 (root/'runs.json').write_text(json.dumps(rows,indent=2)+'\n')
data=dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),numerical_refinement=rows[-1],ground='PEC',nominal='nominal',geometry_assumptions=['Synthetic rectilinear trunk, two branches and homogeneous crown; dimensions in archived input decks','Positive wire clearance; tree touches ideal PEC ground'],material_provenance=['ER 2/3/4 and sigma .005/.01/.02 S/m are arbitrary sensitivity witnesses, not measured wood or foliage data','Constant isotropic nonmagnetic material at 100 MHz; same material throughout this illustrative tree'],excluded_effects=['Anisotropy, spatial moisture variation and measured frequency dispersion','Leaves, bark, roots and unresolved branch structure','Real-ground effects; GN 2 has separate small-fixture acceptance and no FFT support'],scenarios=[r for r in rows if r['name']!='refined'])
(args.output_dir/'m13-05-tree-input.json').write_text(json.dumps(data,indent=2)+'\n')
result=report(data)
(args.output_dir/'m13-05-tree-report.json').write_text(json.dumps(result,indent=2)+'\n')
