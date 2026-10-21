#!/usr/bin/env python3
"""Review the M13-04 control sets separately; never cross-pair reused vacuum jobs."""
import argparse
import csv
import json
from pathlib import Path
from m10_review_results import review


def review_controls(root):
    controls=[]
    for label,sigma in (('lossless',0),('lossy',.02)):
        candidates=[p for p in Path(root).rglob(label) if p.is_dir() and any(p.rglob('spec.json'))]
        if Path(root).name==label and any(Path(root).rglob('spec.json')):candidates.append(Path(root))
        if len(candidates)!=1:raise ValueError('expected exactly one '+label+' result directory')
        control=review(candidates[0])
        if len(control['runs'])!=6 or len(control['pairs'])!=9:raise ValueError('incomplete '+label+' verification set')
        for run in control['runs']:
            if run['ground']!='real' or run['cube_bottom_mm']!=550 or run['epsilon_r'] not in (1,3):
                raise ValueError('unexpected M13-04 fixture: '+str(run))
            if run['epsilon_r']==3 and run['dielectric_sigma_s_m']!=sigma:raise ValueError('wrong control material')
            run['experiment_case']=label
        for pair in control['pairs']:pair['experiment_case']=label
        controls.append(dict(experiment_case=label,**control))
    return dict(experiment='M13-04',controls=controls,
                runs=[r for c in controls for r in c['runs']],pairs=[p for c in controls for p in c['pairs']],
                acceptance='IN REVIEW: DFT reproduction only; explicit decay, soil/cube XML, mesh/domain and native accuracy must still be audited')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();result=review_controls(args.results);args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'review.json').write_text(json.dumps(result,indent=2)+'\n')
    with (args.output_dir/'pairs.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(result['pairs'][0]));writer.writeheader();writer.writerows(result['pairs'])

if __name__=='__main__':main()
