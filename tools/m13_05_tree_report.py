#!/usr/bin/env python3
"""Report explicitly assumed material scenarios without inventing confidence intervals."""
import argparse
import json
import math
from pathlib import Path


def report(data):
    for key in ('geometry_assumptions','material_provenance','excluded_effects'):
        if not isinstance(data.get(key),list) or not data[key] or not all(isinstance(x,str) and x.strip() for x in data[key]):
            raise ValueError('nonempty text list required: '+key)
    if data.get('ground') not in ('free','PEC','GN2-positive-gap-wires'):
        raise ValueError('unsupported ground scope')
    rows=data.get('scenarios',[])
    if len(rows)<2:raise ValueError('at least two material scenarios required')
    names=set();mesh=None
    for row in rows:
        name=row['name']
        if not isinstance(name,str) or not name or name in names:raise ValueError('unique scenario names required')
        names.add(name)
        if not isinstance(row['z_ohm'],list) or len(row['z_ohm'])!=2:raise ValueError('two impedance components required')
        values=[row[k] for k in ('frequency_mhz','epsilon_r','sigma_s_m','pitch_m','true_relative_residual','power_defect','absorbed_power_w')]+row['z_ohm']
        if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in values):raise ValueError('finite numerical results required')
        if row['frequency_mhz']<=0 or row['epsilon_r']<1 or row['sigma_s_m']<0 or row['pitch_m']<=0 or any(row[k]<0 for k in ('true_relative_residual','power_defect','absorbed_power_w')):raise ValueError('passive material and valid diagnostics required')
        if row['true_relative_residual']>1e-8 or row['power_defect']>=.003:raise ValueError('unconverged numerical scenario')
        signature=(row['frequency_mhz'],row['pitch_m'])
        if mesh is not None and mesh!=signature:raise ValueError('material scenarios must share frequency and mesh pitch')
        mesh=signature
    if data.get('nominal') not in names:raise ValueError('nominal scenario missing')
    nominal=next(x for x in rows if x['name']==data['nominal'])
    ranges={key:[min(r[key] for r in rows),max(r[key] for r in rows)] for key in ('epsilon_r','sigma_s_m','absorbed_power_w')}
    ranges['r_ohm']=[min(r['z_ohm'][0] for r in rows),max(r['z_ohm'][0] for r in rows)]
    ranges['x_ohm']=[min(r['z_ohm'][1] for r in rows),max(r['z_ohm'][1] for r in rows)]
    result=dict(schema='nec2c-tree-material-report-1',assumptions={k:data[k] for k in ('geometry_assumptions','material_provenance','excluded_effects','ground')},nominal=nominal,scenarios=rows,sampled_ranges=ranges,
                interpretation='Envelope of supplied scenarios only; not a confidence interval, probability distribution, or guaranteed bound over unsampled materials. Numerical residual and power defect do not bound material or geometry error.')
    if 'binary_sha256' in data:
        digest=data['binary_sha256']
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):raise ValueError('invalid solver binary SHA256')
        result['binary_sha256']=digest
    if 'numerical_refinement' in data:
        refined=data['numerical_refinement']
        if any(refined[k]!=nominal[k] for k in ('epsilon_r','sigma_s_m','frequency_mhz')) or not 0<refined['pitch_m']<nominal['pitch_m']:
            raise ValueError('refinement must preserve nominal material/frequency and reduce pitch')
        validation=dict(data);validation.pop('numerical_refinement')
        validation['scenarios']=[nominal,dict(refined,pitch_m=nominal['pitch_m'])]
        report(validation)
        result['numerical_refinement']=refined
        result['refinement_component_change_ohm']=[abs(a-b) for a,b in zip(refined['z_ohm'],nominal['z_ohm'])]
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=report(json.loads(a.input.read_text()));a.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
