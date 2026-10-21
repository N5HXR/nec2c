#!/usr/bin/env python3
"""Material uncertainty reporting must not hide failed solves or mix meshes."""
import copy
import json
import os
from pathlib import Path
import sys
ROOT=Path(os.environ.get('TEST_SRCDIR','.')).resolve();sys.path.insert(0,str(ROOT/'tools'))
from m13_05_tree_report import report

def main():
    row=dict(name='nominal',epsilon_r=3,sigma_s_m=.01,frequency_mhz=100,pitch_m=.1,z_ohm=[10,20],true_relative_residual=1e-10,power_defect=.001,absorbed_power_w=.01)
    other=dict(row,name='alternate',epsilon_r=4,z_ohm=[8,23])
    data=dict(ground='PEC',nominal='nominal',geometry_assumptions=['Synthetic boxes'],material_provenance=['Assumed values, not measured'],excluded_effects=['Anisotropy'],scenarios=[row,other])
    result=report(data)
    assert result['sampled_ranges']['r_ohm']==[8,10] and result['sampled_ranges']['x_ohm']==[20,23]
    assert 'not a confidence interval' in result['interpretation']
    refined=copy.deepcopy(data);refined['numerical_refinement']=dict(row,name='fine',pitch_m=.05,z_ohm=[10.1,19.8])
    fine_report=report(refined)
    assert fine_report['sampled_ranges']==result['sampled_ranges']
    assert abs(fine_report['refinement_component_change_ohm'][1]-.2)<1e-12
    for key,value in [('epsilon_r',4),('pitch_m',.2),('true_relative_residual',1),('frequency_mhz',200)]:
        bad=copy.deepcopy(refined);bad['numerical_refinement'][key]=value
        try:report(bad)
        except ValueError:pass
        else:raise AssertionError(('refinement',key))
    for key,value in [('pitch_m',.05),('frequency_mhz',200),('sigma_s_m',-.01),('epsilon_r',.5),('true_relative_residual',.1),('power_defect',.1),('absorbed_power_w',float('nan')),('name','nominal'),('z_ohm',[1])]:
        bad=copy.deepcopy(data);bad['scenarios'][1][key]=value
        try:report(bad)
        except ValueError:pass
        else:raise AssertionError((key,value))
    for key in ('geometry_assumptions','material_provenance','excluded_effects'):
        bad=copy.deepcopy(data);bad[key]=[]
        try:report(bad)
        except ValueError:pass
        else:raise AssertionError(key)
    # Retained physical witness must reproduce the published report exactly.
    source=json.loads((ROOT/'tests/m13-05-tree-input.json').read_text())
    expected=json.loads((ROOT/'tests/m13-05-tree-report.json').read_text())
    assert report(source)==expected
    assert len(expected['scenarios'])==5
    print('Material assumptions, sampled envelopes, convergence refusals and retained report reproduction pass')

if __name__=='__main__':main()
