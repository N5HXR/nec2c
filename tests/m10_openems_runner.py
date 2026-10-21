#!/usr/bin/env python3
"""Check the laptop runner's resource controls without OpenEMS installed."""
import importlib.util
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import time

ROOT=Path(os.environ.get('TEST_SRCDIR','.')).resolve()
module_spec=importlib.util.spec_from_file_location('m10_runner',ROOT/'tools/openems_m10_ground.py')
runner=importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(runner)


def main():
    plan=subprocess.run([sys.executable,str(ROOT/'tools/openems_m10_ground.py'),'--mode','plan'],
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    assert plan.returncode==0,plan.stderr
    assert b'unmeasured' in plan.stdout
    bad=subprocess.run([sys.executable,str(ROOT/'tools/openems_m10_ground.py'),'--seconds','-1'],
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    assert bad.returncode!=0
    for bulk in (15,20,30):
        for spacing in (1.25,2.5,5):
            for ground in ('pec','real'):
                axes=runner.build_axes(bulk,spacing,ground)
                metrics=runner.mesh_metrics(axes)
                assert min(metrics['minimum_step_mm'])>=min(spacing,10)-1e-9,metrics
                assert max(metrics['maximum_adjacent_ratio'])<=2+1e-9,metrics
                assert metrics['cells']<2000000,metrics
                assert all(any(abs(v-anchor)<1e-9 for v in axes[2]) for anchor in (0,362.5,550,595,605,650,837.5,1400))
                assert all(abs(a+b)<1e-8 for a,b in zip(axes[1],reversed(axes[1])))
    contact_plan=subprocess.run([sys.executable,str(ROOT/'tools/openems_m10_ground.py'),'--ground','real','--cube-bottom-mm','0'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    assert contact_plan.returncode!=0 and b'contact fixture requires' in contact_plan.stderr
    metrics=runner.mesh_metrics(runner.build_axes(30,5,'real'))
    assert metrics['grid_cells']==[51,50,70],metrics
    assert metrics['minimum_step_mm']==[5,5,10],metrics
    for margin,gap in ((200,10),(0,5)):
        axes=runner.build_axes(20,5,'pec',domain_margin_mm=margin,feed_gap_mm=gap)
        assert axes[0][0]==-650-margin and axes[0][-1]==650+margin
        assert axes[2][-1]==1400+margin
        assert 600-gap/2 in axes[2] and 600+gap/2 in axes[2]
        metrics=runner.mesh_metrics(axes)
        assert min(metrics['minimum_step_mm'])>=5-1e-9,metrics
        assert max(metrics['maximum_adjacent_ratio'])<=2+1e-9,metrics
        assert metrics['cells']<2000000,metrics
    with tempfile.TemporaryDirectory(prefix='m10-mesh-test-') as tmp:
        d=Path(tmp)
        runner.write_nec_reference(d,{'ground':'real','epsilon_r':3,'wire_radius_mm':5})
        reference=(d/'reference.nec').read_text()
        assert '0.005' in reference and 'GN 2 0 0 0 13 .02' in reference
        assert 'DF ' not in (d/'reference-vacuum-legacy.nec').read_text()
        runner.write_nec_reference(d,{'ground':'pec','epsilon_r':3,'wire_radius_mm':5,'cube_bottom_mm':0})
        assert 'BOX 1 .05 -.05 0.0 .15 .05 .65' in (d/'reference.diel').read_text()
        assert 'FR 0 1 0 0 300 0' in (d/'reference.nec').read_text()
        assert 'GN 1' in (d/'reference.nec').read_text()
        runner.write_nec_reference(d,{'ground':'pec','epsilon_r':3,'wire_radius_mm':5,'cube_bottom_mm':0,'dielectric_sigma_s_m':.02})
        assert 'nec2c-dielectric 2' in (d/'reference.diel').read_text()
        assert 'SIGMA 0.02' in (d/'reference.diel').read_text()
        runner.write_nec_reference(d,{'ground':'pec','epsilon_r':3.5,'wire_radius_mm':5,'cube_bottom_mm':0,'dielectric_sigma_s_m':.015,'witness_frequency_mhz':280})
        assert 'FR 0 1 0 0 280 0' in (d/'reference.nec').read_text()
        assert 'ER 3.5 SIGMA 0.015' in (d/'reference.diel').read_text()
        for split in ('x','z'):
            cfg={'ground':'pec','epsilon_r':2,'wire_radius_mm':5,'cube_bottom_mm':0,
                 'dielectric_sigma_s_m':.01,'material_split':split}
            boxes=runner.material_boxes(cfg)
            assert len(boxes)==2 and boxes[0][1:3]==(2,.01) and boxes[1][1:3]==(4,.03)
            axis=0 if split=='x' else 2
            assert boxes[0][4][axis]==boxes[1][3][axis]
            runner.write_nec_reference(d,cfg)
            resource=(d/'reference.diel').read_text()
            assert resource.count('BOX ')==416 and 'ER 4 SIGMA 0.03' in resource
            for bulk,spacing,margin in ((20,5,0),(15,2.5,0),(15,2.5,150)):
                axes=runner.build_axes(bulk,spacing,'pec',domain_margin_mm=margin,
                                       cube_bottom_mm=0,material_split=split)
                assert (100 in axes[0]) if split=='x' else (325 in axes[2])
                metrics=runner.mesh_metrics(axes)
                assert metrics['cells']<2000000 and min(metrics['minimum_step_mm'])>=2.5-1e-9
                assert max(metrics['maximum_adjacent_ratio'])<=2+1e-9
            assert runner.material_boxes(dict(cfg,epsilon_r=1))==[]
        # Returned engine time is authoritative; vacuum CFL alone is not a runtime forecast.
        (d/'worker.log').write_text('Timestep (s) : 1.22394e-12\n'
            'Requested excitation pusle would be 26007 timesteps or 3.18311e-08 s long. Cutting to max number of timesteps!\n'
            'Time for 5000 iterations with 230400.00 cells : 20.29 sec\n')
        (d/'geometry.xml').write_text('<ContinuousStructure><RectilinearGrid DeltaUnit=".001">'+
            ''.join('<'+a+'Lines>0,5,10</'+a+'Lines>' for a in 'XYZ')+
            '</RectilinearGrid></ContinuousStructure>')
        audit=runner.inspect_results(d)[0]
        assert abs(audit['projected_pair_120000_step_seconds']-973.92)<1e-8,audit
        assert audit['source_was_clipped']
    print('Mesh anchors, 2:1 grading, minimum spacing, radius-matched NEC references and pilot-log timing audit pass')
    review_spec=importlib.util.spec_from_file_location('m10_review',ROOT/'tools/m10_review_results.py')
    reviewer=importlib.util.module_from_spec(review_spec);review_spec.loader.exec_module(reviewer)
    import cmath,math,json
    frequency=300.0
    n=400;dt=1/(frequency*1e6*100)
    vphasor=3+4j;iphasor=1-.5j
    voltage=[(i*dt,(vphasor*cmath.exp(2j*math.pi*frequency*1e6*i*dt)).real) for i in range(n)]
    current=[((i+.5)*dt,(iphasor*cmath.exp(2j*math.pi*frequency*1e6*(i+.5)*dt)).real) for i in range(n)]
    assert abs(reviewer.impedance(voltage,current,frequency)-vphasor/iphasor)<1e-10
    with tempfile.TemporaryDirectory(prefix='m10-review-test-') as tmp:
        d=Path(tmp)
        for er in (1,3):
            job=d/str(er);(job/'solver').mkdir(parents=True)
            spec={'action':'verify','ground':'real','epsilon_r':er,'mesh_mm':15,
                  'wire_cell_mm':2.5,'wire_radius_mm':5,'steps':120000}
            (job/'spec.json').write_text(json.dumps(spec))
            (job/'worker.log').write_text('')
            for name,data in (('port_ut_1',voltage),('port_it_1',current)):
                (job/'solver'/name).write_text(''.join(str(t)+' '+str(v)+'\n' for t,v in data))
        report=reviewer.review(d)
        assert len(report['pairs'])==3
        assert all(row['status']=='provisional_incomplete_pair' for row in report['pairs'])
        assert not any(run['normal_early_stop'] for run in report['runs'])
        # Identical bulk grids with different source/domain physics must never
        # be cross-paired, even when imported under one results directory.
        import shutil
        for er in (1,3):
            job=d/('variant-'+str(er))
            shutil.copytree(d/str(er),job)
            spec=json.loads((job/'spec.json').read_text())
            spec['domain_margin_mm']=200
            spec['feed_gap_mm']=5
            (job/'spec.json').write_text(json.dumps(spec))
        report=reviewer.review(d)
        assert len(report['pairs'])==6
        assert {(r['domain_margin_mm'],r['feed_gap_mm']) for r in report['pairs']}=={(0,10),(200,5)}
        for er in (1,3):
            job=d/('contact-'+str(er));shutil.copytree(d/str(er),job)
            spec=json.loads((job/'spec.json').read_text());spec['cube_bottom_mm']=0
            (job/'spec.json').write_text(json.dumps(spec))
        report=reviewer.review(d)
        assert len(report['pairs'])==9
        assert {r['cube_bottom_mm'] for r in report['pairs']}=={0,550}
        for witness,loaded_er,sigma in ((280,3.5,.015),(320,2.5,.025)):
            for source_er in (1,3):
                job=d/('witness-'+str(witness)+'-'+str(source_er))
                shutil.copytree(d/str(source_er),job)
                spec=json.loads((job/'spec.json').read_text())
                spec['epsilon_r']=1 if source_er==1 else loaded_er
                spec['witness_frequency_mhz']=witness
                if source_er!=1: spec['dielectric_sigma_s_m']=sigma
                (job/'spec.json').write_text(json.dumps(spec))
        report=reviewer.review(d)
        assert len(report['pairs'])==15
        for witness,er,sigma in ((280,3.5,.015),(320,2.5,.025)):
            rows=[r for r in report['pairs'] if r['witness_frequency_mhz']==witness]
            assert len(rows)==3 and all(r['loaded_epsilon_r']==er and r['dielectric_sigma_s_m']==sigma for r in rows)
        for split in ('x','z'):
            for er in (1,3):
                job=d/('split-'+split+'-'+str(er));shutil.copytree(d/str(er),job)
                spec=json.loads((job/'spec.json').read_text());spec['material_split']=split
                (job/'spec.json').write_text(json.dumps(spec))
        report=reviewer.review(d)
        assert len(report['pairs'])==21
        for split in ('x','z'):
            assert len([r for r in report['pairs'] if r['material_split']==split])==3
    print('Offline port DFT preserves magnetic half-step timestamps and never promotes incomplete pairs')
    if os.name!='posix':
        print('Plan/argument checks passed; POSIX group termination unavailable')
        return
    with tempfile.TemporaryDirectory(prefix='m10-runner-test-') as tmp:
        d=Path(tmp)
        done=runner.run_worker([sys.executable,'-c','print("finished")'],2,d/'done.log')
        assert done['status']=='returned'
        marker=d/'orphan.txt'
        # A child solver outlives its launcher unless the process group dies.
        child='import time; from pathlib import Path; time.sleep(1); Path('+repr(str(marker))+').write_text("orphan")'
        parent='import subprocess,sys,time; subprocess.Popen([sys.executable,"-c",'+repr(child)+']); time.sleep(10)'
        stopped=runner.run_worker([sys.executable,'-c',parent],.2,d/'stopped.log')
        assert stopped['status']=='timed_out',stopped
        assert stopped['wall_seconds']<3,stopped
        time.sleep(1.1)
        assert not marker.exists(),'timeout left a simulated solver running'
    print('OpenEMS plan, invalid budget, completion and worker/solver group timeout checks passed')

if __name__=='__main__': main()
