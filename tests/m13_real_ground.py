#!/usr/bin/env python3
"""M13-04: passive volume loss plus GN 2 substrate; independent static images."""
import cmath
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import tempfile
import dielectric_ground as g
from dielectric_wire import coupled_impedance,coupled_current

ROOT=Path(os.environ.get('TEST_SRCDIR','.')).resolve()
BINARY=Path(os.environ.get('NEC2C','./nec2c')).resolve()
EPS0=8.8541878128e-12
GN='GN 2 0 0 0 13 .02'
WIRE=g.WIRE.replace('GW 1 11','GW 1 21')
SPEC=importlib.util.spec_from_file_location('static_reference',ROOT/'tools/m11_static_reference.py')
static=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(static)


def solve(d,material='NEC2C-DIELECTRIC 2\nMATERIAL 1 ER 3 SIGMA .02',
          boxes='REGION 1 MATERIAL 1\nBOX 1 .10 -.01 .58 .12 .01 .60',
          geometry='',ground=GN,frequency=300,count=1,step=0,solver='dense',
          excitation='EX 1 1 1 0 37 23 19 0 0 .6',success=True,extra=''):
    (d/'model.diel').write_text(material+'\n'+boxes+'\nMESH MAX_KH .5 MAX_ASPECT 2 MAX_CELLS 100\nEND\n')
    (d/'model.nec').write_text(f'CE\n{geometry}GE\nFR 0 {count} 0 0 {frequency} {step}\n{ground}\n{excitation}\nDF "model.diel"\n{extra}\nXQ\nEN\n')
    r=subprocess.run([str(BINARY),'-s',solver,'-r','1e-11','-i',str(d/'model.nec'),'-o',str(d/'model.out')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    text=(d/'model.out').read_text();assert (r.returncode==0)==success,text+r.stderr.decode()
    return text


def flat(text,n=1,coupled=False):return [v for row in g.currents(text,n,coupled) for v in row]


def check_power(text,coupled=False):
    prefix='RELATIVE POWER DEFECT:' if coupled else 'RELATIVE IMBALANCE:'
    defect=g.scalar(text,prefix);assert defect<.003,defect
    absorbed=g.scalar(text,'DIELECTRIC ABSORBED POWER (W):');assert absorbed>=0
    lower=g.scalar(text,'LOWER-HALF-SPACE POWER:' if 'LOWER-HALF-SPACE POWER:' in text else 'LOWER-HALF-SPACE SCATTERED POWER:')
    assert lower>0
    return absorbed,lower,defect


def test_static_image(d):
    h=.01;height=.01;frequency=1;er=3-.4j;soil=4-.4j
    text=solve(d,material='NEC2C-DIELECTRIC 2\nMATERIAL 1 ER 3 EI .4',
               boxes='REGION 1 MATERIAL 1\nBOX 1 -.005 -.005 .005 .005 .005 .015',
               ground='GN 2 0 0 0 4 -.4',frequency=frequency,
               excitation='EX 1 1 1 0 0 0 0 0 0 0')
    omega=2*math.pi*frequency*1e6;k=omega/299792458
    rho=(soil-1)/(soil+1)
    image=static.diagonal((h,h,h),(0,0,2*height))[0]*(-rho)
    fresnel=(1-cmath.sqrt(soil))/(1+cmath.sqrt(soil))
    background=(cmath.exp(1j*k*height)+fresnel*cmath.exp(-1j*k*height))*math.sin(k*h/2)/(k*h/2)
    expected=1j*omega*EPS0*(er-1)*background/(1+(er-1)/3-(er-1)*image)
    actual=g.currents(text,1)[0][0]
    assert abs(actual/expected-1)<5e-6,(actual,expected,abs(actual/expected-1))
    internal=actual/(1j*omega*EPS0*(er-1))
    joule=.5*omega*EPS0*.4*h**3*abs(internal)**2
    assert abs(g.scalar(text,'DIELECTRIC ABSORBED POWER (W):')/joule-1)<1e-8
    check_power(text)
    print('Independent quasistatic finite-cube image/Fresnel/Joule check: relative current error',abs(actual/expected-1),flush=True)


def test_material_and_ground_limits(d):
    # Matched ground removes reflection while retaining the transmitted-power channel.
    for geometry,ex in [('', 'EX 1 1 1 0 37 23 19 0 0 .6'),(WIRE,'EX 0 1 11 0 1 0')]:
        a=solve(d,geometry=geometry,excitation=ex,ground='GN 2 0 0 0 1 0')
        b=solve(d,geometry=geometry,excitation=ex,ground='')
        g.close(flat(a,coupled=bool(geometry)),flat(b,coupled=bool(geometry)),rtol=1e-9)
        assert abs(g.scalar(a,'DIELECTRIC ABSORBED POWER (W):')/g.scalar(b,'DIELECTRIC ABSORBED POWER (W):')-1)<1e-9
        check_power(a,bool(geometry))
    # Lossless finite-permittivity ground transmits; it is not PEC or free space.
    t=solve(d,ground='GN 2 0 0 0 4 0');check_power(t)
    t=solve(d,material='NEC2C-DIELECTRIC 2\nMATERIAL 1 ER 1 SIGMA .02')
    absorbed,_,_=check_power(t);assert absorbed>0 and len(g.currents(t,1))==1
    # Constant EI and sigma match only when converted at the actual frequency.
    f=300;ei=.02/(2*math.pi*f*1e6*EPS0)
    a=solve(d);b=solve(d,material=f'NEC2C-DIELECTRIC 2\nMATERIAL 1 ER 3 EI {ei:.17g}')
    g.close(flat(a),flat(b),rtol=1e-9);check_power(a)
    a=solve(d,material='NEC2C-DIELECTRIC 2\nMATERIAL 1 ER 3 SIGMA 0')
    b=solve(d,material='NEC2C-DIELECTRIC 1\nMATERIAL 1 ER 3')
    g.close(flat(a),flat(b),rtol=1e-11)
    assert g.scalar(a,'DIELECTRIC ABSORBED POWER (W):')==0
    print('Lossy-volume matched-ground, lossless transmission, EI/sigma and zero-loss legacy limits pass',flush=True)


def test_tables_touching_and_power(d):
    boxes=('REGION 1 MATERIAL 1\nREGION 2 MATERIAL 2\n'
           'BOX 1 .10 -.01 .58 .12 .01 .60\nBOX 2 .12 -.01 .58 .14 .01 .60')
    table=('NEC2C-DIELECTRIC 3\nMATERIAL 1 TABLE\n'
           'SAMPLE 1 FREQ_MHZ 280 ER 2 SIGMA .01\nSAMPLE 1 FREQ_MHZ 320 ER 4 SIGMA .03\n'
           'MATERIAL 2 ER 2 EI .3')
    constant='NEC2C-DIELECTRIC 2\nMATERIAL 1 ER 3 SIGMA .02\nMATERIAL 2 ER 2 EI .3'
    a=solve(d,material=table,boxes=boxes)
    b=solve(d,material=constant,boxes=boxes,solver='gmres')
    g.close(flat(a,2),flat(b,2),rtol=1e-8);check_power(a)
    for ex in ('EX 0 1 11 0 1 0','EX 1 1 1 0 37 23 19 0 0 .6'):
        a=solve(d,material=constant,boxes=boxes,geometry=WIRE,excitation=ex)
        b=solve(d,material=constant,boxes=boxes,geometry=WIRE,excitation=ex,solver='gmres')
        g.close(flat(a,2,True),flat(b,2,True),rtol=1e-8);check_power(a,True)
    # Conductive and fixed-imaginary ground laws must reinitialize on every frequency.
    for gn in (GN,'GN 2 0 0 0 13 -.4'):
        sweep=solve(d,material=table,boxes=boxes,ground=gn,frequency=280,count=2,step=40)
        parts=sweep.split('DIELECTRIC VOLUME INTEGRAL SOLUTION')[1:];assert len(parts)==2
        for f,part in zip((280,320),parts):
            one=solve(d,material=table,boxes=boxes,ground=gn,frequency=f)
            g.close(flat(part,2),flat(one,2),rtol=1e-9)
    print('Touching/table-valued lossy materials over GN 2, solver agreement, driven/plane-wave power and sweep state pass',flush=True)


def test_reciprocity_and_policy(d):
    wires=WIRE+'GW 2 21 .5 .1 .42 .5 .1 .86 .001\n'
    a=solve(d,geometry=wires,excitation='EX 0 1 11 0 1 0')
    b=solve(d,geometry=wires,excitation='EX 0 2 11 0 1 0')
    x=coupled_current(a,32);y=coupled_current(b,11)
    assert abs(x-y)/max(abs(x),abs(y))<2e-4,(x,y)
    for kwargs,message in [({'ground':'GN 0 0 0 0 13 .02'},'GN 0 dielectric coupling'),
                           ({'solver':'fft'},'real-ground FFT is not supported'),
                           ({'boxes':'REGION 1 MATERIAL 1\nBOX 1 .1 -.01 0 .12 .01 .02'},'strictly positive clearance'),
                           ({'geometry':'SM 2 2 .8 -.2 .4 .8 .2 .4\nSC 0 0 .8 .2 .8 .8 -.2 .8\n'},'supports wires only')]:
        text=solve(d,success=False,**kwargs);assert message in text,text
    print('Lossy coupled-port reciprocity and GN 0/patch/contact/FFT scope refusals pass',flush=True)


def test_review_isolation(d):
    import json
    import sys
    sys.path.insert(0,str(ROOT/'tools'))
    from m13_04_review import review_controls
    from m10_review_results import review
    root=d/'review-controls'
    for label,sigma in (('lossless',0),('lossy',.02)):
        for mesh,spacing,margin in ((20,5,0),(15,2.5,0),(15,2.5,150)):
            for er in (1,3):
                job=root/label/('m%d-d%d-er%d'%(mesh,margin,er));(job/'solver').mkdir(parents=True)
                spec=dict(action='verify',ground='real',epsilon_r=er,dielectric_sigma_s_m=sigma if er==3 else 0,
                          mesh_mm=mesh,wire_cell_mm=spacing,wire_radius_mm=5,steps=120000,domain_margin_mm=margin,
                          bandwidth_mhz=240,end_criterion=1e-6,mesh_revision=2,feed_model='finite_box')
                (job/'spec.json').write_text(json.dumps(spec));(job/'worker.log').write_text('')
                dt=1/(300e6*100)
                for name,amplitude,offset in [('port_ut_1',2,0),('port_it_1',1,.5)]:
                    (job/'solver'/name).write_text(''.join('%g %.17g\n'%((i+offset)*dt,amplitude*math.cos(2*math.pi*300e6*(i+offset)*dt)) for i in range(400)))
    result=review_controls(root)
    assert len(result['runs'])==12 and len(result['pairs'])==18
    assert all(p['status']=='provisional_incomplete_pair' for p in result['pairs'])
    for label,sigma in (('lossless',0),('lossy',.02)):
        rows=[p for p in result['pairs'] if p['experiment_case']==label]
        assert len(rows)==9 and all(p['dielectric_sigma_s_m']==sigma for p in rows)
    try:review(root)
    except ValueError as e:assert 'duplicate matching physics' in str(e)
    else:raise AssertionError('generic review should refuse cross-control duplicate vacuum jobs')
    print('Offline control isolation preserves loss labels and never promotes incomplete OpenEMS traces',flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix='m13-04-') as tmp:
        d=Path(tmp);test_static_image(d);test_material_and_ground_limits(d)
        test_tables_touching_and_power(d);test_reciprocity_and_policy(d);test_review_isolation(d)
    print('M13-04 native GN 2 combined-loss checks pass',flush=True)

if __name__=='__main__':main()
