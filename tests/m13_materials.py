#!/usr/bin/env python3
"""M13-03 material interfaces: analytic static blocks and coated-sphere response."""
import cmath
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import dielectric_ground as g

ROOT=Path(os.environ.get('TEST_SRCDIR','.')).resolve()
BINARY=Path(os.environ.get('NEC2C','./nec2c')).resolve()
EPS0=8.8541878128e-12
spec=importlib.util.spec_from_file_location('static',ROOT/'tools/m11_static_reference.py')
static=importlib.util.module_from_spec(spec);spec.loader.exec_module(static)


def solve(d,materials,boxes,frequency=1,mode='dense',ground='',geometry='',excitation='EX 1 1 1 0 0 0 0 0 0 0',success=True):
    resource=materials+'\nREGION 1 MATERIAL 1\nREGION 2 MATERIAL 2\n'
    resource+='\n'.join('BOX %d '%region+' '.join(format(v,'.17g') for v in bounds) for region,bounds in boxes)
    (d/'model.diel').write_text(resource+'\nMESH MAX_KH .5 MAX_ASPECT 2 MAX_CELLS 10000\nEND\n')
    (d/'model.nec').write_text(f'CE\n{geometry}GE\nFR 0 1 0 0 {frequency} 0\n{ground}\n{excitation}\nDF "model.diel"\nXQ\nEN\n')
    result=subprocess.run([str(BINARY),'-s',mode,'-r','1e-11','-i',str(d/'model.nec'),'-o',str(d/'model.out')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    text=(d/'model.out').read_text();assert (result.returncode==0)==success,text+result.stderr.decode()
    return text


def material_text(e1=3-.2j,e2=2-.4j):
    return f'NEC2C-DIELECTRIC 2\nMATERIAL 1 ER {e1.real} EI {-e1.imag}\nMATERIAL 2 ER {e2.real} EI {-e2.imag}'


def flattened(text,n,coupled=False):
    return [v for row in g.currents(text,n,coupled) for v in row]


def static_pair(d):
    h=.02;eps=[3-.2j,2-.4j];omega=2*math.pi*1e6
    boxes=[(1,(-h,-h/2,-h/2,0,h/2,h/2)),(2,(0,-h/2,-h/2,h,h/2,h/2))]
    self_n=static.diagonal((h,h,h),(0,0,0));other_n=static.diagonal((h,h,h),(h,0,0))
    for component,excitation in [(0,'EX 1 1 1 0 0 0 0 0 0 0'),(2,'EX 1 1 1 0 90 0 0 0 0 0')]:
        t=solve(d,material_text(*eps),boxes,excitation=excitation)
        j=g.currents(t,2)
        a=[e-(e-1)*self_n[component] for e in eps]
        b=[-(e-1)*other_n[component] for e in eps]
        kh=omega/299792458*h
        average=math.sin(kh/2)/(kh/2)
        rhs=[1j*omega*EPS0*(e-1)*(1 if component==0 else -1)*average*
             cmath.exp(1j*kh*(i-.5) if component==2 else 0) for i,e in enumerate(eps)]
        det=a[0]*a[1]-b[0]*b[1]
        expected=[(rhs[0]*a[1]-b[0]*rhs[1])/det,(a[0]*rhs[1]-b[1]*rhs[0])/det]
        g.close([row[component] for row in j],expected,rtol=4e-6)
        joule=.5*omega*EPS0*h**3*sum(-e.imag*sum(abs(v/(1j*omega*EPS0*(e-1)))**2 for v in row) for e,row in zip(eps,j))
        assert abs(g.scalar(t,'DIELECTRIC ABSORBED POWER (W):')/joule-1)<1e-9
    print('Independent closed-form static face-contact blocks, normal/tangential polarization and local Joule loss pass',flush=True)


def topology_and_solvers(d):
    h=.02;base=[(1,(0,0,0,h,h,h)),(2,(h,0,0,2*h,h,h))]
    # Face, edge, corner, and a nonconforming partial face share no interior volume.
    for second in [(h,0,0,2*h,h,h),(h,h,0,2*h,2*h,h),(h,h,h,2*h,2*h,2*h),(h,0,0,2*h,h/2,h)]:
        t=solve(d,material_text(),[base[0],(2,second)],frequency=100)
        assert len(g.currents(t,2))==2
    overlapping=solve(d,material_text(),[base[0],(2,(h/2,0,0,1.5*h,h,h))],success=False)
    assert 'physical box interiors overlap' in overlapping
    # Even a tiny overlap must not be snapped into a contact.
    tiny=solve(d,material_text(),[base[0],(2,(h-1e-12,0,0,2*h,h,h))],success=False)
    assert 'physical box interiors overlap' in tiny
    for ground in ('','GN 1'):
        dense=solve(d,material_text(),base,frequency=100,ground=ground)
        for mode in ('gmres','fft'):
            t=solve(d,material_text(),base,frequency=100,ground=ground,mode=mode)
            g.close(flattened(t,2),flattened(dense,2),rtol=2e-8)
            assert abs(g.scalar(t,'DIELECTRIC ABSORBED POWER (W):')/g.scalar(dense,'DIELECTRIC ABSORBED POWER (W):')-1)<1e-8
        assert g.scalar(dense,'RELATIVE IMBALANCE:')<.002
    # Material ID labels and input order carry no physics.
    a=solve(d,material_text(),base,frequency=100)
    b=solve(d,material_text(2-.4j,3-.2j),[(2,base[0][1]),(1,base[1][1])],frequency=100)
    g.close(flattened(a,2),[v for row in reversed(g.currents(b,2)) for v in row],rtol=1e-10)
    # Equal material split must reproduce the same two physical cells.
    a=solve(d,material_text(3-.2j,3-.2j),base,frequency=100)
    b=solve(d,material_text(3-.2j,3-.2j),[(1,x) for _,x in base],frequency=100)
    g.close(flattened(a,2),flattened(b,2),rtol=1e-10)
    # Conductive ER=1 and table-valued materials at a shared interface.
    tables=('NEC2C-DIELECTRIC 3\nMATERIAL 1 TABLE\nSAMPLE 1 FREQ_MHZ 100 ER 1 SIGMA .01\nSAMPLE 1 FREQ_MHZ 300 ER 1 SIGMA .03\n'
            'MATERIAL 2 TABLE\nSAMPLE 2 FREQ_MHZ 100 ER 2 EI .2\nSAMPLE 2 FREQ_MHZ 300 ER 4 EI .4')
    a=solve(d,tables,base,frequency=200,mode='fft')
    e1=1-.02/(2*math.pi*200e6*EPS0)*1j
    b=solve(d,material_text(e1,3-.3j),base,frequency=200)
    g.close(flattened(a,2),flattened(b,2),rtol=2e-8)
    print('Face/edge/corner/partial-face topology, strict overlap refusal, equal-material limit, table identity and dense/GMRES/FFT free/PEC agreement pass',flush=True)
    # Actual wire coupling includes row susceptibility, radiation and absorption.
    boxes=[(1,(.1,-.01,.58,.12,.01,.60)),(2,(.12,-.01,.58,.14,.01,.60))]
    from dielectric_wire import coupled_impedance
    for ground in ('','GN 1'):
        t=solve(d,material_text(),boxes,frequency=300,ground=ground,geometry=g.WIRE,
                excitation='EX 0 1 6 0 1 0')
        for mode in ('gmres','fft'):
            other=solve(d,material_text(),boxes,frequency=300,ground=ground,geometry=g.WIRE,
                        excitation='EX 0 1 6 0 1 0',mode=mode)
            g.close(flattened(t,2,True),flattened(other,2,True),rtol=2e-8)
            assert abs(coupled_impedance(t)-coupled_impedance(other))<1e-7
        assert g.scalar(t,'DIELECTRIC ABSORBED POWER (W):')>0
        assert g.scalar(t,'RELATIVE POWER DEFECT:')<.01
    print('Coupled wire/free-space/PEC touching-material currents and global power pass',flush=True)
    from m10_pec_validation import closed_box
    # Match M13-01 mixed-conductor resolution: 11 wire segments exceed the
    # 0.3% power gate even in a homogeneous-material control.
    boxes=[(1,(.1,-.01,0,.12,.01,.02)),(2,(.12,-.01,0,.14,.01,.02))]
    for ground in ('','GN 1'):
        for geometry,excitation in [(closed_box(),'EX 1 1 1 0 37 23 19 0 0 .6'),
                                    (g.WIRE.replace('GW 1 11','GW 1 21')+closed_box(),'EX 0 1 11 0 1 0')]:
            t=solve(d,material_text(),boxes,frequency=300,ground=ground,
                    geometry=geometry,excitation=excitation)
            f=solve(d,material_text(),boxes,frequency=300,ground=ground,
                    geometry=geometry,excitation=excitation,mode='fft')
            g.close(flattened(t,2,True),flattened(f,2,True),rtol=1e-7)
            assert g.scalar(t,'DIELECTRIC ABSORBED POWER (W):')>0
            assert g.scalar(t,'RELATIVE POWER DEFECT:')<.003,(ground,len(geometry),g.scalar(t,'RELATIVE POWER DEFECT:'))
    print('Heterogeneous patch/mixed-conductor free-space/PEC-contact power and dense/FFT checks pass',flush=True)


def sphere_boxes(n,outer=.03,inner=.018):
    h=2*outer/n;boxes=[]
    for x in range(n):
        for y in range(n):
            for z in range(n):
                indices=(x,y,z);r2=sum(((-n/2+i+.5)*h)**2 for i in indices)
                if r2>=outer**2:continue
                bounds=[(-n/2+i)*h for i in indices]+[(-n/2+i+1)*h for i in indices]
                boxes.append((1 if r2<inner**2 else 2,tuple(bounds)))
    return boxes,h


def coated_alpha(core,shell,outer=.03,inner=.018):
    f=(inner/outer)**3
    return 4*math.pi*outer**3*((shell-1)*(core+2*shell)+f*(core-shell)*(1+2*shell))/((shell+2)*(core+2*shell)+2*f*(shell-1)*(core-shell))


def sphere_checks(d):
    eps=(4-.4j,2-.2j);alpha=coated_alpha(*eps);omega=2*math.pi*1e6;rows=[]
    for n in (8,12,16):
        boxes,h=sphere_boxes(n)
        text=solve(d,material_text(*eps),boxes,mode='fft')
        j=g.currents(text,len(boxes));measured=sum(v[0] for v in j)*h**3/(1j*omega*EPS0)
        error=abs(measured/alpha-1)
        analytic_absorption=-.5*omega*EPS0*alpha.imag
        absorption_error=abs(g.scalar(text,'DIELECTRIC ABSORBED POWER (W):')/analytic_absorption-1)
        rows.append(dict(diameter_cells=n,cells=len(boxes),polarizability_relative_error=error,absorption_relative_error=absorption_error))
        print('Coated sphere '+json.dumps(rows[-1]),flush=True)
        assert error<.08 and absorption_error<.08,rows
    assert rows[-1]['polarizability_relative_error']<.03,rows
    assert rows[-1]['absorption_relative_error']<.03,rows
    assert rows[-1]['polarizability_relative_error']<rows[0]['polarizability_relative_error'],rows
    return rows


def main():
    with tempfile.TemporaryDirectory(prefix='m13-03-') as tmp:
        d=Path(tmp);static_pair(d);topology_and_solvers(d);sphere_checks(d)
    print('M13-03 native multi-material acceptance tests pass',flush=True)

if __name__=='__main__':main()
