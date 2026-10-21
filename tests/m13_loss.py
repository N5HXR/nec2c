#!/usr/bin/env python3
"""M13-01 scalar passive loss: independent polarizability, Joule power and limits."""
import math
import os
import re
from pathlib import Path
import subprocess
import tempfile
import dielectric_ground as g
from dielectric_wire import coupled_impedance
from m10_pec_validation import closed_box

BINARY=Path(os.environ.get('NEC2C','./nec2c')).resolve()
EPS0=8.8541878128e-12
C0=299792458.0


def solve(d, loss='', er=3, ground='', geometry='', mode='dense', frequency=100,
          box='-.01 -.01 -.01 .01 .01 .01', version=2, success=True, sweep=1, wave=False):
    (d/'model.diel').write_text(f'NEC2C-DIELECTRIC {version}\nMATERIAL 1 ER {er} {loss}\n'
        f'REGION 1 MATERIAL 1\nBOX 1 {box}\nMESH MAX_KH .5 MAX_CELLS 100\nEND\n')
    wire=re.search(r'GW\s+1\s+(\d+)',geometry)
    source=(int(wire[1])+1)//2 if wire else 6
    excitation=f'EX 0 1 {source} 0 1 0' if geometry and not wave else 'EX 1 1 1 0 90 0 90 0 0 0'
    if wave: excitation='EX 1 1 1 0 37 23 19 0 0 .6'
    (d/'model.nec').write_text(f'CE\n{geometry}GE\nFR 0 {sweep} 0 0 {frequency} 100\n{ground}\n'
        f'{excitation}\nDF "model.diel"\nXQ\nEN\n')
    r=subprocess.run([str(BINARY),'-s',mode,'-r','1e-12','-i',str(d/'model.nec'),'-o',str(d/'model.out')],
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    text=(d/'model.out').read_text()
    assert (r.returncode==0)==success,text+r.stderr.decode()
    return text


def power(text):
    return g.scalar(text,'DIELECTRIC ABSORBED POWER (W):')


def main():
    with tempfile.TemporaryDirectory(prefix='m13-01-') as tmp:
        d=Path(tmp)
        # Independent electrostatic cube polarizability: depolarization 1/3 by symmetry.
        # Electrically tiny at 1 MHz, so dynamic corrections are negligible.
        t=solve(d,'EI .4',frequency=1)
        j=g.currents(t,1)[0]
        e=3-.4j; omega=2*math.pi*1e6
        expected=1j*omega*EPS0*(e-1)*3/(e+2)
        assert abs(j[1]/expected-1)<3e-6,(j,expected)
        internal=[v/(1j*omega*EPS0*(e-1)) for v in j]
        joule=.5*omega*EPS0*.4*.02**3*sum(abs(v)**2 for v in internal)
        assert abs(power(t)/joule-1)<1e-10
        analytic_absorption=.5*omega*EPS0*.4*.02**3*abs(3/(e+2))**2
        assert abs(power(t)/analytic_absorption-1)<3e-6
        assert power(t)>0
        assert g.scalar(t,'RELATIVE IMBALANCE:')<1e-5
        print('Analytic small-cube polarization, passive sign and independent Joule absorption pass',flush=True)
        # Conductivity and fixed imaginary epsilon agree only at the matching frequency.
        sigma=.02
        for f in (100,200):
            a=solve(d,f'SIGMA {sigma}',frequency=f)
            b=solve(d,f'EI {sigma/(2*math.pi*f*1e6*EPS0):.17g}',frequency=f)
            g.close(g.currents(a,1)[0],g.currents(b,1)[0],rtol=1e-10)
            assert abs(power(a)/power(b)-1)<1e-10
        sweep=solve(d,'SIGMA .02',sweep=2)
        parts=sweep.split('DIELECTRIC VOLUME INTEGRAL SOLUTION')[1:]
        assert len(parts)==2
        for f,part in zip((100,200),parts):
            single=solve(d,'SIGMA .02',frequency=f)
            g.close(g.currents(part,1)[0],g.currents(single,1)[0],rtol=1e-10)
        print('EI/conductivity equivalence and frequency-sweep state pass',flush=True)
        # Test all solver paths and free-space / exact PEC-contact configurations.
        for ground,box in [('', '-.02 -.01 -.01 0 .01 .01\nBOX 1 0 -.01 -.01 .02 .01 .01'),
                           ('GN 1','-.02 -.01 0 0 .01 .02\nBOX 1 0 -.01 0 .02 .01 .02')]:
            dense=solve(d,'EI .4',ground=ground,box=box)
            for mode in ('gmres','fft'):
                other=solve(d,'EI .4',ground=ground,box=box,mode=mode)
                g.close(g.currents(other,2)[0],g.currents(dense,2)[0],rtol=1e-8)
                assert abs(power(other)/power(dense)-1)<1e-8
            assert g.scalar(dense,'RELATIVE IMBALANCE:')<.002
        print('Two-cell lossy dense/GMRES/FFT and PEC-contact power pass',flush=True)
        # Zero loss, old-format compatibility, and smooth approach to zero loss.
        ref=solve(d,version=1)
        for loss in ('','EI 0','SIGMA 0'):
            t=solve(d,loss)
            g.close(g.currents(t,1)[0],g.currents(ref,1)[0],rtol=1e-12)
            assert power(t)==0
        errors=[]
        for ei in (.1,.01,.001):
            t=solve(d,f'EI {ei}')
            errors.append(abs(g.currents(t,1)[0][1]-g.currents(ref,1)[0][1]))
        assert all(b<.11*a for a,b in zip(errors,errors[1:])),errors
        # ER=1 with positive loss is a material, not vacuum.
        t=solve(d,'SIGMA .02',er=1)
        assert len(g.currents(t,1))==1 and power(t)>0
        for loss in ('EI -.1','SIGMA -.1','EI nan','SIGMA inf','TENSOR 1','EI .1 SIGMA .02'):
            solve(d,loss,success=False)
        solve(d,'SIGMA .02',version=1,success=False)
        print('Zero-loss compatibility, continuity, conductive ER=1 and invalid/passive input checks pass',flush=True)
        # Coupled currents and global energy: wire, closed PEC patch, mixed geometry.
        wire=g.WIRE.replace('GW 1 11','GW 1 21')
        for ground in ('','GN 1'):
            for geom,wave in ((wire,False),(closed_box(),True),(wire+closed_box(),False)):
                kw=dict(ground=ground,geometry=geom,wave=wave,box='.1 -.02 0 .14 .02 .04',frequency=300)
                t=solve(d,'SIGMA .02',**kw)
                f=solve(d,'SIGMA .02',mode='fft',**kw)
                g.close(g.currents(f,1,True)[0],g.currents(t,1,True)[0],rtol=1e-7)
                if not wave: assert abs(coupled_impedance(t)-coupled_impedance(f))<1e-6
                assert power(t)>0
                assert g.scalar(t,'RELATIVE POWER DEFECT:')<.003,t
        print('Lossy coupled wire/patch/mixed free-space and PEC-contact energy/FFT checks pass',flush=True)


if __name__=='__main__': main()
