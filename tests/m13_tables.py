#!/usr/bin/env python3
"""M13-02 interpolation, sweep preflight and conservative mesh tests."""
from pathlib import Path
import math
import os
import subprocess
import tempfile
import dielectric_ground as g

BINARY=Path(os.environ.get('NEC2C','./nec2c')).resolve()
EPS0=8.8541878128e-12
TABLE=('NEC2C-DIELECTRIC 3\nMATERIAL 1 TABLE\n'
       'SAMPLE 1 FREQ_MHZ 100 ER 4 SIGMA .01\n'
       'SAMPLE 1 FREQ_MHZ 200 ER 3 SIGMA .02\n'
       'SAMPLE 1 FREQ_MHZ 300 ER 2 SIGMA .03\n')


def constant(er,loss,kind='SIGMA'):
    return f'NEC2C-DIELECTRIC 2\nMATERIAL 1 ER {er:.17g} {kind} {loss:.17g}\n'


def run(d,material=TABLE,frequency=200,count=1,step=0,mode=0,solver='dense',
        box='-.01 -.01 -.01 .01 .01 .01',ground='',geometry='',success=True):
    (d/'model.diel').write_text(material+'REGION 1 MATERIAL 1\nBOX 1 '+box+
        '\nMESH MAX_KH .5 MAX_ASPECT 2 MAX_CELLS 10000\nEND\n')
    ex='EX 0 1 11 0 1 0' if geometry else 'EX 1 1 1 0 37 23 19 0 0 .6'
    (d/'model.nec').write_text(f'CE\n{geometry}GE\nFR {mode} {count} 0 0 {frequency:.17g} {step:.17g}\n{ground}\n{ex}\nDF "model.diel"\nXQ\nEN\n')
    r=subprocess.run([str(BINARY),'-s',solver,'-r','1e-12','-i',str(d/'model.nec'),'-o',str(d/'model.out')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    text=(d/'model.out').read_text()
    assert (r.returncode==0)==success,text+r.stderr.decode()
    return text


def test_nx(d):
    first=run(d,frequency=200)
    (d/'first.diel').write_text((d/'model.diel').read_text())
    deck=(d/'model.nec').read_text().replace('"model.diel"','"first.diel"').replace('EN\n','NX\n')
    second=run(d,constant(2,.01),frequency=200)
    deck+=(d/'model.nec').read_text()
    (d/'joined.nec').write_text(deck)
    result=subprocess.run([str(BINARY),'-i',str(d/'joined.nec'),'-o',str(d/'joined.out')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    text=(d/'joined.out').read_text();assert result.returncode==0,text+result.stderr.decode()
    parts=text.split('DIELECTRIC VOLUME INTEGRAL SOLUTION')[1:]
    assert len(parts)==2
    for part,reference in zip(parts,(first,second)):
        g.close(g.currents(part,1)[0],g.currents(reference,1)[0],rtol=1e-10)
    print('NX releases table data and restores independent material state',flush=True)


def test_fr_changes(d):
    box='-.02 -.02 -.02 .02 .02 .02'
    first=run(d,frequency=100,box=box)
    deck=(d/'model.nec').read_text()
    second=run(d,frequency=300,box=box)
    for f,success in ((300,True),(399,False)):
        joined=deck.replace('EN\n',f'FR 0 1 0 0 {f} 0\nXQ\nEN\n')
        (d/'joined.nec').write_text(joined)
        result=subprocess.run([str(BINARY),'-i',str(d/'joined.nec'),'-o',str(d/'joined.out')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        text=(d/'joined.out').read_text()
        assert (result.returncode==0)==success,text+result.stderr.decode()
        parts=text.split('DIELECTRIC VOLUME INTEGRAL SOLUTION')[1:]
        if success:
            assert len(parts)==2
            for part,ref,n in zip(parts,(first,second),(1,8)):
                g.close([v for row in g.currents(part,n) for v in row],
                        [v for row in g.currents(ref,n) for v in row],rtol=1e-9)
        else: assert len(parts)==1 and 'extrapolation is prohibited' in text
    print('Later FR cards rebuild the mesh and revalidate table coverage',flush=True)


def test_multiple_tables(d):
    geometry='REGION 2 MATERIAL 2\nBOX 2 .05 -.01 -.01 .07 .01 .01\n'
    data=('NEC2C-DIELECTRIC 3\nMATERIAL 2 TABLE\nMATERIAL 1 TABLE\n'
          'SAMPLE 2 FREQ_MHZ 300 ER 4 EI .3\n'
          'SAMPLE 1 FREQ_MHZ 100 ER 4 SIGMA .01\n'
          'SAMPLE 2 FREQ_MHZ 100 ER 2 EI .1\n'
          'SAMPLE 1 FREQ_MHZ 300 ER 2 SIGMA .03\n'+geometry)
    reference=constant(3,.02)+'MATERIAL 2 ER 3 EI .2\n'+geometry
    a=run(d,data);b=run(d,reference)
    for x,y in zip(g.currents(a,2),g.currents(b,2)): g.close(x,y,rtol=1e-9)
    c=run(d,data,solver='gmres')
    for x,y in zip(g.currents(c,2),g.currents(a,2)): g.close(x,y,rtol=1e-9)
    c=run(d,data,solver='fft')
    for x,y in zip(g.currents(c,2),g.currents(a,2)): g.close(x,y,rtol=1e-8)
    print('Independent EI/SIGMA tables retain material identity across all solver paths',flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix='m13-02-') as tmp:
        d=Path(tmp)
        for f in (100,150,200,250,300):
            er=5-f/100; sigma=f/10000
            ref=run(d,constant(er,sigma),frequency=f)
            for solver in ('dense','gmres','fft'):
                table=run(d,frequency=f,solver=solver)
                g.close(g.currents(table,1)[0],g.currents(ref,1)[0],rtol=1e-9)
                assert abs(g.scalar(table,'DIELECTRIC ABSORBED POWER (W):')/g.scalar(ref,'DIELECTRIC ABSORBED POWER (W):')-1)<1e-9
                assert g.scalar(table,'RELATIVE IMBALANCE:')<.002
        print('SIGMA table knots/interpolation and dense/GMRES/FFT absorption match constant-material references',flush=True)
        ei=TABLE.replace('SIGMA .01','EI .2').replace('SIGMA .02','EI .4').replace('SIGMA .03','EI .6')
        a=run(d,ei,frequency=150); b=run(d,constant(3.5,.3,'EI'),frequency=150)
        g.close(g.currents(a,1)[0],g.currents(b,1)[0],rtol=1e-10)
        # Independent static cube polarizability checks actual interpolated complex epsilon.
        small=ei.replace('100 ER','1 ER').replace('200 ER','2 ER').replace('300 ER','3 ER')
        text=run(d,small,frequency=1.5)
        j=g.currents(text,1)[0]
        # Oblique incidence has unit real incident electric vector; compare vector norm.
        expected=2*math.pi*1.5e6*EPS0*abs((3.5-.3j-1)*3/(3.5-.3j+2))
        assert abs(math.sqrt(sum(abs(v)**2 for v in j))/expected-1)<5e-6
        print('EI interpolation and independent small-cube polarizability pass',flush=True)
        kink=TABLE.replace('200 ER 3 SIGMA .02','200 ER 6 SIGMA .03').replace('300 ER 2 SIGMA .03','300 ER 2 SIGMA .02')
        for f,er,sigma in ((150,5,.02),(250,4,.025)):
            a=run(d,kink,frequency=f,solver='fft');b=run(d,constant(er,sigma),frequency=f)
            g.close(g.currents(a,1)[0],g.currents(b,1)[0],rtol=1e-9)
        print('Non-collinear table rows interpolate within each adjacent interval',flush=True)
        conducting=TABLE.replace('ER 4','ER 1').replace('ER 3','ER 1').replace('ER 2','ER 1')
        a=run(d,conducting);b=run(d,constant(1,.02))
        g.close(g.currents(a,1)[0],g.currents(b,1)[0],rtol=1e-9)
        assert g.scalar(a,'DIELECTRIC ABSORBED POWER (W):')>0
        for mode,first,step,count in ((0,100,100,3),(0,300,-100,3),(0,130,40,5),(1,100,math.sqrt(3),3)):
            text=run(d,frequency=first,step=step,count=count,mode=mode,solver='fft')
            parts=text.split('DIELECTRIC VOLUME INTEGRAL SOLUTION')[1:]
            assert len(parts)==count
            f=first
            for part in parts:
                ref=run(d,frequency=f,solver='fft')
                g.close(g.currents(part,1)[0],g.currents(ref,1)[0],rtol=1e-9)
                f=f*step if mode else f+step
        print('Ascending, descending, off-knot and geometric sweeps match separate runs',flush=True)
        for first,step,count in ((99,1,3),(100,101,3),(300,-101,3)):
            t=run(d,frequency=first,step=step,count=count,success=False)
            assert 'extrapolation is prohibited' in t and 'INTEGRAL SOLUTION' not in t
        print('Whole-sweep preflight rejects extrapolation before producing any solution',flush=True)
        # Interior material peak must determine the fixed sweep mesh.
        peak=ei.replace('200 ER 3 EI .4','200 ER 40 EI 0').replace('EI .2','EI 0').replace('EI .6','EI 0')
        kw=dict(frequency=300,box='-.01 -.01 -.01 .01 .01 .01')
        table=run(d,peak,**kw); worst=run(d,constant(40,0,'EI'),**kw); endpoint=run(d,constant(2,0,'EI'),**kw)
        cells=lambda t:int(g.scalar(t,'ACTIVE CELLS:'))
        assert cells(table)==cells(worst)>cells(endpoint)
        assert g.scalar(table,'DIELECTRIC ABSORBED POWER (W):')==0
        print('Interior high-permittivity knot controls mesh; zero-loss table has zero absorption',flush=True)
        # Sort rows and resolve references independently of file order.
        lines=TABLE.splitlines(); shuffled='\n'.join([lines[0],lines[4],lines[3],lines[2],lines[1]])+'\n'
        g.close(g.currents(run(d,shuffled),1)[0],g.currents(run(d),1)[0],rtol=1e-12)
        bad=[TABLE+lines[2]+'\n',TABLE.replace('200 ER 3 SIGMA .02','200 ER 3 EI .02'),
             TABLE.replace('MATERIAL 1 TABLE','MATERIAL 1 ER 3'),TABLE.replace('MATERIAL 1 TABLE','MATERIAL 2 TABLE'),
             '\n'.join(lines[:3])+'\n',TABLE.replace('FREQ_MHZ 100','FREQ_MHZ 0'),
             TABLE.replace('ER 4 SIGMA .01','ER 1 SIGMA 0'),TABLE.replace('SIGMA .01','SIGMA -.01'),
             TABLE.replace('FREQ_MHZ 100','FREQ_MHZ nan'),TABLE.replace('DIELECTRIC 3','DIELECTRIC 2')]
        for value in bad: run(d,value,success=False)
        print('Canonical row ordering, forward references and invalid tables pass',flush=True)
        for ground in ('','GN 1'):
            kw=dict(ground=ground,geometry=g.WIRE.replace('GW 1 11','GW 1 21'),box='.1 -.01 0 .12 .01 .02',frequency=300)
            a=run(d,solver='fft',**kw);b=run(d,constant(2,.03),**kw)
            g.close(g.currents(a,1,True)[0],g.currents(b,1,True)[0],rtol=1e-8)
            assert g.scalar(a,'RELATIVE POWER DEFECT:')<.003
        print('Coupled wire free-space/PEC table currents and power pass',flush=True)
        test_multiple_tables(d)
        test_nx(d)
        test_fr_changes(d)


if __name__=='__main__': main()
