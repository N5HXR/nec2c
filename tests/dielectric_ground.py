#!/usr/bin/env python3
"""M10 PEC integration: explicit images, vacuum compatibility and boundaries."""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

BINARY = Path(os.environ.get('NEC2C','./nec2c')).resolve()
WIRE = 'GW 1 11 0 0 .36 0 0 .84 .001\n'
MIRROR_WIRE = 'GW 2 11 0 0 -.36 0 0 -.84 .001\n'
PATCH = ('SM 2 2 -.2 -.2 .6 .2 -.2 .6\n'
         'SC 0 0 .2 .2 .6 -.2 .2 .6\n')
MIXED_PATCH = ('SM 2 2 .4 -.2 .4 .4 .2 .4\n'
               'SC 0 0 .4 .2 .8 .4 -.2 .8\n')
PATTERN = 'RP 0 3 2 0000 0 0 45 90 0 0'
BOX = (.1,-.02,.58,.14,.02,.62)


def solve(d, geometry='', boxes=(BOX,), er=3, ground='GN 1',
          excitation='EX 1 1 1 0 37 23 19 0 0 .6', output=PATTERN,
          mode='dense', frequency=300, frequencies=1, step=0, success=True):
    resource = 'nec2c-dielectric 1\nMATERIAL 1 ER '+str(er)+'\nREGION 1 MATERIAL 1\n'
    for box in boxes:
        resource += 'BOX 1 '+' '.join(map(str,box))+'\n'
    resource += 'MESH MAX_KH .5 MAX_ASPECT 2 MAX_CELLS 100\nEND\n'
    (d/'model.diel').write_text(resource)
    source = ('CE\n'+geometry+'GE\nFR 0 '+str(frequencies)+' 0 0 '+str(frequency)+' '+str(step)+'\n'+ground+
              '\n'+excitation+'\nDF "model.diel"\n'+output+'\nEN\n')
    (d/'model.nec').write_text(source)
    result=subprocess.run([str(BINARY),'-s',mode,'-r','1e-13','-i',str(d/'model.nec'),
                           '-o',str(d/'model.out')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    text=(d/'model.out').read_text()
    assert (result.returncode==0)==success, text+result.stderr.decode()
    if not success:
        assert 'SOLUTION' not in text, text
    return text


def rows(text,heading,count):
    part=text.split(heading,1)[1]
    result=[]
    for line in part.splitlines():
        try: values=[float(x) for x in line.split()]
        except ValueError: continue
        if values:
            result.append(values)
            if len(result)==count: return result
    raise AssertionError('missing rows '+heading)


def currents(text,count,coupled=False):
    heading='COUPLED VOLUME POLARIZATION CURRENTS' if coupled else 'CELL  REGION MATERIAL'
    return [[complex(*r[i:i+2]) for i in range(len(r)-6,len(r),2)]
            for r in rows(text,heading,count)]


def close(a,b,rtol=2e-7,atol=1e-12):
    assert len(a)==len(b) and len(a)>0,(a,b)
    error=max(abs(x-y) for x,y in zip(a,b))
    scale=max(abs(x) for x in b)
    assert error<=atol+rtol*scale,(a,b,error,scale)


def field(text,magnetic=False):
    r=rows(text,'DIELECTRIC SCATTERED NEAR '+('MAGNETIC' if magnetic else 'ELECTRIC'),1)[0]
    return [complex(*r[i:i+2]) for i in (3,5,7)]


def scalar(text,label):
    return float(re.search(re.escape(label)+r'\s*(\S+)',text)[1])


def test_volume_images(d):
    lower=(BOX[0],BOX[1],-BOX[5],BOX[3],BOX[4],-BOX[2])
    for kind in (1,2,3):
        ex='EX %d 1 1 0 37 23 19 0 0 .6'%kind
        image_kind={1:1,2:3,3:2}[kind]
        image_ex='EX %d 1 1 0 143 23 -19 0 0 .6'%image_kind
        pec=solve(d,excitation=ex)
        direct=solve(d,boxes=(BOX,lower),ground='',excitation=ex)
        reflected=solve(d,boxes=(BOX,lower),ground='',excitation=image_ex)
        jc=currents(pec,1)[0]
        jl,jd=currents(direct,2)
        jlr,jr=currents(reflected,2)
        close(jc,[a+b for a,b in zip(jd,jr)])
        close([(-1 if q<2 else 1)*v for q,v in enumerate(jc)],
              [a+b for a,b in zip(jl,jlr)])
        for pr,dr,rr in zip(rows(pec,'DIELECTRIC SCATTERED FAR FIELD',6),
                            rows(direct,'DIELECTRIC SCATTERED FAR FIELD',6),
                            rows(reflected,'DIELECTRIC SCATTERED FAR FIELD',6)):
            close(pr[2:6],[a+b for a,b in zip(dr[2:6],rr[2:6])])
        iterative=solve(d,excitation=ex,mode='gmres')
        close(currents(iterative,1)[0],jc,rtol=1e-9)
        defect=scalar(pec,'RELATIVE IMBALANCE:')
        assert defect<2e-3,defect
    print('PEC volume: linear/elliptical excitation, explicit image currents/patterns, GMRES and power pass',flush=True)
    # External electric and magnetic image fields, including the ground surface.
    for magnetic in (False,True):
        card=('NH' if magnetic else 'NE')+' 0 1 1 1 .8 .3 0 0 0 0'
        pec=solve(d,output=card)
        direct=solve(d,boxes=(BOX,lower),ground='',output=card)
        reflected=solve(d,boxes=(BOX,lower),ground='',
                        excitation='EX 1 1 1 0 143 23 -19 0 0 .6',output=card)
        got=field(pec,magnetic)
        close(got,[a+b for a,b in zip(field(direct,magnetic),field(reflected,magnetic))])
        if magnetic: assert abs(got[2])<1e-12
        else: assert max(abs(v) for v in got[:2])<1e-12
    print('PEC volume: electric/magnetic images and surface boundary pass',flush=True)


def test_wire_images(d):
    ex='EX 0 1 6 0 1 0'
    pec=solve(d,geometry=WIRE,excitation=ex)
    lower=(BOX[0],BOX[1],-BOX[5],BOX[3],BOX[4],-BOX[2])
    image=solve(d,geometry=WIRE+MIRROR_WIRE,ground='',boxes=(BOX,lower),
                excitation=ex+'\nEX 0 2 6 0 -1 0')
    close(currents(pec,1,True)[0],currents(image,2,True)[1],rtol=3e-7)
    for a,b in zip(rows(pec,'COUPLED WIRE CURRENTS',11),
                   rows(image,'COUPLED WIRE CURRENTS',22)[:11]):
        close(a[-2:],b[-2:],rtol=3e-7)
    for a,b in zip(rows(pec,'COUPLED SCATTERED FAR FIELD',6),
                   rows(image,'COUPLED SCATTERED FAR FIELD',6)):
        close(a[2:],b[2:],rtol=3e-7)
    defect=scalar(pec,'RELATIVE POWER DEFECT:')
    assert defect<.01,defect
    print('PEC wire/volume: explicitly mirrored free-space antenna, pattern and power pass',flush=True)


def test_vacuum_and_patch(d):
    # Vacuum DF must recover the unmodified NEC PEC path. Numeric print precision
    # of legacy wire rows is lower than the coupled diagnostic precision.
    import dielectric_wire as wire
    for name,geometry,ex,box in (
        ('wire',WIRE,'EX 0 1 6 0 1 0',BOX),
        ('patch',PATCH,'EX 1 1 1 0 37 23 19 0 0 .6',(.08,-.02,.95,.12,.02,.99)),
        ('patch-elliptical',PATCH,'EX 2 1 1 0 37 23 19 0 0 .6',(.08,-.02,.95,.12,.02,.99)),
        ('mixed',WIRE+MIXED_PATCH,'EX 0 1 6 0 1 0',(-.14,-.02,.58,-.1,.02,.62))):
        vacuum=solve(d,geometry=geometry,excitation=ex,boxes=(box,),er=1)
        source=(d/'model.nec').read_text().replace('DF "model.diel"\n','')
        (d/'model.nec').write_text(source)
        result=subprocess.run([str(BINARY),'-i',str(d/'model.nec'),'-o',str(d/'legacy.out')],
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        legacy=(d/'legacy.out').read_text()
        assert result.returncode==0,legacy
        if name not in ('patch','patch-elliptical'):
            close([wire.coupled_impedance(vacuum)],[wire.legacy_impedance(legacy)],rtol=6e-5)
            for i in range(1,12):
                close([wire.coupled_current(vacuum,i)],[wire.legacy_current(legacy,i)],rtol=6e-5)
        if name!='wire':
            import dielectric_patch as patch
            got=patch.coupled_patches(vacuum)
            want=patch.legacy_patches(legacy)
            assert len(got)==len(want)==4
            # NEC's legacy Cartesian patch currents are printed with %.2E.
            rounded=[[complex(float(format(v.real,'.2E')),float(format(v.imag,'.2E')))
                      for v in r] for r in got]
            assert rounded==[list(r) for r in want],(rounded,want)
        dense=solve(d,geometry=geometry,excitation=ex,boxes=(box,))
        iterative=solve(d,geometry=geometry,excitation=ex,boxes=(box,),mode='gmres')
        close(currents(dense,1,True)[0],currents(iterative,1,True)[0],rtol=2e-8)
        for a,b in zip(rows(dense,'COUPLED SCATTERED FAR FIELD',6),
                       rows(iterative,'COUPLED SCATTERED FAR FIELD',6)):
            close(a[2:],b[2:],rtol=2e-8)
        print('PEC '+name+': vacuum legacy identity and loaded dense/GMRES pass',flush=True)


def test_patch_images(d):
    import dielectric_patch as patch
    mirror=('SM 2 2 -.2 -.2 -.6 -.2 .2 -.6\n'
            'SC 0 0 .2 .2 -.6 .2 -.2 -.6\n')
    box=(.08,-.02,.95,.12,.02,.99)
    lower=(box[0],box[1],-box[5],box[3],box[4],-box[2])
    pec=solve(d,geometry=PATCH,boxes=(box,))
    direct=solve(d,geometry=PATCH+mirror,boxes=(box,lower),ground='')
    reflected=solve(d,geometry=PATCH+mirror,boxes=(box,lower),ground='',
                    excitation='EX 1 1 1 0 143 23 -19 0 0 .6')
    close(currents(pec,1,True)[0],
          [a+b for a,b in zip(currents(direct,2,True)[1],currents(reflected,2,True)[1])])
    got=patch.coupled_patches(pec)
    a=patch.coupled_patches(direct); b=patch.coupled_patches(reflected)
    close([v for r in got for v in r],
          [u+v for r,s in zip(a[:4],b[:4]) for u,v in zip(r,s)])
    for pr,dr,rr in zip(rows(pec,'COUPLED SCATTERED FAR FIELD',6),
                        rows(direct,'COUPLED SCATTERED FAR FIELD',6),
                        rows(reflected,'COUPLED SCATTERED FAR FIELD',6)):
        close(pr[2:],[a+b for a,b in zip(dr[2:],rr[2:])])
    print('PEC patch/volume: explicit mirrored patches, volume currents and patterns pass',flush=True)


def test_frequency_state(d):
    kwargs={'geometry':WIRE,'excitation':'EX 0 1 6 0 1 0'}
    sweep=solve(d,frequency=290,frequencies=2,step=10,**kwargs)
    parts=sweep.split('---------- COUPLED CONDUCTOR-DIELECTRIC SOLUTION (M7) ----------')[1:]
    assert len(parts)==2,len(parts)
    for frequency,part in zip((290,300),parts):
        independent=solve(d,frequency=frequency,**kwargs)
        close(currents(part,1,True)[0],currents(independent,1,True)[0],rtol=1e-10)
        for a,b in zip(rows(part,'COUPLED WIRE CURRENTS',11),
                       rows(independent,'COUPLED WIRE CURRENTS',11)):
            close(a[-2:],b[-2:],rtol=1e-10)
    print('PEC frequency sweep preserves ground and restores physical geometry between solves',flush=True)


def test_rejections(d):
    for kwargs,fragment in (
        ({'ground':'GN 0 0 0 0 13 .005'},'GN 0'),
        ({'ground':'GN 2 0 0 0 13 .005','geometry':PATCH.replace('.6','1.6')},'wires only'),
        ({'mode':'fft','boxes':((.1,-.02,.581,.14,.02,.621),)},'half-cell lattice'),
        ({'ground':'GN 1 4 0 0 0 0 1 .001'},'without screens'),
        ({'output':'NE 0 1 1 1 .8 .3 -.1 0 0 0'},'z >= 0'),
        ({'output':'NE 1 1 1 1 .8 0 120 0 0 0'},'z >= 0'),
        ({'excitation':'EX 2 1 1 0 37 23 19'},'explicit positive axial ratio'),
        ({'boxes':((.1,0,-.01,.14,.04,.03),)},'z >= 0'),
        ({'output':'RP 0 3 1 0000 0 0 60 0 0 0'},'upper-hemisphere'),
        ({'excitation':'EX 1 1 1 0 120 0 0'},'incident plane wave'),
        ({'geometry':'GW 1 11 0 0 0 0 0 .48 .001\n','excitation':'EX 0 1 6 0 1 0'},'wire endpoint'),
        ({'geometry':'SM 1 1 0 0 .01 .04 0 .01\nSC 0 0 .04 .04 .01 0 .04 .01\n'},'sqrt(area)'),
    ):
        text=solve(d,success=False,**kwargs)
        assert fragment.lower() in text.lower(),text
    print('PEC unsupported ground/backend/domain/contact requests reject explicitly',flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix='nec2c-m10-pec-') as tmp:
        d=Path(tmp)
        test_volume_images(d)
        test_wire_images(d)
        test_vacuum_and_patch(d)
        test_patch_images(d)
        test_frequency_state(d)
        test_rejections(d)

if __name__=='__main__': main()
