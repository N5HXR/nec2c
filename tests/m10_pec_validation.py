#!/usr/bin/env python3
"""PEC physical checks beyond image identities: ports, power and remote ground."""
from pathlib import Path
import tempfile
import dielectric_ground as g
from dielectric_wire import coupled_current, coupled_impedance


def closed_box(n=2):
    # Outward normals: NEC's magnetic-field patch equation requires a closed
    # conductor for a physical scattering/power test, rather than an open sheet.
    faces = [
        [( .3,-.1,.6),(.3,-.1,.8),(.3,.1,.8),(.3,.1,.6)],
        [( .5,-.1,.6),(.5,.1,.6),(.5,.1,.8),(.5,-.1,.8)],
        [( .3,-.1,.6),(.5,-.1,.6),(.5,-.1,.8),(.3,-.1,.8)],
        [( .3, .1,.6),(.3,.1,.8),(.5,.1,.8),(.5,.1,.6)],
        [( .3,-.1,.6),(.3,.1,.6),(.5,.1,.6),(.5,-.1,.6)],
        [( .3,-.1,.8),(.5,-.1,.8),(.5,.1,.8),(.3,.1,.8)],
    ]
    return ''.join(f'SM {n} {n} '+' '.join(str(v) for p in f[:2] for v in p)+
                   '\nSC 0 0 '+' '.join(str(v) for p in f[2:] for v in p)+'\n' for f in faces)


def main():
    with tempfile.TemporaryDirectory(prefix='m10-pec-physical-') as tmp:
        d=Path(tmp)
        # Unequal wire lengths and offsets prevent geometrical exchange
        # symmetry from making the transfer-admittance check tautological.
        geom=g.WIRE+'GW 2 11 .5 .1 .42 .5 .1 .86 .001\n'
        for er in (1,3):
            a=g.solve(d,geometry=geom,er=er,excitation='EX 0 1 6 0 1 0',output='XQ')
            b=g.solve(d,geometry=geom,er=er,excitation='EX 0 2 6 0 1 0',output='XQ')
            x,y=coupled_current(a,17),coupled_current(b,6)
            error=abs(x-y)/max(abs(x),abs(y))
            # Physical ports inherit NEC collocation/basis error (vacuum
            # reference ~7.8e-5); the kernel identity has its own stricter test.
            assert error<1e-4,(er,error)
            print('PEC transfer-admittance reciprocity:',er,error,flush=True)
        free=coupled_impedance(g.solve(d,geometry=g.WIRE,ground='',
                                      excitation='EX 0 1 6 0 1 0',output='XQ'))
        errors=[]
        for h in (0,4,16):
            box=tuple(v+(h if i in (2,5) else 0) for i,v in enumerate(g.BOX))
            geom=f'GW 1 11 0 0 {h+.36} 0 0 {h+.84} .001\n'
            z=coupled_impedance(g.solve(d,geometry=geom,boxes=(box,),
                                       excitation='EX 0 1 6 0 1 0',output='XQ'))
            errors.append(abs(z-free))
        assert errors[1]<.02*errors[0] and errors[2]<.1*errors[1],errors
        assert errors[2]<.005,errors
        print('PEC remote-ground feed errors (ohm):',errors,flush=True)
        for name,geom in [('wire',g.WIRE),('closed patches',closed_box()),
                          ('mixed',g.WIRE+closed_box())]:
            for er in (1,3):
                for kind in (1,2,3):
                    text=g.solve(d,geometry=geom,er=er,output='XQ',
                                 excitation=f'EX {kind} 1 1 0 37 23 19 0 0 .6')
                    defect=g.scalar(text,'RELATIVE POWER DEFECT:')
                    assert g.scalar(text,'EXTINCTION POWER:')>0,text
                    assert defect<.003,(name,er,kind,defect)
                    print('PEC plane-wave power:',name,er,kind,defect,flush=True)
        text=g.solve(d,geometry=closed_box(4),output='XQ')
        defect=g.scalar(text,'RELATIVE POWER DEFECT:')
        assert defect<.003,defect
        print('PEC refined closed-patch power defect:',defect,flush=True)


if __name__=='__main__': main()
