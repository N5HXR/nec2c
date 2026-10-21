#!/usr/bin/env python3
"""PEC contact integration: image equivalence, gap limit and FFT agreement."""
from pathlib import Path
import tempfile
import dielectric_ground as g
from dielectric_wire import coupled_impedance,coupled_current
from m10_pec_validation import closed_box

CONTACT=(.1,-.02,0,.14,.02,.04)


def shifted(gap):
    return tuple(v+(gap if i in (2,5) else 0) for i,v in enumerate(CONTACT))


def test_boundary_reciprocity(d):
    geom=g.WIRE+'GW 2 11 .5 .1 .42 .5 .1 .86 .001\n'
    a=g.solve(d,geometry=geom,boxes=(CONTACT,),excitation='EX 0 1 6 0 1 0',output='XQ')
    b=g.solve(d,geometry=geom,boxes=(CONTACT,),excitation='EX 0 2 6 0 1 0',output='XQ')
    x,y=coupled_current(a,17),coupled_current(b,6)
    error=abs(x-y)/max(abs(x),abs(y))
    # Preserve the same physical port allowance as M10's native-NEC baseline.
    assert error<1e-4,error
    print('PEC contact port reciprocity:',error,flush=True)
    for magnetic in (False,True):
        card=('NH' if magnetic else 'NE')+' 0 1 1 1 .8 .3 0 0 0 0'
        field=g.field(g.solve(d,geometry=g.WIRE+closed_box(),boxes=(CONTACT,),output=card),magnetic)
        residual=abs(field[2]) if magnetic else max(abs(field[0]),abs(field[1]))
        assert residual<1e-14+1e-10*max(abs(v) for v in field),(magnetic,field)
        print('PEC contact mixed boundary:',magnetic,residual,flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix='m11-contact-') as tmp:
        d=Path(tmp)
        lower=(.1,-.02,-.04,.14,.02,0)
        for kind in (1,2,3):
            ex=f'EX {kind} 1 1 0 37 23 19 0 0 .6'
            inverse={1:1,2:3,3:2}[kind]
            contact=g.solve(d,boxes=(CONTACT,),excitation=ex)
            direct=g.solve(d,boxes=(CONTACT,lower),ground='',excitation=ex)
            mirror=g.solve(d,boxes=(CONTACT,lower),ground='',
                           excitation=f'EX {inverse} 1 1 0 143 23 -19 0 0 .6')
            j=g.currents(contact,1)[0]
            g.close(j,[a+b for a,b in zip(g.currents(direct,2)[1],g.currents(mirror,2)[1])])
            fast=g.solve(d,boxes=(CONTACT,),excitation=ex,mode='fft')
            g.close(g.currents(fast,1)[0],j,rtol=2e-7)
            assert g.scalar(contact,'RELATIVE IMBALANCE:')<.002
        print('PEC contact: linear/elliptical explicit images, FFT currents and power pass',flush=True)
        exact=g.currents(g.solve(d,boxes=(CONTACT,)),1)[0]
        errors=[]
        for gap in (.004,.001,.0001,.00001,.000001):
            current=g.currents(g.solve(d,boxes=(shifted(gap),)),1)[0]
            errors.append(max(abs(a-b) for a,b in zip(current,exact)))
        for a,b in zip(errors,errors[1:]): assert b<.5*a,errors
        assert errors[-1]<1e-4*max(abs(v) for v in exact),errors
        print('PEC solved-current approach to contact:',errors,flush=True)
        for name,geom in [('wire',g.WIRE),('patch',closed_box()),('mixed',g.WIRE+closed_box())]:
            dense=g.solve(d,geometry=geom,boxes=(CONTACT,),output='XQ')
            fast=g.solve(d,geometry=geom,boxes=(CONTACT,),output='XQ',mode='fft')
            g.close(g.currents(fast,1,True)[0],g.currents(dense,1,True)[0],rtol=2e-7)
            assert g.scalar(dense,'RELATIVE POWER DEFECT:')<.003,(name,dense)
            print('PEC contact coupled dense/FFT and power:',name,flush=True)
        kwargs=dict(geometry=g.WIRE,excitation='EX 0 1 6 0 1 0',output='XQ')
        exact_text=g.solve(d,boxes=(CONTACT,),**kwargs)
        exact=coupled_impedance(exact_text)
        fast=coupled_impedance(g.solve(d,boxes=(CONTACT,),mode='fft',**kwargs))
        assert abs(fast-exact)<1e-6,(fast,exact)
        errors=[]
        for gap in (.001,.0001,.00001):
            errors.append(abs(coupled_impedance(g.solve(d,boxes=(shifted(gap),),**kwargs))-exact))
        assert errors[1]<.2*errors[0] and errors[2]<.2*errors[1],errors
        print('PEC contact driven feed continuity errors (ohm):',errors,flush=True)
        text=g.solve(d,boxes=(CONTACT,),ground='GN 2 0 0 0 13 .02',success=False)
        assert 'positive clearance' in text,text
        text=g.solve(d,boxes=(shifted(-1e-6),),success=False)
        assert 'z >= 0' in text,text
        print('Real-ground contact and buried volume remain rejected',flush=True)
        test_boundary_reciprocity(d)


if __name__=='__main__': main()
