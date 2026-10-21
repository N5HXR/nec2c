#!/usr/bin/env python3
"""Independent closed-form static PEC self reference for aligned cuboids.

Integrates the Coulomb potential between parallel rectangular faces using
analytic antiderivatives of 1/r, u/r, v/r and uv/r. No production quadrature,
Duffy map, NumPy, or SciPy is used. Prints diagonal N operators (curl curl G).
"""
import math


def primitive(x,y,d):
    r=math.sqrt(x*x+y*y+d*d)
    if r==0: return (0.,0.,0.,0.)
    a=x*math.asinh(y/math.hypot(x,d)) if x else 0.
    b=y*math.asinh(x/math.hypot(y,d)) if y else 0.
    f00=a+b-(d*math.atan(x*y/(d*r)) if d else 0.)
    f10=.5*(y*r+(x*x+d*d)*math.asinh(y/math.hypot(x,d))) if x or d else .5*y*abs(y)
    f01=.5*(x*r+(y*y+d*d)*math.asinh(x/math.hypot(y,d))) if y or d else .5*x*abs(x)
    return f00,f10,f01,r*r*r/3


def segments(length,offset):
    return ((offset-length,offset,length-offset,1),
            (offset,offset+length,length+offset,-1))


def face_integral(a,b,offset_a,offset_b,separation):
    terms=[]
    for lo,hi,A,B in segments(a,offset_a):
        for low,high,C,D in segments(b,offset_b):
            for x,sx in ((lo,-1),(hi,1)):
                for y,sy in ((low,-1),(high,1)):
                    f00,f10,f01,f11=primitive(x,y,abs(separation))
                    terms.append(sx*sy*(A*C*f00+B*C*f10+A*D*f01+B*D*f11))
    return math.fsum(terms)/(4*math.pi)


def diagonal(size,delta):
    out=[0.,0.,0.];volume=size[0]*size[1]*size[2]
    for axis in range(3):
        u,v=[j for j in range(3) if j!=axis]
        terms=[]
        for ts in (-1,1):
            for ss in (-1,1):
                d=delta[axis]+(ts-ss)*size[axis]/2
                terms.append(ts*ss*face_integral(size[u],size[v],delta[u],delta[v],d)/volume)
        value=math.fsum(terms)
        for p in range(3):
            if p!=axis: out[p]+=value
    return out


def reference(size,gap):
    direct=diagonal(size,(0,0,0))
    reflected=diagonal(size,(0,0,size[2]+2*gap))
    return dict(size_m=size,gap_m=gap,direct=direct,reflected=reflected,
                pec=[direct[p]+(-1 if p<2 else 1)*reflected[p] for p in range(3)])


if __name__=='__main__':
    cube=diagonal((1,1,1),(0,0,0))
    assert max(abs(x-2/3) for x in cube)<1e-14,cube
    print('# Independent analytic static PEC self reference; sx,sy,sz,gap,Nxx,Nyy,Nzz')
    for size,gaps in (((.3,.4,.2),(0,.001,.0001,.00001,.000001)),
                      ((.07,.12,.03),(0,.0001,.000001))):
        for gap in gaps:
            row=reference(size,gap)
            print(','.join(format(v,'.17g') for v in tuple(size)+(gap,)+tuple(row['pec'])))
