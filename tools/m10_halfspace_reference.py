#!/usr/bin/env python3
"""Independent axial Sommerfeld oracle: analytic azimuth, adaptive Simpson.

No production code, Bessel routine, NumPy or SciPy is used. exp(+i omega t).
Prints reflected Nxx, Nzz and Hxy for a point current moment at rho=0.
"""
import cmath
import math


def simpson(f, a, b, tolerance=1e-12):
    def refine(a,b,fa,fm,fb,old,tol,depth):
        mid=(a+b)/2
        fl=f((a+mid)/2); fr=f((mid+b)/2)
        left=(mid-a)*(fa+4*fl+fm)/6
        right=(b-mid)*(fm+4*fr+fb)/6
        delta=left+right-old
        if abs(delta)<15*tol:
            return left+right+delta/15
        if depth==0:
            raise RuntimeError('independent oracle failed to converge')
        return (refine(a,mid,fa,fl,fm,left,tol/2,depth-1)+
                refine(mid,b,fm,fr,fb,right,tol/2,depth-1))
    fa,fm,fb=f(a),f((a+b)/2),f(b)
    return refine(a,b,fa,fm,fb,(b-a)*(fa+4*fm+fb)/6,tolerance,24)


def reference(k, height, epsilon):
    def integrand(t, evanescent, component):
        if evanescent:
            v=t/(k*height); u=math.hypot(1,v); w=-1j*v
            jac=1j/(k*height)
        else:
            u=math.sin(t); w=math.cos(t); jac=u
        lower=cmath.sqrt(epsilon-u*u)
        if lower.imag>0: lower=-lower
        rs=(w-lower)/(w+lower)
        rp=(epsilon*w-lower)/(epsilon*w+lower)
        phase=cmath.exp(-1j*k*w*height)
        if component==0: dyad=rs-rp*w*w
        elif component==1: dyad=2*rp*u*u
        else: dyad=w*(rp-rs)/(1j*k)
        return -1j*k**3/(8*math.pi)*jac*phase*dyad
    result=[]
    for component in range(3):
        p=simpson(lambda t: integrand(t,False,component),0,math.pi/2)
        intervals=(0,1,4,16,64)
        e=sum(simpson(lambda t: integrand(t,True,component),a,b)
              for a,b in zip(intervals,intervals[1:]))
        result.append(p+e)
    return result

if __name__=='__main__':
    for label,value in zip(('Nxx','Nzz','Hxy'),reference(2.1,.7,4-.4j)):
        print('{} {:.17g} {:+.17g}j'.format(label,value.real,value.imag))
