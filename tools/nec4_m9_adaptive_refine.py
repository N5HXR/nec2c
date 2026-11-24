#!/usr/bin/env python3
"""Simulate center-certified level-2 to level-3 adaptive F2 refinement."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"tools"))
from nec4_m8_equations import material  # noqa: E402
from nec4_m9_dispatch import dispatch,MODEL_FIT  # noqa: E402
from nec4_m9_grid_experiment import load,values,factor  # noqa: E402


def grid(path):
    meta,rows=load(path);nr,nz,nd=(int(meta[name]) for name in ("nr","nz","nd"))
    f=float(meta["frequency"]);er=float(meta["epsilon"]);sig=float(meta["conductivity"])
    wavelength=float(meta["wavelength"]);lg=float(meta["lg"])
    lower=material(f,er,conductivity=sig);upper=material(f,1)
    k1,k2=complex(lower.wave_number),complex(upper.wave_number)
    data=np.empty((nr,nz,nd,5),complex)
    training=[row for row in rows if row["kind"]=="train"]
    for index,row in enumerate(training):
        ir=index//(nz*nd);iz=(index//nd)%nz;identity=index%nd
        point=tuple(float(row[name]) for name in ("rho","z","depth"))
        scale=factor(point,k1,k2,"f2")
        data[ir,iz,identity]=values(row)*scale if scale else 0
    axes=(np.linspace(0,3*wavelength,nr),np.linspace(0,2*wavelength,nz),
          np.linspace(0,.8*lg,nd))
    return meta,rows,k1,k2,wavelength,lg,data,axes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level2",type=Path,required=True)
    parser.add_argument("--centers",type=Path,required=True)
    parser.add_argument("--level3",type=Path,required=True)
    parser.add_argument("--direct-strip-lambda",type=float,default=.16)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    _,_,k1,k2,wavelength,lg,data2,axes2=grid(args.level2)
    _,rows3,_,_,_,_,data3,axes3=grid(args.level3)
    _,center_rows=load(args.centers)
    interp2=RegularGridInterpolator(axes2,data2,method="linear")
    interp3=RegularGridInterpolator(axes3,data3,method="linear")
    intervals=np.asarray((len(axes2[0])-1,len(axes2[1])-1,len(axes2[2])-1))
    extents=np.asarray((3*wavelength,2*wavelength,.8*lg))
    candidates=[]
    for row in center_rows:
        if row["kind"]!="heldout":continue
        point=np.asarray(tuple(float(row[name]) for name in ("rho","z","depth")))
        route=dispatch(*point[:2],-point[2],k1,k2,wavelength)
        if route.method!=MODEL_FIT or point[1]<args.direct_strip_lambda*wavelength:continue
        expected=values(row);actual=interp2(point)/factor(point,k1,k2,"f2")
        error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
        cell=tuple(np.minimum((point/extents*intervals).astype(int),intervals-1))
        candidates.append((cell,error))
    validations=[row for row in rows3 if row["kind"]=="heldout"]
    reports=[]
    for threshold in (.001,.0025,.005):
        refined={cell for cell,error in candidates if error>threshold}
        new_nodes=set()
        for ir,iz,identity in refined:
            for jr in range(2*ir,2*ir+3):
                for jz in range(2*iz,2*iz+3):
                    for jd in range(2*identity,2*identity+3):
                        if jr%2 or jz%2 or jd%2:new_nodes.add((jr,jz,jd))
        errors=[];local=[];count=0
        for row in validations:
            point=np.asarray(tuple(float(row[name]) for name in ("rho","z","depth")))
            route=dispatch(*point[:2],-point[2],k1,k2,wavelength)
            if route.method!=MODEL_FIT or point[1]<args.direct_strip_lambda*wavelength:continue
            cell=tuple(np.minimum((point/extents*intervals).astype(int),intervals-1))
            interpolator=interp3 if cell in refined else interp2
            expected=values(row);actual=interpolator(point)/factor(point,k1,k2,"f2")
            errors.append(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
            local.extend(np.abs(actual-expected)/np.max(np.abs(expected)));count+=1
        reports.append({"center_threshold":threshold,"refined_cells":len(refined),
          "total_cells":int(np.prod(intervals)),"additional_level3_nodes":len(new_nodes),
          "base_level2_nodes":int(data2.shape[0]*data2.shape[1]*data2.shape[2]),
          "validation_points":count,"rms_vector_error":float(np.sqrt(np.mean(np.square(errors)))),
          "maximum_vector_error":float(np.max(errors)),
          "maximum_component_error_relative_to_local_peak":float(np.max(local))})
    args.output.write_text(json.dumps({"schema":"nec2c-nec4-m9-adaptive-refine",
      "schema_version":1,"direct_strip_lambda0":args.direct_strip_lambda,
      "results":reports},indent=2)+"\n",encoding="utf-8")
    for row in reports:print(f"threshold {row['center_threshold']:.3%}: "
      f"{row['refined_cells']}/{row['total_cells']} cells, "
      f"{row['base_level2_nodes']+row['additional_level3_nodes']} nodes, "
      f"RMS {row['rms_vector_error']:.3%}, max {row['maximum_vector_error']:.3%}, "
      f"component {row['maximum_component_error_relative_to_local_peak']:.3%}")
    print(args.output)


if __name__=="__main__":main()
