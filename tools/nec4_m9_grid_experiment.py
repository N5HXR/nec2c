#!/usr/bin/env python3
"""Measure transformed trilinear interpolation on C direct-oracle grids."""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nec4_m9_tables import storage_factor  # noqa: E402
from nec4_m9_dispatch import dispatch, MODEL_FIT  # noqa: E402
from nec4_m8_equations import material  # noqa: E402


def load(path):
    with path.open(encoding="utf-8") as stream:
        metadata = stream.readline().removeprefix("# ").strip().split()
        metadata = {key: value for key, value in (item.split("=") for item in metadata)}
        rows = list(csv.DictReader(stream))
    return metadata, rows


def values(row):
    return np.asarray([complex(float(row[f"c{i}_real"]),
                               float(row[f"c{i}_imag"])) for i in range(5)])


def factor(point, k1, k2, mode):
    if mode == "raw": return 1 + 0j
    rho, z, depth = point
    if rho == 0 and z == 0 and depth == 0: return 0j
    return storage_factor(int(mode[1]), rho, z, -depth, k1, k2)


def evaluate(path, mode, direct_strip_lambda):
    meta, rows = load(path)
    f=float(meta["frequency"]); er=float(meta["epsilon"]); sig=float(meta["conductivity"])
    wavelength=float(meta["wavelength"]); lg=float(meta["lg"])
    nr,nz,nd=(int(meta[name]) for name in ("nr","nz","nd"))
    lower=material(f,er,conductivity=sig);upper=material(f,1)
    k1,k2=complex(lower.wave_number),complex(upper.wave_number)
    train=[row for row in rows if row["kind"]=="train"]
    heldout=[row for row in rows if row["kind"]=="heldout"]
    grid=np.empty((nr,nz,nd,5),complex)
    for index,row in enumerate(train):
        ir=index//(nz*nd);iz=(index//nd)%nz;identity=index%nd
        point=tuple(float(row[name]) for name in ("rho","z","depth"))
        scale=factor(point,k1,k2,mode)
        grid[ir,iz,identity]=values(row)*scale if scale else 0
    axes=(np.linspace(0,3*wavelength,nr),np.linspace(0,2*wavelength,nz),
          np.linspace(0,.8*lg,nd))
    interp=RegularGridInterpolator(axes,grid,method="linear")
    errors=[];component_errors=[];local_peak_errors=[];locations=[];accepted=0
    for row in heldout:
        point=tuple(float(row[name]) for name in ("rho","z","depth"))
        route=dispatch(point[0],point[1],-point[2],k1,k2,wavelength)
        if route.method!=MODEL_FIT or point[1]<direct_strip_lambda*wavelength: continue
        expected=values(row);scale=factor(point,k1,k2,mode)
        actual=interp(point)/scale
        error=np.linalg.norm(actual-expected)/np.linalg.norm(expected)
        errors.append(error);locations.append(point)
        denominator=np.maximum(np.abs(expected),np.max(np.abs(expected))*1e-11)
        component_errors.extend(np.abs(actual-expected)/denominator)
        local_peak_errors.extend(np.abs(actual-expected)/np.max(np.abs(expected)))
        accepted+=1
    order=np.argsort(errors)
    worst=int(order[-1])
    error_array=np.asarray(errors)
    location_array=np.asarray(locations)
    above={}
    for threshold in (.005,.01,.02,.05):
        selected=location_array[error_array>threshold]
        above[f"{threshold:.3f}"]={
            "count":int(len(selected)),"fraction":float(len(selected)/len(errors)),
            "minimum_normalized":None if not len(selected) else
              list((selected.min(axis=0)/np.asarray((wavelength,wavelength,lg))).astype(float)),
            "maximum_normalized":None if not len(selected) else
              list((selected.max(axis=0)/np.asarray((wavelength,wavelength,lg))).astype(float))}
    return {"material":path.name,"transform":mode,"heldout_points":accepted,
            "direct_strip_lambda0":direct_strip_lambda,
            "rms_vector_error":float(np.sqrt(np.mean(np.square(errors)))),
            "maximum_vector_error":float(np.max(errors)),
            "p95_vector_error":float(np.percentile(errors,95)),
            "p99_vector_error":float(np.percentile(errors,99)),
            "maximum_null_aware_component_error":float(np.max(component_errors)),
            "maximum_component_error_relative_to_local_peak":
              float(np.max(local_peak_errors)),
            "worst_point_m":{"rho":locations[worst][0],"z":locations[worst][1],
                             "source_depth":locations[worst][2]},
            "worst_point_normalized":{"rho_lambda0":locations[worst][0]/wavelength,
              "z_lambda0":locations[worst][1]/wavelength,
              "source_depth_lg":locations[worst][2]/lg},
            "threshold_envelopes_normalized_rho_z_depth":above}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs",type=Path,nargs="+")
    parser.add_argument("--direct-strip-lambda",type=float,default=0.)
    parser.add_argument("--transform",action="append",
                        choices=("raw","f1","f2","f3"))
    parser.add_argument("--output",type=Path,default=ROOT/"tests/nec4-similarity"/"m9-grid-screen.json")
    args=parser.parse_args();results=[]
    for path in args.inputs:
        for mode in (args.transform or ("raw","f1","f2","f3")):
            results.append(evaluate(path,mode,args.direct_strip_lambda))
    args.output.write_text(json.dumps({"schema":"nec2c-nec4-m9-grid-screen",
        "schema_version":1,"results":results},indent=2)+"\n",encoding="utf-8")
    for row in results: print(f"{row['material']} {row['transform']}: "
        f"RMS {row['rms_vector_error']:.3%}, max {row['maximum_vector_error']:.3%}, "
        f"component {row['maximum_null_aware_component_error']:.3%} ({row['heldout_points']} points)")
    print(args.output)


if __name__=="__main__": main()
