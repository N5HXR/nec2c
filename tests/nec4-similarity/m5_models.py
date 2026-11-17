#!/usr/bin/env python3
"""M5 free-end, source-gap, load-gap, and reciprocity acceptance corpus."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def run(deck, engine, directory, stem, jobs=1):
    model = directory / f"{stem}.nec"
    result = directory / f"{stem}-{engine}-{jobs}.json"
    report = directory / f"{stem}-{engine}-{jobs}.out"
    model.write_text(deck, encoding="ascii")
    command = [str(BINARY), "--wire-engine", engine, "--ground-engine", "nec2",
               "--backend", "builtin", "--jobs", str(jobs), "--threads", "1",
               "--result-json", str(result), "-i", str(model), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise AssertionError((command, completed.stdout, completed.stderr,
                              report.read_text(errors="replace")[-2000:]
                              if report.exists() else ""))
    return json.loads(result.read_text(encoding="utf-8")), report.read_text(
        encoding="utf-8", errors="replace")


def monopole(count):
    return ("CM M5 END CAP AND SOURCE GAP REFINEMENT\nCE\n"
            f"GW 1 {count} 0 0 0 0 0 .25 .01\nGE 1\nGN 1\n"
            "EX 0 1 1 0 1 0\nFR 0 1 0 0 299.792458 0\nXQ\nEN\n")


def two_dipoles(source_tag, count=31):
    center=(count+1)//2
    return ("CM M5 RECIPROCITY\nCE\n"
            f"GW 1 {count} -.1 0 -.24 -.1 0 .24 .001\n"
            f"GW 2 {count} .1 0 -.24 .1 0 .24 .001\nGE 0\n"
            "FR 0 1 0 0 299.792458 0\n"
            f"EX 0 {source_tag} {center} 0 1 0\n"
            "RP 0 5 8 1000 0 0 22.5 45 1000 0\nEN\n")


def mismatched_source():
    return ("CM M5 SOURCE LENGTH VALIDITY\nCE\n"
            "GW 1 1 0 0 -.24 0 0 -.02 .001\n"
            "GW 2 1 0 0 -.02 0 0 .02 .001\n"
            "GW 3 1 0 0 .02 0 0 .24 .001\nGE 0\n"
            "FR 0 1 0 0 299.792458 0\nEX 0 2 1 0 1 0\nXQ\nEN\n")


def impedance(document):
    value=document["jobs"][0]["executions"][0]["feeds"][0]["impedance_ohm"]
    return complex(value["re"],value["im"])


def currents(document):
    rows=document["jobs"][0]["executions"][0]["wire_currents"]
    return [complex(row["current_a"]["re"],row["current_a"]["im"])
            for row in rows]


def ripple(values):
    scale=max(map(abs,values))
    return max(abs(values[i]-.5*(values[i-1]+values[i+1]))
               for i in range(1,len(values)-1))/scale


def numeric_delta(left,right):
    values=[]
    if isinstance(left,dict) and isinstance(right,dict):
        for key in left.keys() & right.keys():
            if key not in {"elapsed_seconds","threads","requested","used"}:
                values.extend(numeric_delta(left[key],right[key]))
    elif isinstance(left,list) and isinstance(right,list):
        if len(left)!=len(right): raise AssertionError("result length mismatch")
        for a,b in zip(left,right): values.extend(numeric_delta(a,b))
    elif (isinstance(left,(int,float)) and not isinstance(left,bool) and
          isinstance(right,(int,float)) and not isinstance(right,bool)):
        values.append(abs(left-right)/max(1.,abs(left),abs(right)))
    return values


def main():
    evidence={"schema":"nec2c-nec4-m5-acceptance","schema_version":1}
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m5-") as name:
        directory=Path(name)
        sweep={engine:[] for engine in ("nec2","nec4")}
        for count in (10,20,40,80,120):
            for engine in sweep:
                document,_=run(monopole(count),engine,directory,
                               f"monopole-{count}-{engine}")
                sweep[engine].append({"segments":count,
                                      "length_radius":25/count,
                                      "impedance_ohm":{
                                          "re":impedance(document).real,
                                          "im":impedance(document).imag},
                                      "current_ripple":ripple(currents(document))})
        if sweep["nec4"][3]["current_ripple"] >= .2*sweep["nec2"][3]["current_ripple"]:
            raise AssertionError(("cap ripple reduction",sweep))
        evidence["source_refinement"]=sweep

        _,warning=run(mismatched_source(),"nec4",directory,"source-mismatch")
        if "1 SOURCE LENGTH MISMATCHES" not in warning:
            raise AssertionError("source-length mismatch was not diagnosed")
        model=directory/"source-mismatch-strict.nec"
        model.write_text(mismatched_source(),encoding="ascii")
        strict=subprocess.run([str(BINARY),"--wire-engine","nec4",
            "--model-check","strict","-i",str(model),"-o",
            str(directory/"source-mismatch-strict.out")],cwd=ROOT,text=True,
            capture_output=True)
        if strict.returncode==0: raise AssertionError("strict mismatch accepted")
        evidence["source_length_mismatch"]={"warned":True,"strict_rejected":True}

        ex5=monopole(10).replace("EX 0 1 1", "EX 5 1 1")
        _,ex5_report=run(ex5,"nec4",directory,"legacy-ex5")
        if "1 LEGACY EX5 SOURCES" not in ex5_report:
            raise AssertionError("EX5 audit marker missing")
        evidence["legacy_ex5"]={"retained":True,"gap_cap_correction":False}

        transfers=[]
        center=15
        for source,target in ((1,2),(2,1)):
            document,_=run(two_dipoles(source),"nec4",directory,
                           f"reciprocity-{source}")
            values=currents(document)
            transfers.append(values[(target-1)*31+center])
        reciprocity=abs(transfers[0]-transfers[1])/max(1.,*map(abs,transfers))
        if reciprocity>1e-11: raise AssertionError(("reciprocity",transfers,reciprocity))
        evidence["reciprocity"]={"transfer_currents_a":[
            {"re":v.real,"im":v.imag} for v in transfers],"scaled_defect":reciprocity}

        loaded=(ROOT/"tests/performance/m6-loaded-network-402.nec").read_text()
        serial,_=run(loaded,"nec4",directory,"loaded-network",1)
        parallel,_=run(loaded,"nec4",directory,"loaded-network",8)
        delta=max(numeric_delta(serial,parallel),default=0.)
        if delta>1e-12: raise AssertionError(("parallel load/network",delta))
        power=serial["jobs"][0]["executions"][0]["power"]
        accounted=sum(power[key] or 0. for key in
                      ("air_radiated_w","structure_loss_w","network_loss_w"))
        defect=abs(power["accepted_w"]-accounted)/power["accepted_w"]
        if defect>1e-6: raise AssertionError(("power",power,defect))
        evidence["loaded_network"]={"serial_parallel_scaled_delta":delta,
                                     "power_balance_relative":defect,
                                     "power":power}

    output=ROOT/"tests/nec4-similarity/m5-acceptance.json"
    output.write_text(json.dumps(evidence,indent=2)+"\n",encoding="utf-8")
    try:
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,4.2))
        for engine,style in (("nec2","o--"),("nec4","o-")):
            rows=sweep[engine]; x=[r["length_radius"] for r in rows]
            axes[0].plot(x,[r["current_ripple"] for r in rows],style,label=engine)
            axes[1].plot(x,[r["impedance_ohm"]["re"] for r in rows],style,
                         label=f"{engine} R")
            axes[1].plot(x,[r["impedance_ohm"]["im"] for r in rows],style,
                         label=f"{engine} X")
        axes[0].set_yscale("log");axes[0].set_xlabel("segment length / radius")
        axes[0].set_ylabel("normalized current ripple")
        axes[1].set_xlabel("segment length / radius");axes[1].set_ylabel("ohms")
        for axis in axes: axis.grid(True,alpha=.3);axis.legend()
        fig.tight_layout();fig.savefig(ROOT/"tests/nec4-similarity/m5-cap-convergence.png",
                                      dpi=150);plt.close(fig)
    except ImportError:
        pass
    print("M5 source refinement:")
    for engine,rows in sweep.items():
        print(engine+": "+", ".join(
            f"L/a={r['length_radius']:.4g} ripple={r['current_ripple']:.3e} "
            f"Z={r['impedance_ohm']['re']:.4f}{r['impedance_ohm']['im']:+.4f}j"
            for r in rows))
    print(f"M5 reciprocity defect {reciprocity:.3e}; loaded/network power defect "
          f"{defect:.3e}; serial/parallel delta {delta:.3e}")


if __name__=="__main__":
    main()
