#!/bin/bash
set -eu
cd -- "$(dirname -- "$0")"
python_bin="${PYTHON:-python3}"
pack_results() {
    if [ -d results ]; then
        tar -czf m13-02-results.tar.gz results
        echo "Results: $PWD/m13-02-results.tar.gz"
        echo 'Return with: scp m13-02-results.tar.gz rxtx:/tmp/'
    fi
}
trap pack_results EXIT
run_witness() {
    witness="$1"; er="$2"; sigma="$3"
    run_stage() {
        "$python_bin" openems_m10_ground.py --ground pec --cube-bottom-mm 0 \
          --loaded-epsilon-r "$er" --dielectric-sigma "$sigma" \
          --witness-frequency-mhz "$witness" --workdir "results/witness-$witness" "$@"
    }
    echo "M13-02 witness at $witness MHz: ER=$er, sigma=$sigma S/m."
    run_stage --mode pilot --pilot-epsilon-r "$er" --mesh-mm 20 --wire-cell-mm 5 --seconds 120
    run_stage --mode verify --mesh-mm 20 --wire-cell-mm 5 --seconds 600
    run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --seconds 900
    run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --domain-margin-mm 150 --seconds 1200
}
run_witness 280 3.5 .015
run_witness 320 2.5 .025
