#!/bin/bash
set -eu
cd -- "$(dirname -- "$0")"
python_bin="${PYTHON:-python3}"
pack_results() {
    if [ -d results ]; then
        tar -czf m13-03-results.tar.gz results
        echo "Results: $PWD/m13-03-results.tar.gz"
        echo 'Return with: scp m13-03-results.tar.gz rxtx:/tmp/'
    fi
}
trap pack_results EXIT
for split in x z; do
    run_stage() {
        "$python_bin" openems_m10_ground.py --ground pec --cube-bottom-mm 0 \
          --material-split "$split" --loaded-epsilon-r 2 --dielectric-sigma .01 \
          --workdir "results/split-$split" "$@"
    }
    echo "M13-03 split-$split: ER=2/sigma=.01 and ER=4/sigma=.03 S/m."
    run_stage --mode pilot --pilot-epsilon-r 2 --mesh-mm 20 --wire-cell-mm 5 --seconds 180
    run_stage --mode verify --mesh-mm 20 --wire-cell-mm 5 --seconds 900
    run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --seconds 1200
    run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --domain-margin-mm 150 --seconds 1500
done
