#!/bin/bash
set -eu
cd -- "$(dirname -- "$0")"
python_bin="${PYTHON:-python3}"
pack_results() {
    if [ -d results ]; then
        tar -czf m13-04-results.tar.gz results
        echo "Results: $PWD/m13-04-results.tar.gz"
        echo 'Return with: scp m13-04-results.tar.gz rxtx:/tmp/'
    fi
}
trap pack_results EXIT
for loss in lossless lossy; do
    sigma=0
    if [ "$loss" = lossy ]; then sigma=.02; fi
    run_stage() {
        "$python_bin" openems_m10_ground.py --ground real --cube-bottom-mm 550 \
          --loaded-epsilon-r 3 --dielectric-sigma "$sigma" \
          --workdir "results/$loss" "$@"
    }
    echo "M13-04 $loss cube: ER=3, sigma=$sigma S/m; soil ER=13, sigma=.02 S/m."
    run_stage --mode pilot --pilot-epsilon-r 3 --mesh-mm 20 --wire-cell-mm 5 --seconds 240
    run_stage --mode verify --mesh-mm 20 --wire-cell-mm 5 --seconds 1200
    run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --seconds 1800
    run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --domain-margin-mm 150 --seconds 2400
done
