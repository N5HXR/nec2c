#!/bin/bash
# Run this script in the existing laptop OpenEMS Python environment.
set -eu
cd -- "$(dirname -- "$0")"
python_bin="${PYTHON:-python3}"
pack_results() {
    if [ -d results ]; then
        tar -czf m13-01-results.tar.gz results
        echo "Results: $PWD/m13-01-results.tar.gz"
        echo 'Return with: scp m13-01-results.tar.gz rxtx:/tmp/'
    fi
}
trap pack_results EXIT
run_stage() {
    "$python_bin" openems_m10_ground.py --ground pec --cube-bottom-mm 0 --dielectric-sigma .02 --workdir results "$@"
}
echo 'M13-01: isotropic epsilon_r=3, conductivity=0.02 S/m, PEC-contact pillar.'
run_stage --mode pilot --pilot-epsilon-r 3 --mesh-mm 20 --wire-cell-mm 5 --seconds 120
run_stage --mode verify --mesh-mm 20 --wire-cell-mm 5 --seconds 600
run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --seconds 900
# Independent domain enlargement at the fine mesh; same physical geometry.
run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --domain-margin-mm 150 --seconds 1200
