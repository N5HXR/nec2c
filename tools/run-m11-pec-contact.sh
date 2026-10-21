#!/bin/bash
# Run in the existing OpenEMS Python environment on the laptop.
set -u
cd -- "$(dirname -- "$0")" || exit 1
python_bin="${PYTHON:-python3}"
pack_results() {
    if [ -d results ]; then
        tar -czf m11-pec-contact-results.tar.gz results || return 1
        echo "Results packaged: $PWD/m11-pec-contact-results.tar.gz"
        echo 'Return them with: scp m11-pec-contact-results.tar.gz rxtx:/tmp/'
    fi
}
run_stage() {
    "$python_bin" openems_m10_ground.py --ground pec --cube-bottom-mm 0 --workdir results "$@"
    status=$?
    if [ "$status" -ne 0 ]; then
        pack_results
        echo 'Stage failed or exhausted its budget. Return the archive for review, or rerun this script to reuse completed configurations.' >&2
        exit "$status"
    fi
}
echo 'M11 grounded dielectric pillar: loaded pilot (60-second cap).'
run_stage --mode pilot --pilot-epsilon-r 3 --mesh-mm 20 --wire-cell-mm 5 --seconds 60
echo 'Coarse vacuum/loaded comparison (180-second pair cap).'
run_stage --mode verify --mesh-mm 20 --wire-cell-mm 5 --seconds 180
echo 'Fine vacuum/loaded comparison (360-second pair cap).'
run_stage --mode verify --mesh-mm 15 --wire-cell-mm 2.5 --seconds 360
pack_results
