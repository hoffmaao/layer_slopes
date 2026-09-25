#!/usr/bin/env bash
# Run the regression suite with whichever MATLAB this machine has.
#
#   bash tests/run_tests.sh              # finds MATLAB itself
#   MATLAB=/path/to/bin/matlab bash tests/run_tests.sh
#
# MATLAB is rarely on PATH, so after PATH this looks where it is usually
# installed: /Applications on a Mac, /opt/sw/matlab on the CReSIS servers,
# /usr/local/MATLAB elsewhere. The newest release found wins. Exits with the
# suite's status: 0 when every test passes.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
m="${MATLAB:-}"
if [ -z "$m" ]; then
    if command -v matlab >/dev/null 2>&1; then
        m="$(command -v matlab)"
    else
        for c in /Applications/MATLAB_R20*.app/bin/matlab \
                 /opt/sw/matlab/*/bin/matlab /usr/local/MATLAB/*/bin/matlab; do
            if [ -x "$c" ]; then m="$c"; fi      # globs sort ascending: newest last
        done
    fi
fi
if [ -z "$m" ]; then
    echo "run_tests.sh: no MATLAB found; set MATLAB=/path/to/bin/matlab" >&2
    exit 2
fi
echo "run_tests.sh: $m"
cd "$here"
exec "$m" -batch "run_tests"
