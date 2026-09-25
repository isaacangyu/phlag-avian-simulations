#!/usr/bin/env bash
# Run `simulate <branch>` for every branch at once, one core per branch.
#
# Usage:
#   scripts/script.sh                              # run the BRANCHES list below
#   scripts/script.sh 10X/up/N554 10X/down/N109    # run the branches given
#
# Needs the phlag-avian-simulations env active (so `simulate` is on PATH).
# Output goes to the terminal, each line prefixed with its branch.

set -euo pipefail

BRANCHES=(
    10X/up/N118
    10X/up/N159
    10X/up/N220
    10X/up/N228
    10X/up/N233
)

# Extra flags passed to every `simulate` call, e.g. (-g --gt-threads 1).
# Keep --gt-threads at 1 so each branch uses one core.
EXTRA_ARGS=()

if [ "$#" -gt 0 ]; then
    BRANCHES=("$@")
fi

run_branch() {
    local branch="$1"
    echo "[start] $branch"
    if simulate "$branch" "${EXTRA_ARGS[@]}" 2>&1 | sed -u "s#^#[$branch] #"; then
        echo "[done]  $branch"
    else
        echo "[FAIL]  $branch"
        return 1
    fi
}

echo "Running ${#BRANCHES[@]} branch(es) on ${#BRANCHES[@]} core(s)"

for branch in "${BRANCHES[@]}"; do
    run_branch "$branch" &
done

failed=0
for _ in "${BRANCHES[@]}"; do
    wait -n || failed=1
done

exit "$failed"
