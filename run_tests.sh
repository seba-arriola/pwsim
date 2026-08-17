#!/usr/bin/env bash
# Validation script for pwsim.
#
# Usage:
#   ./run_tests.sh            quick validation suite (fast, reduced inputs)
#   ./run_tests.sh full       quick suite + full params under ASan (slower)
#   ./run_tests.sh valgrind   quick suite + valgrind on reduced inputs
#
# The reduced inputs exercise every code path (applyTF=0/1, stat 0/1/2,
# calcfs 1/2, radpat 0/1) while keeping runtimes short enough for
# sanitizers and valgrind.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

MODE="${1:-quick}"
FAIL=0

run_asan() {
    local params="$1"
    rm -rf output
    echo "  [ASan] $params"
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ./pwsim_san "$params" >/dev/null 2>/tmp/opencode/pwsim_asan.log
    if [ $? -ne 0 ]; then
        echo "    FAILED (see /tmp/opencode/pwsim_asan.log)"; cat /tmp/opencode/pwsim_asan.log; FAIL=1; return
    fi
    if grep -qE "ERROR: AddressSanitizer|runtime error:|LeakSanitizer" /tmp/opencode/pwsim_asan.log; then
        echo "    FAILED (sanitizer report)"; cat /tmp/opencode/pwsim_asan.log; FAIL=1; return
    fi
    # basic sanity on the output
    if ! ls output/synt_*.dat >/dev/null 2>&1; then
        echo "    FAILED: no output files"; FAIL=1; return
    fi
    if grep -lE "nan|inf" output/synt_*.dat >/dev/null 2>&1; then
        echo "    FAILED: NaN/Inf in output"; FAIL=1; return
    fi
    echo "    ok"
}

run_valgrind() {
    local params="$1"
    rm -rf output
    echo "  [valgrind] $params"
    local out
    out=$(valgrind --leak-check=full --track-origins=yes --error-exitcode=99 ./pwsim "$params" 2>&1)
    if [ $? -ne 0 ]; then
        echo "    FAILED"; echo "$out" | tail -30; FAIL=1; return
    fi
    local errs
    errs=$(echo "$out" | grep -E "ERROR SUMMARY: [1-9]|definitely lost: [1-9]|indirectly lost: [1-9]|possibly lost: [1-9]")
    if [ -n "$errs" ]; then
        echo "    FAILED"; echo "$errs"; FAIL=1; return
    fi
    echo "    ok"
}

echo "==> Building normal executable"
make clean >/dev/null 2>&1
make pwsim >/dev/null 2>&1

echo "==> Building AddressSanitizer + UBSan executable"
make sanitize >/dev/null 2>&1

echo "==> Quick validation suite (ASan)"
for t in test_base test_satf test_calcfs2 test_radpat0; do
    run_asan "test/$t.dat"
done

if [ "$MODE" = "full" ] || [ "$MODE" = "all" ]; then
    echo "==> Full params under ASan"
    run_asan "params/params_nsf_full.dat"
fi

if [ "$MODE" = "valgrind" ] || [ "$MODE" = "all" ]; then
    echo "==> valgrind on reduced inputs"
    for t in test_base test_satf test_calcfs2 test_radpat0; do
        run_valgrind "test/$t.dat"
    done
fi

echo
if [ "$FAIL" -eq 0 ]; then
    echo "==> All validation steps passed"
else
    echo "==> Validation FAILED"
    exit 1
fi
