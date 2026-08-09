#!/usr/bin/env bash
# Assert the licensing claim PUBLISHING.md and the fork READMEs make: every
# source file this fork CHANGES relative to its upstream base carries an
# in-file, DATED modification notice naming the CORRECT licence clause.
#
# History, because it is the argument for this script existing:
#   - the claim was made universally and was untrue in all three forks;
#   - the first version of this checker was marker-only, so it passed on zero
#     arguments, on a failed `git diff`, on an undated notice, on changed files
#     whose extension it did not look at, and -- most usefully -- on five MPACK
#     files whose notices cited GPLv2 §2a while the files are LGPL-3-only.
# Every one of those is now a failure, and each has a self-test below.
#
#   bash patches/check_notices.sh <clone-dir> [<clone-dir> ...]
#   bash patches/check_notices.sh --self-test
#
# Needs FULL history (it derives the upstream base), so a shallow clone fails
# loudly rather than passing vacuously. CI runs it after `git fetch --unshallow`.
set -uo pipefail

# Resolve a REAL grep. This environment exports a shell function that shadows
# grep with ugrep, whose regex handling differs enough to make this checker
# report the opposite answer -- found while writing the self-tests below. A
# verification tool must not depend on the caller's shell.
GREP=$(command -p -v grep 2>/dev/null || true)
[ -x "${GREP:-}" ] || GREP=/usr/bin/grep
[ -x "$GREP" ] || { echo "FAIL: no grep binary found" >&2; exit 2; }

MARKER='MODIFIED from upstream|MODIFICATION NOTICE|NEW FILE'
# Any source-ish file that can carry a notice. The old .cpp/.h-only filter meant
# a changed .c, .hpp, .cc or .inc silently needed nothing.
SRC_RE='\.(c|cc|cpp|cxx|h|hh|hpp|hxx)$'
# patch fragments and vendored third-party trees are not this fork's source
SKIP_RE='(^|/)(patches|download|work)/'
DATE_RE='[0-9]{4}-[0-9]{2}-[0-9]{2}'

check_clone() {
    local clone=$1 name rc=0
    name=$(basename "$clone")

    local first base
    first=$(git -C "$clone" log --format='%H %an' 2>/dev/null | "$GREP" Zechuan | tail -1 | cut -d' ' -f1)
    if [ -z "$first" ]; then
        echo "FAIL: $name -- no fork commits found (not a clone of this fork, or shallow)" >&2
        return 1
    fi
    if ! base=$(git -C "$clone" rev-parse "$first^" 2>/dev/null) || [ -z "$base" ]; then
        echo "FAIL: $name -- upstream base unreachable (shallow clone: use a full clone or a bundle)" >&2
        return 1
    fi

    # The diff must SUCCEED. Previously its failure produced an empty list, and
    # an empty list read as "nothing to check" -- the empty-equals-clean trap.
    local changed
    if ! changed=$(git -C "$clone" diff --name-only "$base..HEAD" 2>/dev/null); then
        echo "FAIL: $name -- git diff against the upstream base failed; refusing to report a pass" >&2
        return 1
    fi
    local files
    files=$(printf '%s\n' "$changed" | "$GREP" -E "$SRC_RE" | "$GREP" -vE "$SKIP_RE" || true)
    if [ -z "$files" ]; then
        echo "FAIL: $name -- no changed source files found at all; a fork that changes nothing is not what this checks" >&2
        return 1
    fi

    local n=0
    while read -r p; do
        [ -n "$p" ] || continue
        [ -f "$clone/$p" ] || continue          # deleted upstream files carry nothing
        n=$((n + 1))
        local f="$clone/$p"
        if ! "$GREP" -qE "($MARKER)" "$f"; then
            echo "FAIL: $name/$p is changed from upstream but carries no notice" >&2
            rc=1; continue
        fi
        # the notice must be DATED: "in-file, dated" is the documented policy
        # NOT a pipeline: `grep -q` exits at the first match, the producer gets
        # SIGPIPE, and `set -o pipefail` turns that into a failed test -- which
        # made this report "no ISO date" for a correctly dated file whose notice
        # happens to be followed by a very long line. Capture, then test.
        local ctx
        ctx=$("$GREP" -E "($MARKER)" -A6 "$f")
        if ! printf '%s' "$ctx" | "$GREP" -qE "$DATE_RE"; then
            echo "FAIL: $name/$p has a notice with no ISO date" >&2
            rc=1
        fi
        # and must name the licence the FILE is actually under
        local want got
        if "$GREP" -q "Lesser General Public License" "$f"; then want=LGPL; else want=GPL; fi
        if "$GREP" -q "LGPL-3" "$f"; then got=LGPL; else got=GPL; fi
        if [ "$want" = LGPL ] && [ "$got" != LGPL ]; then
            echo "FAIL: $name/$p is LGPL-3 (its own header says so) but its notice cites GPLv2" >&2
            rc=1
        fi
        if [ "$want" = GPL ] && "$GREP" -qE "($MARKER).*LGPL-3" "$f"; then
            echo "FAIL: $name/$p is GPLv2 but its notice cites LGPL-3" >&2
            rc=1
        fi
    done <<< "$files"

    [ $rc -eq 0 ] && echo "ok   $name -- all $n changed source files carry a dated, correctly-licensed notice"
    return $rc
}

self_test() {
    local t rc=0
    t=$(mktemp -d)
    # a clone with one changed file, varied to break each rule in turn
    git init -q "$t/r" && cd "$t/r"
    git config user.email a@b; git config user.name Upstream
    printf 'GNU General Public License\nint x;\n' > a.cpp
    git add a.cpp && git commit -qm upstream
    git config user.name "Zechuan Zheng"
    printf 'GNU General Public License\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\nint x=1;\n' > a.cpp
    git add a.cpp && git commit -qm fork
    cd - >/dev/null
    check_clone "$t/r" >/dev/null 2>&1 || { echo "SELF-TEST FAIL: a well-formed clone did not pass" >&2; rc=1; }
    # undated
    printf 'GNU General Public License\n/* MODIFIED from upstream (GPLv2 2a notice): x. */\nint x=1;\n' > "$t/r/a.cpp"
    check_clone "$t/r" >/dev/null 2>&1 && { echo "SELF-TEST FAIL: an undated notice passed" >&2; rc=1; }
    # wrong licence label
    printf 'GNU Lesser General Public License\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\nint x=1;\n' > "$t/r/a.cpp"
    check_clone "$t/r" >/dev/null 2>&1 && { echo "SELF-TEST FAIL: an LGPL file citing GPLv2 passed" >&2; rc=1; }
    # missing notice
    printf 'GNU General Public License\nint x=1;\n' > "$t/r/a.cpp"
    check_clone "$t/r" >/dev/null 2>&1 && { echo "SELF-TEST FAIL: a missing notice passed" >&2; rc=1; }
    # not a repository
    check_clone "$t" >/dev/null 2>&1 && { echo "SELF-TEST FAIL: a non-repository passed" >&2; rc=1; }
    rm -rf "$t"
    [ $rc -eq 0 ] && echo "ok   self-test: well-formed passes; undated, mislabelled, missing and non-repository all fail"
    return $rc
}

if [ "${1:-}" = --self-test ]; then self_test; exit $?; fi

# Zero arguments used to exit 0 -- a check that ran on nothing and reported success.
[ $# -ge 1 ] || { echo "FAIL: no clone directories given; refusing to report a pass on nothing" >&2
                  echo "usage: check_notices.sh <clone-dir> [<clone-dir> ...] | --self-test" >&2; exit 2; }

rc=0
for clone in "$@"; do check_clone "$clone" || rc=1; done
exit $rc
