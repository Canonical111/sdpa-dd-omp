#!/usr/bin/env bash
# Regression guard for KNOWN-STALE published claims.
#
# WHY THIS EXISTS. Four separate times in this campaign a claim was corrected in one copy and left
# standing in another: SDPA_SPCHOL_LOG's `0` semantics fixed at one call site and not the pattern;
# the source's dependency comment left contradicting the code it described; the 30-significant-
# digit claim fixed in the review documents and left in the shipped ones; the release recorded as
# "pending" after it shipped. Every one of them would have been caught by exactly this.
#
# WHAT IT IS NOT. It is not a substitute for tying each published table to its TSV metadata --
# that is what the .meta files are for. It checks a short list of claims KNOWN to have been wrong,
# in the CURRENT-CLAIM regions of the MAINTAINED documents, and nothing else.
#
# DELIBERATELY NARROW. A global grep for "7.17" or "46.6 s" becomes a nuisance the moment a
# labelled historical passage legitimately mentions those values -- and this project keeps such
# passages on purpose. So each rule below either targets an exact sentence, or excludes lines that
# sit inside an explicitly-labelled historical/erratum block.
set -u
LC_ALL=C; export LC_ALL          # the Unicode multiplication sign must compare identically on
                                 # macOS and Linux; a locale-dependent collation makes this flaky
cd "$(dirname "$0")/.." || exit 2 # resolve from the repository root, not the caller's cwd

fail=0
note() { printf '  %-6s %s\n' "$1" "$2"; }
bad()  { note FAIL "$1"; fail=$((fail+1)); }

# Lines inside a blockquote (`> `) are labelled historical/erratum context in these documents.
# Stripping them is what lets the guard coexist with deliberately-preserved wrong numbers.
current_claims() { grep -vE '^[[:space:]]*>' "$1"; }

# --- 1. the release must not describe itself as unreleased -------------------------------------
for f in README.md BENCHMARKS.md INSTALL.md; do
  if grep -qE 'release is pending|pending and not yet tagged' "$f"; then
    bad "$f still says the v7.1.3-omp.3 release is pending"
  fi
done

# --- 2. the dE4 headline ratio ------------------------------------------------------------------
# 7.17x/7.2x came from dividing the CURRENT 24-thread median into the HISTORICAL 1-thread median.
# Allowed inside an erratum or a labelled historical block; never as a current claim.
for f in README.md INSTALL.md BENCHMARKS.md; do
  if current_claims "$f" | grep -qE '7\.17|7\.2×.*(speedup|whole-solver)'; then
    bad "$f states 7.17x/7.2x as a current claim (supported figure: 7.13x)"
  fi
done

# --- 3. the 1-thread baseline must not be presented as a current-release measurement -------------
# 46.6 s is the HISTORICAL unmodified fork at one thread. It may appear, but not in a row that
# also claims to be this release's 1-thread cell.
if current_claims INSTALL.md | grep -E '46\.6 s' | grep -qE 'sparse route \(its default\)'; then
  if ! grep -qE 'historical|unmodified fork, 1 thread' INSTALL.md; then
    bad "INSTALL.md presents 46.6 s as a current-release 1-thread cell without a historical label"
  fi
fi

# --- 4. the convergence claim must name the tolerances actually tested ---------------------------
# Markup-insensitive: strip * and _ before matching, or a bolded "**any**" slips past a pattern
# that allows a single asterisk. (That exact hole was found by this script's own negative control.)
for f in README.md INSTALL.md BENCHMARKS.md; do
  if tr -d '*_' < "$f" | grep -qiE 'under any (reachable )?tolerance|any tolerance (dd|it) can reach'; then
    bad "$f claims non-convergence under ANY tolerance; only two settings were tested"
  fi
done

# --- 5. superseded pre-race-fix tables must carry a historical label -----------------------------
# The v.2 upstream/scaling tables are kept on purpose. They must announce which build they are.
if grep -qE '^\| \*\*this fork\*\* \| 47\.08' BENCHMARKS.md; then
  grep -qE 'Historical: `v7\.1\.3-omp\.2` campaign' BENCHMARKS.md \
    || bad "BENCHMARKS.md keeps the v.2 scaling table without its historical label"
fi
if grep -qE '6\.44 s — 7\.31×' README.md; then
  grep -qE 'Historical: `v7\.1\.3-omp\.2` campaign' README.md \
    || bad "README.md keeps the v.2 upstream table without its historical label"
fi

# --- 6. the jittery dE3-fill cell must not be quoted to three significant figures ---------------
# Six runs across two campaigns span 4.34-4.90 s (12%). "5.06x" and "4.92x" have both been
# published from it; neither is supportable. The published wording must stay approximate.
for f in README.md INSTALL.md BENCHMARKS.md; do
  # Bare ratio anywhere in a current claim: 5.06 and 4.92 are specific enough not to collide
  # with anything else this project publishes, and requiring a nearby "fill"/"default route"
  # missed the sentence that states the ratio with neither word on the line.
  if current_claims "$f" | tr -d '*_' | grep -qE '(5\.06|4\.92)×'; then
    bad "$f quotes the dE3-fill advantage to 3 significant figures; its cell spreads 12% (say ~5x)"
  fi
done

if [ "${1:-}" != "--self-test" ]; then
  if [ "$fail" -eq 0 ]; then
    note ok "documentation claims match the archived evidence"
    exit 0
  fi
  echo "$fail stale published claim(s); see review/DD-OPEN-ITEMS-AND-PLAN-2026-08-24.md" >&2
  exit 1
fi

# ------------------------------------------------------------------------------------------------
# --self-test: NEGATIVE CONTROL. Reintroduce each stale claim into a scratch copy of the docs and
# require the corresponding rule to fire.
#
# This exists because the first version of rule 4 was DEAD: its pattern allowed a single asterisk
# while the document used `**any**`, so the check passed on text it was written to reject. It was
# found by running exactly this. A guard whose rules are never shown to fail is the same trap the
# guard is guarding against.
# ------------------------------------------------------------------------------------------------
[ "$fail" -eq 0 ] || { echo "self-test needs a clean starting state; fix the $fail failure(s) first" >&2; exit 1; }

tmp=$(mktemp -d) || exit 2
trap 'rm -rf "$tmp"' EXIT
cp README.md INSTALL.md BENCHMARKS.md tests/check_doc_claims.sh "$tmp/" 2>/dev/null
mkdir -p "$tmp/tests" && cp tests/check_doc_claims.sh "$tmp/tests/"

probes=0
dead=0
probe() { # label file sed-expression
  probes=$((probes+1))
  cp README.md INSTALL.md BENCHMARKS.md "$tmp/" 2>/dev/null
  ( cd "$tmp" && sed -i.bak "$3" "$2" && bash tests/check_doc_claims.sh >/dev/null 2>&1 )
  if [ $? -eq 0 ]; then
    note DEAD "$1 -- rule did not fire"
    dead=$((dead+1))
  else
    note fires "$1"
  fi
}

probe "release described as pending"        README.md    's|Fixed in \[`v7.1.3-omp.3`\]|the `v7.1.3-omp.3` release is pending and [|'
probe "dE4 headline back to 7.17x"          INSTALL.md   's/7\.28×/7.17×/'
probe "46.6 s as a current 1-thread cell"   INSTALL.md   's/46\.612 s/46.6 s/'
probe "non-convergence under ANY tolerance" INSTALL.md   's/either of the two tolerance settings/**any** tolerance dd can reach/'
probe "v.2 scaling table unlabelled"        BENCHMARKS.md 's/Historical: `v7.1.3-omp.2` campaign/Current campaign/'
probe "v.2 upstream table unlabelled"       README.md    's/Historical: `v7.1.3-omp.2` campaign/Current campaign/'
probe "dE3-fill quoted to 3 sig figs"       INSTALL.md   's/— about 5× faster on 2.3× less/— 4.92× faster on 2.3× less/'

echo
if [ "$dead" -eq 0 ]; then
  note ok "$probes/$probes rules demonstrated able to fail"
  exit 0
fi
echo "$dead of $probes rules could not be made to fire -- they prove nothing" >&2
exit 1
