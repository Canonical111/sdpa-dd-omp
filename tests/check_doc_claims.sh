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
# RETIRED and folded into rule 8. It used to read "supported figure: 7.13x" -- and 7.13x is now
# itself the over-precise value rule 8 refuses, because nine runs of that cell span 8.2%. Two rules
# asserting different "correct" headlines is worse than one, so rule 8 is the only one.
#
# --- 3. the 1-thread baseline must not pose as a current measurement ----------------------------
# RETIRED, and it had become ACTIVELY WRONG. It fired when INSTALL's "sparse route (its default)"
# row contained "46.6 s", on the assumption that 46.6 could only be the HISTORICAL unmodified
# fork's 1-thread median (46.587 s). Since the re-campaign that row holds a genuine CURRENT
# 1-thread measurement of 46.555 s -- which rounds to 46.6. The rule would now reject a true
# statement.
#
# The invariant it was really protecting -- a table must not mix builds -- is served better and
# positively by rule 3b: the scaling table has to name where it came from.
if grep -qE '^\| dE4, sparse route \(its default\)' INSTALL.md; then
  grep -qE 'Measured on current `main`' INSTALL.md \
    || bad "INSTALL.md's scaling table does not say which build it was measured on"
fi

# --- 4. the convergence claim must name the tolerances actually tested ---------------------------
# Markup-insensitive: strip * and _ before matching, or a bolded "**any**" slips past a pattern
# that allows a single asterisk. (That exact hole was found by this script's own negative control.)
for f in README.md INSTALL.md BENCHMARKS.md; do
  if tr -d '*_' < "$f" | grep -qiE 'under any (reachable )?tolerance|any tolerance (dd|it) can reach'; then
    bad "$f claims non-convergence under ANY tolerance; only two settings were tested"
  fi
done

# --- 5. the comparison tables must be v.3 numbers, and must name the build -----------------------
# The v.2 campaign's fork rows were replaced by a v.3 rerun. Two things can go wrong: a v.2 value
# creeping back in, or a table losing the statement of which build it measured. Both are checked.
# (The earlier form of this rule guarded a "historical" label on the v.2 tables. When those tables
# were replaced the rule silently stopped guarding anything -- the --self-test below caught that,
# which is what it is for.)
for pat in '47\.08' '6\.44 s' '7\.31×' '7\.29×' '13\.5×' '2\.65×' '60\.84' '434\.62'; do
  for f in README.md BENCHMARKS.md; do
    if current_claims "$f" | grep -qE "$pat"; then
      bad "$f still carries the superseded v7.1.3-omp.2 value $pat"
    fi
  done
done
if grep -qE '^\| upstream `6eaad8d`' BENCHMARKS.md; then
  grep -qE 'both builds `v7\.1\.3-omp\.3` and' BENCHMARKS.md \
    || bad "BENCHMARKS.md's upstream tables do not name the fork build they measured"
fi
if grep -qE 'upstream, 24 threads' README.md; then
  grep -qE "Measured 2026-08-24 on \`v7\.1\.3-omp\.3\`" README.md \
    || bad "README.md's upstream table does not name the fork build and date"
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

# --- 7. main is ahead of the latest release: say so wherever a main-only figure is quoted -------
# The two-bit overlap map is on main and NOT in v7.1.3-omp.3, so the memory figures in these
# documents describe something a reader cannot download. A v7.1.3-omp.4 tag is deferred pending
# review, and until it exists the divergence must be stated rather than left for the reader to
# discover. When v.4 IS cut, this rule should be deleted along with the notices it guards.
if current_claims INSTALL.md | grep -qE '382\.0 MB|277\.5 MB'; then
  grep -qF "not the latest release" INSTALL.md \
    || bad "INSTALL.md quotes main-only memory figures without saying they are not in the release"
fi
if grep -qE 'v7\.1\.3-omp\.3' README.md; then
  grep -qE 'These documents describe .main., which is ahead of the latest release' README.md \
    || bad "README.md points at the release without noting that main is ahead of it"
fi

# --- 8. the dE4 headline must not be quoted to three significant figures ------------------------
# Nine runs across three campaigns span 8.2% on that cell, giving ratios of 7.13x, 7.24x and 7.25x.
# Every earlier revision picked one campaign's median and published it to three figures, and the
# number moved on every re-measurement. "About 7.2x" is what the data supports.
# Narrow to the HEADLINE ratio only. 7.24x also appears legitimately as the fork-versus-UPSTREAM
# advantage, which is a different quantity measured against a different denominator, so the rule
# keys on 7.13x -- the specific over-precise headline that kept recurring -- plus a positive
# requirement that the approximation is stated.
for f in README.md INSTALL.md BENCHMARKS.md; do
  if current_claims "$f" | tr -d '*_' | grep -qE '7\.13×'; then
    bad "$f quotes the dE4 headline as 7.13x; its cell spreads 8.2% across 9 runs (say ~7.2x)"
  fi
done
grep -qE 'about a \*\*7\.2×\*\* improvement' README.md \
  || bad "README.md no longer states the dE4 headline as an approximation"

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
probe "scaling table stops naming build"    INSTALL.md   's/Measured on current `main`/Measured/'
probe "non-convergence under ANY tolerance" INSTALL.md   's/either of the two tolerance settings/**any** tolerance dd can reach/'
probe "superseded v.2 value returns"        BENCHMARKS.md 's/| 46\.62 | 47\.06 |/| 47.08 | 47.29 |/'
probe "upstream table stops naming build"   BENCHMARKS.md 's/both builds `v7.1.3-omp.3` and/both builds and/'
# Anchored on "about 5×", which is the wording the rule protects, rather than on a neighbouring
# memory ratio that legitimately changes when the memory is re-measured. (The previous probe was
# tied to "2.3× less" and went dead the moment that became 2.4×.)
probe "dE3-fill quoted to 3 sig figs"       INSTALL.md   's/about 5× faster/4.92× faster/'
probe "dE4 headline to 3 sig figs"          README.md    's/about a \*\*7.2×\*\* improvement/a **7.13×** improvement/'
probe "main-ahead notice dropped (INSTALL)" INSTALL.md   's/not the latest release/not the newest thing/'
probe "main-ahead notice dropped (README)"  README.md    's/These documents describe `main`, which is ahead of the latest release/Current figures/'

echo
if [ "$dead" -eq 0 ]; then
  note ok "$probes/$probes rules demonstrated able to fail"
  exit 0
fi
echo "$dead of $probes rules could not be made to fire -- they prove nothing" >&2
exit 1
