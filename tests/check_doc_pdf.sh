#!/usr/bin/env bash
# doc/technical.pdf is a COMMITTED build product of doc/technical.tex, so it is a second copy of
# every claim in the source -- the same two-copies-of-one-claim shape this repository has been
# bitten by four times, except that here the second copy is the one a reader actually opens.
#
# WHAT THIS ESTABLISHES, precisely. doc/build.sh writes the SHA-256 of technical.tex into the PDF's
# own Info dictionary (\pdfinfo, raw ASCII, greppable). Rule 1 recomputes that hash and looks for it
# INSIDE THE COMMITTED PDF. So the check reads the artefact it is making a claim about, and a source
# edited without a rebuild fails.
#
# WHAT IT DOES NOT ESTABLISH. It detects STALENESS, not tampering: someone who edited the PDF body
# while leaving the Info dictionary alone would still pass. Defending against that would mean
# rebuilding and comparing rendered output, which needs a TeX installation CI does not have. The
# message below says "records the hash of", not "is a faithful rendering of", for that reason.
#
# HISTORY, because the first version was weaker than it read. It compared the source against a
# SIDECAR stamp file and never touched the PDF at all -- so a stale PDF next to a freshly written
# hash passed both halves while the guard printed "was built from". An independent review caught it.
# The version before THAT compared which file git committed later, which needs no TeX and would
# have been vacuous in CI, because actions/checkout is shallow and `git log -1 -- <path>` returns
# nothing there. Two designs that could not mean what they said, before one that can.
set -u
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.." || exit 2

tex=doc/technical.tex
pdf=doc/technical.pdf
fail=0
note() { printf '  %-6s %s\n' "$1" "$2"; }

for f in "$tex" "$pdf"; do
  [ -f "$f" ] || { note FAIL "missing $f (run doc/build.sh)"; exit 1; }
done

# --- 1. the committed PDF must record the hash of the committed source --------------------------
if command -v sha256sum >/dev/null 2>&1; then
  now=$(sha256sum "$tex" | awk '{print $1}')
else
  now=$(shasum -a 256 "$tex" | awk '{print $1}')
fi

if grep -a -q "/SourceSHA256 ($now)" "$pdf"; then
  note ok "$pdf records the SHA-256 of the current $tex"
else
  got=$(grep -a -o '/SourceSHA256 ([0-9a-f]*)' "$pdf" | head -1)
  note FAIL "$tex has changed since $pdf was built -- run doc/build.sh and commit both"
  note ""   "  source now  $now"
  note ""   "  PDF records ${got:-<no /SourceSHA256 field at all>}"
  fail=$((fail+1))
fi

# --- 2. if a TeX toolchain is present, the source must still compile cleanly ---------------------
if command -v latexmk >/dev/null 2>&1; then
  tmp=$(mktemp -d) || exit 2
  trap 'rm -rf "$tmp"' EXIT
  cp "$tex" "$tmp/" || exit 2
  # srchash.tex is deliberately NOT copied: the source has to build without it (\IfFileExists), so
  # a fresh clone can compile the document with no build script.
  if ( cd "$tmp" && latexmk -pdf -interaction=nonstopmode -halt-on-error technical.tex ) >"$tmp/out" 2>&1; then
    note ok "$tex compiles"
    # Both directions. The earlier version checked \hbox only while its message said "overfull
    # boxes", which is the same overclaim as the sidecar stamp, in miniature.
    for kind in hbox vbox; do
      if grep -q "Overfull \\\\$kind" "$tmp/technical.log" 2>/dev/null; then
        note FAIL "$tex compiles but produces overfull ${kind}es -- content runs off the page"
        grep "Overfull \\\\$kind" "$tmp/technical.log" | sed 's/^/         /'
        fail=$((fail+1))
      fi
    done
  else
    note FAIL "$tex no longer compiles"
    tail -30 "$tmp/out" | sed 's/^/         /'
    fail=$((fail+1))
  fi
else
  note skip "no latexmk here; ran the source-binding check only"
fi

[ "$fail" -eq 0 ] && exit 0
echo "$fail problem(s) with the published technical document" >&2
exit 1
