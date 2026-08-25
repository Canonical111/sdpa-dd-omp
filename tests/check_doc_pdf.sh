#!/usr/bin/env bash
# doc/technical.pdf is a COMMITTED BUILD PRODUCT of doc/technical.tex. That is deliberate -- the
# document has to be readable in a browser, without a TeX installation -- but it recreates the
# exact failure mode this repository has hit four times: a claim corrected in one copy and left
# standing in another. Here the two copies are the source and the PDF a reader actually opens.
#
# The link between them is doc/technical.pdf.stamp, which doc/build.sh writes and nothing else
# does: it holds the SHA-256 of the .tex the committed PDF was built from. Comparing hashes needs
# no TeX and no git history, which matters because CI has neither (actions/checkout is shallow, so
# a "which file was committed later" test would pass vacuously there -- that was this guard's first
# design, and it would have been dead in the only place it runs).
set -u
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.." || exit 2

tex=doc/technical.tex
pdf=doc/technical.pdf
stamp=doc/technical.pdf.stamp
fail=0
note() { printf '  %-6s %s\n' "$1" "$2"; }

for f in "$tex" "$pdf" "$stamp"; do
  [ -f "$f" ] || { note FAIL "missing $f (run doc/build.sh)"; exit 1; }
done

# --- 1. the committed PDF must have been built from the committed source ------------------------
if command -v sha256sum >/dev/null 2>&1; then
  now=$(sha256sum "$tex" | awk '{print $1}')
else
  now=$(shasum -a 256 "$tex" | awk '{print $1}')
fi
want=$(tr -d '[:space:]' < "$stamp")

if [ "$now" = "$want" ]; then
  note ok "$pdf was built from the current $tex"
else
  note FAIL "$tex has changed since $pdf was built -- run doc/build.sh and commit both"
  note ""   "  source now $now"
  note ""   "  stamp says $want"
  fail=$((fail+1))
fi

# --- 2. if a TeX toolchain is present, the source must still compile cleanly ---------------------
if command -v latexmk >/dev/null 2>&1; then
  tmp=$(mktemp -d) || exit 2
  trap 'rm -rf "$tmp"' EXIT
  cp "$tex" "$tmp/" || exit 2
  if ( cd "$tmp" && latexmk -pdf -interaction=nonstopmode -halt-on-error technical.tex ) >"$tmp/out" 2>&1; then
    note ok "$tex compiles"
    # An overfull box is a real defect in a document whose tables carry the reference material:
    # a column that runs off the page loses text a reader needs.
    if grep -q 'Overfull \\hbox' "$tmp/technical.log" 2>/dev/null; then
      note FAIL "$tex compiles but produces overfull boxes -- content runs off the page"
      grep 'Overfull \\hbox' "$tmp/technical.log" | sed 's/^/         /'
      fail=$((fail+1))
    fi
  else
    note FAIL "$tex no longer compiles"
    tail -30 "$tmp/out" | sed 's/^/         /'
    fail=$((fail+1))
  fi
else
  note skip "no latexmk here; ran the stamp check only"
fi

[ "$fail" -eq 0 ] && exit 0
echo "$fail problem(s) with the published technical document" >&2
exit 1
