#!/usr/bin/env bash
# Rebuild doc/technical.pdf from doc/technical.tex, binding the two together.
#
# The binding is the point. technical.pdf is a COMMITTED build product, so it is a second copy of
# every claim in the source -- the exact shape this repository has been bitten by four times. So
# the SHA-256 of the source is written into the PDF itself (Info dictionary, via \pdfinfo, which
# pdftex emits as raw ASCII), and tests/check_doc_pdf.sh verifies it with nothing but grep.
#
# An earlier design wrote the hash to a SIDECAR file instead. That was weaker than it read: the
# guard compared the source against the sidecar and never touched the PDF at all, so a stale PDF
# beside a freshly written hash passed both halves of the check while the guard printed "built
# from". Caught in review. The hash now lives where it is actually load-bearing.
set -eu
cd "$(dirname "$0")"

command -v latexmk >/dev/null 2>&1 || {
  echo "latexmk not found (on macOS: /Library/TeX/texbin is not on PATH)" >&2; exit 2; }

sha() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
  else shasum -a 256 "$1" | awk '{print $1}'; fi
}

# The hash covers technical.tex ONLY. srchash.tex is generated and \input, so the source's own text
# never contains the hash -- no fixed point to chase.
H=$(sha technical.tex)
printf '\\newcommand{\\SourceHash}{%s}\n' "$H" > srchash.tex

latexmk -pdf -interaction=nonstopmode -halt-on-error technical.tex >/dev/null

for kind in hbox vbox; do
  if grep -q "Overfull \\\\$kind" technical.log 2>/dev/null; then
    echo "refusing to publish: the document has overfull ${kind}es -- content runs off the page" >&2
    grep "Overfull \\\\$kind" technical.log >&2
    latexmk -c >/dev/null 2>&1 || true
    exit 1
  fi
done

latexmk -c >/dev/null 2>&1 || true

# Refuse to leave behind a PDF that does not carry the hash: if \pdfinfo silently did nothing --
# a different engine, a future hyperref -- the guard downstream would fail with a confusing
# message. Fail here instead, where the cause is visible.
grep -a -q "/SourceSHA256 ($H)" technical.pdf || {
  echo "refusing to publish: technical.pdf does not carry /SourceSHA256 ($H)." >&2
  echo "The \\pdfinfo stamp did not reach the PDF -- check the engine is pdflatex." >&2
  exit 1; }

echo "technical.pdf rebuilt and stamped with the source hash $H"
