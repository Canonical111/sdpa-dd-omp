#!/usr/bin/env bash
# Rebuild doc/technical.pdf from doc/technical.tex, and stamp the source hash beside it.
#
# The stamp is what makes the committed PDF checkable. `tests/check_doc_pdf.sh` compares the
# current source's SHA-256 against it, so a .tex edited without a rebuild fails CI instead of
# shipping a PDF that disagrees with its own source. Nothing else writes the stamp.
set -eu
cd "$(dirname "$0")"

command -v latexmk >/dev/null 2>&1 || {
  echo "latexmk not found (on macOS: /Library/TeX/texbin is not on PATH)" >&2; exit 2; }

latexmk -pdf -interaction=nonstopmode -halt-on-error technical.tex >/dev/null

if grep -q 'Overfull \\hbox' technical.log 2>/dev/null; then
  echo "refusing to stamp: the document has overfull boxes -- content runs off the page" >&2
  grep 'Overfull \\hbox' technical.log >&2
  latexmk -c >/dev/null 2>&1 || true
  exit 1
fi

latexmk -c >/dev/null 2>&1 || true

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum technical.tex | awk '{print $1}' > technical.pdf.stamp
else
  shasum -a 256 technical.tex | awk '{print $1}' > technical.pdf.stamp
fi

echo "technical.pdf rebuilt; stamp $(cat technical.pdf.stamp)"
