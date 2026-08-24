#!/usr/bin/env bash
# Every relative Markdown link in the maintained documents must resolve in a fresh clone.
#
# WHY. README/INSTALL/BENCHMARKS cited five benchmark files under `review/artifacts/...` -- a
# directory that exists in the private evidence repository and NOT in this one. A reader arriving
# through GitHub could not follow any of the provenance the documents promised. INSTALL also linked
# `doc/technical.pdf`, which has never existed here. Both are the same failure: a path that is real
# on the author's machine and nowhere else.
set -u
LC_ALL=C; export LC_ALL
cd "$(dirname "$0")/.." || exit 2

docs="README.md INSTALL.md BENCHMARKS.md RUNTIME.md bench/README.md bench/probes/README.md"
missing=0
checked=0

for d in $docs; do
  [ -f "$d" ] || { echo "  FAIL  missing document: $d"; missing=$((missing+1)); continue; }
  base=$(dirname "$d")
  # [text](target) -- skip URLs, anchors and mailto
  grep -oE '\]\([^)]+\)' "$d" | sed 's/^](//; s/)$//' | while read -r t; do
    case "$t" in http*|\#*|mailto:*) continue;; esac
    printf '%s\t%s\n' "$base" "${t%%#*}"
  done
done > /tmp/_links.$$

while IFS=$'\t' read -r base target; do
  [ -n "$target" ] || continue
  checked=$((checked+1))
  if [ ! -e "$base/$target" ]; then
    echo "  FAIL  $base/$target does not exist"
    missing=$((missing+1))
  fi
done < /tmp/_links.$$
total=$(wc -l < /tmp/_links.$$ | tr -d ' ')
rm -f /tmp/_links.$$

# Also refuse any reference to the private evidence repository's layout.
if grep -rn 'review/artifacts' $docs 2>/dev/null; then
  echo "  FAIL  a document cites review/artifacts/, which is not in this repository"
  missing=$((missing+1))
fi

if [ "$missing" -eq 0 ]; then
  echo "  ok    $total relative links resolve"
  exit 0
fi
echo "$missing broken link(s)" >&2
exit 1
