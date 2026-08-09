#!/usr/bin/env bash
# Assert the licensing claim PUBLISHING.md and the fork READMEs make: every
# source file this fork CHANGES relative to its upstream base carries an
# in-file, DATED modification notice naming the licence clause of the licence
# THAT FILE is under.
#
#   bash check_notices.sh <clone-dir> [<clone-dir> ...]
#   bash check_notices.sh --self-test
#
# Needs FULL history (it derives the upstream base), so a shallow clone fails
# loudly rather than passing vacuously. CI runs it after `git fetch --unshallow`.
#
# WHY THIS IS PYTHON INSIDE A .sh: the logic needs the notice's own COMMENT
# BLOCK, not "the marker plus six lines" -- an undated notice followed by an
# unrelated date passed the line-window version, and a licence token anywhere
# in the file satisfied the licence test. Both were reproduced. Two earlier
# shell versions also failed for reasons unrelated to the policy (an
# unparenthesised alternation; a shell function shadowing grep with ugrep; and
# SIGPIPE under pipefail when grep -q exited early). The name and CLI are kept
# because CI and the release gate reference them.
set -uo pipefail
exec python3 - "$@" <<'PYEOF'
import os
import re
import subprocess
import sys
import tempfile

MARKER = re.compile(r"MODIFIED from upstream|MODIFICATION NOTICE|NEW FILE")
DATE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")
# The exact extensions this policy covers. Named explicitly rather than
# described as "every source extension", which was never true (review2 §24.2).
SRC_EXT = (".c", ".cc", ".cpp", ".cxx", ".inc", ".inl", ".ipp", ".tpp",
           ".h", ".hh", ".hpp", ".hxx")
# Vendored third-party trees and patch fragments are not this fork's source.
SKIP_DIR = re.compile(r"(^|/)(patches|download|work)/")

# Licence classification. Three families actually occur in these trees:
#   GPL   -- the SDPA sources themselves
#   LGPL  -- qd's bundled mpack
#   BSD   -- dd's and gmp's bundled mplapack (2-clause, no GNU licence at all)
# An unrecognised header is a FAILURE, not a default: treating "not Lesser" as
# GPL is what let 13 BSD-headered files carry GPLv2-labelled notices and be
# reported correctly-licensed (review2 §24.2).
CLAUSE = {"LGPL": re.compile(r"LGPL-3"),
          "GPL": re.compile(r"GPLv2"),
          "BSD": re.compile(r"BSD")}


def classify(text):
    if "Lesser General Public License" in text:
        return "LGPL"
    if "GNU General Public License" in text:
        return "GPL"
    if "Redistribution and use in source and binary" in text:
        return "BSD"
    return None


def notice_blocks(text):
    """The comment block(s) containing a marker -- the /* ... */ the marker sits
    in, or the contiguous run of // lines around it. Checks bind to THIS, so a
    date or licence token elsewhere in the file cannot satisfy them."""
    out = []
    for m in MARKER.finditer(text):
        i = m.start()
        open_c = text.rfind("/*", 0, i)
        close_c = text.find("*/", i)
        if open_c != -1 and (close_c != -1) and text.find("*/", open_c, i) == -1:
            out.append(text[open_c:close_c + 2])
            continue
        # line comment: walk the contiguous // run
        ls = text[:i].count("\n")
        lines = text.split("\n")
        a = b = ls
        while a > 0 and lines[a - 1].lstrip().startswith("//"):
            a -= 1
        while b + 1 < len(lines) and lines[b + 1].lstrip().startswith("//"):
            b += 1
        out.append("\n".join(lines[a:b + 1]))
    return out


def check_clone(clone):
    name = os.path.basename(os.path.normpath(clone))
    bad = []

    def git(*a):
        return subprocess.run(("git", "-C", clone) + a, capture_output=True, text=True)

    log = git("log", "--format=%H %an")
    if log.returncode != 0:
        print("FAIL: %s -- not a git repository" % name, file=sys.stderr)
        return 1
    forks = [l.split(" ", 1)[0] for l in log.stdout.splitlines() if "Zechuan" in l]
    if not forks:
        print("FAIL: %s -- no fork commits found (not a clone of this fork, or shallow)" % name, file=sys.stderr)
        return 1
    base = git("rev-parse", forks[-1] + "^")
    if base.returncode != 0 or not base.stdout.strip():
        print("FAIL: %s -- upstream base unreachable (shallow clone: use a full clone or a bundle)" % name, file=sys.stderr)
        return 1
    base = base.stdout.strip()

    diff = git("diff", "--name-only", base + "..HEAD")
    if diff.returncode != 0:
        print("FAIL: %s -- git diff against the upstream base failed; refusing to report a pass" % name, file=sys.stderr)
        return 1
    files = [p for p in diff.stdout.splitlines()
             if p.endswith(SRC_EXT) and not SKIP_DIR.search(p)]
    if not files:
        print("FAIL: %s -- no changed source files found at all; a fork that changes nothing is not what this checks" % name, file=sys.stderr)
        return 1

    n = 0
    for p in files:
        f = os.path.join(clone, p)
        if not os.path.isfile(f):
            continue                       # deleted upstream files carry nothing
        n += 1
        try:
            text = open(f, encoding="utf-8", errors="replace").read()
        except OSError as e:
            bad.append("%s/%s unreadable (%s)" % (name, p, e)); continue

        blocks = notice_blocks(text)
        if not blocks:
            bad.append("%s/%s is changed from upstream but carries no notice" % (name, p)); continue
        if not any(DATE.search(b) for b in blocks):
            bad.append("%s/%s has a notice with no ISO date INSIDE the notice block" % (name, p))

        want = classify(text)
        if want is None:
            bad.append("%s/%s has an unrecognised licence header -- classify it rather than "
                       "assuming a default" % (name, p)); continue
        if not any(CLAUSE[want].search(b) for b in blocks):
            named = [k for k in CLAUSE if any(CLAUSE[k].search(b) for b in blocks)]
            bad.append("%s/%s is %s (its own header says so) but its notice names %s"
                       % (name, p, want, "/".join(named) if named else "no licence"))

    for b in bad:
        print("FAIL: " + b, file=sys.stderr)
    if not bad:
        print("ok   %s -- all %d changed source files carry a dated notice naming their own licence" % (name, n))
    return 1 if bad else 0


# ------------------------------------------------------------------ self-tests
HEADERS = {"GPL": "GNU General Public License",
           "LGPL": "GNU Lesser General Public License",
           "BSD": "Redistribution and use in source and binary forms"}


def _mkclone(d, body, fname="a.cpp"):
    subprocess.run(["git", "init", "-q", d], check=True)
    for k, v in (("user.email", "a@b"), ("user.name", "Upstream")):
        subprocess.run(["git", "-C", d, "config", k, v], check=True)
    open(os.path.join(d, fname), "w").write("/* upstream */\nint x;\n")
    subprocess.run(["git", "-C", d, "add", fname], check=True)
    subprocess.run(["git", "-C", d, "commit", "-qm", "up"], check=True)
    subprocess.run(["git", "-C", d, "config", "user.name", "Zechuan Zheng"], check=True)
    open(os.path.join(d, fname), "w").write(body)
    subprocess.run(["git", "-C", d, "add", fname], check=True)
    subprocess.run(["git", "-C", d, "commit", "-qm", "fork"], check=True)


def self_test():
    ok_body = "/* %s */\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\nint x=1;\n" % HEADERS["GPL"]
    cases = [
        # (name, body, filename, expect_pass)
        ("well-formed GPL", ok_body, "a.cpp", True),
        ("well-formed LGPL", "/* %s */\n/* MODIFICATION NOTICE (LGPL-3 4a), 2026-08-09: x. */\nint x=1;\n" % HEADERS["LGPL"], "a.cpp", True),
        ("well-formed BSD", "/* %s */\n/* MODIFIED from upstream (BSD 2-clause), 2026-08-09: x. */\nint x=1;\n" % HEADERS["BSD"], "a.cpp", True),
        ("undated notice", "/* %s */\n/* MODIFIED from upstream (GPLv2 2a notice): x. */\nint x=1;\n" % HEADERS["GPL"], "a.cpp", False),
        ("undated notice, unrelated date nearby",
         "/* %s */\n/* MODIFIED from upstream (GPLv2 2a notice): x. */\n// unrelated 2026-01-01\nint x=1;\n" % HEADERS["GPL"], "a.cpp", False),
        ("LGPL file, GPLv2 notice", "/* %s */\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\nint x=1;\n" % HEADERS["LGPL"], "a.cpp", False),
        ("LGPL file, GPLv2 notice, stray LGPL-3 elsewhere",
         "/* %s */\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\n// see LGPL-3 elsewhere\nint x=1;\n" % HEADERS["LGPL"], "a.cpp", False),
        ("GPL file, LGPL-3 notice (reverse)", "/* %s */\n/* MODIFIED from upstream\n   (LGPL-3 4a), 2026-08-09: x. */\nint x=1;\n" % HEADERS["GPL"], "a.cpp", False),
        ("BSD file, GPLv2 notice", "/* %s */\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\nint x=1;\n" % HEADERS["BSD"], "a.cpp", False),
        ("unknown licence header", "/* Proprietary, all rights reserved */\n/* MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x. */\nint x=1;\n", "a.cpp", False),
        ("missing notice", "/* %s */\nint x=1;\n" % HEADERS["GPL"], "a.cpp", False),
        ("alternate extension .inc is covered", "/* %s */\nint x=1;\n" % HEADERS["GPL"], "a.inc", False),
        ("line-comment notice", "// %s\n// MODIFIED from upstream (GPLv2 2a notice), 2026-08-09: x.\nint x=1;\n" % HEADERS["GPL"], "a.cpp", True),
        # the SIGPIPE regression: a correct notice followed by a very long line
        ("dated notice followed by a 40 KB line", ok_body + "// " + ("z" * 40000) + "\n", "a.cpp", True),
    ]
    rc = 0
    for label, body, fname, expect in cases:
        d = tempfile.mkdtemp()
        try:
            _mkclone(os.path.join(d, "r"), body, fname)
            got = check_clone(os.path.join(d, "r")) == 0
            if got != expect:
                print("SELF-TEST FAIL: %s -- expected %s, got %s"
                      % (label, "pass" if expect else "fail", "pass" if got else "fail"), file=sys.stderr)
                rc = 1
        finally:
            subprocess.run(["rm", "-rf", d])
    # structural guards
    d = tempfile.mkdtemp()
    if check_clone(d) == 0:
        print("SELF-TEST FAIL: a non-repository passed", file=sys.stderr); rc = 1
    # a repository with fork commits but no changed SOURCE files
    e = tempfile.mkdtemp()
    _mkclone(os.path.join(e, "r"), "data\n", "notes.txt")
    if check_clone(os.path.join(e, "r")) == 0:
        print("SELF-TEST FAIL: a clone with no changed source files passed", file=sys.stderr); rc = 1
    subprocess.run(["rm", "-rf", d, e])
    if rc == 0:
        print("ok   self-test: %d notice cases plus non-repository and empty-changed-set guards, "
              "all behaving as documented" % len(cases))
    return rc


args = sys.argv[1:]
if args[:1] == ["--self-test"]:
    sys.exit(self_test())
if not args:
    # Zero arguments used to exit 0 -- a check that ran on nothing, and said ok.
    print("FAIL: no clone directories given; refusing to report a pass on nothing", file=sys.stderr)
    print("usage: check_notices.sh <clone-dir> [<clone-dir> ...] | --self-test", file=sys.stderr)
    sys.exit(2)
sys.exit(max(check_clone(a) for a in args))
PYEOF
