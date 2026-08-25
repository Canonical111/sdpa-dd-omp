#!/usr/bin/env python3
"""Check the published benchmark rows against the current headline claims made from them.

The point is not that the numbers are *good* -- it is that a reader can confirm the tables in
README/INSTALL/BENCHMARKS were actually derived from these rows, and that the rows are complete
enough to support them. Four classes of check:

  structure    every TSV parses, has a header, and has a .meta beside it
  completeness every cell has the repeat count its metadata claims
  health       no row hides a nonzero return code or a missing measurement
  derivation   the medians recomputed here equal the figures the documents publish

The last one is the reason this exists. A benchmark archive nobody recomputes from is decoration.
Run: python3 bench/validate.py   (from the repository root; CI runs it on every push)
"""
import csv, os, statistics as st, sys

HERE = os.path.dirname(os.path.abspath(__file__))
D    = os.path.join(HERE, "dd-port3-2026-08-24")
fails, checks = [], 0

def ck(ok, what, detail=""):
    global checks
    checks += 1
    if not ok:
        fails.append(f"{what}{(' :: ' + detail) if detail else ''}")
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}" + (f"   {detail}" if detail and not ok else ""))

def rows(name):
    with open(os.path.join(D, name)) as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")]

def med(rs, field, **eq):
    v = [float(r[field]) for r in rs if all(r[k] == str(w) for k, w in eq.items()) and r[field] not in ("", "NA")]
    return st.median(v) if v else None

def n_of(rs, **eq):
    return sum(1 for r in rs if all(r[k] == str(w) for k, w in eq.items()))

# ---- structure -------------------------------------------------------------------------------
print("structure")
tsvs = sorted(f for f in os.listdir(D) if f.endswith(".tsv"))
ck(len(tsvs) >= 12, f"{len(tsvs)} TSVs published")
for t in tsvs:
    ck(os.path.exists(os.path.join(D, t + ".meta")), f"{t} has provenance metadata")
    ck(len(rows(t)) > 0, f"{t} parses and is non-empty")

# ---- health ----------------------------------------------------------------------------------
print("\nhealth")
for t in tsvs:
    rs = rows(t)
    if "rc" in (rs[0].keys() if rs else {}):
        bad = [r for r in rs if r["rc"] not in ("0", "", "NA")]
        ck(not bad, f"{t}: every row returned 0", f"{len(bad)} nonzero rc")
    blank = [r for r in rs if any(v in ("", None) for v in r.values())]
    ck(not blank, f"{t}: no blank fields", f"{len(blank)} rows with a gap")

# ---- completeness ----------------------------------------------------------------------------
print("\ncompleteness")
h = rows("dd_final_headline_postfix.tsv")
for p, m in [("dE4", "auto"), ("dE3", "auto"), ("dE3", "fill")]:
    ck(n_of(h, problem=p, mode=m) == 3, f"headline {p}/{m}: 3 repeats")
s = rows("dd_v3_scaling_1t_24t.tsv")
for p, m in [("dE4", "auto"), ("dE3", "auto"), ("dE3", "fill")]:
    for th in ("1", "24"):
        ck(n_of(s, problem=p, mode=m, threads=th) == 3, f"scaling {p}/{m} at {th}t: 3 repeats")
u = rows("dd_upstream_vs_v3.tsv")
for b in ("upstream-6eaad8d", "v7.1.3-omp.3"):
    for p in ("dE4", "dE3"):
        for th in ("1", "2", "4", "8", "16", "24"):
            ck(n_of(u, build=b, problem=p, threads=th) == 3, f"upstream {b}/{p}/{th}t: 3 repeats")
f7 = rows("dd_fill_seven_structures.tsv")
for m in ("2439", "4489", "5278", "6067", "8359", "10614", "11227"):
    for mode in ("auto", "fill"):
        ck(n_of(f7, m=m, route_mode=mode) == 3, f"fill m={m}/{mode}: 3 repeats")

# ---- derivation: recompute what the documents publish -----------------------------------------
print("\nderivation -- recomputing the published figures from these rows")
# BOTH sides of the 7.13x come from published rows. This used to be a hard-coded literal sourced
# from a TSV that exists only in the private evidence repository, so the validator derived one side
# of a public claim and trusted the other -- caught in review.
hist = rows("dd_historical_baseline.tsv")
HIST_DE4_24T = med(hist, "mainloop_s", problem="dE4", threads="24", mode="-")
ck(HIST_DE4_24T is not None, "historical dE4 baseline is present in the published archive")

def near(a, b, tol=0.005):
    return a is not None and abs(a - b) / b <= tol

d4 = med(h, "mainloop_s", problem="dE4", mode="auto")
ck(near(d4, 6.495, 0.01), "headline campaign 1: dE4 6.495 s", f"got {d4}")

# The published headline is POOLED across all three campaigns of this cell, because a single
# campaign's median is not stable to three figures -- nine runs span 8.2%. Recomputed here rather
# than asserted, and the ~7.2x in README must bracket it.
cur = rows("dd_current_scaling.tsv")
pool = ([float(r["mainloop_s"]) for r in h  if r["problem"] == "dE4" and r["mode"] == "auto"] +
        [float(r["mainloop_s"]) for r in s  if r["problem"] == "dE4" and r["mode"] == "auto"
                                            and r["threads"] == "24"] +
        [float(r["mainloop_s"]) for r in cur if r["problem"] == "dE4" and r["mode"] == "default"
                                            and r["threads"] == "24"])
ck(len(pool) == 9, "dE4 headline cell: 9 published runs across 3 campaigns", str(len(pool)))
spread = (max(pool) - min(pool)) / st.median(pool) * 100
ck(spread > 5.0, "that cell really is jittery, which is why ~7.2x is an approximation",
   f"spread {spread:.1f}%")
ratio = HIST_DE4_24T / st.median(pool)
ck(7.15 <= ratio <= 7.25, "pooled headline lands in the stated ~7.2x", f"got {ratio:.3f}")

# current-build cells, the ones INSTALL now publishes
for p_, m_, th, want in [("dE4","default","1",46.555), ("dE4","default","24",6.390),
                         ("dE3","default","1",31.361), ("dE3","default","24",4.417),
                         ("dE3","legacy","1",176.287), ("dE3","legacy","24",22.526)]:
    got = med(cur, "mainloop_s", problem=p_, mode=m_, threads=th)
    ck(near(got, want, 0.01), f"current build: {p_}/{m_} at {th}t = {want} s", f"got {got}")
ck(len({r["objValPrimal"] for r in cur if r["problem"] == "dE3"}) == 1,
   "current build: dE3 gives ONE objective across both routes and both thread counts")

for p, m, th, want in [("dE4","auto","1",46.612), ("dE4","auto","24",6.399),
                       ("dE3","auto","1",176.357), ("dE3","auto","24",22.645),
                       ("dE3","fill","1",31.370), ("dE3","fill","24",4.475)]:
    got = med(s, "mainloop_s", problem=p, mode=m, threads=th)
    ck(near(got, want, 0.01), f"INSTALL scaling: {p}/{m} at {th}t = {want} s", f"got {got}")

for b, p, th, want in [("upstream-6eaad8d","dE4","1",49.24), ("upstream-6eaad8d","dE4","2",49.14),
                       ("upstream-6eaad8d","dE4","4",47.80), ("upstream-6eaad8d","dE4","8",47.14),
                       ("upstream-6eaad8d","dE4","16",46.90), ("upstream-6eaad8d","dE4","24",46.95),
                       ("v7.1.3-omp.3","dE4","1",46.62), ("v7.1.3-omp.3","dE4","2",47.06),
                       ("v7.1.3-omp.3","dE4","4",24.26), ("v7.1.3-omp.3","dE4","8",14.14),
                       ("v7.1.3-omp.3","dE4","16",8.29), ("v7.1.3-omp.3","dE4","24",6.49),
                       ("upstream-6eaad8d","dE3","1",434.85), ("upstream-6eaad8d","dE3","2",409.14),
                       ("upstream-6eaad8d","dE3","4",211.32), ("upstream-6eaad8d","dE3","8",120.69),
                       ("upstream-6eaad8d","dE3","16",71.08), ("upstream-6eaad8d","dE3","24",60.92),
                       ("v7.1.3-omp.3","dE3","1",176.26), ("v7.1.3-omp.3","dE3","2",168.11),
                       ("v7.1.3-omp.3","dE3","4",88.20), ("v7.1.3-omp.3","dE3","8",49.63),
                       ("v7.1.3-omp.3","dE3","16",27.86), ("v7.1.3-omp.3","dE3","24",22.81)]:
    got = med(u, "mainloop_s", build=b, problem=p, threads=th)
    ck(near(got, want, 0.01), f"BENCHMARKS upstream: {b}/{p} at {th}t = {want} s", f"got {got}")

for m, want_a, want_f in [("2439",1.924,0.606), ("4489",9.219,3.593), ("5278",14.502,3.877),
                          ("6067",22.974,4.498), ("8359",66.203,15.131),
                          ("10614",150.473,22.362), ("11227",175.807,47.944)]:
    for mode, want in (("auto",want_a), ("fill",want_f)):
        got = med(f7, "mainloop_s", m=m, route_mode=mode)
        ck(near(got, want, 0.01), f"INSTALL fill: m={m}/{mode} = {want} s", f"got {got}")

b = rows("dd_bitset_memory.tsv")
for p, mo, want in [("dE4","auto",382.0), ("dE3","fill",277.5)]:
    got = med(b, "rss_kb", build="bitset", problem=p, mode=mo)
    ck(near(got/1024, want, 0.01), f"INSTALL memory: {p}/{mo} = {want} MB", f"got {got/1024:.1f}")

# ---- the property the whole release is about --------------------------------------------------
print("\ndeterminism recorded in the rows")
objs = {r["objValPrimal"] for r in s if r["problem"] == "dE4"}
ck(len(objs) == 1, "scaling: dE4 gives ONE objective across 1 and 24 threads", f"{len(objs)} distinct")
objs = {r["objValPrimal"] for r in s if r["problem"] == "dE3"}
ck(len(objs) == 1, "scaling: dE3 gives ONE objective across both routes and both thread counts",
   f"{len(objs)} distinct")
t6 = rows("dd_race_truss6.tsv")
racing = {r["objValPrimal"] for r in t6 if r["build"] == "v7.1.3-omp.2"}
fixed  = {r["objValPrimal"] for r in t6 if r["build"] == "fixed"}
ck(len(racing) == 5, "truss6: the withdrawn release gives 5 distinct objectives in 5 runs", f"{len(racing)}")
ck(len(fixed) == 1, "truss6: the fixed build gives one", f"{len(fixed)}")

print(f"\n{checks - len(fails)}/{checks} checks passed")
if fails:
    print("\nFAILURES:", file=sys.stderr)
    for f in fails:
        print("  " + f, file=sys.stderr)
    sys.exit(1)
