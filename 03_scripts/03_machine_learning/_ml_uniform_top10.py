# NOTE ON INPUTS: this script reads its inputs from ../../02_derived_data/, which is
# NOT distributed with this repository (see README section 1). Rebuild it by running
# the upstream step, or obtain the tables from the manuscript's supplementary files.
# -*- coding: utf-8 -*-
"""
Consensus of top-10 feature rankings from 10 learners; identifies 5 hub genes.
"""
import os, csv, io
from collections import Counter, OrderedDict
import os
_HERE = os.path.dirname(os.path.abspath(__file__))

BASE = os.path.join(_HERE, "..", "..", "02_derived_data", "per_learner_rankings")
# learner -> (file, ascending? True=higher is better)
FILES = OrderedDict([
    ("RF",       ("RF.csv",        True)),
    ("XGBoost",  ("XGBoost.csv",   True)),
    ("GBM",      ("GBM.csv",       True)),
    ("KNN",      ("KNN.csv",       True)),
    ("DecisionTree", ("DecisionTree.csv", True)),
    ("C5.0",     ("C5.0.csv",      True)),
    ("GLM",      ("GLM.csv",       True)),
    ("NNET",     ("NNET.csv",      True)),
    ("LASSO",    ("LASSO.csv",     True)),
    ("SVM",      ("SVM_RFE.csv",       True)),
])

def load(path):
    """Load rankings [(gene, value)]. Handles 2-col (var,imp) or 3-col (idx,var,imp)."""
    full = os.path.join(BASE, path)
    rows = []
    with io.open(full, "r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        next(r, None)  # header
        for row in r:
            if len(row) < 2:
                continue
            try:
                v = float(row[-1])
            except ValueError:
                continue
            g = row[-2].strip()
            if g:
                rows.append((g, v))
    rows.sort(key=lambda x: x[1], reverse=True)
    return rows

print("=" * 96)
print("Per-learner top-1 rankings")
print("=" * 96)
ranks = OrderedDict()
for algo, (p, _) in FILES.items():
    r = load(p)
    ranks[algo] = r
    print(f"  {algo:<14}  = {len(r):>4}   1 = {r[0][0]:<12} ={r[0][1]:.4f}   10 = {r[9][0]:<12}" if len(r) >= 10 else f"  {algo:<14} ={len(r)}")

print()
print("=" * 96)
print("Top-10 per learner")
print("=" * 96)
top10 = OrderedDict()
for algo, r in ranks.items():
    t = [g for g, _ in r[:10]]
    top10[algo] = t
    print(f"  {algo:<14} {' | '.join(t)}")

print()
print("=" * 96)
print("Intersection: RF ∩ XGBoost ∩ GBM top-10")
print("=" * 96)
s_rf, s_xgb, s_gbm = (set(top10[a]) for a in ("RF", "XGBoost", "GBM"))
inter3 = s_rf & s_xgb & s_gbm
print(f"   = {sorted(inter3)}   (n={len(inter3)})")
print(f"   5 hub genes {{ZEB2, NAP1L2, MST1, PLCB3, ASNS}}  match? {inter3 == {'ZEB2','NAP1L2','MST1','PLCB3','ASNS'}}")

print()
print("=" * 96)
print("Consensus count: top-10 across 10 learners")
print("=" * 96)
cnt = Counter()
for algo, t in top10.items():
    for g in set(t):
        cnt[g] += 1
total = sum(cnt.values())
print(f"   = {len(cnt)} genes    = {total}")
print()
print("               ")
print("  " + "-" * 88)
for i, (g, c) in enumerate(sorted(cnt.items(), key=lambda x: (-x[1], x[0])), 1):
    hit = [a for a in top10 if g in set(top10[a])]
    star = "  ← hub" if g in inter3 else ""
    print(f"  {i:>3}   {g:<12} {c:>4}/10   {', '.join(hit)}{star}")
    if i >= 25:
        break

print()
print("=" * 96)
print("5 hub genes consensus count")
print("=" * 96)
ordered = sorted(cnt.items(), key=lambda x: (-x[1], x[0]))
pos = {g: i for i, (g, _) in enumerate(ordered, 1)}
for g in ["ZEB2", "ASNS", "NAP1L2", "MST1", "PLCB3"]:
    print(f"  {g:<10}  {cnt[g]}/10    rank {pos[g]} ")
print()
print("  Other genes with >=5 consensus:")
for g, c in ordered:
    if g not in inter3 and c >= 5:
        print(f"    {g:<12} {c}/10    rank {pos[g]} ")

print()
print("=" * 96)
print("Top-6 consensus genes and hub count")
print("=" * 96)
top6 = [g for g, _ in ordered[:6]]
print(f"   top 6 = {top6}")
print(f"   hub genes in top 6: {len([g for g in top6 if g in inter3])}/5")

print()
print("=" * 96)
print("Intersection with LASSO / SVM added")
print("=" * 96)
for combo in [("RF","XGBoost","GBM","LASSO","SVM"),
              ("RF","XGBoost","GBM","LASSO"),
              ("RF","XGBoost","GBM","SVM"),
              ("RF","XGBoost","GBM","KNN"),
              ("RF","XGBoost","GBM","KNN","DecisionTree","C5.0")]:
    s = None
    for a in combo:
        s = set(top10[a]) if s is None else (s & set(top10[a]))
    tag = " = 5 hub genes" if s == {'ZEB2','NAP1L2','MST1','PLCB3','ASNS'} else ""
    print(f"  {'∩'.join(combo):<42} -> {sorted(s)}  (n={len(s)}){tag}")
