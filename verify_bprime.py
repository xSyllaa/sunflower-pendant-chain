"""Verification of the B' bound  mu_3 >= (1127 + sqrt(925536))^(1/17) = 1.5678015201...
Python 3 standard library only. Usage: python3 verify_bprime.py   (expected last line: ALL CHECKS PASS)

1. Rebuild F (75 sets on the edges of K5) and H (121 sets on Z_11) from their definitions [Cla2026] and compare with
   data/F.txt, data/H.txt.
2. Rebuild B = DEGKM glue of F at point 2 and H at point 0 (4349 sets, 19 points) and check B is contained in B'.
3. Check that the 15 extra sets of B' lie on the F side (points 0..8).
4. Exhaustively check B' (data/Bprime.txt, 4364 sets): no duplicates, Sperner, no 3-sunflower.
   Sunflower test (bucket method): for every member S_k, group all other members S_i by S_i & S_k; members S_i, S_j of
   one group with S_i & S_j equal to the group key form a sunflower with S_k. Negative controls are run first.
5. Port counts N_ab at ports 0 (in) and 1 (out), transfer matrix, and exact integer comparisons with the previous
   records (Cla2026 1.5666969787..., same 17 points per period; Syl2026 step 1.5675517557..., exponent 1/60)."""
import sys, time, hashlib, os
from itertools import combinations
from math import isqrt
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
ok_all = True
def check(name, cond):
    global ok_all
    print(("PASS  " if cond else "FAIL  ") + name); ok_all &= bool(cond)
def read(fn):
    out = []
    for l in open(os.path.join(HERE, fn)):
        l = l.strip()
        if l and not l.startswith("#"):
            m = 0
            for t in l.split(): m |= 1 << int(t)
            out.append(m)
    return out
def mask(it): return sum(1 << i for i in it)
# 1. F and H from definitions
E = list(combinations(range(5), 2)); idx = {e: i for i, e in enumerate(E)}
def es(edges): return mask(idx[tuple(sorted(e))] for e in edges)
Fd = set()
for c in range(5): Fd.add(es([(c, w) for w in range(5) if w != c]))                       # stars K_{1,4}
for P in combinations(range(5), 2): Fd.add(es([(a, b) for a in P for b in range(5) if b not in P]))  # K_{2,3}
for c in range(5):                                                                       # complements of chairs
    for z in range(5):
        for w in range(5):
            if len({c, z, w}) < 3: continue
            x, y = [t for t in range(5) if t not in (c, z, w)]
            Fd.add(1023 ^ es([(c, x), (c, y), (c, z), (z, w)]))
Hd = set()
for base in ({0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 6, 7}, {0, 1, 2, 6, 9}):
    for a in (1, 3, 4, 5, 9):
        for b in range(11): Hd.add(mask((a * s + b) % 11 for s in base))
F = read("data/F.txt"); H = read("data/H.txt"); Bp = read("data/Bprime.txt")
check("F rebuilt from definition: 75 sets, equals data/F.txt", len(Fd) == 75 and sorted(Fd) == sorted(F))
check("H rebuilt from definition: 121 sets, equals data/H.txt", len(Hd) == 121 and sorted(Hd) == sorted(H))
# 2. B = glue(F at 2, H at 0)
fm = [0, 1, None, 2, 3, 4, 5, 6, 7, 8]
def fpart(C): return mask(fm[p] for p in range(10) if C >> p & 1 and p != 2)
def hpart(D): return mask(8 + q for q in range(1, 11) if D >> q & 1)
B = {fpart(C) | hpart(D) for C in Fd for D in Hd if (C >> 2 & 1) != (D & 1)}
SB = set(Bp)
check("B (glue) has 4349 sets and is contained in B'", len(B) == 4349 and B <= SB)
extra = [S for S in Bp if S not in B]
check("B' = B + 15 extra sets, all on points 0..8", len(extra) == 15 and all(S < (1 << 9) for S in extra))
# 4. exhaustive check
def sperner_sunflower(fam):
    n = len(fam); dup = n - len(set(fam)); sp = 0; sf = 0
    for i in range(n):
        a = fam[i]
        for j in range(n):
            if i != j and a != fam[j] and a & fam[j] == a: sp += 1
    for k in range(n):
        Sk = fam[k]; groups = {}
        for i in range(n):
            if i != k: groups.setdefault(fam[i] & Sk, []).append(fam[i])
        for key, g in groups.items():
            if len(g) < 2: continue
            for x in range(len(g)):
                gx = g[x]
                for y in range(x + 1, len(g)):
                    if gx & g[y] == key: sf += 1
    return dup, sp, sf // 3
# negative controls
check("control: planted sunflower {0,1},{0,2},{0,3} detected", sperner_sunflower([3, 5, 9])[2] == 1)
check("control: planted Sperner pair detected", sperner_sunflower([3, 7])[1] == 1)
check("control: B' + {19},{20},{21} (planted sunflower, empty core) detected at full size",
      sperner_sunflower(Bp + [1 << 19, 1 << 20, 1 << 21])[2] >= 1)
t = time.time(); dup, sp, sf = sperner_sunflower(Bp)
check("B' exhaustive: 4364 sets, 0 duplicates, 0 Sperner pairs, 0 3-sunflowers (%.0f s)" % (time.time() - t),
      len(Bp) == 4364 and dup == 0 and sp == 0 and sf == 0)
from functools import reduce
check("B' uses exactly the 19 points 0..18", reduce(lambda x, y: x | y, Bp) == (1 << 19) - 1)
# 5. matrix and inequalities
N = {(a, b): 0 for a in (0, 1) for b in (0, 1)}
for S in Bp: N[(S & 1, S >> 1 & 1)] += 1
T = [[N[(1, 0)], N[(1, 1)]], [N[(0, 0)], N[(0, 1)]]]
check("T = [[1127,1488],[622,1127]]", T == [[1127, 1488], [622, 1127]])
a, disc = 1127, 1488 * 622
def gt(q, a, disc, L):  # a + sqrt(disc) > q^L
    x = q ** L - a; return x < 0 or x * x < disc
def lt(q, a, disc, L):  # a + sqrt(disc) < q^L
    x = q ** L - a; return x > 0 and x * x > disc
check("rho^(1/17) in (1.5678015201, 1.5678015203)", gt(Fr(15678015201, 10**10), a, disc, 17) and lt(Fr(15678015203, 10**10), a, disc, 17))
# vs Cla2026: 1127+sqrt(925536) > 1042+sqrt(1044828)  <=>  85 > sqrt(1044828)-sqrt(925536)
check("exact: 1127+sqrt(925536) > 1042+2*sqrt(261207) [Cla2026, same 17 points]",
      (isqrt(1044828) + 1) - isqrt(925536) < 85)
q = Fr(156770, 10**5)
check("exact: rate > 1.56770 > rho_60^(1/60) [Syl2026 1.5675517557...]",
      gt(q, a, disc, 17) and lt(q, 279351870092, 378843386627 * 148811965456, 60))
print("sha256 data/Bprime.txt:", hashlib.sha256(open(os.path.join(HERE, "data/Bprime.txt"), "rb").read()).hexdigest())
print("ALL CHECKS PASS" if ok_all else "SOME CHECK FAILED"); sys.exit(0 if ok_all else 1)
