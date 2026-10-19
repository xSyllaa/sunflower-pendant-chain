#!/usr/bin/env python3
"""verify.py: standalone checker (Python 3 standard library only, no external dependencies) for
mu_3^S >= (1123 + sqrt(921182))^(1/17) = 1.5675245... (see README.md). Reads data/F.txt, data/H.txt, data/B.txt.
Checks: (1) F.txt and H.txt rebuilt from their one-line descriptions match the files; (2) F, H, B are Sperner and
3-sunflower-free (three distinct sets with A&B == A&C == B&C), exhaustively; (3) B.txt equals the DEGKM gluing of F
(at point 2) with H (at point 0); (4) the transfer matrix of B at ports x=0, z=1; (5) exact comparison of
rho(T_B)^(1/17) with the 49a record (1042 + 2*sqrt(261207))^(1/17); (6) explicit two-copy chain size vs the recursion."""
import sys, os, time
from fractions import Fraction
from math import isqrt
from itertools import combinations
here = os.path.dirname(os.path.abspath(__file__))
def read(fn):
    out = []
    for line in open(os.path.join(here, "data", fn)):
        if line.startswith("#") or not line.strip(): continue
        out.append(sum(1 << int(v) for v in line.split()))
    return out
def check_family(name, fam):
    assert len(set(fam)) == len(fam), name + ": duplicate member"
    for a in fam:
        for b in fam:
            if a != b and a & b == a: raise SystemExit(f"{name}: not Sperner")
    n = len(fam)
    for i in range(n):
        A = fam[i]; buck = {}
        for j in range(n):
            if j != i: buck.setdefault(fam[j] & A, []).append(fam[j])
        for j in range(n):
            if j == i: continue
            Bm = fam[j]; I = A & Bm
            for C in buck.get(I, ()):
                if C != Bm and C & Bm == I: raise SystemExit(f"{name}: 3-sunflower {A:b} {Bm:b} {C:b}")
    print(f"{name}: {n} distinct members, Sperner, 3-sunflower-free (exhaustive)")
# (1) rebuild F and H
E5 = list(combinations(range(5), 2)); eid = {e: i for i, e in enumerate(E5)}
def em(es): return sum(1 << eid[tuple(sorted(e))] for e in es)
Fr = set(em([e for e in E5 if v in e]) for v in range(5))
for S in combinations(range(5), 2):
    Fr.add(em([(a, b) for a in S for b in range(5) if b not in S]))
for es in combinations(E5, 4):   # complement is a chair: a tree with degree sequence 3,2,1,1,1
    deg = [sum(v in e for e in es) for v in range(5)]
    par = list(range(5))
    def fnd(x):
        while par[x] != x: x = par[x]
        return x
    acyclic = True
    for a, b in es:
        ra, rb = fnd(a), fnd(b)
        if ra == rb: acyclic = False
        par[ra] = rb
    if acyclic and sorted(deg) == [1, 1, 1, 2, 3]:
        Fr.add(em([e for e in E5 if e not in es]))
Hr = set()
for base in ({0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 6, 7}, {0, 1, 2, 6, 9}):
    for a in (1, 3, 4, 5, 9):
        for b in range(11): Hr.add(sum(1 << ((a * x + b) % 11) for x in base))
F = read("F.txt"); H = read("H.txt"); B = read("B.txt")
assert set(F) == Fr and len(F) == 75, "F.txt does not match its description"
assert set(H) == Hr and len(H) == 121, "H.txt does not match its description"
print("F.txt, H.txt match their descriptions (75 and 121 sets)")
t0 = time.time()
check_family("F", F); check_family("H", H); check_family("B", B)
print(f"  ({time.time()-t0:.0f}s)")
# (3) gluing
def drop(M, p): return (M & ((1 << p) - 1)) | ((M >> (p + 1)) << p)
G = set()
for C in F:
    for D in H:
        if (C >> 2 & 1) + (D & 1) == 1: G.add(drop(C, 2) | (drop(D, 0) << 9))
assert G == set(B) and len(G) == len(B) == 4349, "B is not the gluing of F and H"
assert max(B) < 1 << 19
print("B.txt = DEGKM gluing of F (point 2) and H (point 0), 4349 sets on 19 points")
# (4) transfer matrix at ports x=0 (in), z=1 (out)
N = [[0, 0], [0, 0]]
for M in B: N[M & 1][M >> 1 & 1] += 1
T = [[N[1][0], N[1][1]], [N[0][0], N[0][1]]]
print("N[a][b] = #{D: [0 in D]=a, [1 in D]=b} =", N, " T_B = [[N10,N11],[N00,N01]] =", T)
# (5) exact comparison: rho(T)^17 > rho0^17 where rho0 = 1042 + 2 sqrt(261207); i.e. rho(T) > rho0 (same exponent 17)
tr = T[0][0] + T[1][1]; det = T[0][0] * T[1][1] - T[0][1] * T[1][0]; disc = tr * tr - 4 * det
S = 10 ** 50
rho_lo = (Fraction(tr) + Fraction(isqrt(disc * S * S), S)) / 2
rho0_hi = 1042 + 2 * Fraction(isqrt(261207 * S * S) + 1, S)
print(f"rho(T_B) = ({tr} + sqrt({disc}))/2 >= {float(rho_lo):.10f};  record rho0 = 1042+2*sqrt(261207) <= {float(rho0_hi):.10f}")
assert rho_lo > rho0_hi
# exact symbolic: rho(T_B) = 1123 + sqrt(921182) > 1042 + sqrt(1044828)  <=>  81 + sqrt(921182) > sqrt(1044828)
assert T[0][0] == T[1][1] == 1123 and T[0][1] * T[1][0] == 921182
lhs2 = 81 * 81 + 921182; rest = 1044828 - lhs2   # (81+s)^2 = 6561 + 921182 + 162 s > 1044828 <=> 162 s > rest
assert rest < 0 or (162 * 162 * 921182 > rest * rest)
print("exact: 1123 + sqrt(921182) > 1042 + 2*sqrt(261207); both chains add 17 points per period")
print("rate rho(T_B)^(1/17) ~ %.10f  vs record %.10f" % (float(rho_lo) ** (1 / 17), float(rho0_hi) ** (1 / 17)))
# (6) two-copy chain: explicit size equals (u,w)·T·1
u = sum(1 for M in B if not M >> 1 & 1); w = len(B) - u
two = sum(1 for C in B for D in B if (C >> 1 & 1) + (D & 1) == 1)
pred = u * (T[0][0] + T[0][1]) + w * (T[1][0] + T[1][1])
assert two == pred
print(f"two-copy chain: {two} members = (u,w)·T·1 with (u,w)=({u},{w})")
print("ALL CHECKS PASS")
