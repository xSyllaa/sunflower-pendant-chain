#!/usr/bin/env python3
"""verify.py: standalone checker (Python 3 standard library only, no external dependencies) for
Step 1:  mu_3^S >= (1123 + sqrt(921182))^(1/17) = 1.5675245... (B-chain), and
Theorem: mu_3^S >= rho^(1/60) = 1.5675517557..., rho = 279351870092 + sqrt(378843386627*148811965456) (unit chain).
See README.md. Reads data/F.txt, data/H.txt, data/B.txt; builds the 28-point gadget with gadget.py.
Part I (step 1):
Checks: (1) F.txt and H.txt rebuilt from their one-line descriptions match the files; (2) F, H, B are Sperner and
3-sunflower-free (three distinct sets with A&B == A&C == B&C), exhaustively; (3) B.txt equals the DEGKM gluing of F
(at point 2) with H (at point 0); (4) the transfer matrix of B at ports x=0, z=1; (5) exact comparison of
rho(T_B)^(1/17) with the 49a record (1042 + 2*sqrt(261207))^(1/17); (6) explicit two-copy chain size vs the recursion."""
import sys, os, time, random
from decimal import Decimal, getcontext
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
t0 = time.time(); t_start = t0
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

print("Part I (step 1) done  (%.0fs)" % (time.time() - t_start))
# ---------------------------------------------------------------------------------------------------------------
# Part II: the 60-point unit
# Block method. Split the ground set into blocks P0 | P1 | P2. Three sets A, B, C form a 3-sunflower iff they are
# pairwise distinct and on every block their projections a, b, c satisfy a&b == a&c == b&c ("weak" triple, repeats
# allowed). A is a proper subset of B iff a <= b on every block and A != B. Group the family by its P0-projection p;
# the members with projection p are grouped into rectangles R1 x R2 (all members with the same set of P2-projections
# for their P1-projection form one rectangle), so every member lies in exactly one rectangle and every element of
# R1 x R2 is a member. For a triple of rectangles the achievable (block1, block2) choices are independent, so a
# sunflower exists iff for some weak triple of P0-projections, some rectangles, and some achievable equality patterns
# e0, e1, e2 on the three blocks, the three sets are pairwise distinct (no pair equal on all blocks). Exhaustive.
def block_check(name, fam, blocks, verbose=True):
    """blocks = [(shift, width)] * 3; returns (sperner_ok, sunflower_free); raises nothing."""
    assert len(set(fam)) == len(fam), name + ": duplicate member"
    (s0, w0), (s1, w1), (s2, w2) = blocks
    m0, m1, m2 = (1 << w0) - 1, (1 << w1) - 1, (1 << w2) - 1
    allbits = ((m0 << s0) | (m1 << s1) | (m2 << s2))
    assert all(M & ~allbits == 0 for M in fam), name + ": point outside the blocks"
    by0 = {}
    for M in fam:
        by0.setdefault(M >> s0 & m0, {}).setdefault(M >> s1 & m1, set()).add(M >> s2 & m2)
    rects = {}
    for p, d in by0.items():
        g = {}
        for p1, p2s in d.items(): g.setdefault(frozenset(p2s), set()).add(p1)
        rects[p] = [(frozenset(r1), r2) for r2, r1 in g.items()]
    assert sum(len(r1) * len(r2) for rs in rects.values() for r1, r2 in rs) == len(fam)
    # equality pattern of a triple: bit0 a==b, bit1 a==c, bit2 b==c
    def pat(a, b, c): return (a == b) | (a == c) << 1 | (b == c) << 2
    cache3 = {}
    def pats3(X, Y, Z):
        k = (X, Y, Z)
        if k not in cache3:
            out = set()
            for a in X:
                for b in Y:
                    ab = a & b
                    for c in Z:
                        if a & c == ab and b & c == ab: out.add(pat(a, b, c))
            cache3[k] = out
        return cache3[k]
    cache2 = {}
    def pats2(X, Y):   # {True if equal, False if proper subset} over a in X, b in Y with a <= b
        k = (X, Y)
        if k not in cache2:
            cache2[k] = {a == b for a in X for b in Y if a & b == a}
        return cache2[k]
    P = list(rects)
    sperner = True
    for pa in P:
        for pb in P:
            if pa & pb != pa: continue
            for ra in rects[pa]:
                for rb in rects[pb]:
                    E1 = pats2(ra[0], rb[0]); E2 = pats2(ra[1], rb[1])
                    if not E1 or not E2: continue
                    if any(not (pa == pb and e1 and e2) for e1 in E1 for e2 in E2):
                        sperner = False
    free = True
    for pa in P:
        for pb in P:
            ab = pa & pb
            for pc in P:
                if pa & pc != ab or pb & pc != ab: continue
                e0 = pat(pa, pb, pc)
                for ra in rects[pa]:
                    for rb in rects[pb]:
                        for rc in rects[pc]:
                            E1 = pats3(ra[0], rb[0], rc[0])
                            if not E1: continue
                            E2 = pats3(ra[1], rb[1], rc[1])
                            if any(e0 & e1 & e2 == 0 for e1 in E1 for e2 in E2):
                                free = False
    if verbose:
        print(f"{name}: {len(fam)} distinct members; {len(P)} P0-projections, "
              f"{sum(len(r) for r in rects.values())} rectangles; Sperner={sperner}, 3-sunflower-free={free} (exhaustive, block method)")
    return sperner, free

def check_quiet(fam):
    """True iff fam is 3-sunflower-free (fam assumed Sperner); brute force"""
    return brute(fam)[1]

def brute(fam):
    S = set(fam); sp = all(not (a != b and a & b == a) for a in fam for b in fam)
    sf = True
    for i, a in enumerate(fam):
        for j in range(i + 1, len(fam)):
            b = fam[j]; I = a & b
            for c in fam[j + 1:]:
                if a & c == I and b & c == I: sf = False
    return sp, sf

# (7) self-test of block_check against the brute-force checker on random small families with block structure
t1 = time.time()
rng = random.Random(20261019)
Bblocks = [(0, 9), (9, 10), (19, 0)]          # B: F-part (9 points) | H-part (10 points) | empty block
agree = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
for trial in range(60):
    fam = rng.sample(B, 60)
    for _ in range(rng.randrange(0, 3)):
        fam.append(rng.randrange(1, 1 << 19))
    fam = list(dict.fromkeys(fam))
    r1 = block_check("t", fam, Bblocks, verbose=False); r2 = brute(fam)
    assert r1 == r2, ("block checker disagrees with brute force", trial, r1, r2)
    agree[r1] += 1
# 3-block test families: random subfamilies of the gadget-like glue of three small pieces, plus random extra sets
import gadget as gd
X = gd.glue(gd.glue(F, 10, 4, H, 0), 19, 5, H, 0)
Gblocks = [(0, 8), (8, 10), (18, 10)]
for trial in range(60):
    fam = rng.sample(X, 50)
    for _ in range(rng.randrange(0, 3)):
        fam.append(rng.randrange(1, 1 << 28))
    fam = list(dict.fromkeys(fam))
    r1 = block_check("t", fam, Gblocks, verbose=False); r2 = brute(fam)
    assert r1 == r2, ("block checker disagrees with brute force", trial, r1, r2)
    agree[r1] += 1
# uniform test families (always Sperner) on 3 small blocks, so that sunflowers are frequent
pts = [0, 1, 2, 3, 8, 9, 10, 18, 19, 20]
for trial in range(60):
    k = rng.choice([3, 4]); fam = set()
    while len(fam) < 22:
        fam.add(sum(1 << p for p in rng.sample(pts, k)))
    fam = list(fam)
    r1 = block_check("t", fam, Gblocks, verbose=False); r2 = brute(fam)
    assert r1 == r2, ("block checker disagrees with brute force", trial, r1, r2)
    agree[r1] += 1
assert agree[(True, False)] > 0 and agree[(True, True)] > 0 and agree[(False, True)] > 0
assert block_check("B", B, Bblocks, verbose=False) == (True, True)
print("(7) block checker = brute force on 180 random test families (outcomes (Sperner, free): %s); on B: (True, True)  (%.0fs)"
      % (dict((str(k), v) for k, v in agree.items()), time.time() - t1))

# (8) the explicit gadget
t1 = time.time()
Gd = gd.build(F, H)
GADGET_SHA256 = "4ada967cda0b447dd71a541b5d0bb6705e9ee2ab0f1111ff532eb60ec2c01b81"
assert len(Gd) == len(set(Gd)) == 250659 and max(Gd) < 1 << 28
assert gd.sha256(Gd) == GADGET_SHA256, "gadget differs from the recorded one"
print("(8) gadget built by gadget.py: 250659 distinct sets on 28 points, sha256(canonical text) =", GADGET_SHA256)
assert block_check("gadget", Gd, Gblocks) == (True, True), "gadget is not Sperner and 3-sunflower-free"
print("    (%.0fs)" % (time.time() - t1))

# positive control at full size: H plus one extra set s (H+s Sperner, not sunflower-free) gives a gadget-size
# family that the block checker must reject
H2 = None
for sx in range(1, 1 << 11):
    if sx & 1 or any(a & sx == a or a & sx == sx for a in H): continue
    if not check_quiet(H + [sx]):
        H2 = H + [sx]; break
assert H2 is not None and brute(H2) == (True, False)
bad = gd.build(F, H2)
r = block_check("control: gadget built from H + {%s}" % " ".join(str(p) for p in range(11) if sx >> p & 1), bad, Gblocks)
assert r == (True, False), "positive control not detected"
# (9) root counts of the gadget, explicit and by the marginal formula
h1 = sum(D & 1 for D in H); h0 = len(H) - h1
def wt(bit, avoid, contain): return avoid if bit else contain   # glued point in C -> partner sets avoid their point
g_contain = sum(M & 1 for M in Gd); g_avoid = len(Gd) - g_contain
f_avoid = sum(wt(C >> 4 & 1, h0, h1) * wt(C >> 6 & 1, h0, h1) for C in F if not C & 1)
f_contain = sum(wt(C >> 4 & 1, h0, h1) * wt(C >> 6 & 1, h0, h1) for C in F if C & 1)
assert (g_avoid, g_contain) == (f_avoid, f_contain) == (98308, 152351)
print(f"(9) H at point 0: (avoid, contain) = ({h0}, {h1}); gadget root: (avoid, contain) = ({g_avoid}, {g_contain}), explicit = formula")

# (10) the unit: core F, in-port 0, out-port 1, gadget roots glued at core points 2 and 3 (62 points, 60 per period)
U = [[0, 0], [0, 0]]
for C in F:
    U[C & 1][C >> 1 & 1] += wt(C >> 2 & 1, g_avoid, g_contain) * wt(C >> 3 & 1, g_avoid, g_contain)
N00, N01, N10, N11 = U[0][0], U[0][1], U[1][0], U[1][1]
assert (N00, N01, N10, N11) == (148811965456, 279351870092, 279351870092, 378843386627)
assert sum(map(sum, U)) == sum(wt(C >> 2 & 1, g_avoid, g_contain) * wt(C >> 3 & 1, g_avoid, g_contain) for C in F)
print(f"(10) unit: N00={N00} N01={N01} N10={N10} N11={N11};  T = [[N10,N11],[N00,N01]],"
      f" rho = {N01} + sqrt({N11}*{N00})")

# (11) exact comparison rho^(1/60) > (1123+sqrt(921182))^(1/17):  rho^17 > rhoB^60, with rho >= lo, rhoB <= hi
S = 10 ** 40
lo = N01 * S + isqrt(N11 * N00 * S * S)              # rho  >= lo / S
hi = 1123 * S + isqrt(921182 * S * S) + 1            # rhoB <= hi / S
assert lo ** 17 * S ** 43 > hi ** 60
print("(11) exact: rho^17 > (1123+sqrt(921182))^60, hence rho^(1/60) > (1123+sqrt(921182))^(1/17) > (1042+2*sqrt(261207))^(1/17)")
getcontext().prec = 60
rho = Decimal(N01) + (Decimal(N11) * Decimal(N00)).sqrt()
rate = (rho.ln() / 60).exp(); rB = ((Decimal(1123) + Decimal(921182).sqrt()).ln() / 17).exp()
r0 = ((Decimal(1042) + 2 * Decimal(261207).sqrt()).ln() / 17).exp()
print(f"     rho = {rho:.6f};  rate rho^(1/60) = {rate:.13f}  vs step 1 {rB:.13f}  vs record {r0:.13f}")
print("total %.0fs" % (time.time() - t_start))
print("ALL CHECKS PASS")
