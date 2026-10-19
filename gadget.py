#!/usr/bin/env python3
"""gadget.py: deterministic generator for the 28-point gadget family (Python 3 standard library only).

Gadget = copy of F (75 sets on the edges of K5) with a copy of H (121 sets on Z_11) glued, by the DEGKM 1997
Lemma 4.4 construction, at F-point 4 = (1,2) and at F-point 6 = (1,4) (H is glued at its point 0 both times).
Points of the gadget (0..27):
   0..7   F-points 0,1,2,3,5,7,8,9 = edges (0,1),(0,2),(0,3),(0,4),(1,3),(2,3),(2,4),(3,4); point 0 = (0,1) is the root
   8..17  H-points 1..10 of the copy glued at F-point 4
  18..27  H-points 1..10 of the copy glued at F-point 6
Usage: python3 gadget.py            -> prints the member count, root counts and sha256 of the canonical text
       python3 gadget.py --write    -> also writes data/gadget.txt (one set per line, points ascending,
                                       lines in ascending order of the bitmask sum(2^p))
"""
import os, sys, hashlib
here = os.path.dirname(os.path.abspath(__file__))

def read(fn):
    out = []
    for line in open(os.path.join(here, "data", fn)):
        if line.startswith("#") or not line.strip(): continue
        out.append(sum(1 << int(v) for v in line.split()))
    return out

def drop(M, p):
    """remove point p from bitmask M, shifting higher points down by one"""
    return (M & ((1 << p) - 1)) | ((M >> (p + 1)) << p)

def glue(X, nX, a, Y, b):
    """DEGKM gluing of X (on nX points) at point a with Y at point b:
    {(C-a) u D : a in C, b not in D}  u  {C u (D-b) : a not in C, b in D}; X's points first, then Y's."""
    out = []
    for C in X:
        ca = C >> a & 1
        for D in Y:
            if ca + (D >> b & 1) == 1:
                out.append(drop(C, a) | (drop(D, b) << (nX - 1)))
    return out

def build(F, H):
    X1 = glue(F, 10, 4, H, 0)          # 19 points: F-points 0,1,2,3,5,6,7,8,9 then H1-points 1..10
    return glue(X1, 19, 5, H, 0)       # point 5 of X1 is F-point 6; result has 28 points

def canonical_text(G):
    lines = []
    for M in sorted(G):
        lines.append(" ".join(str(p) for p in range(28) if M >> p & 1))
    return "\n".join(lines) + "\n"

def sha256(G):
    return hashlib.sha256(canonical_text(G).encode()).hexdigest()

if __name__ == "__main__":
    G = build(read("F.txt"), read("H.txt"))
    contain = sum(M & 1 for M in G)
    print(f"gadget: {len(G)} members ({len(set(G))} distinct) on 28 points; root 0: avoid {len(G)-contain}, contain {contain}")
    print("sha256(canonical data/gadget.txt) =", sha256(G))
    if "--write" in sys.argv:
        open(os.path.join(here, "data", "gadget.txt"), "w").write(canonical_text(G))
        print("wrote data/gadget.txt")
