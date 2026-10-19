# A pendant chain for the 3-sunflower-free capacity

**Theorem.** Let f(n) be the largest size of a family of subsets of an n-element set with no three distinct members
A, B, C such that A∩B = A∩C = B∩C, and let μ₃ = lim f(n)^{1/n} (the Erdős–Szemerédi 3-sunflower-free capacity). Then

  μ₃ ≥ (1123 + √921182)^{1/17} = 1.5675245285…

This is above the previous best lower bound listed in T. Tao's optimization-constants table (constant 49a),
(1042 + 2√261207)^{1/17} = 1.5666969787… [Cla2026]. Like the other asterisked rows there, it is a limit obtained from
an infinite chain; no single finite uniform family of this rate is claimed.

Run `python3 verify.py` (Python 3, standard library only, about 10 seconds). The expected last line is
`ALL CHECKS PASS`; the output we got is in `verify_output.txt`. Optional second checker (naive triple loop, C):
`cc -O2 -o sf_tri extra/sf_tri.c && ./sf_tri data/B.txt` (output in `extra/sf_tri_output.txt`).

## The construction

Points of the 10-point family F are the edges of K₅ (vertices 0–4, edges in lexicographic order
(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4) = points 0–9).

* **F** (75 sets on 10 points, from [Cla2026]): the 5 stars K₁,₄, the 10 complete bipartite graphs K₂,₃, and the 60
  graphs whose complement in K₅ is a chair (degree sequence 3,2,1,1,1). `data/F.txt`.
* **H** (121 sets on Z₁₁, from [Cla2026]): the images of {0,…,6}, {0,1,2,3,4,6,7} and {0,1,2,6,9} under x ↦ ax+b
  with a a nonzero square mod 11. `data/H.txt`.
* **B** (4349 sets on 19 points, new): glue F at the edge (0,3) (point 2) to H at the point 0 with the DEGKM lemma
  below:
  B = {(C∖{2}) ∪ D : 2 ∈ C ∈ F, 0 ∉ D ∈ H} ∪ {C ∪ (D∖{0}) : 2 ∉ C ∈ F, 0 ∈ D ∈ H}.
  In `data/B.txt` the points of F other than 2 are renumbered 0–8 in order and the points 1–10 of H become 9–18.
  The chain ports of B are its points 0 and 1, i.e. the edges (0,1) and (0,2) of K₅: together with the glued edge
  (0,3) these are three edges at one vertex.

Then chain copies of B: glue the out-port (point 1) of each copy to the in-port (point 0) of the next one.
Equivalently: a chain of copies of F (glued (0,2) → (0,1)) in which every copy of F carries a pendant copy of H at its
edge (0,3). [Cla2026] alternated F and H *along* the chain; here H hangs off it.

## Gluing lemma (Deuber, Erdős, Gunderson, Kostochka, Meyer 1997, Lemma 4.4)

**Lemma.** Let 𝒜 and ℬ be Sperner, 3-sunflower-free families on disjoint ground sets X and Y, a ∈ X, b ∈ Y. Put
𝒜₁ = {C∖{a} : a ∈ C ∈ 𝒜}, 𝒜₀ = {C ∈ 𝒜 : a ∉ C}, and ℬ₁, ℬ₀ likewise with b. Then
𝒢 = {P ∪ Q : P ∈ 𝒜₁, Q ∈ ℬ₀} ∪ {P ∪ Q : P ∈ 𝒜₀, Q ∈ ℬ₁} is a Sperner, 3-sunflower-free family on
(X ∪ Y)∖{a,b}, and |𝒢| = |𝒜₁||ℬ₀| + |𝒜₀||ℬ₁|.

*Proof.* Every member of 𝒢 is written uniquely as P ∪ Q with P ⊆ X∖{a}, Q ⊆ Y∖{b}. 𝒜₁ ∩ 𝒜₀ = ∅, since P and
P ∪ {a} cannot both lie in the Sperner family 𝒜; likewise ℬ₁ ∩ ℬ₀ = ∅. So the two parts of 𝒢 are disjoint and the
count follows.

*Sperner.* Let P ∪ Q ⊆ P′ ∪ Q′ be members, so P ⊆ P′ and Q ⊆ Q′. If both lie in 𝒜₁⊗ℬ₀, then P∪{a} ⊆ P′∪{a} in 𝒜
and Q ⊆ Q′ in ℬ, so P = P′ and Q = Q′. Same for 𝒜₀⊗ℬ₁. If P∪Q ∈ 𝒜₁⊗ℬ₀ and P′∪Q′ ∈ 𝒜₀⊗ℬ₁, then
Q ⊊ Q′ ∪ {b} with both in ℬ, which is impossible. The other mixed case is symmetric (using 𝒜).

*No 3-sunflower.* Let three distinct members Pᵢ ∪ Qᵢ have equal pairwise intersections. Intersecting with X∖{a} and
with Y∖{b}, the triples (Pᵢ) and (Qᵢ) each have equal pairwise intersections.
* All three in 𝒜₁⊗ℬ₀: the Qᵢ lie in ℬ. If Q₁ = Q₂ then Q₁ = Q₁∩Q₂ = Q₁∩Q₃ ⊆ Q₃, so Q₁ = Q₃ (Sperner). Otherwise the
  Qᵢ would be pairwise distinct (by the same argument) and form a 3-sunflower in ℬ. So Q₁ = Q₂ = Q₃, the Pᵢ are
  distinct, and the sets Pᵢ ∪ {a} form a 3-sunflower in 𝒜. Contradiction. Symmetric for 𝒜₀⊗ℬ₁.
* Two in 𝒜₁⊗ℬ₀, one in 𝒜₀⊗ℬ₁: Q₁, Q₂, Q₃∪{b} lie in ℬ, and since b ∉ Q₁, Q₂ they have equal pairwise
  intersections. If Q₁ ≠ Q₂ they are three distinct sets (Q₃∪{b} contains b) forming a sunflower in ℬ. If Q₁ = Q₂,
  then Q₁ ⊆ Q₃ ⊊ Q₃ ∪ {b}, contradicting Sperner. The case of two in 𝒜₀⊗ℬ₁ is symmetric. ∎

B is obtained from F and H by one application of the lemma (both are Sperner and sunflower-free; `verify.py` also
checks B itself exhaustively).

## The chain and its count (transfer matrix as in [IRSI2026], [Cla2026])

Let N_ab = #{D ∈ B : [0 ∈ D] = a, [1 ∈ D] = b}. `verify.py` computes N₀₀ = 622, N₀₁ = N₁₀ = 1123, N₁₁ = 1481.
Let R₁ = B and R_{k+1} = 𝒢(R_k, out-port of its last copy; a fresh copy of B, its point 0). Each R_k is Sperner and
sunflower-free by the lemma (the two families glued are always on disjoint ground sets), on 19k − 2(k−1) = 17k + 2
points. If (u_k, w_k) counts members of R_k avoiding / containing the current out-port, then members containing it
combine with copy members avoiding point 0 and vice versa, so

  (u_{k+1}, w_{k+1}) = (u_k, w_k) · T,  T = [[N₁₀, N₁₁], [N₀₀, N₀₁]] = [[1123, 1481], [622, 1123]].

T is positive with Perron root ρ = 1123 + √(1481·622) = 1123 + √921182 ≈ 2082.7823, and (u₁, w₁) = (1745, 2604) is
positive, so |R_k| ≥ c ρ^k for some c > 0 (e.g. expand (u₁,w₁) in the positive Perron eigenvector). Since f is
nondecreasing and the limit μ₃ exists, μ₃ = lim f(17k+2)^{1/(17k+2)} ≥ ρ^{1/17} = 1.5675245285…

**Comparison with the previous record.** Both bounds add 17 points per period, so it suffices that
1123 + √921182 > 1042 + 2√261207 = 1042 + √1044828, i.e. 81 + √921182 > √1044828. Squaring: 6561 + 921182 + 162√921182
> 1044828 ⇔ 162√921182 > 117085 ⇔ 162²·921182 = 24175500408 > 117085² = 13708897225. True.

## What verify.py checks

1. `data/F.txt` and `data/H.txt` equal the families rebuilt from the descriptions above.
2. F, H and B have distinct members, are Sperner and 3-sunflower-free (exhaustive check of all pairs and third sets).
3. `data/B.txt` equals the gluing of F (point 2) and H (point 0).
4. The counts N_ab of B at ports 0, 1 and the matrix T.
5. The exact inequality above, and rigorous rational bounds on both Perron roots.
6. The explicit two-copy chain has 9,087,960 members, as the recursion predicts.

Independent reproductions: three separate fresh-context AI reviewers checked the claim without this repository's
code. Their notes record: F, H and B rebuilt from the descriptions; T and the exact inequality recomputed; the lemma
re-proved and its statement checked against the DEGKM paper; B glued to F (157,504 sets on 27 points) and the two-copy
chain (9,087,960 sets on 36 points) checked exhaustively for Sperner and sunflower-freeness.

## Credits

* W. A. Deuber, P. Erdős, D. S. Gunderson, A. V. Kostochka, A. G. Meyer, *Intersection statements for systems of
  sets*, J. Combin. Theory Ser. A 79 (1997) 118–132: the gluing lemma (Lemma 4.4).
* [IRSI2026] Institute for Recursive Self Improvement, https://github.com/Aspect5/sunflower-free-capacity: the exact
  2×2 transfer-matrix count of the DEGKM chain.
* [Cla2026] https://github.com/joshua17-create/sunflower-lower-bound: the families F and H and the previous record.
* T. Tao et al., optimization constants table, https://github.com/teorth/optimizationproblems (constant 49a).
* Found and checked with the V3R1T4S kit (https://github.com/elder-plinius/V3R1T4S) by an AI agent (Grok Bot) for
  Sylla (GitHub: xSyllaa).

## Status

Not peer-reviewed. The construction and the checks were produced by an AI agent and re-checked by three independent
AI reviewers; no human mathematician has checked it yet. Corrections are welcome as issues.

License: MIT (see LICENSE).
