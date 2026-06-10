# TODO — Mathlib Graph Theory PRs

Updated: 2026-05-22

## PR Stack

```
#39530 (graph-action) ← #39548 (coset-graph) ← #39550 (sabidussi) ← #39551 (lorimer)
                                                                          ↑
#39653 (css-surface-codes) ← #39654 (surface-instances) ← #39649 (langer)
                                                                          ↑
                             #39650 (named-graphs) ← #39651 (voltage-graphs)
                                                                          ↑
                                                     #39684 (meinhold-quotients)
```

---

## HIGH PRIORITY

### Worktree A: `worktree-meinhold` — Meinhold Sabidussi + Lorimer
**Branch**: `feat/meinhold-quotients` (PR #39684)
**File**: `Archive/MeinholdGraphs.lean`

1. **Sabidussi proof for Meinhold-72**: `meinhold72SabidussiIso`
   - Meinhold-72 = Sab(G₂(2), PSL(3,2), D), 72 = 12096/168 vertices
   - Need G₂(2) generators as permutations of Fin 72 (compute in GAP)
   - BFS witness words for 72 vertices (small, fast)
   - Use SabidussiWitness machinery from PR #39654

2. **Lorimer quotient**: prove Clayworth → Meinhold-72 is the coset projection
   - `meinholdBlockMap = cosetProjection C₃ PSL(3,2)`
   - Requires: Clayworth Sabidussi iso, Meinhold-72 Sabidussi iso, show they commute with block map
   - First concrete use of `quotient_cosetGraph_iso` from PR #39551

### Worktree B: `worktree-zhou-lorimer` — Zhou quotient via Lorimer
**Branch**: `feat/surface-instances` (PR #39654)
**File**: `Archive/ZhouGraph.lean`

1. **Lorimer quotient**: prove Zhou-3 → Zhou-6 is the coset projection
   - `zhouBlockMap = cosetProjection S₃ D₁₂` via S₃ ≤ D₁₂
   - Simplest Lorimer instantiation (182 → 91 vertices)
   - Both graphs already have Sabidussi proofs

2. **Arc-transitivity**: prove Zhou-3 is arc-transitive via `lorimer_forward`
   - Connection set is single double coset S₃·a·S₃ where a² ∈ S₃

### Worktree C: `worktree-dodec-lorimer` — Named graphs via Lorimer
**Branch**: `feat/named-graphs-archive` (PR #39650)
**File**: `Archive/NamedGraphs.lean`

1. **Lorimer quotient**: prove dodecahedron → Petersen is the coset projection
   - `dodecAntipodalMap = cosetProjection C₃ D₅` via C₃ ≤ D₅ in A₅
   - Both graphs already have Sabidussi proofs

2. **Arc-transitivity** for dodecahedron and Petersen via `lorimer_forward`

3. **Petersen verification**: prove `petersenGraph` matches GP(5,2) or Kneser K(5,2)
   - Currently only defined as quotient; no independent verification

4. **Girth**: `petersenGraph_girth = 5`, `dodecahedronGraph_girth = 5`

---

## MEDIUM PRIORITY

### Worktree D: `worktree-css-endtoend` — CSS k = 2g instantiation
**Branch**: `feat/surface-instances` (PR #39654)
**File**: `Archive/HeawoodSurface.lean`

1. **Connectivity**: prove `heawoodSurface.toSimpleGraph.Connected`
2. **Dual connectivity**: prove dual graph of Heawood surface is connected
3. **Instantiate `css_k_eq_2g_from_surface`**: Heawood → [[21, 2]]
   - First formally verified CSS code parameters
4. If easy, repeat for Bolza (genus 2, [[24, 4]])

### Worktree E: `worktree-voltage-lorimer` — Voltage graphs arc-transitivity
**Branch**: `feat/voltage-graphs` (PR #39651)
**File**: `Archive/VoltageGraphs.lean`

1. Prove Heawood arc-transitive via `lorimer_forward`
2. Prove Möbius-Kantor arc-transitive via `lorimer_forward`
3. Prove `voltageGraphK2` is always bipartite
4. Make explicit: voltageGraphK2 = Sab(Z_m ⋊ Z₂, ⟨1⟩, D)

### Girth proofs (can be added to any Archive branch)

- Heawood (6), Möbius-Kantor (6), Klein (6)
- All checkable by native_decide on these small graphs
- Needs `SimpleGraph.girth` — check if already in Mathlib

---

## LOW PRIORITY

- Klein quartic Z₇ quotient structure (8-vertex quotient)
- Remaining Meinhold family: 504, 576, 1344, 2016 vertices
- `GraphAction.toMonoidHom : G →* Aut(Gamma)`
- `IsEdgeTransitive` for semisymmetric graphs
- `cosetGraph.degree` and `cosetGraph.Connected` characterization
- Cayley graph specialization (H = ⊥)

---

## PROBABLY NOT NECESSARY RIGHT NOW

- Distance-regular intersection arrays ({6,4,4;1,1,3} for Langer)
- Diameter proofs
- CSS code distance bounds
- |G₂(2)| = 12096 formalized
- Q₆₂ ↔ CellularSurface connection
- Meinhold-1344 and Meinhold-2016

---

## Parallel Worktree Plan

Five worktrees, each touching different files:

| Worktree | Branch | Key file(s) | Blocked by |
|----------|--------|-------------|------------|
| A: `worktree-meinhold` | `feat/meinhold-quotients` | `MeinholdGraphs.lean` | PRs 1-4 in Mathlib |
| B: `worktree-zhou-lorimer` | `feat/surface-instances` | `ZhouGraph.lean` | PRs 1-4 in Mathlib |
| C: `worktree-dodec-lorimer` | `feat/named-graphs-archive` | `NamedGraphs.lean` | PRs 1-4 in Mathlib |
| D: `worktree-css-endtoend` | `feat/surface-instances` | `HeawoodSurface.lean` | PR 8 in Mathlib |
| E: `worktree-voltage-lorimer` | `feat/voltage-graphs` | `VoltageGraphs.lean` | PRs 1-4 in Mathlib |

**A, B, C** are the Lorimer instantiation family — same pattern, different graphs.
**D** is independent (CSS chain, not Sabidussi chain).
**E** is lower priority, same pattern as A/B/C.

**B and D** share a branch but touch different files — can run in parallel.

**Critical path**: PRs #39530-#39551 (the four Mathlib library files) must merge first. All Lorimer work depends on them. CSS end-to-end (D) only needs PR #39653.
