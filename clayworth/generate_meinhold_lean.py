#!/usr/bin/env python3
"""Generate Lean 4 definition of Meinhold-72 and block map from Clayworth."""

import subprocess
import tempfile
import os

# GAP script to extract Meinhold-72 adj and block map
gap_script = r"""
LoadPackage("grape");;

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
cwC3 := First(
    Filtered(ConjugacyClassesSubgroups(G), c -> Size(Representative(c)) = 3),
    function(cl)
        return Number([1..63], function(p)
            return p^GeneratorsOfGroup(Representative(cl))[1] = p;
        end) = 0;
    end
);;
cwC3 := Representative(cwC3);;
cwAct := FactorCosetAction(G, cwC3);;
cwPerm := Image(cwAct);;
n := 4032;;

# Connected self-paired orbital
stab_cw1 := Stabilizer(cwPerm, 1);;
orbsCW := Orbits(stab_cw1, [2..n]);;
gamma := fail;;
for orb in Filtered(orbsCW, function(o) return Length(o) = 3; end) do
    v := orb[1];;
    stabv := Stabilizer(cwPerm, v);;
    orbOf1 := First(Orbits(stabv, [1..n]), function(o) return 1 in o; end);;
    if orbOf1 = fail or Length(orbOf1) <> 3 then continue; fi;
    testG := EdgeOrbitsGraph(cwPerm, List(orb, function(vv) return [1,vv]; end), n);;
    if IsConnectedGraph(testG) then gamma := testG; break; fi;
od;

# K = PSL(3,2) for Meinhold-72
ints := IntermediateSubgroups(G, cwC3);;
K := fail;;
for Ki in ints.subgroups do
    if Size(Ki) <> 168 then continue; fi;
    Kimage := Image(cwAct, Ki);;
    block := Set(Orbit(Kimage, 1));;
    blocks := Orbit(cwPerm, block, OnSets);;
    blockMap := ListWithIdenticalEntries(n, 0);;
    for bi in [1..72] do
        for v in blocks[bi] do blockMap[v] := bi; od;
    od;
    # Check quotient degree
    qAdj := List([1..72], function(i) return []; end);;
    for v in [1..n] do
        bv := blockMap[v];
        for nbr in Adjacency(gamma, v) do
            bnbr := blockMap[nbr];
            if bv <> bnbr then AddSet(qAdj[bv], bnbr); fi;
        od;
    od;
    if Length(qAdj[1]) = 7 then
        K := Ki; break;
    fi;
od;

if K = fail then Print("ERROR: K not found\n"); QuitGap(1); fi;

# Output Meinhold-72 neighbor data (0-indexed, 7 neighbors per vertex)
Print("MEINHOLD72_NEIGHBORS:[\n");
for i in [1..72] do
    for j in qAdj[i] do
        Print(j-1, ",");
    od;
od;
Print("]\n");

# Output block map (0-indexed: Clayworth vertex → Meinhold-72 vertex)
Print("BLOCK_MAP:[\n");
for v in [1..n] do
    Print(blockMap[v]-1, ",");
od;
Print("]\n");

# Output Clayworth neighbor data for verification
# (first 10 vertices only, to cross-check)
Print("CLAYWORTH_SAMPLE:");
for v in [1..10] do
    Print(Adjacency(gamma, v), ",");
od;
Print("\n");

QuitGap(0);
"""

with tempfile.NamedTemporaryFile(mode='w', suffix='.g', delete=False) as f:
    f.write(gap_script)
    gap_file = f.name

print("Computing Meinhold-72 data via GAP...")
result = subprocess.run(['gap', '-q', gap_file], capture_output=True, text=True, timeout=300)
os.unlink(gap_file)

# Parse neighbor data — output may span multiple lines
m72_neighbors = []
block_map = []

output = result.stdout
# Extract MEINHOLD72_NEIGHBORS
start = output.find('MEINHOLD72_NEIGHBORS:[')
if start >= 0:
    end = output.find(']', start)
    data = output[start + len('MEINHOLD72_NEIGHBORS:['):end]
    m72_neighbors = [int(x) for x in data.replace('\n', '').split(',') if x.strip()]

# Extract BLOCK_MAP
start = output.find('BLOCK_MAP:[')
if start >= 0:
    end = output.find(']', start)
    data = output[start + len('BLOCK_MAP:['):end]
    block_map = [int(x) for x in data.replace('\n', '').split(',') if x.strip()]

if not m72_neighbors or not block_map:
    print("ERROR: Failed to parse GAP output")
    print(result.stdout[:2000])
    print(result.stderr[:1000])
    exit(1)

print(f"Meinhold-72 neighbors: {len(m72_neighbors)} entries ({len(m72_neighbors)//7} vertices × 7)")
print(f"Block map: {len(block_map)} entries")

# Generate Lean file
lean_lines = []
lean_lines.append("""/-
Copyright (c) 2026 Robin Langer. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Robin Langer
-/
import Archive.ClayworthSabidussi

/-!
# Meinhold-72: quotient of the Clayworth graph via PSL(3,2)

The **Meinhold-72** graph (72 vertices, 7-regular, girth 4, diameter 5) is the
block quotient of the Clayworth graph (4032 vertices, cubic) by the intermediate
subgroup K ≅ PSL(3,2) of order 168.

`Aut(Meinhold-72) = G₂(2)`, order 12096.

Named after Nick Meinhold, co-founder of Imagineering Melbourne.

## The Langer graph is NOT a quotient of Clayworth

The Langer graph (63 vertices, 6-regular) is `Sab(G₂(2), H₁₉₂)` where
H₁₉₂ = Stab_G(p) is the point stabiliser of GH(2,2). The Clayworth C₃
generator has **zero fixed points** on the 63 GH(2,2) points, so C₃ is not
contained in any conjugate of H₁₉₂. No coset projection exists.

Furthermore, the three order-192 overgroups of C₃ have structure
(((C₄×C₄):C₂):C₃):C₂, NOT conjugate to H₁₉₂ ≅ (SL(2,3):C₄):C₂.

## The Meinhold family

Five quotients of Clayworth named for their sparsity and symmetry:

| Graph | V | Deg | Girth | Aut | Via K |
|-------|-----|-----|-------|-----|-------|
| Meinhold-72 | 72 | 7 | 4 | G₂(2) | PSL(3,2) |
| Meinhold-504 | 504 | 6 | 6 | G₂(2) | S₄ |
| Meinhold-576 | 576 | 7 | 6 | G₂(2) | C₇:C₃ |
| Meinhold-1344 | 1344 | 9 | 6 | 4·G₂(2) | C₃×C₃ |
| Meinhold-2016 | 2016 | 6 | 8 | G₂(2) | C₆ |

## References

* Robin Langer, *Symmetric Graphs and their Quotients*, arXiv:1306.4798
-/

set_option linter.style.longLine false
set_option linter.style.nativeDecide false

/-! ## Meinhold-72 neighbor data -/
""")

# Neighbor data (72 × 7 = 504 entries)
lean_lines.append(f"private def meinhold72NeighborData : Array (Fin 72) := #[")
for i in range(0, len(m72_neighbors), 14):
    chunk = m72_neighbors[i:i+14]
    lean_lines.append("  " + ", ".join(str(x) for x in chunk) + ",")
lean_lines.append("]")

lean_lines.append("""
/-- The Meinhold-72 graph: 72 vertices, 7-regular, girth 4, diameter 5. -/
def meinhold72Graph : SimpleGraph (Fin 72) where
  Adj v w := v ≠ w ∧ w ∈ (List.range 7).map
    (fun i => meinhold72NeighborData[v.val * 7 + i])
  symm := by
    intro v w ⟨hne, hmem⟩
    constructor
    · exact hne.symm
    · native_decide
  loopless v h := h.1 rfl
""")

# Block map
lean_lines.append("/-! ## Block map: Clayworth → Meinhold-72 -/\n")
lean_lines.append(f"private def meinholdBlockMapData : Array (Fin 72) := #[")
for i in range(0, len(block_map), 30):
    chunk = block_map[i:i+30]
    lean_lines.append("  " + ", ".join(str(x) for x in chunk) + ",")
lean_lines.append("]")

lean_lines.append("""
/-- The block map from Clayworth (4032 vertices) to Meinhold-72 (72 vertices). -/
def meinholdBlockMap : Fin 4032 → Fin 72 :=
  fun v => meinholdBlockMapData[v.val]

/-- The block map is surjective. -/
theorem meinholdBlockMap_surjective : Function.Surjective meinholdBlockMap := by
  native_decide

/-- The block map is a graph homomorphism: adjacent Clayworth vertices map to
    adjacent or equal Meinhold-72 vertices. -/
theorem meinholdBlockMap_adj {v w : Fin 4032}
    (h : clayworthGraph.Adj v w) :
    meinholdBlockMap v = meinholdBlockMap w ∨ meinhold72Graph.Adj (meinholdBlockMap v) (meinholdBlockMap w) := by
  native_decide
""")

lean_path = os.path.expanduser("~/git/mathlib4-fork/Archive/MeinholdGraphs.lean")
with open(lean_path, 'w') as f:
    f.write('\n'.join(lean_lines))

print(f"\nGenerated {lean_path}")
print(f"  Neighbor data: {len(m72_neighbors)} entries")
print(f"  Block map: {len(block_map)} entries")
