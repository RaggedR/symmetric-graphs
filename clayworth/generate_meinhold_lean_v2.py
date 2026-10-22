#!/usr/bin/env python3
"""Generate Lean 4 Meinhold-72 data using the EXACT Clayworth orbital from ClayworthGraph.lean.

The Lean file defines vertex 0's neighbors as {65, 1480, 2592} (0-indexed).
In GAP (1-indexed), vertex 1 has neighbors {66, 1481, 2593}.
We must use this specific orbital to ensure consistency.
"""

import subprocess, tempfile, os

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

# Find the SPECIFIC orbital matching ClayworthGraph.lean:
# Vertex 1 (GAP 1-indexed) has neighbors {66, 1481, 2593}
stab_cw1 := Stabilizer(cwPerm, 1);;
orbsCW := Orbits(stab_cw1, [2..n]);;

# The orbital from vertex 1 containing 66
targetOrb := First(Filtered(orbsCW, function(o) return Length(o) = 3; end),
    function(o) return 66 in o; end);;

if targetOrb = fail then
    Print("ERROR: Could not find orbital containing vertex 66 from vertex 1\n");
    QuitGap(1);
fi;
Print("Orbital: ", targetOrb, "\n");

# Build graph with this orbital using GRAPE
gamma := EdgeOrbitsGraph(cwPerm,
    List(targetOrb, function(v) return [1, v]; end), n);;
Print("Connected: ", IsConnectedGraph(gamma), "\n");
Print("Degree: ", Length(Adjacency(gamma, 1)), "\n");
Print("Neighbors of 1: ", Adjacency(gamma, 1), "\n");

# Verify matches Lean data (0-indexed: {65, 1480, 2592})
Assert(0, 66 in Adjacency(gamma, 1));
Assert(0, 1481 in Adjacency(gamma, 1));
Assert(0, 2593 in Adjacency(gamma, 1));

# K = PSL(3,2) for Meinhold-72
ints := IntermediateSubgroups(G, cwC3);;
K := fail;; chosenAdj := fail;; chosenBlockMap := fail;;
for Ki in ints.subgroups do
    if Size(Ki) <> 168 then continue; fi;
    Kimage := Image(cwAct, Ki);;
    block := Set(Orbit(Kimage, 1));;
    blocks := Orbit(cwPerm, block, OnSets);;
    blockMap := ListWithIdenticalEntries(n, 0);;
    for bi in [1..72] do
        for v in blocks[bi] do blockMap[v] := bi; od;
    od;
    qAdj := List([1..72], function(i) return []; end);;
    for v in [1..n] do
        bv := blockMap[v];
        for nbr in Adjacency(gamma, v) do
            bnbr := blockMap[nbr];
            if bv <> bnbr then AddSet(qAdj[bv], bnbr); fi;
        od;
    od;
    if Length(qAdj[1]) = 7 then
        K := Ki; chosenAdj := qAdj; chosenBlockMap := blockMap;
        break;
    fi;
od;

if K = fail then
    # Try all PSL(3,2) subgroups with different degrees
    Print("No 7-regular quotient found. Degrees: ");
    for Ki in ints.subgroups do
        if Size(Ki) <> 168 then continue; fi;
        Kimage := Image(cwAct, Ki);;
        block := Set(Orbit(Kimage, 1));;
        blocks := Orbit(cwPerm, block, OnSets);;
        blockMap := ListWithIdenticalEntries(n, 0);;
        for bi in [1..72] do
            for v in blocks[bi] do blockMap[v] := bi; od;
        od;
        qAdj := List([1..72], function(i) return []; end);;
        for v in [1..n] do
            bv := blockMap[v];
            for nbr in Adjacency(gamma, v) do
                bnbr := blockMap[nbr];
                if bv <> bnbr then AddSet(qAdj[bv], bnbr); fi;
            od;
        od;
        Print(Length(qAdj[1]), " ");
    od;
    Print("\n");
    # Pick any available one
    for Ki in ints.subgroups do
        if Size(Ki) <> 168 then continue; fi;
        Kimage := Image(cwAct, Ki);;
        block := Set(Orbit(Kimage, 1));;
        blocks := Orbit(cwPerm, block, OnSets);;
        blockMap := ListWithIdenticalEntries(n, 0);;
        for bi in [1..72] do
            for v in blocks[bi] do blockMap[v] := bi; od;
        od;
        qAdj := List([1..72], function(i) return []; end);;
        for v in [1..n] do
            bv := blockMap[v];
            for nbr in Adjacency(gamma, v) do
                bnbr := blockMap[nbr];
                if bv <> bnbr then AddSet(qAdj[bv], bnbr); fi;
            od;
        od;
        K := Ki; chosenAdj := qAdj; chosenBlockMap := blockMap;
        Print("Using degree ", Length(qAdj[1]), " quotient\n");
        break;
    od;
fi;

# Output neighbor data (0-indexed)
Print("MEINHOLD72_DEG:", Length(chosenAdj[1]), "\n");
Print("MEINHOLD72_NEIGHBORS:[");
for i in [1..72] do
    for j in chosenAdj[i] do
        Print(j-1, ",");
    od;
od;
Print("]\n");

# Output block map (0-indexed)
Print("BLOCK_MAP:[");
for v in [1..n] do
    Print(chosenBlockMap[v]-1, ",");
od;
Print("]\n");

# Verify: for each Clayworth edge, blocks are adjacent or equal
badCount := 0;;
for v in [1..n] do
    bv := chosenBlockMap[v];
    for nbr in Adjacency(gamma, v) do
        bnbr := chosenBlockMap[nbr];
        if bv <> bnbr and not (bnbr in chosenAdj[bv]) then
            badCount := badCount + 1;
        fi;
    od;
od;
Print("VERIFY_HOMOMORPHISM:", badCount, "\n");

QuitGap(0);
"""

with tempfile.NamedTemporaryFile(mode='w', suffix='.g', delete=False) as f:
    f.write(gap_script)
    gap_file = f.name

print("Computing Meinhold-72 with EXACT Clayworth orbital...")
result = subprocess.run(['gap', '-q', gap_file], capture_output=True, text=True, timeout=300)
os.unlink(gap_file)

print(result.stdout[:500])

# Parse
output = result.stdout
m72_neighbors = []
block_map = []
deg = 0

for line in output.split('\n'):
    if line.startswith('MEINHOLD72_DEG:'):
        deg = int(line.split(':')[1])
    elif line.startswith('VERIFY_HOMOMORPHISM:'):
        bad = int(line.split(':')[1])
        print(f"Homomorphism verification: {bad} bad edges")

start = output.find('MEINHOLD72_NEIGHBORS:[')
if start >= 0:
    end = output.find(']', start)
    data = output[start + len('MEINHOLD72_NEIGHBORS:['):end]
    m72_neighbors = [int(x) for x in data.replace('\n', '').split(',') if x.strip()]

start = output.find('BLOCK_MAP:[')
if start >= 0:
    end = output.find(']', start)
    data = output[start + len('BLOCK_MAP:['):end]
    block_map = [int(x) for x in data.replace('\n', '').split(',') if x.strip()]

print(f"Degree: {deg}")
print(f"Neighbors: {len(m72_neighbors)} entries ({len(m72_neighbors)//deg} × {deg})")
print(f"Block map: {len(block_map)} entries")

# Generate Lean neighbor array
print(f"\n--- Lean neighbor data ({deg}-regular) ---")
print(f"private def meinhold72NeighborData : Array (Fin 72) := #[")
for i in range(0, len(m72_neighbors), 14):
    chunk = m72_neighbors[i:i+14]
    print("  " + ", ".join(str(x) for x in chunk) + ",")
print("]")

# Generate adjacency function
print(f"\nprivate def meinhold72AdjBool (u v : Fin 72) : Bool :=")
print(f"  u ≠ v &&")
print(f"  (let b := {deg} * u.val")
checks = [f"   meinhold72NeighborData[b + {i}]'(by omega) == v" for i in range(deg)]
print(" ||\n".join(checks) + ")")

# Generate block map array
print(f"\nprivate def meinholdBlockMapData : Array (Fin 72) := #[")
for i in range(0, len(block_map), 30):
    chunk = block_map[i:i+30]
    print("  " + ", ".join(str(x) for x in chunk) + ",")
print("]")

# Write to file
lean_path = os.path.expanduser("~/git/mathlib4-fork/Archive/MeinholdGraphs_data.txt")
with open(lean_path, 'w') as f:
    f.write(f"-- Degree: {deg}\n")
    f.write(f"-- Neighbor count: {len(m72_neighbors)}\n")
    f.write(f"-- Block map count: {len(block_map)}\n\n")

    f.write("private def meinhold72NeighborData : Array (Fin 72) := #[\n")
    for i in range(0, len(m72_neighbors), 14):
        chunk = m72_neighbors[i:i+14]
        f.write("  " + ", ".join(str(x) for x in chunk) + ",\n")
    f.write("]\n\n")

    f.write("private def meinholdBlockMapData : Array (Fin 72) := #[\n")
    for i in range(0, len(block_map), 30):
        chunk = block_map[i:i+30]
        f.write("  " + ", ".join(str(x) for x in chunk) + ",\n")
    f.write("]\n")

print(f"\nData written to {lean_path}")
