# Compute the Zhou-3 → Zhou-6 covering map
# using the actual Lean vertex numbering.
#
# Strategy: the generators in Zhou3Sabidussi.lean generate PSL(2,13) ≤ Sym(182).
# Find the unique block system with blocks of size 2, then match the quotient
# to the Zhou-6 graph from Zhou6Sabidussi.lean.

# Zhou-3 generators (copied from Zhou3Sabidussi.lean, 0-indexed → 1-indexed)
zhou3g1 := PermList([
  3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 29, 30, 33, 34, 35, 36, 18, 17, 37, 38, 39, 40, 41, 42, 43, 44,
  47, 48, 49, 50, 10, 9, 51, 52, 55, 56, 57, 58, 1, 2, 63, 64, 46, 45, 69, 70, 73, 74, 77, 78, 81, 82, 85, 86, 89, 90,
  83, 84, 93, 94, 97, 98, 99, 100, 103, 104, 107, 108, 5, 6, 113, 114, 117, 118, 121, 122, 123, 124, 127, 128, 101, 102, 61, 62, 131, 132,
  135, 136, 137, 138, 141, 142, 143, 144, 145, 146, 13, 14, 22, 21, 60, 59, 106, 105, 151, 152, 68, 67, 155, 156, 159, 160, 161, 162, 163, 164,
  167, 168, 171, 172, 169, 170, 126, 125, 96, 95, 71, 72, 120, 119, 134, 133, 32, 31, 53, 54, 130, 129, 115, 116, 173, 174, 80, 79, 150, 149,
  148, 147, 66, 65, 165, 166, 111, 112, 154, 153, 178, 177, 176, 175, 181, 182, 109, 110, 88, 87, 25, 26, 158, 157, 92, 91, 140, 139, 75, 76,
  179, 180
]);;

zhou3g2 := PermList([
  5, 6, 9, 10, 13, 14, 17, 18, 21, 22, 25, 26, 1, 2, 31, 32, 19, 20, 7, 8, 3, 4, 27, 28, 36, 35, 34, 33, 45, 46,
  40, 39, 24, 23, 12, 11, 53, 54, 16, 15, 59, 60, 61, 62, 65, 66, 67, 68, 71, 72, 75, 76, 79, 80, 83, 84, 87, 88, 91, 92,
  86, 85, 95, 96, 29, 30, 101, 102, 105, 106, 109, 110, 111, 112, 115, 116, 119, 120, 37, 38, 125, 126, 99, 100, 44, 43, 129, 130, 133, 134,
  41, 42, 139, 140, 98, 97, 64, 63, 55, 56, 47, 48, 138, 137, 147, 148, 149, 150, 49, 50, 153, 154, 157, 158, 51, 52, 163, 164, 165, 166,
  169, 170, 132, 131, 155, 156, 173, 174, 57, 58, 175, 176, 144, 143, 141, 142, 171, 172, 168, 167, 177, 178, 90, 89, 113, 114, 69, 70, 151, 152,
  107, 108, 73, 74, 81, 82, 145, 146, 127, 128, 121, 122, 179, 180, 77, 78, 94, 93, 161, 162, 104, 103, 159, 160, 124, 123, 135, 136, 117, 118,
  181, 182
]);;

G := Group(zhou3g1, zhou3g2);;
Print("Group generated, order = ", Size(G), "\n");
Print("Transitive: ", IsTransitive(G, [1..182]), "\n");

# Find block systems
# Use Blocks() which returns a list of blocks
blockSys := Blocks(G, [1..182]);;
Print("Default block system: ", Length(blockSys), " blocks of size ", Length(blockSys[1]), "\n");

# If default gives a non-trivial system, check if it's size 2
if Length(blockSys[1]) <> 2 then
  # Try to find size-2 blocks by looking for block containing 1 of size 2
  # The subgroup lattice: S_3 < D_12 < PSL(2,13)
  # D_12/S_3 = Z_2, so blocks should have size 2
  # Try AllBlocks from GRAPE or manual approach
  Print("Looking for size-2 blocks...\n");
  # A block of size 2 through point 1: find which point pairs with 1
  # under the D_12 action. The kernel of the projection S_3\G → D_12\G
  # acts within fibres.
  stab1 := Stabilizer(G, 1);;
  Print("Stabilizer(1) order: ", Size(stab1), "\n");
  # stab1 has order 6 (= |S_3|). Elements of D_12\S_3 permute within fibres.
  # Look for normalizer
  norm := Normalizer(G, stab1);;
  Print("Normalizer(Stab(1)) order: ", Size(norm), "\n");
  # norm should have order 12 (= |D_12|)
  # The extra element of D_12 \ S_3 maps 1 to its fibre partner
  for g in norm do
    if not g in stab1 then
      partner := 1^g;
      if partner <> 1 then
        Print("Fibre partner of 1: ", partner, "\n");
        blockSys := Blocks(G, [1..182], [1, partner]);;
        Print("Block system: ", Length(blockSys), " blocks of size ", Length(blockSys[1]), "\n");
        break;
      fi;
    fi;
  od;
fi;

# Build covering map: vertex v maps to the index of its block
coverMap := ListWithIdenticalEntries(182, 0);;
for i in [1..Length(blockSys)] do
  for v in blockSys[i] do
    coverMap[v] := i;
  od;
od;

# Build the quotient graph: blocks i,j adjacent iff some vertex in block i
# is adjacent to some vertex in block j
# First we need the Zhou-3 edge list to check adjacency
# Read edges from the Lean file format
zhou3Edges := [
  [1, 98], [1, 120], [1, 132],
  [2, 109], [2, 127], [2, 170],
  [3, 98], [3, 68], [3, 105],
  [4, 18], [4, 142], [4, 164],
  [5, 43], [5, 140], [5, 151],
  [6, 103], [6, 149], [6, 170],
  [7, 126], [7, 134], [7, 161],
  [8, 56], [8, 83], [8, 122],
  [9, 93], [9, 119], [9, 155],
  [10, 48], [10, 67], [10, 173]
];
# ... this approach is getting unwieldy. Let me use the graph from the group action instead.

# Build Zhou-3 adjacency from the orbital
stab1 := Stabilizer(G, 1);;
orbs := Orbits(stab1, [1..182]);;
Print("Stabilizer orbits: ", List(orbs, Length), "\n");
cubicOrb := First(orbs, o -> Length(o) = 3);;
Print("Cubic orbital: ", cubicOrb, "\n");

# The quotient graph on the block system
# Block i ~ Block j iff ∃ u ∈ block_i, v ∈ block_j with u ~ v (in orbital)
quotEdges := Set([]);;
for i in [1..Length(blockSys)] do
  for u in blockSys[i] do
    # neighbors of u via the orbital
    for n in cubicOrb do
      g := RepresentativeAction(G, 1, u);
      if g <> fail then
        w := n^g;
        j := coverMap[w];
        if i <> j then
          AddSet(quotEdges, Set([i, j]));
        fi;
      fi;
    od;
  od;
od;

Print("Quotient edges: ", Length(quotEdges), "\n");
Print("Expected for Zhou-6 (91 vertices, 6-regular): 273 edges\n");

# Now we need to find the isomorphism between quotient graph (on blocks 1..91)
# and zhou6Graph (on Fin 91 in Lean, vertices 0..90)
#
# Build the quotient graph as a GRAPE graph
quotAdj := function(x, y) return Set([x, y]) in quotEdges; end;;
quotGrp := Action(G, blockSys, OnSets);;
quotGraph := Graph(quotGrp, [1..Length(blockSys)], OnPoints, quotAdj, true);;

# Build Zhou-6 graph from Lean edge data
# Zhou-6 edges (0-indexed in Lean, 1-indexed here)
zhou6EdgeList := [
  [1, 9], [1, 45], [1, 49], [1, 64], [1, 79], [1, 85],
  [2, 6], [2, 17], [2, 28], [2, 43], [2, 44], [2, 72],
  [3, 13], [3, 25], [3, 58], [3, 68], [3, 75], [3, 86],
  [4, 20], [4, 55], [4, 62], [4, 69], [4, 76], [4, 84],
  [5, 10], [5, 38], [5, 58], [5, 72], [5, 78], [5, 88],
  [6, 18], [6, 29], [6, 46], [6, 52], [6, 73],
  [7, 15], [7, 47], [7, 55], [7, 78], [7, 83], [7, 87],
  [8, 10], [8, 26], [8, 35], [8, 47], [8, 51], [8, 60],
  [9, 14], [9, 31], [9, 41], [9, 82], [9, 90],
  [10, 21], [10, 53], [10, 70],
  [11, 14], [11, 18], [11, 47], [11, 59], [11, 72], [11, 89],
  [12, 32], [12, 36], [12, 48], [12, 67], [12, 76], [12, 80],
  [13, 19], [13, 52], [13, 55], [13, 60], [13, 79],
  [14, 46], [14, 71], [14, 88],
  [15, 22], [15, 39], [15, 62], [15, 84], [15, 88],
  [16, 22], [16, 31], [16, 69], [16, 70], [16, 78], [16, 85],
  [17, 23], [17, 46], [17, 67], [17, 79], [17, 90],
  [18, 64], [18, 71], [18, 87],
  [19, 30], [19, 44], [19, 61], [19, 81], [19, 90],
  [20, 35], [20, 37], [20, 58], [20, 75], [20, 79],
  [21, 24], [21, 25], [21, 54], [21, 59], [21, 68],
  [22, 41], [22, 56], [22, 63], [22, 86],
  [23, 33], [23, 49], [23, 51], [23, 61], [23, 84],
  [24, 28], [24, 48], [24, 62], [24, 82], [24, 91],
  [25, 30], [25, 39], [25, 64], [25, 87],
  [26, 30], [26, 37], [26, 40], [26, 66], [26, 89],
  [27, 29], [27, 42], [27, 45], [27, 53], [27, 57], [27, 63],
  [28, 56], [28, 63], [28, 82], [28, 90],
  [29, 69], [29, 71], [29, 77], [29, 87],
  [30, 57], [30, 65], [30, 80],
  [31, 38], [31, 42], [31, 50], [31, 75], [31, 83],
  [32, 43], [32, 50], [32, 54], [32, 86], [32, 91],
  [33, 39], [33, 65], [33, 73], [33, 76], [33, 82],
  [34, 46], [34, 48], [34, 50], [34, 57], [34, 81], [34, 85],
  [35, 44], [35, 52], [35, 61], [35, 74],
  [36, 43], [36, 45], [36, 66], [36, 77], [36, 83],
  [37, 43], [37, 53], [37, 56], [37, 67], [37, 74],
  [38, 48], [38, 60], [38, 64], [38, 71], [38, 73],
  [39, 70], [39, 72], [39, 89],
  [40, 42], [40, 51], [40, 68], [40, 74], [40, 78],
  [41, 54], [41, 59], [41, 65], [41, 66], [41, 83],
  [42, 44], [42, 73], [42, 85],
  [43, 89],
  [44, 61], [44, 76],
  [45, 53], [45, 60], [45, 70], [45, 80],
  [46, 77], [46, 88],
  [47, 57], [47, 66], [47, 69], [47, 81],
  [48, 68], [48, 74],
  [49, 51], [49, 58], [49, 62], [49, 91],
  [50, 52], [50, 63], [50, 65], [50, 80],
  [51, 67], [51, 77], [51, 84],
  [52, 59], [52, 85], [52, 90],
  [53, 75], [53, 81],
  [54, 56], [54, 86], [54, 88],
  [55, 57], [55, 68], [55, 75], [55, 91],
  [56, 73], [56, 80], [56, 91],
  [57, 70], [57, 77],
  [58, 67], [58, 81], [58, 83],
  [59, 63], [59, 74], [59, 78],
  [60, 66], [60, 82], [60, 86],
  [61, 75], [61, 87], [61, 89],
  [62, 65], [62, 71], [62, 79],
  [63, 69], [63, 76],
  [64, 70], [64, 72], [64, 77],
  [65, 74], [65, 90],
  [66, 84], [66, 87],
  [67, 71], [67, 85],
  [68, 81],
  [69, 86],
  [70, 91],
  [71, 82],
  [72, 74],
  [73, 76],
  [74, 79],
  [75, 81],
  [75, 85],
  [76, 84],
  [76, 85],
  [77, 81],
  [78, 89],
  [83, 86]
];;

zhou6Adj := function(x, y)
  return [x,y] in zhou6EdgeList or [y,x] in zhou6EdgeList;
end;;

zhou6GraphGRAPE := Graph(Group(()), [1..91], OnPoints, zhou6Adj, true);;

Print("Quotient graph: ", OrderGraph(quotGraph), " vertices, ",
  Size(UndirectedEdges(quotGraph)), " edges\n");
Print("Zhou-6 graph: ", OrderGraph(zhou6GraphGRAPE), " vertices, ",
  Size(UndirectedEdges(zhou6GraphGRAPE)), " edges\n");

# Find isomorphism between quotient and Zhou-6
iso := GraphIsomorphism(quotGraph, zhou6GraphGRAPE);;
if iso <> fail then
  Print("Isomorphism found!\n");

  # Build the final covering map: vertex v (1-indexed) in Zhou-3
  # maps to iso(block_index(v)) (1-indexed) in Zhou-6
  finalCoverMap := List([1..182], v -> v^iso);;
  # Wait, iso acts on block indices. Let me fix this.
  finalCoverMap := List([1..182], v -> coverMap[v]^iso);;

  # Print 0-indexed for Lean
  Print("\n-- Final covering map (0-indexed for Lean) --\n");
  Print("#[");
  for i in [1..182] do
    if i > 1 then Print(", "); fi;
    Print(finalCoverMap[i] - 1);
  od;
  Print("]\n");

  # Verify homomorphism using the orbital
  isHom := true;
  for v in [1..182] do
    g := RepresentativeAction(G, 1, v);
    if g <> fail then
      for n in cubicOrb do
        w := n^g;
        fv := finalCoverMap[v];
        fw := finalCoverMap[w];
        if fv <> fw and not zhou6Adj(fv, fw) then
          Print("FAILED at edge ", v, " ~ ", w, "\n");
          isHom := false;
        fi;
      od;
    fi;
  od;
  Print("Graph homomorphism verified: ", isHom, "\n");

  # Verify fibre sizes
  fibres := List([1..91], j -> Length(Filtered([1..182], i -> finalCoverMap[i] = j)));
  Print("All fibres size 2: ", ForAll(fibres, f -> f = 2), "\n");
else
  Print("No isomorphism found!\n");
  # Try matching manually by checking degrees
  Print("Quotient degrees: ", List([1..5], v -> Length(Adjacency(quotGraph, v))), "\n");
  Print("Zhou-6 degrees: ", List([1..5], v -> Length(Adjacency(zhou6GraphGRAPE, v))), "\n");
fi;
