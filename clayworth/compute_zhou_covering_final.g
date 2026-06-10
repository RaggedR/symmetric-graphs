# Compute the Zhou-3 → Zhou-6 covering map using Lean vertex numbering.
# Uses edge data extracted directly from the Lean files.

LoadPackage("GRAPE");;

# Zhou-3 generators (from Zhou3Sabidussi.lean, 0-indexed → 1-indexed)
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
Print("Group order = ", Size(G), "\n");

# Block system of size 2
blockSys := Blocks(G, [1..182]);;
Print("Blocks: ", Length(blockSys), " of size ", Length(blockSys[1]), "\n");

# Covering map: vertex → block index
coverMap := ListWithIdenticalEntries(182, 0);;
for i in [1..Length(blockSys)] do
  for v in blockSys[i] do
    coverMap[v] := i;
  od;
od;

# Cubic orbital
stab1 := Stabilizer(G, 1);;
cubicOrb := First(Orbits(stab1, [1..182]), o -> Length(o) = 3);;

# Quotient edges
quotEdgeSet := Set([]);;
for v in [1..182] do
  g := RepresentativeAction(G, 1, v);
  if g <> fail then
    for n in cubicOrb do
      w := n^g;
      if coverMap[v] <> coverMap[w] then
        AddSet(quotEdgeSet, Set([coverMap[v], coverMap[w]]));
      fi;
    od;
  fi;
od;
Print("Quotient edges: ", Length(quotEdgeSet), "\n");

# Zhou-6 edges (extracted from Lean, 1-indexed)
Read("/tmp/zhou6_edges.g");
zhou6AdjSet := Set(List(zhou6EdgeList, p -> Set(p)));;
Print("Zhou-6 edges: ", Length(zhou6AdjSet), "\n");

# Build GRAPE graphs
quotAdj := function(x, y) return Set([x, y]) in quotEdgeSet; end;;
quotGrp := Action(G, blockSys, OnSets);;
quotGrapeGraph := Graph(quotGrp, [1..91], OnPoints, quotAdj, true);;

zhou6Adj := function(x, y) return Set([x, y]) in zhou6AdjSet; end;;
zhou6GrapeGraph := Graph(Group(()), [1..91], OnPoints, zhou6Adj, true);;

Print("Quotient: ", OrderGraph(quotGrapeGraph), "v, ",
  Size(UndirectedEdges(quotGrapeGraph)), "e\n");
Print("Zhou-6:   ", OrderGraph(zhou6GrapeGraph), "v, ",
  Size(UndirectedEdges(zhou6GrapeGraph)), "e\n");

# Check isomorphism
Print("Isomorphic: ", IsIsomorphicGraph(quotGrapeGraph, zhou6GrapeGraph), "\n");

iso := GraphIsomorphism(quotGrapeGraph, zhou6GrapeGraph);;
if iso <> fail then
  Print("Isomorphism found!\n");
  finalMap := List([1..182], v -> coverMap[v]^iso);;

  # Verify homomorphism
  isHom := true;
  for v in [1..182] do
    g := RepresentativeAction(G, 1, v);
    if g <> fail then
      for n in cubicOrb do
        w := n^g;
        fv := finalMap[v];
        fw := finalMap[w];
        if fv <> fw and not (Set([fv, fw]) in zhou6AdjSet) then
          Print("FAILED: ", v, "~", w, " → ", fv, "~", fw, "\n");
          isHom := false;
        fi;
      od;
    fi;
  od;
  Print("Homomorphism verified: ", isHom, "\n");
  Print("All fibres size 2: ",
    ForAll([1..91], j -> Length(Filtered([1..182], i -> finalMap[i] = j)) = 2), "\n");

  # Output for Lean (0-indexed)
  Print("\n-- Lean Array (0-indexed) --\n");
  Print("private def zhou3ToZhou6Data : Array (Fin 91) := #[\n  ");
  for i in [1..182] do
    if i > 1 then Print(", "); fi;
    if i > 1 and (i-1) mod 20 = 0 then Print("\n  "); fi;
    Print(finalMap[i] - 1);
  od;
  Print("\n]\n");
else
  Print("No isomorphism — graphs are not isomorphic.\n");
  # Debug: check degree sequences
  quotDegs := SortedList(List([1..91], v -> Length(Adjacency(quotGrapeGraph, v))));;
  zhou6Degs := SortedList(List([1..91], v -> Length(Adjacency(zhou6GrapeGraph, v))));;
  Print("Quotient degree seq: ", Set(quotDegs), "\n");
  Print("Zhou-6 degree seq:   ", Set(zhou6Degs), "\n");
fi;
