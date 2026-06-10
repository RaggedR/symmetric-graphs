# Investigate: are the Zhou-3 quotient and Zhou-6 truly non-isomorphic?
# If so, we need to regenerate Zhou-6 from the correct D_12.

LoadPackage("GRAPE");;

# Zhou-3 generators (from Lean)
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

G3 := Group(zhou3g1, zhou3g2);;

# Zhou-6 generators (from Lean)
zhou6g1 := PermList([
  2, 4, 6, 8, 9, 11, 13, 14, 16, 18, 19, 21, 7, 22, 24, 26, 28, 27, 31, 33, 35, 1, 38, 40, 42, 43, 45, 17, 48, 49,
  50, 52, 54, 32, 25, 56, 57, 59, 61, 62, 63, 47, 66, 65, 58, 70, 12, 41, 67, 3, 72, 74, 76, 77, 78, 80, 39, 71, 82, 83,
  60, 84, 55, 86, 46, 5, 36, 64, 44, 87, 10, 34, 23, 51, 88, 85, 89, 29, 15, 30, 20, 90, 37, 79, 53, 75, 69, 91, 81, 73,
  68
]);;

zhou6g2 := PermList([
  3, 5, 7, 6, 10, 12, 1, 15, 17, 2, 20, 4, 18, 23, 25, 27, 29, 30, 32, 34, 36, 37, 39, 41, 8, 44, 46, 47, 9, 13,
  51, 53, 55, 11, 56, 54, 58, 60, 14, 52, 64, 65, 67, 68, 69, 16, 71, 57, 48, 42, 73, 75, 19, 21, 79, 81, 49, 22, 59, 70,
  45, 74, 85, 24, 50, 80, 72, 26, 61, 38, 28, 43, 31, 82, 40, 86, 77, 66, 33, 78, 35, 62, 87, 63, 84, 90, 91, 88, 89, 76,
  83
]);;

G6 := Group(zhou6g1, zhou6g2);;

Print("G3 order: ", Size(G3), "  G6 order: ", Size(G6), "\n");

# Check: are G3 and G6 both PSL(2,13)?
Print("G3 isomorphic to G6: ", IdGroup(G3) = IdGroup(G6), "\n");

# The block system of G3 acting on [1..182]
blockSys := Blocks(G3, [1..182]);;
Print("G3 block system: ", Length(blockSys), " blocks of size ", Length(blockSys[1]), "\n");

# Stabilizer subgroups
stab3 := Stabilizer(G3, 1);;
Print("Stab_G3(1) order: ", Size(stab3), " structure: ", StructureDescription(stab3), "\n");

stab6 := Stabilizer(G6, 1);;
Print("Stab_G6(1) order: ", Size(stab6), " structure: ", StructureDescription(stab6), "\n");

# Find which block contains vertex 1
block_of_1 := First(blockSys, b -> 1 in b);;
Print("Block containing 1: ", block_of_1, "\n");

# The setwise stabilizer of the block = the D_12 for the quotient
blockStab := Stabilizer(G3, Set(block_of_1), OnSets);;
Print("Block stabilizer order: ", Size(blockStab), " structure: ", StructureDescription(blockStab), "\n");

# Check: is stab3 ≤ blockStab?
Print("Stab(1) ≤ Block_stab: ", IsSubgroup(blockStab, stab3), "\n");

# Now let's check the girth of both 6-regular graphs
# Build the quotient graph
coverMap := ListWithIdenticalEntries(182, 0);;
for i in [1..Length(blockSys)] do
  for v in blockSys[i] do coverMap[v] := i; od;
od;

cubicOrb := First(Orbits(stab3, [1..182]), o -> Length(o) = 3);;
quotEdgeSet := Set([]);;
for v in [1..182] do
  g := RepresentativeAction(G3, 1, v);
  if g <> fail then
    for n in cubicOrb do
      w := n^g;
      if coverMap[v] <> coverMap[w] then
        AddSet(quotEdgeSet, Set([coverMap[v], coverMap[w]]));
      fi;
    od;
  fi;
od;

# Zhou-6 edge set
Read("/tmp/zhou6_edges.g");
zhou6AdjSet := Set(List(zhou6EdgeList, p -> Set(p)));;

# Compute girths
quotAdj := function(x, y) return Set([x, y]) in quotEdgeSet; end;;
zhou6Adj := function(x, y) return Set([x, y]) in zhou6AdjSet; end;;

# Use GRAPE for girth
quotGrape := Graph(Group(()), [1..91], OnPoints, quotAdj, true);;
zhou6Grape := Graph(Group(()), [1..91], OnPoints, zhou6Adj, true);;

Print("\nQuotient graph girth: ", Girth(quotGrape), "\n");
Print("Zhou-6 graph girth:   ", Girth(zhou6Grape), "\n");

# Also check if they're REALLY non-isomorphic
Print("Isomorphic check: ", IsIsomorphicGraph(quotGrape, zhou6Grape), "\n");

# Check if there are multiple conjugacy classes of D_12 in PSL(2,13)
psl := G3;  # This is PSL(2,13) acting on 182 points
# Conjugacy classes of subgroups of order 12
subs12 := Filtered(ConjugacyClassesSubgroups(psl), c -> Size(Representative(c)) = 12);;
Print("\nConjugacy classes of order-12 subgroups: ", Length(subs12), "\n");
for c in subs12 do
  Print("  ", StructureDescription(Representative(c)), " (", Length(c), " conjugates)\n");
od;

# Check which D_12 class our block stabilizer belongs to
for i in [1..Length(subs12)] do
  if blockStab in subs12[i] then
    Print("Block stabilizer is in class ", i, "\n");
  fi;
od;
