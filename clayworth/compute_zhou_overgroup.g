LoadPackage("GRAPE");;

# Zhou-3 generators from Lean (0-indexed), converted to 1-indexed GAP permutations
# g1 forward map
g1f := [120,121,114,109,118,45,13,116,115,101,110,76,104,49,12,117,105,111,112,15,103,11,113,108,46,106,48,72,73,8,74,75,44,47,102,100,119,14,107,165,174,95,149,156,34,179,97,158,178,136,85,161,129,31,147,130,170,171,138,70,36,40,82,7,1,71,57,122,140,153,181,38,172,33,163,150,135,59,26,87,157,61,164,146,166,173,148,65,96,32,134,94,77,81,5,66,64,132,133,24,155,58,144,90,128,168,167,145,92,176,68,175,98,9,67,52,55,137,177,20,35,151,154,141,88,27,159,143,126,160,86,89,56,4,60,51,127,17,162,142,25,91,28,169,131,79,41,93,63,62,50,124,29,30,84,80,69,0,139,37,152,180,22,43,39,42,2,10,78,6,99,83,3,53,54,19,125,18,16,123,21,23];

# g2 forward map
g2f := [150,163,172,135,115,33,94,166,157,64,133,156,59,26,87,110,132,101,134,77,65,96,148,173,32,116,13,114,149,95,174,165,24,5,61,66,146,81,164,175,167,178,177,137,90,155,158,128,144,58,136,104,52,98,176,68,117,168,49,12,76,34,179,97,9,20,35,67,55,151,145,92,142,89,159,162,60,19,78,140,125,37,169,180,152,113,119,14,120,73,44,153,71,181,6,29,21,63,53,99,126,17,127,129,51,111,143,160,161,112,15,105,109,85,27,4,25,56,141,86,88,154,131,139,138,80,100,102,47,103,130,122,16,10,18,3,50,43,124,123,79,118,72,106,48,70,36,170,22,28,0,69,84,91,121,45,11,8,46,74,107,108,75,1,38,31,7,40,57,82,147,171,2,23,30,39,54,42,41,62,83,93];

# Convert to 1-indexed
g1f_gap := List(g1f, x -> x + 1);
g2f_gap := List(g2f, x -> x + 1);

g1 := PermList(g1f_gap);
g2 := PermList(g2f_gap);

if g1 = fail then Print("g1 is fail!\n"); fi;
if g2 = fail then Print("g2 is fail!\n"); fi;

Print("g1 order: ", Order(g1), "\n");
Print("g2 order: ", Order(g2), "\n");

G := Group(g1, g2);
Print("Group order: ", Size(G), "\n");

# The block containing vertex 1 (= vertex 0 in Lean) is {1, 138}
# because zhouBlockData[0] = 0 and zhouBlockData[137] = 0
setstab := Stabilizer(G, Set([1, 138]), OnSets);
Print("Setwise stab order: ", Size(setstab), "\n");
Print("Structure: ", StructureDescription(setstab), "\n");

stab := Stabilizer(G, 1);
Print("Point stab order: ", Size(stab), "\n");
Print("Point stab <= Setwise stab? ", IsSubgroup(setstab, stab), "\n");

# Find an involution in setstab that swaps 1 <-> 138
swappers := Filtered(Elements(setstab), g -> 1^g = 138 and Order(g) = 2);
Print("Number of involutions swapping 1<->138: ", Length(swappers), "\n");

if Length(swappers) > 0 then
  s := swappers[1];
  Print("Swapper: ", s, "\n");
  Print("Swapper order: ", Order(s), "\n");

  # Output 0-indexed for Lean
  imgs := List([1..182], i -> i^s - 1);
  Print("SWAPPER := ", imgs, ";\n");
fi;

# Find 2 generators for setstab
sg := SmallGeneratingSet(setstab);
Print("\nSetstab generators (", Length(sg), "):\n");
for i in [1..Length(sg)] do
  Print("  gen", i, " order=", Order(sg[i]), " maps 1->", 1^sg[i], " 138->", 138^sg[i], "\n");
  imgs := List([1..182], j -> j^sg[i] - 1);
  Print("KGEN", i, " := ", imgs, ";\n");
od;

# Verify: connection set D for Zhou-3
# D = {g in G : zhouGraph.Adj 0 (g(0))}
# Neighbors of vertex 1 (=0) in Zhou-3 are: 4 (=3), 38 (=37), 69 (=68)
# So connectionSet elements send 1 to 4, 38, or 69
connElts := Filtered(Elements(G), g -> 1^g in [4, 38, 69]);
Print("\n|connectionSet| = ", Length(connElts), "\n");

# Check D ∩ K = ∅
dCapK := Intersection(Set(connElts), Set(Elements(setstab)));
Print("|D ∩ K| = ", Length(dCapK), "\n");

quit;
