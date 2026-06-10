# Generate correct Zhou-6 = quotient of Zhou-3 by S_3 ≤ D_12 block system
# The existing Zhou6Sabidussi.lean has Sab(PSL(2,13), A_4), not D_12!
# This script computes the correct Zhou-6 = Sab(PSL(2,13), D_12).

LoadPackage("GRAPE");;

# Zhou-3 generators (from Lean, 0→1 indexed)
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
Print("PSL(2,13) order: ", Size(G), "\n");

# Block system
blockSys := Blocks(G, [1..182]);;
Print("Block system: ", Length(blockSys), " blocks of size ", Length(blockSys[1]), "\n");

# Covering map
coverMap := ListWithIdenticalEntries(182, 0);;
for i in [1..Length(blockSys)] do
  for v in blockSys[i] do coverMap[v] := i; od;
od;

# The quotient action: G acts on blocks by OnSets
quotAction := ActionHomomorphism(G, blockSys, OnSets);;
quotGroup := Image(quotAction);;
Print("Quotient group order: ", Size(quotGroup), " degree: ", NrMovedPoints(quotGroup), "\n");

# Generators of quotient action (images of zhou3g1, zhou3g2)
quotGen1 := Image(quotAction, zhou3g1);;
quotGen2 := Image(quotAction, zhou3g2);;

# Verify transitivity
Print("Quotient transitive: ", IsTransitive(quotGroup, [1..91]), "\n");

# Compute the quotient graph edges
stab1 := Stabilizer(G, 1);;
cubicOrb := First(Orbits(stab1, [1..182]), o -> Length(o) = 3);;

quotEdges := [];;
quotEdgeSet := Set([]);;
for v in [1..182] do
  g := RepresentativeAction(G, 1, v);
  if g <> fail then
    for n in cubicOrb do
      w := n^g;
      bv := coverMap[v];
      bw := coverMap[w];
      if bv < bw then
        if not Set([bv, bw]) in quotEdgeSet then
          AddSet(quotEdgeSet, Set([bv, bw]));
          Add(quotEdges, [bv, bw]);
        fi;
      fi;
    od;
  fi;
od;

Print("Quotient edges: ", Length(quotEdges), "\n");

# Verify 6-regular
degs := List([1..91], v -> Length(Filtered(quotEdges, e -> v in e)));;
Print("Degree sequence: ", Set(degs), "\n");

# BFS witness words for the quotient
BFS := function(gens, n, start)
  local words, queue, visited, v, w, i, gen;
  words := ListWithIdenticalEntries(n, fail);
  words[start] := [];
  queue := [start];
  visited := BlistList([1..n], [start]);
  while Length(queue) > 0 do
    v := Remove(queue, 1);
    for i in [1..Length(gens)] do
      for gen in [gens[i], gens[i]^-1] do
        w := v^gen;
        if not visited[w] then
          visited[w] := true;
          if gen = gens[i] then
            words[w] := Concatenation(words[v], [2*(i-1)]);
          else
            words[w] := Concatenation(words[v], [2*i-1]);
          fi;
          Add(queue, w);
        fi;
      od;
    od;
  od;
  return words;
end;;

quotGens := [quotGen1, quotGen2];;
quotWords := BFS(quotGens, 91, 1);;
Print("BFS complete: ", ForAll(quotWords, w -> w <> fail), "\n");

# === OUTPUT FOR LEAN ===

Print("\n\n========== LEAN OUTPUT ==========\n\n");

# 1. Covering map (0-indexed)
Print("-- Covering map: Zhou-3 (Fin 182) → Zhou-6 (Fin 91) --\n");
Print("private def zhou3ToZhou6Data : Array (Fin 91) := #[\n  ");
for i in [1..182] do
  if i > 1 then Print(", "); fi;
  if i > 1 and (i-1) mod 20 = 0 then Print("\n  "); fi;
  Print(coverMap[i] - 1);
od;
Print("\n]\n\n");

# 2. Quotient edges (0-indexed)
Sort(quotEdges);;
Print("-- Zhou-6 edges (correct, from quotient) --\n");
Print("private def zhou6Edges : Array (Fin 91 × Fin 91) := #[\n");
for i in [1..Length(quotEdges)] do
  Print("  (", quotEdges[i][1]-1, ", ", quotEdges[i][2]-1, ")");
  if i < Length(quotEdges) then Print(","); fi;
  Print("\n");
od;
Print("]\n\n");

# 3. Quotient generators (0-indexed)
Print("-- Zhou-6 gen1 forward --\n");
Print("private def zhou6Gen1Fwd : Array (Fin 91) := #[\n  ");
for i in [1..91] do
  if i > 1 then Print(", "); fi;
  if i > 1 and (i-1) mod 20 = 0 then Print("\n  "); fi;
  Print(i^quotGen1 - 1);
od;
Print("\n]\n");

Print("-- Zhou-6 gen1 inverse --\n");
Print("private def zhou6Gen1Inv : Array (Fin 91) := #[\n  ");
for i in [1..91] do
  if i > 1 then Print(", "); fi;
  if i > 1 and (i-1) mod 20 = 0 then Print("\n  "); fi;
  Print(i^(quotGen1^-1) - 1);
od;
Print("\n]\n");

Print("-- Zhou-6 gen2 forward --\n");
Print("private def zhou6Gen2Fwd : Array (Fin 91) := #[\n  ");
for i in [1..91] do
  if i > 1 then Print(", "); fi;
  if i > 1 and (i-1) mod 20 = 0 then Print("\n  "); fi;
  Print(i^quotGen2 - 1);
od;
Print("\n]\n");

Print("-- Zhou-6 gen2 inverse --\n");
Print("private def zhou6Gen2Inv : Array (Fin 91) := #[\n  ");
for i in [1..91] do
  if i > 1 then Print(", "); fi;
  if i > 1 and (i-1) mod 20 = 0 then Print("\n  "); fi;
  Print(i^(quotGen2^-1) - 1);
od;
Print("\n]\n");

# 4. BFS witness words (0-indexed generators)
Print("\n-- Zhou-6 BFS witness words --\n");
Print("private def zhou6WitnessWords : Fin 91 → List (Fin 4) := fun v =>\n");
Print("  match v with\n");
for i in [1..91] do
  Print("  | ", i-1, " => ", quotWords[i], "\n");
od;
Print("\n");

# 5. Verify covering map properties
Print("\n-- Verification --\n");
# Homomorphism check
isHom := true;
for v in [1..182] do
  g := RepresentativeAction(G, 1, v);
  if g <> fail then
    for n in cubicOrb do
      w := n^g;
      fv := coverMap[v];
      fw := coverMap[w];
      if fv <> fw and not (Set([fv, fw]) in quotEdgeSet) then
        isHom := false;
      fi;
    od;
  fi;
od;
Print("Homomorphism: ", isHom, "\n");

# Fibre check
Print("All fibres size 2: ",
  ForAll([1..91], j -> Length(Filtered([1..182], i -> coverMap[i] = j)) = 2), "\n");

# Equivariance: coverMap commutes with generators
equivar := true;
for v in [1..182] do
  if coverMap[v^zhou3g1] <> coverMap[v]^quotGen1 then equivar := false; fi;
  if coverMap[v^zhou3g2] <> coverMap[v]^quotGen2 then equivar := false; fi;
od;
Print("Equivariant: ", equivar, "\n");
