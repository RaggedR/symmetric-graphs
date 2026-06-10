# Compute Sabidussi witness data for named graphs
# Output: generators as permutations + BFS witness words

# Helper: BFS from vertex `base` using generators and their inverses
BFSWitnesses := function(gens, n, base)
  local words, visited, queue, front, i, g, gi, v, nv;
  words := List([1..n], i -> fail);
  words[base] := [];
  visited := BlistList([1..n], [base]);
  queue := [base];
  front := 1;
  while front <= Length(queue) do
    v := queue[front];
    front := front + 1;
    for gi in [1..Length(gens)] do
      for g in [gens[gi], gens[gi]^(-1)] do
        nv := nv^g;  # wrong, should be image
      od;
    od;
  od;
  return words;
end;

# Actually, let me write this properly for GAP

# -- CUBE GRAPH --
# Vertices 1..8 (GAP is 1-indexed, Lean is 0-indexed)
# Cube edges (Hamming distance 1 on 3-bit strings):
# 000=1, 001=2, 010=3, 011=4, 100=5, 101=6, 110=7, 111=8
cubeEdges := [[1,2],[1,3],[1,5],[2,4],[2,6],[3,4],[3,7],[4,8],[5,6],[5,7],[6,8],[7,8]];
cubeGamma := Graph(Group(()), [1..8], OnPoints, function(x,y) return [x,y] in cubeEdges or [y,x] in cubeEdges; end, true);

cubeAut := AutGroupGraph(cubeGamma);
Print("Cube: |Aut| = ", Size(cubeAut), "\n");

# Get two generators
cubeGens := GeneratorsOfGroup(cubeAut);
Print("Cube generators:\n");
for g in cubeGens do
  Print("  ", ListPerm(g, 8), "\n");
od;

# BFS from vertex 1
BFSFromBase := function(gens, n, base)
  local words, visited, queue, v, gi, sign, perm, nv, word;
  words := [];
  for v in [1..n] do words[v] := fail; od;
  words[base] := [];
  visited := [base];
  queue := [base];
  while Length(queue) > 0 do
    v := Remove(queue, 1);
    for gi in [1..Length(gens)] do
      # Forward: generator gi (encoded as 2*(gi-1))
      nv := v ^ gens[gi];
      if words[nv] = fail then
        words[nv] := Concatenation(words[v], [2*(gi-1)]);
        Add(visited, nv);
        Add(queue, nv);
      fi;
      # Inverse: generator gi^-1 (encoded as 2*(gi-1)+1)
      nv := v ^ (gens[gi]^(-1));
      if words[nv] = fail then
        words[nv] := Concatenation(words[v], [2*(gi-1)+1]);
        Add(visited, nv);
        Add(queue, nv);
      fi;
    od;
  od;
  return words;
end;

cubeWords := BFSFromBase(cubeGens, 8, 1);
Print("Cube BFS words (1-indexed, gen indices 0-based):\n");
for v in [1..8] do
  Print("  vertex ", v, ": ", cubeWords[v], "\n");
od;

# Print as Lean arrays (0-indexed)
Print("\n-- Lean output (0-indexed) --\n");
Print("-- Generators as forward arrays:\n");
for i in [1..Length(cubeGens)] do
  Print("gen", i, "Fwd := #[");
  for v in [1..8] do
    Print(v^cubeGens[i] - 1);
    if v < 8 then Print(", "); fi;
  od;
  Print("]\n");
  Print("gen", i, "Inv := #[");
  for v in [1..8] do
    Print(v^(cubeGens[i]^(-1)) - 1);
    if v < 8 then Print(", "); fi;
  od;
  Print("]\n");
od;
Print("-- BFS witness words:\n");
for v in [1..8] do
  Print("  ", v-1, " -> ", cubeWords[v], "\n");
od;
