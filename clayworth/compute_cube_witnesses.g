# Cube Sabidussi witnesses
# Vertices: 1..8 (GAP), binary: 1=000, 2=001, 3=010, 4=011, 5=100, 6=101, 7=110, 8=111

CubeAdj := function(u, v)
  return Number([0,1,2], b -> QuoInt((u-1), 2^b) mod 2 <> QuoInt((v-1), 2^b) mod 2) = 1;
end;

# Rotation group S_4 generators:
# g1: 90° rotation around z-axis: 1->3->7->5 and 2->4->8->6
g1 := (1,3,7,5)(2,4,8,6);
# g2: 120° rotation around (1,1,1) diagonal: fixes 1,8, cycles 2->3->5 and 4->7->6
g2 := (2,3,5)(4,7,6);

# Verify adjacency preservation
for g in [g1, g1^-1, g2, g2^-1] do
  for u in [1..8] do
    for v in [u+1..8] do
      if CubeAdj(u,v) <> CubeAdj(u^g, v^g) then
        Print("FAIL!\n"); QuitGap(0);
      fi;
    od;
  od;
od;
Print("All generators preserve adjacency.\n");

G := Group(g1, g2);
Print("|G| = ", Size(G), "\n");
Print("Transitive: ", IsTransitive(G, [1..8]), "\n");

# BFS: word encoding 0=g1, 1=g1^-1, 2=g2, 3=g2^-1
allGens := [g1, g1^(-1), g2, g2^(-1)];
words := List([1..8], i -> fail);
words[1] := [];
queue := [1];
while Length(queue) > 0 do
  v := Remove(queue, 1);
  for gi in [1..4] do
    nv := v ^ allGens[gi];
    if words[nv] = fail then
      words[nv] := Concatenation(words[v], [gi-1]);
      Add(queue, nv);
    fi;
  od;
od;

# Print Lean data (0-indexed)
Print("\n-- LEAN DATA (0-indexed) --\n");
Print("g1Fwd: "); for v in [1..8] do Print(v^g1-1, " "); od; Print("\n");
Print("g1Inv: "); for v in [1..8] do Print(v^(g1^-1)-1, " "); od; Print("\n");
Print("g2Fwd: "); for v in [1..8] do Print(v^g2-1, " "); od; Print("\n");
Print("g2Inv: "); for v in [1..8] do Print(v^(g2^-1)-1, " "); od; Print("\n");
Print("Words:\n");
for v in [1..8] do Print("  ", v-1, ": ", words[v], "\n"); od;

QuitGap(0);
