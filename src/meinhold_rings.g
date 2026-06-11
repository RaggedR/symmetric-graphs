# Find the best free cyclic subgroup of Aut(Meinhold-72) for a ring layout.
# A cyclic <g> acts freely iff no nontrivial power of g fixes a vertex;
# then the 72 vertices fall into 72/|g| generator-ordered rings of size |g|.
# Run from src/: gap --quitonbreak meinhold_rings.g < /dev/null

LoadPackage("grape");;
Read("symmetric-layout.g");;

rec72 := GraphLibrary.meinhold72();;
gamma := Graph(Group(()), [1..rec72.n], OnPoints,
               function(x,y) return y in rec72.adj[x]; end, true);;
A := AutomorphismGroup(gamma);;
Print("Aut order = ", Size(A), "\n");

# Collect one representative per conjugacy class, check free action
best := fail;;
bestOrd := 1;;
for c in ConjugacyClasses(A) do
  g := Representative(c);
  m := Order(g);
  if m > 1 and 72 mod m = 0 then
    # free <=> every nontrivial power is fixed-point-free
    free := true;
    for k in [1..m-1] do
      if NrMovedPoints(g^k) < 72 then free := false; break; fi;
    od;
    if free and m > bestOrd then
      bestOrd := m; best := g;
    fi;
    if free then
      Print("free cyclic of order ", m, " (", 72/m, " rings)\n");
    fi;
  fi;
od;

if best = fail then
  Print("NO free cyclic subgroup found\n");
else
  Print("\nBEST: order ", bestOrd, " -> ", 72/bestOrd, " rings of ", bestOrd, "\n");
  # Export generator-ordered orbits
  seen := BlistList([1..72], []);;
  Print("ORBITS := [\n");
  for v in [1..72] do
    if not seen[v] then
      orb := [v];
      seen[v] := true;
      w := v^best;
      while w <> v do
        Add(orb, w); seen[w] := true; w := w^best;
      od;
      Print("  ", orb, ",\n");
    fi;
  od;
  Print("];\n");
fi;
