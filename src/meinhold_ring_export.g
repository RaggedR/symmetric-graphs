# Export generator-ordered cyclic orbits of Meinhold-72 for ring layouts.
# Candidates: C7 (10 rings of 7 + 2 fixed), C8 (6x8 + 6x4), C12 (4x12 + 4x6).
# Output: /tmp/meinhold_orbits.json
# Run: gap --quitonbreak -q meinhold_ring_export.g < /dev/null

LoadPackage("grape");;
Read("/tmp/meinhold_adj.g");;

gamma := Graph(Group(()), [1..72], OnPoints,
               function(x,y) return y in MEINHOLD_ADJ[x]; end, true);;
A := AutomorphismGroup(gamma);;

# Pick, per target order, the class rep whose orbit profile has NO fixed
# points except for order 7 (where 2 fixed points are unavoidable).
exportOrbits := function(g)
  local seen, out, v, orb, w;
  seen := BlistList([1..72], []);
  out := [];
  for v in [1..72] do
    if not seen[v] then
      orb := [v]; seen[v] := true; w := v^g;
      while w <> v do Add(orb, w); seen[w] := true; w := w^g; od;
      Add(out, orb);
    fi;
  od;
  return out;
end;;

best := rec();;
for c in ConjugacyClasses(A) do
  g := Representative(c);
  m := Order(g);
  if m in [7, 8, 12] then
    sizes := Collected(List(OrbitsDomain(Group(g), [1..72]), Length));
    # prefer profiles whose smallest orbit is largest (fewest degenerate rings)
    key := Concatenation("c", String(m));
    if not IsBound(best.(key)) or
       Minimum(List(sizes, s -> s[1])) > best.(key).minsize then
      best.(key) := rec(g := g, minsize := Minimum(List(sizes, s -> s[1])),
                        sizes := sizes);
    fi;
  fi;
od;

out := OutputTextFile("/tmp/meinhold_orbits.json", false);;
SetPrintFormattingStatus(out, false);;
PrintTo(out, "{\n");
first := true;;
for key in RecNames(best) do
  if not first then AppendTo(out, ",\n"); fi;
  first := false;
  AppendTo(out, "\"", key, "\": ", exportOrbits(best.(key).g));
  Print(key, ": profile ", best.(key).sizes, "\n");
od;
AppendTo(out, "\n}\n");
CloseStream(out);;
Print("written /tmp/meinhold_orbits.json\n");
