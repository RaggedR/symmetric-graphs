# Compute the Zhou-3 → Zhou-6 covering map
# Zhou-3 = PSL(2,13)/S_3 (182 vertices, cubic)
# Zhou-6 = PSL(2,13)/D_12 (91 vertices, 6-regular)
# The covering map is the natural projection S_3\G → D_12\G

LoadPackage("GRAPE");;

psl213 := PSL(2,13);;
Print("PSL(2,13) order = ", Size(psl213), "\n");

# Find S_3 and D_12 subgroups with S_3 ≤ D_12
subs := ConjugacyClassesSubgroups(psl213);;
s3subs := Filtered(subs, c -> Size(Representative(c)) = 6 and
  IdGroup(Representative(c)) = [6, 1]);;  # S_3
d12subs := Filtered(subs, c -> Size(Representative(c)) = 12 and
  IdGroup(Representative(c)) = [12, 4]);;  # D_12

s3sub := Representative(s3subs[1]);;
Print("S_3 subgroup found, order ", Size(s3sub), "\n");

# Find D_12 containing our S_3
d12sub := fail;
for c in d12subs do
  for g in Elements(c) do
    if IsSubgroup(g, s3sub) then
      d12sub := g;
      break;
    fi;
  od;
  if d12sub <> fail then break; fi;
od;
Print("D_12 subgroup found, order ", Size(d12sub), "\n");
Print("S_3 ≤ D_12: ", IsSubgroup(d12sub, s3sub), "\n");

# Compute the coset actions
zhou3Action := FactorCosetAction(psl213, s3sub);;
zhou3Perm := Image(zhou3Action);;
Print("Zhou-3 degree = ", NrMovedPoints(zhou3Perm), "\n");

zhou6Action := FactorCosetAction(psl213, d12sub);;
zhou6Perm := Image(zhou6Action);;
Print("Zhou-6 degree = ", NrMovedPoints(zhou6Perm), "\n");

# The covering map: for each S_3-coset gS_3, find the D_12-coset gD_12
# In terms of coset actions: coset i in Zhou-3 corresponds to gS_3 for some g.
# The image is the D_12-coset containing g, which is the orbit of the D_12
# action containing this g.

# Get coset representatives for S_3
s3cosets := RightCosets(psl213, s3sub);;
Print("Number of S_3-cosets: ", Length(s3cosets), "\n");

d12cosets := RightCosets(psl213, d12sub);;
Print("Number of D_12-cosets: ", Length(d12cosets), "\n");

# Build covering map: for S_3-coset i, find which D_12-coset contains it
coveringMap := [];
for i in [1..182] do
  # coset i in Zhou-3 is s3cosets[i] = s3sub * g_i
  rep := Representative(s3cosets[i]);
  # Find which D_12-coset this belongs to
  for j in [1..91] do
    if rep in d12cosets[j] then
      Add(coveringMap, j);
      break;
    fi;
  od;
od;

Print("Covering map computed, length = ", Length(coveringMap), "\n");

# Verify: each D_12-coset has exactly 2 preimages (since [D_12:S_3] = 2)
fibres := List([1..91], j -> Length(Filtered([1..182], i -> coveringMap[i] = j)));
Print("All fibres size 2: ", ForAll(fibres, f -> f = 2), "\n");

# Now we need to match vertex numbering with the Lean files.
# The Lean files use specific edge lists, so we need to verify the map
# is consistent with those graphs.

# Build Zhou-3 graph edges from the coset action
zhou3Gens := SmallGeneratingSet(zhou3Perm){[1..2]};;
stab1 := Stabilizer(zhou3Perm, 1);;
orbs := Orbits(stab1, [1..182]);;
cubicOrb := First(orbs, o -> Length(o) = 3);;
Print("Zhou-3 cubic orbital at vertex 1: ", cubicOrb, "\n");

# Generate all Zhou-3 edges from the orbital
zhou3Edges := Set([]);;
for v in [1..182] do
  for u in cubicOrb do
    # Map orbital to vertex v
    # Need a group element g taking 1 to v
    g := RepresentativeAction(zhou3Perm, 1, v);
    if g <> fail then
      w := u^g;
      AddSet(zhou3Edges, Set([v, w]));
    fi;
  od;
od;
Print("Zhou-3 edges computed: ", Length(zhou3Edges), "\n");

# Print the covering map (0-indexed for Lean)
Print("\n-- Covering map (0-indexed) --\n");
Print("coveringMap := [");
for i in [1..182] do
  if i > 1 then Print(", "); fi;
  Print(coveringMap[i] - 1);
od;
Print("]\n");

# Verify the map is a graph homomorphism:
# If {u,v} is an edge of Zhou-3, then {f(u), f(v)} should be an edge of Zhou-6
zhou6Gens := SmallGeneratingSet(zhou6Perm){[1..2]};;
stab1_6 := Stabilizer(zhou6Perm, 1);;
orbs6 := Orbits(stab1_6, [1..91]);;
sixOrb := First(orbs6, o -> Length(o) = 6);;
Print("Zhou-6 six-orbital at vertex 1: ", sixOrb, "\n");

zhou6Edges := Set([]);;
for v in [1..91] do
  for u in sixOrb do
    g := RepresentativeAction(zhou6Perm, 1, v);
    if g <> fail then
      w := u^g;
      AddSet(zhou6Edges, Set([v, w]));
    fi;
  od;
od;
Print("Zhou-6 edges computed: ", Length(zhou6Edges), "\n");

# Check homomorphism
isHom := true;
for e in zhou3Edges do
  el := Elements(e);
  img := Set([coveringMap[el[1]], coveringMap[el[2]]]);
  if img[1] <> img[2] and not (img in zhou6Edges) then
    Print("FAILED: edge ", e, " maps to non-edge ", img, "\n");
    isHom := false;
  fi;
od;
Print("Is graph homomorphism: ", isHom, "\n");
