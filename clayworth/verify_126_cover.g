LoadPackage("grape");;

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
stab192 := Stabilizer(G, 1);;

# Build Langer
orbs63 := Orbits(stab192, [1..63]);;
orb6_langer := First(orbs63, o -> Length(o) = 6);;
langer := EdgeOrbitsGraph(G, [[1, orb6_langer[1]]], 63);;
Print("Langer: V=63, 6-reg, girth=", Girth(langer), ", diam=", Diameter(langer), "\n");;

# K = SL(2,3):C4 inside H_192
subsOfStab := ConjugacyClassesSubgroups(stab192);;
K96 := First(List(subsOfStab, Representative),
    k -> Size(k) = 96 and StructureDescription(k) = "SL(2,3) : C4");;

act := FactorCosetAction(G, K96);;
im := Image(act);;
stab1 := Stabilizer(im, 1);;
orbs := Orbits(stab1, [1..126]);;

orb6list := Filtered(orbs, o -> Length(o) = 6 and not (1 in o));;

# Orbital 1 is disconnected (two copies of 63). Focus on orbital 2.
gamma := EdgeOrbitsGraph(im, [[1, orb6list[2][1]]], 126);;
gr := Girth(gamma);;
diam := Diameter(gamma);;
Print("\n=== 126-VERTEX COVER (orbital 2) ===\n");;
Print("  Connected, girth=", gr, ", diam=", diam, "\n");;

autGrp := AutGroupGraph(gamma);;
Print("  |Aut| = ", Size(autGrp), ", ", StructureDescription(autGrp), "\n");;

# Distance distribution
Print("  Dist: ");;
for d in [1..diam] do
    Print("d", d, "=", Length(DistanceSet(gamma,d,[1])), " ");;
od;
Print("\n");;

# Distance-regularity
isDistReg := true;;
dists1 := List([1..diam], d -> Length(DistanceSet(gamma,d,[1])));;
for v in [2..126] do
    dv := List([1..diam], d -> Length(DistanceSet(gamma,d,[v])));;
    if dv <> dists1 then isDistReg := false; break; fi;
od;
Print("  Distance-regular: ", isDistReg, "\n");;

if isDistReg then
    b := []; c := [];;
    for d in [0..diam] do
        w := First([1..126], u -> Distance(gamma,1,u) = d);;
        nbrs := Adjacency(gamma, w);;
        if d < diam then Add(b, Number(nbrs, u -> Distance(gamma,1,u) = d+1)); fi;
        if d > 0 then Add(c, Number(nbrs, u -> Distance(gamma,1,u) = d-1)); fi;
    od;
    Print("  Intersection array: {", b, "; ", c, "}\n");;
fi;

# Block quotient → Langer
Print("\n  Block analysis...\n");;
allBl := AllBlocks(im);;
Print("  Block sizes: ", List(allBl, Length), "\n");;

block2 := First(allBl, bs -> Length(bs) = 2 and 1 in bs);;
if block2 <> fail then
    allBlocks := Orbit(im, Set(block2), OnSets);;
    Print("  Block system: ", Length(allBlocks), " blocks of size ",
          Length(allBlocks[1]), "\n");;

    blockOf := function(v)
        local i;
        for i in [1..Length(allBlocks)] do
            if v in allBlocks[i] then return i; fi;
        od;
        return fail;
    end;;

    qAdj := List([1..Length(allBlocks)], i -> []);;
    for i in [1..Length(allBlocks)] do
        for v in allBlocks[i] do
            for w in Adjacency(gamma, v) do
                j := blockOf(w);;
                if j <> i then AddSet(qAdj[i], j); fi;
            od;
        od;
    od;

    degs := List(qAdj, Length);;
    Print("  Quotient degrees: ", Set(degs), "\n");;
    if Length(Set(degs)) = 1 and degs[1] = 6 then
        qGamma := Graph(Group(()), [1..63], OnPoints,
            function(x,y) return y in qAdj[x]; end, true);;
        Print("  *** QUOTIENT ISO TO LANGER: ",
              IsIsomorphicGraph(qGamma, langer), " ***\n");;
    fi;
else
    Print("  No block of size 2 containing vertex 1\n");;
fi;

# SECOND 63-VERTEX G_2(2) GRAPH
Print("\n\n=== SECOND 63-VERTEX G_2(2) GRAPH ===\n");;
maxclasses := MaximalSubgroupClassReps(G);;
M4 := First(maxclasses, m -> Size(m) = 192 and not IsConjugate(G, m, stab192));;
Print("M4 = ", StructureDescription(M4), "\n");;
act4 := FactorCosetAction(G, M4);;
im4 := Image(act4);;
stab4 := Stabilizer(im4, 1);;
orbs4 := Orbits(stab4, [1..63]);;
orb6_other := First(orbs4, o -> Length(o) = 6);;
graph2 := EdgeOrbitsGraph(im4, [[1, orb6_other[1]]], 63);;
Print("  girth=", Girth(graph2), ", diam=", Diameter(graph2),
      ", iso_to_langer=", IsIsomorphicGraph(graph2, langer), "\n");;

aut2 := AutGroupGraph(graph2);;
Print("  |Aut|=", Size(aut2), ", Aut=", StructureDescription(aut2), "\n");;

Print("  Dist: ");;
for d in [1..Diameter(graph2)] do
    Print("d", d, "=", Length(DistanceSet(graph2,d,[1])), " ");;
od;
Print("\n");;

isDistReg := true;;
dists1 := List([1..Diameter(graph2)], d -> Length(DistanceSet(graph2,d,[1])));;
for v in [2..63] do
    dv := List([1..Diameter(graph2)], d -> Length(DistanceSet(graph2,d,[v])));;
    if dv <> dists1 then isDistReg := false; break; fi;
od;
Print("  Distance-regular: ", isDistReg, "\n");;

if isDistReg then
    d2 := Diameter(graph2);;
    b := []; c := [];;
    for d in [0..d2] do
        w := First([1..63], u -> Distance(graph2,1,u) = d);;
        nbrs := Adjacency(graph2, w);;
        if d < d2 then Add(b, Number(nbrs, u -> Distance(graph2,1,u) = d+1)); fi;
        if d > 0 then Add(c, Number(nbrs, u -> Distance(graph2,1,u) = d-1)); fi;
    od;
    Print("  Intersection array: {", b, "; ", c, "}\n");;
fi;

Print("\n  Is graph2 the complement of Langer?\n");;
langerComp := ComplementGraph(langer);;
Print("  Complement of Langer: deg=", Length(Adjacency(langerComp, 1)),
      ", iso_to_graph2=", IsIsomorphicGraph(langerComp, graph2), "\n");;

# What IS graph2? Check if it's a known combinatorial object
Print("\n  Graph2 properties:\n");;
Print("  V=63, deg=6, girth=", Girth(graph2), ", diam=", Diameter(graph2), "\n");;
Print("  Edges: ", Length(UndirectedEdges(graph2)), "\n");;
# Check lambda, mu parameters
nbrs1 := Adjacency(graph2, 1);;
for v in nbrs1 do
    lambda := Length(Intersection(Adjacency(graph2, v), nbrs1));;
    Print("  lambda(1,", v, ")=", lambda, " ");;
od;
Print("\n");;

# mu values for non-adjacent pairs
non_nbrs := Difference([1..63], Union(nbrs1, [1]));;
mus := [];;
for v in non_nbrs do
    mu := Length(Intersection(Adjacency(graph2, v), nbrs1));;
    AddSet(mus, mu);;
od;
Print("  mu values: ", mus, "\n");;

Print("\n=== DONE ===\n");;
QuitGap(0);
