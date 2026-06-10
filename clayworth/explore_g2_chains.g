LoadPackage("grape");;

# =============================================================
# Explore ALL G_2(2) coset graphs
# Question: Is the Langer graph at the bottom of its own "chain"
# of G_2(2) graphs, or is it truly singular?
# =============================================================

# Build G_2(2) as Normalizer of PSU(3,3) in S_63
psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
Print("G = G_2(2), |G| = ", Size(G), "\n\n");;

# Part 1: ALL maximal subgroups of G_2(2)
Print("=== MAXIMAL SUBGROUPS OF G_2(2) ===\n");;
maxclasses := MaximalSubgroupClassReps(G);;
Print("Number of conjugacy classes of maximal subgroups: ", Length(maxclasses), "\n\n");;
for i in [1..Length(maxclasses)] do
    M := maxclasses[i];;
    Print("Class ", i, ": |M| = ", Size(M),
          ", [G:M] = ", Index(G, M),
          ", ", StructureDescription(M), "\n");;
od;

# Part 2: Subgroup lattice above H_192 (the Langer stabiliser)
Print("\n=== SUBGROUP LATTICE ABOVE H_192 (Langer stabiliser) ===\n");;
stab192 := Stabilizer(G, 1);;  # This IS H_192
Print("H_192 = Stab(1), |H| = ", Size(stab192),
      ", ", StructureDescription(stab192), "\n");;
ints := IntermediateSubgroups(G, stab192);;
Print("Intermediate subgroups: ", Length(ints.subgroups), "\n");;
if Length(ints.subgroups) = 0 then
    Print("NONE — H_192 is maximal. Langer has NO non-trivial quotients.\n");;
fi;

# Part 3: Maximal subgroup analysis — which give the Langer stabiliser?
Print("\n=== WHICH MAXIMAL CLASS IS H_192? ===\n");;
for i in [1..Length(maxclasses)] do
    M := maxclasses[i];;
    if Size(M) = 192 then
        if IsConjugate(G, M, stab192) then
            Print("Class ", i, " is conjugate to H_192 (Langer stabiliser)\n");;
        else
            Print("Class ", i, " has order 192 but is NOT conjugate to H_192\n");;
            Print("  Structure: ", StructureDescription(M), "\n");;
            # What graph does THIS give?
            act := FactorCosetAction(G, M);;
            im := Image(act);;
            stab1 := Stabilizer(im, 1);;
            orbs := Orbits(stab1, [1..63]);;
            subdeg := SortedList(List(orbs, Length));;
            Print("  Subdegrees: ", subdeg, "\n");;
        fi;
    fi;
od;

# Part 4: ALL cubic connected G_2(2) graphs
# Use GRAPE's EdgeOrbitsGraph for correct orbital graph construction
Print("\n\n=== ALL CUBIC CONNECTED G_2(2) GRAPHS ===\n");;
Print("(Checking all conjugacy classes of subgroups)\n\n");;

ccls := ConjugacyClassesSubgroups(G);;
Print("Total conjugacy classes of subgroups: ", Length(ccls), "\n\n");;

cubicGraphs := [];;
for ci in [1..Length(ccls)] do
    H := Representative(ccls[ci]);;
    idx := Index(G, H);;

    # Skip trivial cases
    if idx < 4 or idx > 10000 then continue; fi;

    # Build coset action
    act := FactorCosetAction(G, H);;
    im := Image(act);;

    # Check for cubic orbital
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..idx]);;
    cubOrbs := Filtered(orbs, o -> Length(o) = 3 and not (1 in o));;

    if Length(cubOrbs) > 0 then
        for co in cubOrbs do
            # Check if the orbital is self-paired (needed for undirected graph)
            # An orbital O from vertex 1 is self-paired if
            # for some v in O, 1 is in the orbit of v under stab_v
            v := co[1];;
            stabv := Stabilizer(im, v);;
            orbFromV := Orbit(stabv, 1);;
            selfPaired := (Length(orbFromV) <= 3);;  # crude check

            # Better: use EdgeOrbitsGraph
            edge := [1, co[1]];;
            gamma := EdgeOrbitsGraph(im, [edge], idx);;

            if not IsSimpleGraph(gamma) then
                continue;
            fi;

            if IsConnectedGraph(gamma) then
                gr := Girth(gamma);;
                diam := Diameter(gamma);;
                prim := IsPrimitive(im, [1..idx]);;
                Print("CUBIC GRAPH: |V|=", idx,
                      ", |H|=", Size(H),
                      ", H=", StructureDescription(H),
                      ", girth=", gr,
                      ", diam=", diam,
                      ", prim=", prim,
                      "\n");;
                Add(cubicGraphs, rec(
                    vertices := idx,
                    stabOrder := Size(H),
                    stabDesc := StructureDescription(H),
                    girth := gr,
                    diameter := diam,
                    primitive := prim,
                    subgroup := H
                ));
            fi;
        od;
    fi;
od;

Print("\n=== SUMMARY: ALL CUBIC CONNECTED G_2(2) GRAPHS ===\n");;
Print("Found ", Length(cubicGraphs), " cubic connected graphs:\n\n");;
for r in cubicGraphs do
    Print("  V=", r.vertices,
          ", H=", r.stabDesc,
          " (|H|=", r.stabOrder, ")",
          ", girth=", r.girth,
          ", diam=", r.diameter,
          ", primitive=", r.primitive, "\n");;
od;

# Part 5: Covering relationships among cubic graphs
Print("\n\n=== COVERING RELATIONSHIPS ===\n");;
Print("Sab(G,H1) covers Sab(G,H2) iff H1 ≤ H2^g for some g\n\n");;

for i in [1..Length(cubicGraphs)] do
    for j in [1..Length(cubicGraphs)] do
        if i = j then continue; fi;
        if cubicGraphs[i].vertices > cubicGraphs[j].vertices then
            # Graph i is larger (smaller stabiliser) — potential cover of j
            Hi := cubicGraphs[i].subgroup;;
            Hj := cubicGraphs[j].subgroup;;
            # Check if some conjugate of Hi is contained in Hj
            found := false;;
            for g in RightTransversal(G, Normalizer(G, Hi)) do
                if IsSubgroup(Hj, Hi^g) then
                    found := true;;
                    break;;
                fi;
            od;
            if found then
                Print("  Sab(G,", cubicGraphs[i].stabDesc,
                      ") [V=", cubicGraphs[i].vertices,
                      "] COVERS Sab(G,", cubicGraphs[j].stabDesc,
                      ") [V=", cubicGraphs[j].vertices, "]\n");;
            fi;
        fi;
    od;
od;

# Part 6: Subgroups of H_192 that give cubic covers of Langer
Print("\n\n=== SUBGROUPS OF H_192 GIVING CUBIC COVERS OF LANGER ===\n");;
Print("These are graphs that map onto Langer via block quotient\n\n");;

subsOfStab := ConjugacyClassesSubgroups(stab192);;
Print("Conjugacy classes of subgroups of H_192: ", Length(subsOfStab), "\n\n");;

langerCovers := [];;
for ci in [1..Length(subsOfStab)] do
    H := Representative(subsOfStab[ci]);;
    idx := Index(G, H);;
    if idx <= 63 or idx > 10000 then continue; fi;

    act := FactorCosetAction(G, H);;
    im := Image(act);;
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..idx]);;
    cubOrbs := Filtered(orbs, o -> Length(o) = 3 and not (1 in o));;

    if Length(cubOrbs) > 0 then
        for co in cubOrbs do
            gamma := EdgeOrbitsGraph(im, [[1, co[1]]], idx);;
            if not IsSimpleGraph(gamma) then continue; fi;
            if IsConnectedGraph(gamma) then
                gr := Girth(gamma);;
                diam := Diameter(gamma);;
                fib := Index(stab192, H);;
                Print("LANGER COVER: V=", idx,
                      ", |H|=", Size(H),
                      ", H=", StructureDescription(H),
                      ", girth=", gr,
                      ", diam=", diam,
                      ", fibre=", fib,
                      " (", fib, "-fold cover)\n");;
                Add(langerCovers, rec(
                    vertices := idx,
                    stabOrder := Size(H),
                    stabDesc := StructureDescription(H),
                    girth := gr,
                    diameter := diam,
                    fibreSize := fib
                ));
            fi;
        od;
    fi;
od;

if Length(langerCovers) = 0 then
    Print("NONE — no cubic connected covers of Langer exist within G_2(2)\n");;
fi;

# Part 7: Check if Clayworth is in the list (it should be: H=C_3, V=4032)
Print("\n\n=== IS CLAYWORTH IN A CHAIN WITH LANGER? ===\n");;
c3cls := Filtered(ccls, c -> Size(Representative(c)) = 3);;
Print("C_3 conjugacy classes: ", Length(c3cls), "\n");;
for ci in [1..Length(c3cls)] do
    H := Representative(c3cls[ci]);;
    nfixed := Number([1..63], p -> p^GeneratorsOfGroup(H)[1] = p);;
    Print("  Class ", ci, ": ", Size(c3cls[ci]),
          " conjugates, fixed points = ", nfixed, "\n");;

    # Is this C_3 contained in (a conjugate of) H_192?
    contained := false;;
    for g in RightTransversal(G, Normalizer(G, H)) do
        if IsSubgroup(stab192, H^g) then
            contained := true;;
            break;;
        fi;
    od;
    Print("    C_3 ≤ H_192 (up to conj): ", contained, "\n");;
od;

Print("\n=== DONE ===\n");;
QuitGap(0);
