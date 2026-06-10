LoadPackage("grape");;

# =============================================================
# The Langer graph is 6-REGULAR, not cubic!
# Its orbital has size 6 from the subdegree [1, 6, 24, 32].
# We need to check ALL orbital graphs, not just cubic ones.
# =============================================================

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
Print("G = G_2(2), |G| = ", Size(G), "\n\n");;

stab192 := Stabilizer(G, 1);;
Print("H_192 = ", StructureDescription(stab192), "\n");;

# Confirm Langer's subdegrees
orbs := Orbits(stab192, [1..63]);;
Print("Langer subdegrees: ", SortedList(List(orbs, Length)), "\n\n");;

# The orbital of size 6 gives the Langer graph
orb6 := First(orbs, o -> Length(o) = 6);;
Print("The 6-orbital: ", orb6, "\n");;
langer := EdgeOrbitsGraph(G, [[1, orb6[1]]], 63);;
Print("Langer: simple=", IsSimpleGraph(langer),
      ", connected=", IsConnectedGraph(langer),
      ", girth=", Girth(langer), "\n\n");;

# Now systematically check ALL subgroups for graphs of ANY regularity
# that are connected and have the same automorphism group
Print("=== ALL CONNECTED G_2(2)-ORBITAL GRAPHS ===\n");;
Print("(Checking orbital graphs for all subgroup actions)\n\n");;

ccls := ConjugacyClassesSubgroups(G);;
allGraphs := [];;

for ci in [1..Length(ccls)] do
    H := Representative(ccls[ci]);;
    idx := Index(G, H);;
    if idx < 4 or idx > 10000 then continue; fi;

    act := FactorCosetAction(G, H);;
    im := Image(act);;
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..idx]);;

    # Check each non-trivial orbital
    for oi in [1..Length(orbs)] do
        o := orbs[oi];;
        if 1 in o then continue; fi;  # skip fixed-point orbit
        deg := Length(o);;
        if deg > 20 then continue; fi;  # skip very high degree

        gamma := EdgeOrbitsGraph(im, [[1, o[1]]], idx);;
        if not IsSimpleGraph(gamma) then continue; fi;
        if not IsConnectedGraph(gamma) then continue; fi;

        gr := Girth(gamma);;
        diam := Diameter(gamma);;
        prim := IsPrimitive(im, [1..idx]);;

        Print("GRAPH: V=", idx,
              ", deg=", deg,
              ", |H|=", Size(H),
              ", H=", StructureDescription(H),
              ", girth=", gr,
              ", diam=", diam,
              ", prim=", prim, "\n");;

        Add(allGraphs, rec(
            vertices := idx,
            degree := deg,
            stabOrder := Size(H),
            stabDesc := StructureDescription(H),
            girth := gr,
            diameter := diam,
            primitive := prim,
            subgroup := H
        ));
    od;
od;

Print("\n=== SUMMARY: ALL CONNECTED G_2(2)-ORBITAL GRAPHS (deg ≤ 20) ===\n");;
Print("Found ", Length(allGraphs), " graphs:\n\n");;

# Sort by vertex count then degree
SortBy(allGraphs, r -> [r.vertices, r.degree]);;
for r in allGraphs do
    Print("  V=", r.vertices,
          ", deg=", r.degree,
          ", H=", r.stabDesc,
          " (|H|=", r.stabOrder, ")",
          ", girth=", r.girth,
          ", diam=", r.diameter,
          ", prim=", r.primitive, "\n");;
od;

# Part 2: Check specifically — is ANYTHING a cover of the Langer graph?
# A graph Sab(G,K) with K < H_192 covers Langer
Print("\n\n=== COVERS OF LANGER (subgroups K < H_192) ===\n");;
Print("Looking for ANY connected orbital graph from subgroups of H_192\n\n");;

subsOfStab := ConjugacyClassesSubgroups(stab192);;
langerCovers := [];;

for ci in [1..Length(subsOfStab)] do
    K := Representative(subsOfStab[ci]);;
    idx := Index(G, K);;
    if idx <= 63 or idx > 10000 then continue; fi;

    act := FactorCosetAction(G, K);;
    im := Image(act);;
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..idx]);;

    for o in orbs do
        if 1 in o then continue; fi;
        deg := Length(o);;
        if deg > 20 then continue; fi;

        gamma := EdgeOrbitsGraph(im, [[1, o[1]]], idx);;
        if not IsSimpleGraph(gamma) then continue; fi;
        if not IsConnectedGraph(gamma) then continue; fi;

        gr := Girth(gamma);;
        fib := Index(stab192, K);;
        Print("LANGER COVER: V=", idx,
              ", deg=", deg,
              ", |K|=", Size(K),
              ", K=", StructureDescription(K),
              ", girth=", gr,
              ", fibre=", fib, "\n");;
        Add(langerCovers, rec(v := idx, deg := deg, fib := fib));;
    od;
od;

if Length(langerCovers) = 0 then
    Print("NONE — Langer has NO covers of degree ≤ 20 within G_2(2)!\n");;
fi;

# Part 3: The OTHER order-192 maximal subgroup
# It gives a DIFFERENT 63-vertex action. What graph does it produce?
Print("\n\n=== THE OTHER ORDER-192 MAXIMAL SUBGROUP ===\n");;
maxclasses := MaximalSubgroupClassReps(G);;
for i in [1..Length(maxclasses)] do
    M := maxclasses[i];;
    if Size(M) = 192 and not IsConjugate(G, M, stab192) then
        Print("Class ", i, ": ", StructureDescription(M), "\n");;
        act := FactorCosetAction(G, M);;
        im := Image(act);;
        stab1 := Stabilizer(im, 1);;
        orbs := Orbits(stab1, [1..63]);;
        Print("  Subdegrees: ", SortedList(List(orbs, Length)), "\n");;

        # Build the 6-orbital graph
        orb6 := First(orbs, o -> Length(o) = 6);;
        if orb6 <> fail then
            gamma := EdgeOrbitsGraph(im, [[1, orb6[1]]], 63);;
            Print("  6-orbital: simple=", IsSimpleGraph(gamma));;
            if IsSimpleGraph(gamma) then
                Print(", connected=", IsConnectedGraph(gamma));;
                if IsConnectedGraph(gamma) then
                    Print(", girth=", Girth(gamma),
                          ", diam=", Diameter(gamma));;
                    # Is it isomorphic to Langer?
                    Print(", iso_to_langer=",
                          IsIsomorphicGraph(gamma, langer));;
                fi;
            fi;
            Print("\n");;
        fi;

        # Check all small orbitals
        for o in orbs do
            if 1 in o then continue; fi;
            deg := Length(o);;
            if deg > 20 then continue; fi;
            gamma2 := EdgeOrbitsGraph(im, [[1, o[1]]], 63);;
            if IsSimpleGraph(gamma2) and IsConnectedGraph(gamma2) then
                Print("  Orbital deg ", deg, ": girth=", Girth(gamma2),
                      ", diam=", Diameter(gamma2),
                      ", iso_to_langer=", IsIsomorphicGraph(gamma2, langer),
                      "\n");;
            fi;
        od;
    fi;
od;

Print("\n=== DONE ===\n");;
QuitGap(0);
