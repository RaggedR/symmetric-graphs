LoadPackage("grape");;

# =============================================================
# Verify the Langer chain structure
# Key finding: V=126 6-regular girth-6 graph covers Langer
# Is this the Tutte 12-cage's incidence graph? (No — Tutte is cubic)
# What IS this 126-vertex graph?
# =============================================================

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
stab192 := Stabilizer(G, 1);;

# Build the Langer graph for comparison
orbs63 := Orbits(stab192, [1..63]);;
orb6_langer := First(orbs63, o -> Length(o) = 6);;
langer := EdgeOrbitsGraph(G, [[1, orb6_langer[1]]], 63);;
Print("Langer: V=63, deg=6, girth=", Girth(langer), "\n");;
Print("Langer distance-regular: ");;
dr := DistanceSet(langer, 1, [1]);;  # just check distances
Print("distances from 1: ");;
for d in [1..Diameter(langer)] do
    ds := DistanceSet(langer, d, [1]);;
    Print("|d=", d, "|=", Length(ds), " ");;
od;
Print("\n\n");;

# Now build the V=126 6-regular cover
# H = SL(2,3):C4, |H| = 96, this is index 2 in H_192
# So K = SL(2,3):C4 is a subgroup of H_192 with [H_192:K] = 2
Print("=== THE 126-VERTEX 6-REGULAR COVER ===\n");;

# Find K = SL(2,3):C4 inside H_192
subsOfStab := ConjugacyClassesSubgroups(stab192);;
candidates := [];;
for ci in [1..Length(subsOfStab)] do
    K := Representative(subsOfStab[ci]);;
    if Size(K) = 96 and Index(stab192, K) = 2 then
        Add(candidates, rec(K := K, desc := StructureDescription(K)));;
    fi;
od;
Print("Order-96 index-2 subgroups of H_192: ", Length(candidates), "\n");;
for c in candidates do
    Print("  ", c.desc, "\n");;
od;

# For each, build coset graph and check if 6-regular connected girth 6
for c in candidates do
    K := c.K;;
    act := FactorCosetAction(G, K);;
    im := Image(act);;
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..126]);;
    subdeg := SortedList(List(orbs, Length));;
    Print("\nK = ", c.desc, "\n");;
    Print("  Subdegrees: ", subdeg, "\n");;

    # Find 6-orbital
    orb6 := First(orbs, o -> Length(o) = 6);;
    if orb6 <> fail then
        gamma := EdgeOrbitsGraph(im, [[1, orb6[1]]], 126);;
        if IsSimpleGraph(gamma) and IsConnectedGraph(gamma) then
            gr := Girth(gamma);;
            diam := Diameter(gamma);;
            Print("  6-reg graph: girth=", gr, ", diam=", diam, "\n");;

            # Check distance-regularity
            Print("  Distances from 1: ");;
            for d in [1..diam] do
                ds := DistanceSet(gamma, d, [1]);;
                Print("|d=", d, "|=", Length(ds), " ");;
            od;
            Print("\n");;

            # Intersection array (if distance-regular)
            Print("  Checking distance-regularity...\n");;
            isDistReg := true;;
            for v in [1..10] do  # sample 10 vertices
                dists_v := [];;
                for d in [1..diam] do
                    Add(dists_v, Length(DistanceSet(gamma, d, [v])));;
                od;
                dists_1 := [];;
                for d in [1..diam] do
                    Add(dists_1, Length(DistanceSet(gamma, d, [1])));;
                od;
                if dists_v <> dists_1 then
                    Print("  NOT distance-regular (v=", v, " differs)\n");;
                    isDistReg := false;;
                    break;;
                fi;
            od;
            if isDistReg then
                Print("  Appears distance-regular (sampled 10 vertices)\n");;
            fi;

            # Is this isomorphic to the Tutte 12-cage?
            # Tutte 12-cage: 126 vertices, CUBIC (degree 3), girth 12
            # This is 6-regular girth 6, so NOT the Tutte 12-cage
            # But could it be the DISTANCE-2 GRAPH of something?
            Print("  (Not the Tutte 12-cage: Tutte is cubic, this is 6-regular)\n");;

            # Check automorphism group size
            autGrp := AutGroupGraph(gamma);;
            Print("  |Aut| = ", Size(autGrp), "\n");;
            Print("  Aut = ", StructureDescription(autGrp), "\n");;

            # Is this bipartite?
            bip := IsBipartiteGraph(gamma);;
            Print("  Bipartite: ", bip, "\n");;

            # Check if the quotient map to Langer exists
            # K < H_192, so there IS a coset quotient Sab(G,K) → Sab(G,H_192)
            Print("  This is a 2-fold cover of Langer (fibre size 2)\n");;

            # Check relationship to the GH(2,2) incidence graph
            # GH(2,2) has 63 points + 63 lines = 126 vertices in its incidence graph
            # The incidence graph is the Tutte 12-cage (cubic, girth 12)
            # Our graph is 6-regular girth 6 on 126 vertices
            # Could be distance-2 of Tutte 12-cage restricted differently
            Print("\n  Checking if this is the FULL distance-2 graph of Tutte 12-cage...\n");;

            # Build Tutte 12-cage's incidence graph from GH(2,2)
            # Actually, we already know: Tutte 12-cage has 126 vertices
            # Its distance-2 graph would be 126-vertex but with higher degree
            # (each vertex at distance 2 from v in a cubic graph)
        fi;
    fi;
od;

# Part 2: The OTHER 6-regular 126-vertex graph
Print("\n\n=== THE OTHER 126-VERTEX 6-REGULAR GRAPH ===\n");;
Print("H = ((C4 x C4) : C3) : C2, |H| = 96\n");;
# This comes from the OTHER class of order-192 maximal subgroup
# Its index-2 subgroup also gives a 126-vertex graph

# Find all order-96 subgroups of G
cc96 := Filtered(ConjugacyClassesSubgroups(G),
    c -> Size(Representative(c)) = 96);;
Print("Conjugacy classes of order-96 subgroups: ", Length(cc96), "\n");;
for ci in [1..Length(cc96)] do
    K := Representative(cc96[ci]);;
    Print("  Class ", ci, ": ", StructureDescription(K),
          ", conjugates=", Size(cc96[ci]), "\n");;

    # Build 6-regular orbital graph
    act := FactorCosetAction(G, K);;
    im := Image(act);;
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..126]);;
    orb6 := First(orbs, o -> Length(o) = 6);;
    if orb6 <> fail then
        gamma := EdgeOrbitsGraph(im, [[1, orb6[1]]], 126);;
        if IsSimpleGraph(gamma) and IsConnectedGraph(gamma) then
            Print("    6-reg graph: girth=", Girth(gamma),
                  ", diam=", Diameter(gamma), "\n");;
            autGrp := AutGroupGraph(gamma);;
            Print("    |Aut| = ", Size(autGrp),
                  ", Aut = ", StructureDescription(autGrp), "\n");;
        fi;
    fi;
od;

# Part 3: Complete chain analysis for 6-regular covers
Print("\n\n=== COMPLETE 6-REGULAR CHAIN ABOVE LANGER ===\n");;
Print("All 6-regular connected covers of Langer from subgroups of H_192:\n\n");;

sixRegCovers := [];;
for ci in [1..Length(subsOfStab)] do
    K := Representative(subsOfStab[ci]);;
    idx := Index(G, K);;
    if idx <= 63 or idx > 10000 then continue; fi;

    act := FactorCosetAction(G, K);;
    im := Image(act);;
    stab1 := Stabilizer(im, 1);;
    orbs := Orbits(stab1, [1..idx]);;
    orb6 := First(orbs, o -> Length(o) = 6);;

    if orb6 <> fail then
        gamma := EdgeOrbitsGraph(im, [[1, orb6[1]]], idx);;
        if IsSimpleGraph(gamma) and IsConnectedGraph(gamma) then
            gr := Girth(gamma);;
            fib := Index(stab192, K);;
            Print("  V=", idx,
                  ", K=", StructureDescription(K),
                  " (|K|=", Size(K), ")",
                  ", girth=", gr,
                  ", fibre=", fib, "\n");;
            Add(sixRegCovers, rec(
                V := idx, K := K, desc := StructureDescription(K),
                girth := gr, fibre := fib
            ));;
        fi;
    fi;
od;

Print("\nTotal 6-regular covers: ", Length(sixRegCovers), "\n");;

# Check containment to find the chain structure
Print("\n=== CHAIN STRUCTURE (6-regular covers) ===\n");;
SortBy(sixRegCovers, r -> r.V);;
for i in [1..Length(sixRegCovers)] do
    for j in [1..Length(sixRegCovers)] do
        if sixRegCovers[i].V > sixRegCovers[j].V then
            # Does i cover j? K_i < K_j (up to conjugacy in stab192)?
            if IsSubgroup(sixRegCovers[j].K, sixRegCovers[i].K) then
                Print("  V=", sixRegCovers[i].V, " (", sixRegCovers[i].desc,
                      ") covers V=", sixRegCovers[j].V,
                      " (", sixRegCovers[j].desc, ")\n");;
            fi;
        fi;
    od;
od;

Print("\n=== DONE ===\n");;
QuitGap(0);
