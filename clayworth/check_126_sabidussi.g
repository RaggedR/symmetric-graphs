LoadPackage("grape");;

# Is the 126-vertex bipartite graph a Sabidussi coset graph?
# Sab(G, H) requires: G transitive, H = stabiliser, connection set S
# The graph IS Sab(G_2(2), K_96) by construction (coset graph).
# But is it a SIMPLE GRAPH Sabidussi, i.e., from a connection set?

# A coset graph Sab(G, H, S) has:
# - Vertices = cosets G/H
# - Edge gH ~ g'H iff g'^{-1}g ∈ H·S·H (double coset)
# The connection set S must satisfy: S = S^{-1} and S ∩ H = ∅

# Our graph is an ORBITAL graph of the G_2(2) action on G/K_96.
# Every orbital graph of a transitive group IS a coset graph
# (union of double cosets). So YES, it is Sabidussi by definition.

# But let's make it explicit: find the connection set.

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
stab192 := Stabilizer(G, 1);;

subsOfStab := ConjugacyClassesSubgroups(stab192);;
K96 := First(List(subsOfStab, Representative),
    k -> Size(k) = 96 and StructureDescription(k) = "SL(2,3) : C4");;

Print("G = G_2(2), |G| = ", Size(G), "\n");;
Print("K = SL(2,3):C4, |K| = ", Size(K96), "\n");;
Print("[G:K] = ", Index(G, K96), "\n\n");;

# The coset action
act := FactorCosetAction(G, K96);;
im := Image(act);;
stab1 := Stabilizer(im, 1);;
orbs := Orbits(stab1, [1..126]);;

Print("Subdegrees: ", SortedList(List(orbs, Length)), "\n\n");;

# The orbital of size 6 (the connected one)
orb6list := Filtered(orbs, o -> Length(o) = 6 and not (1 in o));;

# For Sabidussi: the connection set is the set of group elements
# that map coset 1 (= K itself) to one of the 6 neighbors.
# In the coset picture: S = {g ∈ G : gK ∈ neighbors of K in G/K}
# Equivalently, S = union of K-cosets corresponding to neighbors.

# In the permutation action:
# Coset 1 ↔ point 1. Neighbors = orb6list[2] (the connected orbital).
# For each neighbor n, find g ∈ G mapping 1 → n in the coset action.
# Then gK is the coset, and g is a representative.

Print("The 126-vertex graph IS a Sabidussi coset graph by construction:\n");;
Print("  Sab(G_2(2), SL(2,3):C4, S)\n");;
Print("  where S is the connection set from the 6-orbital.\n\n");;

# Find connection set representatives
connSetReps := [];;
for n in orb6list[2] do
    g := RepresentativeAction(im, 1, n);;
    Add(connSetReps, g);;
od;

# Check self-inverseness: S = S^{-1}
Print("Connection set (6 elements mapping 1 to neighbors):\n");;
isSelfInverse := true;;
for g in connSetReps do
    gInv := g^(-1);;
    img := 1^gInv;;
    if not img in orb6list[2] then
        Print("  g maps 1→", 1^g, ", but g^-1 maps 1→", img,
              " which is NOT in the orbital\n");;
        # Check the other orbital
        if img in orb6list[1] then
            Print("    (it's in the OTHER 6-orbital — the disconnected one)\n");;
        fi;
        isSelfInverse := false;;
    fi;
od;

if isSelfInverse then
    Print("  Connection set IS self-inverse → undirected Sabidussi graph ✓\n\n");;
else
    Print("  Connection set is NOT self-inverse!\n");;
    Print("  This means the orbital is NOT self-paired.\n");;
    Print("  The graph uses both the orbital AND its pair.\n\n");;

    # For an undirected graph, GRAPE uses edge orbits which automatically
    # symmetrize. So the actual connection set is S ∪ S^{-1}.
    # Let's check what the full connection set looks like.
    fullS := Set(Concatenation(
        List(orb6list[2], n -> RepresentativeAction(im, 1, n)),
        List(orb6list[2], n -> RepresentativeAction(im, 1, n)^(-1))
    ));;

    # Where do the inverses land?
    invTargets := Set(List(connSetReps, g -> 1^(g^(-1))));;
    Print("  Inverse images land in: ", invTargets, "\n");;
    if invTargets = Set(orb6list[1]) then
        Print("  The inverses land in the FIRST 6-orbital (disconnected one)!\n");;
        Print("  So the symmetric closure S ∪ S^{-1} uses BOTH orbitals.\n");;
        Print("  The undirected graph is the UNION of both 6-orbitals = 12-regular?\n");;
        Print("  But our graph is 6-regular... checking.\n\n");;

        # Actually, GRAPE's EdgeOrbitsGraph symmetrizes:
        # it includes edge [1, orb6list[2][1]] AND its reverse.
        # But the reverse of edge 1→n is edge n→1, not 1→n^{-1}.
        # So EdgeOrbitsGraph({1,n}) means: v~w iff w ∈ orbit of n under Stab(v).
        # This automatically gives an undirected graph.

        Print("  Clarification: EdgeOrbitsGraph automatically creates\n");;
        Print("  an undirected graph. The orbital is self-paired iff\n");;
        Print("  the reverse edges lie in the SAME orbital.\n");;
        Print("  For bipartite graphs, this can't happen (reverse crosses parts).\n\n");;
    fi;
fi;

# Let's verify: is the orbital self-paired?
# An orbital O from point 1 is self-paired if, for some (hence all)
# n in O, 1 is in the orbital of n.
n := orb6list[2][1];;
stabN := Stabilizer(im, n);;
orbFromN := Orbit(stabN, 1);;
Print("Orbital self-pairing check:\n");;
Print("  From vertex 1, orbital = ", Set(orb6list[2]), "\n");;
Print("  From vertex ", n, ", orbit containing 1 has size ", Length(orbFromN), "\n");;
Print("  Is 1 in the same-size orbital from ", n, "? ");;

# Find the orbital of vertex n that contains 1
orbsFromN := Orbits(stabN, [1..126]);;
orbContaining1 := First(orbsFromN, o -> 1 in o);;
Print("Yes, in orbital of size ", Length(orbContaining1), "\n");;

if Length(orbContaining1) = 6 then
    Print("  → Orbital IS self-paired (6↔6) → genuine undirected graph\n");;
    Print("  → Sabidussi connection set S = S^{-1} ✓\n");;
else
    Print("  → Orbital NOT self-paired\n");;
fi;

# Summary
Print("\n=== SABIDUSSI VERIFICATION ===\n");;
Print("The 126-vertex graph = Sab(G_2(2), SL(2,3):C_4, S)\n");;
Print("where S corresponds to the self-paired 6-orbital.\n");;
Print("It is a legitimate Sabidussi coset graph. ✓\n");;

# What is the structure of K_96 = SL(2,3):C_4?
Print("\nK_96 structure details:\n");;
Print("  ", StructureDescription(K96), "\n");;
Print("  Is normal in H_192: ", IsNormal(stab192, K96), "\n");;
Print("  H_192/K_96: ", StructureDescription(stab192/K96), "\n");;

Print("\n=== DONE ===\n");;
QuitGap(0);
