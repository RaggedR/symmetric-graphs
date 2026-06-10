# build_triality_hexagon.g
# Construct the collinearity graph of GH(2,8) from 3D4(2)
# and export adjacency + orbit data for the Python renderer.
#
# The graph: 819 vertices, valency 18, distance-regular, diameter 3.
# Automorphism group: 3D4(2), order 211,341,312.

LoadPackage("grape");

Print("Constructing 3D4(2)...\n");
G := SimpleGroup("3D4(2)");
Print("Order: ", Size(G), "\n");

# Get the maximal parabolic P1 (stabilizer of a point of GH(2,8))
maxes := MaximalSubgroupClassReps(G);
P1 := maxes[1];  # index 819
Print("P1 order: ", Size(P1), ", index: ", Size(G)/Size(P1), "\n");

# Build the coset graph: vertices = cosets of P1
# Edge between gP1 and hP1 iff h^{-1}g in P1*a*P1 (the valency-18 suborbit)
Print("Computing double cosets...\n");
dc := DoubleCosets(G, P1, P1);
lens := List(dc, d -> Size(d)/Size(P1));
Print("Suborbit lengths: ", lens, "\n");

# Find the double coset of length 18 (collinearity)
collinear_dc := First(dc, d -> Size(d)/Size(P1) = 18);
rep := Representative(collinear_dc);
Print("Collinearity representative found, order ", Order(rep), "\n");

# Construct the orbital graph using GRAPE
Print("Building GRAPE graph (819 vertices, valency 18)...\n");
gamma := Graph(G, RightCosets(G, P1), OnRight, function(x, y)
    return RepresentativeAction(G, x, y) <> fail and
           Representative(y) * Representative(x)^(-1) in collinear_dc;
end, true);

# Actually, GRAPE's Orbital approach is better:
# Use the action on cosets directly
Print("Using orbital graph construction...\n");

# Get a permutation representation on 819 points
cos := RightCosets(G, P1);
act := Action(G, cos, OnRight);
Print("Permutation group on ", NrMovedPoints(act), " points\n");

# Find the suborbit of length 18
stab := Stabilizer(act, 1);
orbs := Orbits(stab, [1..819]);
Print("Point stabilizer orbits: ", List(orbs, Length), "\n");

# The orbit of length 18 gives us the neighbors of vertex 1
orb18 := First(orbs, o -> Length(o) = 18);
Print("Neighbors of vertex 1: ", orb18, "\n");

# Build full adjacency list using the group action
Print("Computing full adjacency (819 vertices)...\n");
adj := List([1..819], i -> []);
gens := GeneratorsOfGroup(act);

# Use the orbital: for each g in G, vertex g(1) is adjacent to g(v) for v in orb18
# More efficient: use the GRAPE package
gamma2 := EdgeOrbitsGraph(act, List(orb18, v -> [1, v]));
Print("Graph constructed. Vertices: ", gamma2.order, "\n");
Print("Verifying regularity...\n");
degs := List([1..819], v -> Length(Adjacency(gamma2, v)));
Print("Min degree: ", Minimum(degs), ", Max degree: ", Maximum(degs), "\n");

# Find a good cyclic subgroup for the layout
# 819 = 9 * 91 = 9 * 7 * 13
# Look for elements of order 7, 9, 13, or 91 for rotational symmetry
Print("\nLooking for cyclic subgroups for layout...\n");
Print("819 = 9 * 7 * 13\n");

# Try order 9: gives 91 vertices per orbit (ring)
# Try order 7: gives 117 vertices per orbit
# Try order 13: gives 63 vertices per orbit

# Find element of order 13 in the permutation group
found := false;
for g in List([1..100], i -> PseudoRandom(act)) do
    if Order(g) mod 13 = 0 then
        h := g^(Order(g)/13);
        if Order(h) = 13 then
            C13 := Group(h);
            orbs_c := Orbits(C13, [1..819]);
            sizes_c := Collected(List(orbs_c, Length));
            if Length(sizes_c) = 1 and sizes_c[1] = [13, 63] then
                Print("Found C13 with 63 orbits of size 13\n");
                found := true;
                break;
            fi;
        fi;
    fi;
od;

if not found then
    # Try order 7
    for g in List([1..100], i -> PseudoRandom(act)) do
        if Order(g) mod 7 = 0 then
            h := g^(Order(g)/7);
            if Order(h) = 7 then
                C7 := Group(h);
                orbs_c := Orbits(C7, [1..819]);
                sizes_c := Collected(List(orbs_c, Length));
                Print("C7 orbit structure: ", sizes_c, "\n");
                if ForAll(orbs_c, o -> Length(o) = 7) then
                    Print("Found C7 with 117 orbits of size 7\n");
                    C13 := C7;  # reuse variable name
                    found := true;
                    break;
                fi;
            fi;
        fi;
    od;
fi;

if not found then
    # Try order 9
    for g in List([1..200], i -> PseudoRandom(act)) do
        if Order(g) mod 9 = 0 then
            h := g^(Order(g)/9);
            if Order(h) = 9 then
                C9 := Group(h);
                orbs_c := Orbits(C9, [1..819]);
                sizes_c := Collected(List(orbs_c, Length));
                Print("C9 orbit structure: ", sizes_c, "\n");
                if ForAll(orbs_c, o -> Length(o) = 9) then
                    Print("Found C9 with 91 orbits of size 9\n");
                    C13 := C9;
                    found := true;
                    break;
                fi;
            fi;
        fi;
    od;
fi;

# Export adjacency list
Print("\n=== ADJACENCY OUTPUT ===\n");
for v in [1..819] do
    Print(v, ":", JoinStringsWithSeparator(List(Adjacency(gamma2, v), String), ","), "\n");
od;
Print("=== END ADJACENCY ===\n");

# Export orbit structure if found
if found then
    Print("\n=== ORBIT OUTPUT ===\n");
    for i in [1..Length(orbs_c)] do
        Print(i, ":", JoinStringsWithSeparator(List(orbs_c[i], String), ","), "\n");
    od;
    Print("=== END ORBIT ===\n");
fi;

Print("\nDone.\n");
