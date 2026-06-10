# Verify Langer graph primitivity
# G_2(2) = PSU(3,3).2, order 12096
# Langer = Sab(G_2(2), H_192) on 63 vertices

# G_2(2)' = PSU(3,3), order 6048
G2prime := PSU(3,3);;
Print("PSU(3,3) order: ", Size(G2prime), "\n");;

# G_2(2) = Aut(PSU(3,3)) has order 12096
# But we need the permutation representation on 63 points
# The Langer graph has Aut = G_2(2) acting on 63 points

# Use the Langer graph directly — read from the edge list
# and compute its automorphism group

# Actually, let's use the split Cayley hexagon construction
# The 63 points of Q(6,2) carry a G_2(2) action

# Alternative approach: check primitivity of the existing graph
# by loading the edge data and using nauty/bliss

# Simplest: construct G_2(2) as a matrix group over GF(2)
# G_2(2) preserves the quadratic form Q(x) = x0x4 + x1x5 + x2x6 + x3^2

# Let's just check: does PSU(3,3) have a transitive action on 63 points?
# PSU(3,3) has order 6048 = 63 * 96. So stabiliser would have order 96.
# G_2(2) has order 12096 = 63 * 192. Stabiliser order 192.

# Build G_2(2) as a permutation group via its action on GH(2,2) points
# Use AllPrimitiveGroups to check if there's a primitive group of degree 63
# and order 12096

prims := AllPrimitiveGroups(DegreeOperation, 63);;
Print("Primitive groups of degree 63:\n");;
for p in prims do
    Print("  Order ", Size(p), ", ", StructureDescription(p), "\n");;
od;

# Check if any has order 12096
g2_prim := Filtered(prims, p -> Size(p) = 12096);;
Print("\nPrimitive groups of degree 63 and order 12096: ", Length(g2_prim), "\n");;

if Length(g2_prim) > 0 then
    Print("G_2(2) IS primitive on 63 points\n");;

    G2 := g2_prim[1];;
    stab := Stabilizer(G2, 1);;
    Print("Stabiliser order: ", Size(stab), "\n");;
    Print("Stabiliser structure: ", StructureDescription(stab), "\n");;
    Print("Is maximal: ", stab in MaximalSubgroupClassReps(G2), "\n");;

    # Check intermediate subgroups
    ints := IntermediateSubgroups(G2, stab);;
    Print("Intermediate subgroups: ", Length(ints.subgroups), "\n");;
    if Length(ints.subgroups) = 0 then
        Print("Stabiliser is MAXIMAL — Langer is PRIMITIVE\n");;
    else
        for K in ints.subgroups do
            Print("  |K| = ", Size(K), ", [G:K] = ", Index(G2, K), "\n");;
        od;
    fi;
else
    Print("G_2(2) is NOT in the primitive groups database for degree 63\n");;
    Print("Checking transitive groups...\n");;

    # Try TransitiveGroup
    # or just construct from PSU(3,3)
    G2prime := PSU(3,3);;

    # Find a subgroup of index 63 in PSU(3,3)
    # 6048/63 = 96
    subs96 := Filtered(ConjugacyClassesSubgroups(G2prime),
        c -> Size(Representative(c)) = 96);;
    Print("Subgroups of PSU(3,3) of order 96: ", Length(subs96), "\n");;

    for si in [1..Length(subs96)] do
        H := Representative(subs96[si]);;
        cosets := RightCosets(G2prime, H);;
        if Length(cosets) = 63 then
            act := ActionHomomorphism(G2prime, cosets, OnRight);;
            P := Image(act);;
            is_prim := IsPrimitive(P);;
            Print("  Class ", si, ": primitive = ", is_prim, "\n");;
        fi;
    od;
fi;
