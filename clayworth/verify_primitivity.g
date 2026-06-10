# Verify primitivity of Klein and Langer graphs

# === KLEIN: Sab(PSL(2,7), C_3) ===
Print("=== KLEIN GRAPH ===\n");;
G1 := PSL(2,7);;
Print("PSL(2,7) order: ", Size(G1), "\n");;

# Find C_3 subgroups
c3_classes := Filtered(ConjugacyClassesSubgroups(G1),
    c -> Size(Representative(c)) = 3);;
Print("C_3 conjugacy classes: ", Length(c3_classes), "\n");;

for ci in [1..Length(c3_classes)] do
    H := Representative(c3_classes[ci]);;
    Print("\n  Class ", ci, ": index ", Index(G1, H), "\n");;
    if Index(G1, H) = 56 then
        # This is the Klein graph stabiliser
        cosets := RightCosets(G1, H);;
        act := ActionHomomorphism(G1, cosets, OnRight);;
        P := Image(act);;
        is_prim := IsPrimitive(P);;
        Print("  Primitive: ", is_prim, "\n");;

        if not is_prim then
            blocks := AllBlocks(P);;
            for b in blocks do
                Print("    Block size ", Length(b), "\n");;
            od;
        fi;

        # Intermediate subgroups
        ints := IntermediateSubgroups(G1, H);;
        Print("  Intermediate subgroups:\n");;
        for K in ints.subgroups do
            Print("    |K| = ", Size(K), ", [G:K] = ", Index(G1, K),
                  ", ", StructureDescription(K), "\n");;
        od;
        if Length(ints.subgroups) = 0 then
            Print("    NONE — H is maximal, graph is PRIMITIVE\n");;
        fi;
    fi;
od;

# === LANGER: Sab(G_2(2), H_192) ===
Print("\n\n=== LANGER GRAPH ===\n");;

# G_2(2) as automorphism group of GH(2,2)
# We can build it from the Chevalley group
# G2(2) has a faithful permutation representation on 63 points
# Use the existing computation from the project

# Alternative: use AtlasGroup
G2 := AtlasGroup("U3(3).2");;  # G_2(2) ~ U(3,3).2 ? No...

# Actually G_2(2) is not in ATLAS under that name. Let me construct it.
# G_2(2) has order 12096 and acts on 63 points
# It's the automorphism group of the split Cayley hexagon GH(2,2)

# Use the Chevalley group construction
G2 := ChevalleyG(2, 2);;
if G2 = fail then
    # Try alternative
    G2 := SimpleGroup("G2(2)'");;
fi;

if IsBound(G2) and G2 <> fail then
    Print("G_2(2) order: ", Size(G2), "\n");;
else
    Print("Could not construct G_2(2), trying alternative...\n");;
    # G_2(2)' (derived subgroup) has order 6048, G_2(2) has order 12096
    # PSU(3,3) has order 6048, PSU(3,3).2 has order 12096
    G2 := PSU(3,3);;
    Print("PSU(3,3) order: ", Size(G2), "\n");;
    # G_2(2) = PSU(3,3).2
    G2full := AutomorphismGroup(G2);;
    Print("Aut(PSU(3,3)) order: ", Size(G2full), "\n");;
fi;
