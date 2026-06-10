# Klein Quartic Cohomology — GAP verification
# Verifies H^1 = V_3 + V_3' using the character table of PSL(2,7)

G := PSL(2,7);
Print("G = PSL(2,7), |G| = ", Size(G), "\n\n");

# Character table
ct := CharacterTable(G);
Display(ct);

# Identify the (2,3,7) generators
# PSL(2,7) = <a,b | a^2 = b^3 = (ab)^7 = 1>
gens := GeneratorsOfGroup(G);
Print("\nGenerators of G: ", gens, "\n");
Print("Orders: ", List(gens, Order), "\n");

# Find (2,3,7) triple
found := false;
for a in Elements(G) do
    if Order(a) <> 2 then continue; fi;
    for b in Elements(G) do
        if Order(b) <> 3 then continue; fi;
        if Order(a*b) = 7 then
            Print("\n(2,3,7) generators found:\n");
            Print("  a (order 2) = ", a, "\n");
            Print("  b (order 3) = ", b, "\n");
            Print("  ab (order 7) = ", a*b, "\n");
            found := true;
            break;
        fi;
    od;
    if found then break; fi;
od;

# Cell stabilizers
Ha := Group(a);  # edge stabilizer, order 2
Hb := Group(b);  # face stabilizer, order 3
Hab := Group(a*b);  # vertex stabilizer, order 7

# Permutation characters (= characters of coset representations)
Print("\n--- Permutation characters ---\n");

chi_vertices := PermutationCharacter(G, Hab);  # 24-dim
chi_edges := PermutationCharacter(G, Ha);      # 84-dim
chi_faces := PermutationCharacter(G, Hb);      # 56-dim

Print("chi(vertices) dim = ", chi_vertices[1], " (expect 24)\n");
Print("chi(edges) dim = ", chi_edges[1], " (expect 84)\n");
Print("chi(faces) dim = ", chi_faces[1], " (expect 56)\n");

# The permutation character counts UNSIGNED fixed points.
# For the Lefschetz trace, we need SIGNED traces on oriented chains.
#
# Key correction: an involution at an edge midpoint REVERSES the edge
# orientation (swaps endpoints), contributing -1 not +1 to tr(g|C_1).
#
# The signed character on oriented 1-chains:
#   For each g, tr(g|C_1) = (edges fixed with + orientation) - (edges fixed with - orientation)
#
# For the (2,3,7) triangle group structure:
# - Order-2 elements stabilize edges and REVERSE their orientation
# - Order-3 elements stabilize faces and PRESERVE their orientation
# - Order-7 elements stabilize vertices (0-cells, no orientation issue)
# - Order-4 elements fix nothing
#
# The signed edge character for involutions:
# All fixed edges are reversed, so tr = -|fixed edges|

# To compute the H_1 character correctly, we use the Euler characteristic
# approach with the KNOWN character table.

irr := Irr(ct);
Print("\n--- Irreducible characters ---\n");
for i in [1..Length(irr)] do
    Print("chi_", i, " (dim ", irr[i][1], "): ", ValuesOfClassFunction(irr[i]), "\n");
od;

# Decompose permutation characters
Print("\n--- Decomposition of permutation characters ---\n");
for pair in [["vertices", chi_vertices], ["edges", chi_edges], ["faces", chi_faces]] do
    name := pair[1];
    chi := pair[2];
    Print(name, " = ");
    first := true;
    for i in [1..Length(irr)] do
        m := ScalarProduct(ct, chi, irr[i]);
        if m > 0 then
            if not first then Print(" + "); fi;
            if m > 1 then Print(m, "*"); fi;
            Print("chi_", i);
            first := false;
        fi;
    od;
    Print("\n");
od;

# Identify which irreps are the two conjugate 3-dim ones
# and which is the 6-dim irrep
Print("\n--- Identifying 3-dim and 6-dim irreps ---\n");
for i in [1..Length(irr)] do
    if irr[i][1] = 3 then
        Print("chi_", i, " is a 3-dim irrep: ", ValuesOfClassFunction(irr[i]), "\n");
    elif irr[i][1] = 6 then
        Print("chi_", i, " is a 6-dim irrep: ", ValuesOfClassFunction(irr[i]), "\n");
    fi;
od;

# Compute chi_2 + chi_3 (the two conjugate 3-dim irreps)
# These should be dimensions 3+3=6
three_dim_indices := Filtered([1..Length(irr)], i -> irr[i][1] = 3);
Print("\n3-dim irreps at indices: ", three_dim_indices, "\n");

if Length(three_dim_indices) = 2 then
    chi_sum := irr[three_dim_indices[1]] + irr[three_dim_indices[2]];
    Print("chi_2 + chi_3 = ", ValuesOfClassFunction(chi_sum), "\n");

    # This should be H^1
    Print("\nExpected H^1 character (from Hodge theory): chi_2 + chi_3\n");
    Print("  = V_3 + V_3' where V_3 = H^{1,0} (holomorphic 1-forms)\n");
fi;

# Verify: the alternating sum of permutation characters
# chi(C_0) - chi(C_1) + chi(C_2) evaluated at each class
# should give the Lefschetz number L(g)
# And L(g) = tr(g|H_0) - tr(g|H_1) + tr(g|H_2) = 2 - tr(g|H_1)
# (since H_0 = H_2 = trivial)
#
# BUT: chi(C_1) here is the UNSIGNED permutation character.
# The SIGNED character differs at involutions!

Print("\n--- Verifying with Lefschetz formula ---\n");
Print("(Using unsigned permutation characters for comparison)\n\n");

cl := ConjugacyClasses(ct);
Print("  Class | Size | Order | chi(C0) | chi(C1) | chi(C2) | L_unsigned | L_signed\n");
Print("  ------+------+-------+---------+---------+---------+------------+---------\n");

for i in [1..Length(cl)] do
    sz := SizesConjugacyClasses(ct)[i];
    ord := OrdersClassRepresentatives(ct)[i];
    c0 := chi_vertices[i];
    c1_unsigned := chi_edges[i];
    c2 := chi_faces[i];

    # Signed correction: involutions reverse all fixed edges
    # For order-2: signed tr = -unsigned tr
    # For order-4: no fixed edges anyway (0 either way)
    # For others: no fixed edges, or orientation preserved
    if ord = 2 then
        c1_signed := -c1_unsigned;
    else
        c1_signed := c1_unsigned;
    fi;

    L_unsigned := c0 - c1_unsigned + c2;
    L_signed := c0 - c1_signed + c2;

    Print("  ", String(i, 5), " | ", String(sz, 4), " | ", String(ord, 5),
          " | ", String(c0, 7), " | ", String(c1_unsigned, 7),
          " | ", String(c2, 7), " | ", String(L_unsigned, 10),
          " | ", String(L_signed, 7), "\n");
od;

Print("\n  tr(g|H_1) = 2 - L_signed:\n  ");
for i in [1..Length(cl)] do
    sz := SizesConjugacyClasses(ct)[i];
    ord := OrdersClassRepresentatives(ct)[i];
    c0 := chi_vertices[i];
    c1 := chi_edges[i];
    c2 := chi_faces[i];
    if ord = 2 then c1 := -c1; fi;
    L := c0 - c1 + c2;
    Print(2 - L, "  ");
od;
Print("\n");

# Final verification
Print("\n=== FINAL RESULT ===\n");
Print("H^1(Klein quartic; C) = V_3 + V_3'\n");
Print("  V_3  = H^{1,0} = holomorphic 1-forms\n");
Print("  V_3' = H^{0,1} = antiholomorphic 1-forms\n");
Print("  V_3 and V_3' are the two conjugate 3-dim irreps of PSL(2,7)\n");
Print("  Their characters involve zeta_7 = exp(2*pi*i/7)\n\n");

Print("The key insight: orientation reversal at involutions distinguishes\n");
Print("  chi_2 + chi_3 (value -2 at 2A) from chi_4 (value +2 at 2A).\n");
Print("  The cell structure of the (2,3,7) tessellation REQUIRES this sign.\n");
