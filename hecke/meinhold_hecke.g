# meinhold_hecke.g — Hecke algebras H(G₂(2), K) for each Meinhold graph
#
# For each intermediate subgroup C₃ ≤ K ≤ G₂(2), computes:
#   - Number of double cosets |K\G/K| = dim H(G₂(2), K)
#   - Valencies of each orbital
#   - Whether the Hecke algebra is commutative (Gelfand pair)
#   - Structure constants T_i * T_j = Σ c_{ij}^k T_k
#
# Usage: gap -q meinhold_hecke.g

# ─────────────────────────────────────────────
# Construct G₂(2) = PSU(3,3):C₂, order 12096
# ─────────────────────────────────────────────

G := PrimitiveGroup(63, 2);;
Print("G = ", StructureDescription(G), "\n");;
Print("|G| = ", Size(G), "\n\n");;

# ─────────────────────────────────────────────
# Find conjugacy classes of subgroups
# ─────────────────────────────────────────────

Print("Computing conjugacy classes of subgroups...\n");;
subs := ConjugacyClassesSubgroups(G);;
Print("Found ", Length(subs), " conjugacy classes\n\n");;

# ─────────────────────────────────────────────
# Hecke algebra computation
# ─────────────────────────────────────────────

ComputeHeckeAlgebra := function(G, K, name)
    local dcosets, r, i, j, k, valencies, n, struct_const,
          rep_k, z, prod, is_commutative, is_first, dc_index_map,
          elt_list, dc_of_elt, el;

    Print("\n╔═══════════════════════════════════════════════════╗\n");;
    Print("║  H(G₂(2), ", name, ")\n");;
    Print("╚═══════════════════════════════════════════════════╝\n\n");;

    Print("K = ", StructureDescription(K), ", |K| = ", Size(K), "\n");;
    n := Size(G) / Size(K);;
    Print("|G/K| = ", n, " vertices\n");;

    # Double cosets
    dcosets := DoubleCosets(G, K, K);;
    r := Length(dcosets);;
    Print("rank = |K\\G/K| = ", r, "\n\n");;

    # Valencies
    valencies := List(dcosets, d -> Size(d) / Size(K));;
    Print("Orbital valencies:\n");;
    for i in [1..r] do
        if Size(dcosets[i]) = Size(K) then
            Print("  T_", i, ": val ", valencies[i], " (identity)\n");;
        else
            Print("  T_", i, ": val ", valencies[i], "\n");;
        fi;
    od;
    Print("  sum = ", Sum(valencies), " (should be ", n, ")\n\n");;

    # ─────────────────────────────────────────
    # Symmetry check: KgK = Kg⁻¹K for all g?
    # If yes, the algebra is commutative (Gelfand pair).
    # ─────────────────────────────────────────

    is_commutative := true;;
    for i in [1..r] do
        if not Representative(dcosets[i])^(-1) in dcosets[i] then
            is_commutative := false;;
            Print("Asymmetric double coset: T_", i, "\n");;
        fi;
    od;
    if is_commutative then
        Print("All KgK = Kg⁻¹K → symmetric pair → GELFAND PAIR (commutative)\n\n");;
    else
        Print("Some KgK ≠ Kg⁻¹K → NOT necessarily a Gelfand pair\n");;
        Print("(Need full structure constants to confirm)\n\n");;
    fi;

    # ─────────────────────────────────────────
    # Structure constants (exact, via element enumeration)
    #
    # c_{ij}^k = |{z ∈ KaᵢK : z⁻¹·rₖ ∈ KaⱼK}| / |K|
    #
    # We precompute a lookup: for each g ∈ G, which double coset?
    # ─────────────────────────────────────────

    if n <= 600 then
        Print("Computing structure constants (|G| = ", Size(G), ")...\n");;

        # Build element → double coset index map
        elt_list := AsSortedList(G);;
        dc_of_elt := List([1..Size(G)], x -> 0);;
        for i in [1..r] do
            for el in dcosets[i] do
                j := PositionSorted(elt_list, el);;
                dc_of_elt[j] := i;;
            od;
        od;

        # Verify all elements are mapped
        if 0 in dc_of_elt then
            Print("ERROR: some elements not in any double coset!\n");;
            return fail;
        fi;

        struct_const := List([1..r], i -> List([1..r], j -> List([1..r], k -> 0)));;

        for k in [1..r] do
            rep_k := Representative(dcosets[k]);;
            for i in [1..r] do
                for z in dcosets[i] do
                    prod := z^(-1) * rep_k;;
                    j := dc_of_elt[PositionSorted(elt_list, prod)];;
                    struct_const[i][j][k] := struct_const[i][j][k] + 1;;
                od;
            od;
        od;

        # Normalize by |K|
        for i in [1..r] do
            for j in [1..r] do
                for k in [1..r] do
                    struct_const[i][j][k] := struct_const[i][j][k] / Size(K);;
                od;
            od;
        od;

        # Print multiplication table
        Print("\nMultiplication table:\n\n");;
        for i in [1..r] do
            for j in [1..r] do
                Print("  T_", i, " · T_", j, " = ");;
                is_first := true;;
                for k in [1..r] do
                    if struct_const[i][j][k] <> 0 then
                        if not is_first then Print(" + "); fi;;
                        if struct_const[i][j][k] <> 1 then
                            Print(struct_const[i][j][k]);;
                            Print("*");;
                        fi;
                        Print("T_", k);;
                        is_first := false;;
                    fi;
                od;
                if is_first then Print("0"); fi;;
                Print("\n");;
            od;
        od;

        # Verify commutativity from structure constants
        is_commutative := true;;
        for i in [1..r] do
            for j in [i+1..r] do
                for k in [1..r] do
                    if struct_const[i][j][k] <> struct_const[j][i][k] then
                        is_commutative := false;;
                    fi;
                od;
            od;
        od;
        Print("\nCommutative (verified): ", is_commutative, "\n");;

    else
        Print("(n = ", n, " > 600: skipping structure constants)\n");;
    fi;

    Print("\n");;
    return rec(rank := r, valencies := valencies, commutative := is_commutative);
end;;


# ─────────────────────────────────────────────
# Identify the Meinhold subgroups and run
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Meinhold family: Hecke algebras of G₂(2)\n");;
Print("═══════════════════════════════════════════════════\n");;

# Helper: find all subgroup classes of given order, print info
FindSubs := function(order)
    local result, sc, H;
    result := [];;
    for sc in subs do
        H := Representative(sc);;
        if Size(H) = order then
            Add(result, H);;
        fi;
    od;
    return result;
end;;

# Meinhold-72:   K = PSL(3,2), |K| = 168
k168 := FindSubs(168);;
Print("\nOrder 168: ", Length(k168), " class(es): ");;
for K in k168 do Print(StructureDescription(K), " "); od; Print("\n");;

# Meinhold-504:  K = S₄, |K| = 24  (need to find the right S₄)
k24 := FindSubs(24);;
Print("Order 24: ", Length(k24), " class(es): ");;
for K in k24 do Print(StructureDescription(K), " "); od; Print("\n");;

# Meinhold-576:  K = C₇:C₃, |K| = 21
k21 := FindSubs(21);;
Print("Order 21: ", Length(k21), " class(es): ");;
for K in k21 do Print(StructureDescription(K), " "); od; Print("\n");;

# Meinhold-1344: K = C₃×C₃, |K| = 9
k9 := FindSubs(9);;
Print("Order 9: ", Length(k9), " class(es): ");;
for K in k9 do Print(StructureDescription(K), " "); od; Print("\n");;

# Meinhold-2016: K = C₆, |K| = 6
k6 := FindSubs(6);;
Print("Order 6: ", Length(k6), " class(es): ");;
for K in k6 do Print(StructureDescription(K), " "); od; Print("\n");;

# ─── Run Hecke algebra for each ───

# PSL(3,2) — unique conjugacy class
if Length(k168) > 0 then
    ComputeHeckeAlgebra(G, k168[1], "PSL(3,2)");;
fi;

# S₄ — multiple classes, try each that is actually S₄
for i in [1..Length(k24)] do
    if StructureDescription(k24[i]) = "S4" then
        ComputeHeckeAlgebra(G, k24[i],
            Concatenation("S4-class", String(i)));;
        break;  # just the first S₄ class (conjugate classes give iso algebras)
    fi;
od;

# C₇:C₃
if Length(k21) > 0 then
    ComputeHeckeAlgebra(G, k21[1], "C7:C3");;
fi;

# C₃×C₃
for i in [1..Length(k9)] do
    if StructureDescription(k9[i]) = "C3 x C3" then
        ComputeHeckeAlgebra(G, k9[i], "C3xC3");;
        break;
    fi;
od;

# C₆
for i in [1..Length(k6)] do
    if StructureDescription(k6[i]) = "C6" then
        ComputeHeckeAlgebra(G, k6[i], "C6");;
        break;
    fi;
od;

Print("\n═══ DONE ═══\n");;
