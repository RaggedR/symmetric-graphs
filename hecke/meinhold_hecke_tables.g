# meinhold_hecke_tables.g — Compute and save multiplication tables for all
# five Meinhold Hecke algebras H(G₂(2), K).
#
# Outputs one file per algebra: meinhold_mult_NNN.txt
#
# Usage: gap -q meinhold_hecke_tables.g

# ─────────────────────────────────────────────
# G₂(2) = PSU(3,3):C₂, order 12096
# ─────────────────────────────────────────────

G := PrimitiveGroup(63, 2);;
Print("G = ", StructureDescription(G), ", |G| = ", Size(G), "\n\n");;

subs := ConjugacyClassesSubgroups(G);;

# ─────────────────────────────────────────────
# Core computation
# ─────────────────────────────────────────────

ComputeAndSave := function(G, K, graphname)
    local dcosets, r, i, j, k, valencies, n, struct_const,
          rep_k, z, prod, elt_list, dc_of_elt, el,
          filename, is_first, line, is_comm;

    n := Size(G) / Size(K);;
    dcosets := DoubleCosets(G, K, K);;
    r := Length(dcosets);;
    valencies := List(dcosets, d -> Size(d) / Size(K));;

    Print("── ", graphname, " ──\n");;
    Print("  K = ", StructureDescription(K), ", |K| = ", Size(K),
          ", n = ", n, ", rank = ", r, "\n");;

    # Build element → double coset index map
    Print("  Building element lookup...\n");;
    elt_list := AsSortedList(G);;
    dc_of_elt := List([1..Size(G)], x -> 0);;
    for i in [1..r] do
        for el in dcosets[i] do
            j := PositionSorted(elt_list, el);;
            dc_of_elt[j] := i;;
        od;
    od;

    # Compute structure constants
    Print("  Computing ", r, "×", r, " structure constants...\n");;
    struct_const := List([1..r], i -> List([1..r], j -> List([1..r], k -> 0)));;

    for k in [1..r] do
        if k mod 20 = 0 then Print("    k = ", k, "/", r, "\n"); fi;;
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

    # Check commutativity
    is_comm := true;;
    for i in [1..r] do
        for j in [i+1..r] do
            for k in [1..r] do
                if struct_const[i][j][k] <> struct_const[j][i][k] then
                    is_comm := false;;
                fi;
            od;
        od;
    od;
    Print("  Commutative: ", is_comm, "\n");;

    # Write to file
    filename := Concatenation("meinhold_mult_", String(n), ".txt");;
    Print("  Writing to ", filename, "...\n");;
    PrintTo(filename, "");;

    AppendTo(filename, "# Hecke algebra H(G₂(2), ", StructureDescription(K), ")\n");;
    AppendTo(filename, "# Graph: ", graphname, "\n");;
    AppendTo(filename, "# |G| = ", Size(G), ", |K| = ", Size(K),
             ", n = |G/K| = ", n, "\n");;
    AppendTo(filename, "# rank = ", r, " (number of double cosets)\n");;
    AppendTo(filename, "# Commutative: ", is_comm, "\n");;
    AppendTo(filename, "#\n");;

    # Valencies
    AppendTo(filename, "# Orbital valencies:\n");;
    for i in [1..r] do
        if valencies[i] = 1 then
            AppendTo(filename, "#   T_", i, ": val 1 (identity)\n");;
        else
            AppendTo(filename, "#   T_", i, ": val ", valencies[i], "\n");;
        fi;
    od;
    AppendTo(filename, "#\n");;
    AppendTo(filename, "# Multiplication table: T_i * T_j = sum_k c_{ij}^k T_k\n");;
    AppendTo(filename, "#\n\n");;

    for i in [1..r] do
        for j in [1..r] do
            is_first := true;;
            line := Concatenation("T_", String(i), " * T_", String(j), " = ");;
            for k in [1..r] do
                if struct_const[i][j][k] <> 0 then
                    if not is_first then
                        Append(line, " + ");;
                    fi;
                    if struct_const[i][j][k] <> 1 then
                        Append(line, Concatenation(String(struct_const[i][j][k]), "*"));;
                    fi;
                    Append(line, Concatenation("T_", String(k)));;
                    is_first := false;;
                fi;
            od;
            if is_first then Append(line, "0"); fi;;
            AppendTo(filename, line, "\n");;
        od;
    od;

    # Also write raw structure constants as GAP-readable data
    AppendTo(filename, "\n# Raw structure constants (GAP format):\n");;
    AppendTo(filename, "# struct_const[i][j][k] = coefficient of T_k in T_i * T_j\n");;
    AppendTo(filename, "struct_const_", String(n), " := ", struct_const, ";\n");;
    AppendTo(filename, "valencies_", String(n), " := ", valencies, ";\n");;

    Print("  Done.\n\n");;
    return rec(rank := r, commutative := is_comm);
end;;

# ─────────────────────────────────────────────
# Find the right subgroups and compute
# ─────────────────────────────────────────────

FindFirst := function(order, desc)
    local sc, H;
    for sc in subs do
        H := Representative(sc);;
        if Size(H) = order and StructureDescription(H) = desc then
            return H;
        fi;
    od;
    # Fallback: just match order
    for sc in subs do
        H := Representative(sc);;
        if Size(H) = order then return H; fi;
    od;
    return fail;
end;;

Print("═══════════════════════════════════════════════════\n");;
Print("Computing all five Meinhold Hecke algebra tables\n");;
Print("═══════════════════════════════════════════════════\n\n");;

ComputeAndSave(G, FindFirst(168, "PSL(3,2)"), "Meinhold-72");;
ComputeAndSave(G, FindFirst(24, "S4"), "Meinhold-504");;
ComputeAndSave(G, FindFirst(21, "C7 : C3"), "Meinhold-576");;
ComputeAndSave(G, FindFirst(9, "C3 x C3"), "Meinhold-1344");;
ComputeAndSave(G, FindFirst(6, "C6"), "Meinhold-2016");;

Print("═══ ALL DONE ═══\n");;
