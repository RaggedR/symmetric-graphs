# psl2_graphs.g
# Construct cubic arc-transitive graphs from PSL(2,q) via Sabidussi construction.
#
# Two families:
#   Split torus:     q ≡ 1 mod 3 → H = N(T_s), |V| = q(q+1)/2
#   Non-split torus: q ≡ 2 mod 3 → H = N(T_ns), |V| = q(q-1)/2
#
# Usage: gap psl2_graphs.g

LoadPackage("grape");

# ============================================================
# Core: build the Sabidussi graph Sab(G, H, HaH)
# ============================================================

# Given a group G, subgroup H, and involution a (with a^2 in H, a not in H),
# build the coset graph with vertices = right cosets of H,
# and edges Hx ~ Hy iff xy^{-1} in HaH.
SabidussiGraph := function(G, H, a)
    local cosets, n, adj, i, j, dcoset, x, y;

    cosets := RightCosets(G, H);
    n := Length(cosets);
    adj := List([1..n], i -> []);

    # Precompute the double coset HaH
    dcoset := DoubleCoset(H, a, H);

    for i in [1..n] do
        x := Representative(cosets[i]);
        for j in [i+1..n] do
            y := Representative(cosets[j]);
            if x * y^(-1) in dcoset then
                Add(adj[i], j);
                Add(adj[j], i);
            fi;
        od;
    od;

    for i in [1..n] do adj[i] := Set(adj[i]); od;

    return rec(
        n := n,
        adj := adj,
        cosets := cosets,
        valency := Length(adj[1])
    );
end;


# ============================================================
# PSL(2,q) graph construction
# ============================================================

# Build cubic graph from PSL(2,q) using the normalizer of a maximal torus.
# For q ≡ 1 mod 3: uses split torus (stabilizer of 2 points on P^1)
# For q ≡ 2 mod 3: uses non-split torus (from F_{q^2})
# Returns: rec(n, adj, valency, name, aut_order, aut_structure, ...)
PSL2Graph := function(q)
    local G, H, a, result, graph, aut, torus_type,
          field, one, mat, gen_t, gen_w, T, NH,
          involutions, inv, best_a, best_val, val, dcoset,
          psl, candidates, c, n_vertices;

    if not IsPrimeInt(q) or q < 5 then
        Print("Need q >= 5 prime\n");
        return fail;
    fi;

    if q mod 3 = 0 then
        Print("q = ", q, " is divisible by 3, neither family applies cleanly\n");
        return fail;
    fi;

    # Build PSL(2,q) as a permutation group acting on P^1(F_q) = {1,...,q+1}
    # GAP's PSL function gives a permutation group
    G := PSL(2, q);

    Print("PSL(2,", q, "): order ", Size(G), "\n");

    if q mod 3 = 1 then
        torus_type := "split";
        n_vertices := q * (q + 1) / 2;
        Print("Split torus family (q ≡ 1 mod 3)\n");
        Print("Expected vertices: ", n_vertices, "\n");

        # The split torus in PSL(2,q) acting on P^1 = {1,...,q+1}:
        # It's the stabilizer of a pair of points, e.g., {1, 2} (first two points).
        # Its normalizer is the setwise stabilizer of {1, 2}.
        H := Stabilizer(G, Set([1, 2]), OnSets);

    else  # q mod 3 = 2
        torus_type := "non-split";
        n_vertices := q * (q - 1) / 2;
        Print("Non-split torus family (q ≡ 2 mod 3)\n");
        Print("Expected vertices: ", n_vertices, "\n");

        # For the non-split torus, we need the stabilizer of an "imaginary pair"
        # — a pair of points in P^1(F_{q^2}) not in P^1(F_q).
        #
        # In permutation group terms: the non-split torus acts fixed-point-freely
        # on P^1(F_q). Its normalizer is a dihedral group D_{(q+1)/2} of order q+1.
        #
        # Strategy: find a cyclic subgroup of order (q+1)/2 in G, then take normalizer.
        # Elements of order dividing (q+1)/2 come from the non-split torus.

        # Find an element of order (q+1)/2 (or (q+1) if gcd issues)
        T := fail;
        for c in ConjugacyClasses(G) do
            if Order(Representative(c)) = (q + 1) / 2 then
                T := Group(Representative(c));
                break;
            fi;
        od;

        if T = fail then
            # Try order q+1 (might happen for small q)
            for c in ConjugacyClasses(G) do
                if Order(Representative(c)) = (q + 1) then
                    T := Group(Representative(c));
                    break;
                fi;
            od;
        fi;

        if T = fail then
            Print("Could not find non-split torus element for q=", q, "\n");
            # Fall back: search for elements whose order divides (q+1)/2
            for c in ConjugacyClasses(G) do
                if (q + 1) / 2 mod Order(Representative(c)) = 0 and
                   Order(Representative(c)) > 1 and
                   Length(Orbits(Group(Representative(c)), [1..q+1])) = 2 then
                    # Fixed-point-free with 2 orbits = non-split torus element
                    T := Group(Representative(c));
                    break;
                fi;
            od;
        fi;

        if T = fail then
            Print("FAILED to find non-split torus for q=", q, "\n");
            return fail;
        fi;

        H := Normalizer(G, T);
        Print("Non-split torus normalizer found: |H| = ", Size(H), "\n");
    fi;

    Print("|H| = ", Size(H), " (", StructureDescription(H), ")\n");
    Print("|G/H| = ", Size(G) / Size(H), " vertices\n");

    # Now find involutions a not in H with [H : a^{-1}Ha ∩ H] = 3
    # All involutions in PSL(2,q) are conjugate (for q odd), so we just
    # need to find one outside H with the right intersection index.

    best_a := fail;
    best_val := 0;

    # Get all involutions (elements of order 2)
    involutions := [];
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 2 then
            Append(involutions, AsList(c));
            break;  # all involutions are conjugate in PSL(2,q) for q odd
        fi;
    od;

    Print("Number of involutions: ", Length(involutions), "\n");

    # Test involutions outside H
    for inv in involutions do
        if not inv in H then
            dcoset := DoubleCoset(H, inv, H);
            val := Size(dcoset) / Size(H);
            if val = 3 then
                best_a := inv;
                best_val := 3;
                Print("Found cubic involution! |HaH|/|H| = 3\n");
                break;
            fi;
        fi;
    od;

    if best_a = fail then
        Print("No cubic involution found. Trying all valencies:\n");
        for inv in involutions do
            if not inv in H then
                dcoset := DoubleCoset(H, inv, H);
                val := Size(dcoset) / Size(H);
                Print("  valency = ", val, "\n");
                if best_a = fail then
                    best_a := inv;
                    best_val := val;
                fi;
                break;
            fi;
        od;
    fi;

    if best_a = fail then
        Print("No suitable involution found!\n");
        return fail;
    fi;

    Print("\nBuilding Sabidussi graph with valency ", best_val, "...\n");
    result := SabidussiGraph(G, H, best_a);

    Print("Built graph: ", result.n, " vertices, valency ", result.valency, "\n");

    # Compute automorphism group via GRAPE
    Print("Computing automorphism group...\n");
    graph := Graph(Group(()), [1..result.n],
        function(x,g) return x; end,
        function(x,y) return y in result.adj[x]; end, true);
    aut := AutGroupGraph(graph);

    result.name := Concatenation("PSL2_", String(q));
    result.q := q;
    result.torus_type := torus_type;
    result.aut_order := Size(aut);
    result.aut_structure := StructureDescription(aut);
    result.psl_order := Size(G);
    result.stabilizer_order := Size(H);
    result.stabilizer_structure := StructureDescription(H);

    Print("\n=== RESULT ===\n");
    Print("PSL(2,", q, ") / ", StructureDescription(H), "\n");
    Print("Torus type: ", torus_type, "\n");
    Print("Vertices: ", result.n, "\n");
    Print("Valency: ", result.valency, "\n");
    Print("|Aut(Γ)| = ", result.aut_order, "  ", result.aut_structure, "\n");
    Print("|PSL(2,", q, ")| = ", result.psl_order, "\n");
    Print("Arc-transitive: ", result.psl_order >= result.n * result.valency, "\n");
    Print("==============\n\n");

    return result;
end;


# ============================================================
# Main: scan over small primes
# ============================================================

Print("========================================\n");
Print("PSL(2,q) cubic arc-transitive graphs\n");
Print("========================================\n\n");

PSL2Results := [];

for q in [5, 7, 11, 13, 17, 19, 23, 29, 31] do
    if q mod 3 <> 0 then
        Print("────────────────────────────────\n");
        Print("q = ", q, "\n");
        Print("────────────────────────────────\n");
        Add(PSL2Results, PSL2Graph(q));
        Print("\n");
    fi;
od;

# Summary table
Print("\n\n");
Print("╔══════╦═══════════╦═════════╦═════════╦═══════════════════════════════╗\n");
Print("║  q   ║ |V|       ║ valency ║ torus   ║ Aut(Γ)                        ║\n");
Print("╠══════╬═══════════╬═════════╬═════════╬═══════════════════════════════╣\n");
for result in PSL2Results do
    if result <> fail then
        Print("║ ", String(result.q, 4), " ║ ",
              String(result.n, 9), " ║ ",
              String(result.valency, 7), " ║ ",
              String(result.torus_type, 7), " ║ ",
              result.aut_structure, "\n");
    fi;
od;
Print("╚══════╩═══════════╩═════════╩═════════╩═══════════════════════════════╝\n");
