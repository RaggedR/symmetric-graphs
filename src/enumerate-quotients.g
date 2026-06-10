# enumerate-quotients.g
# Exhaustive enumeration of ALL quotients of cubic arc-transitive graphs (n <= 30)
#
# For each graph Γ in the Foster census, computes Aut(Γ), enumerates ALL
# conjugacy classes of subgroups H ≤ Aut(Γ), and for each computes the
# quotient graph Γ/H.
#
# Two conjugate subgroups give isomorphic quotients (they just relabel vertices),
# so one representative per conjugacy class suffices.
#
# Usage: gap -q enumerate-quotients.g
#        gap -q enumerate-quotients.g -- petersen    (single graph)

skip_main := true;
Read("symmetric-layout.g");

# ============================================================
# Quotient graph identification
# ============================================================

# Count edges in an adjacency list
CountEdges := function(qadj)
    local count, i, j;
    count := 0;
    for i in [1..Length(qadj)] do
        for j in qadj[i] do
            if j > i then count := count + 1; fi;
        od;
    od;
    return count;
end;

# Identify the quotient graph type
IdentifyQuotient := function(n_orbs, qadj, valencies)
    local n_edges, max_val, min_val, is_connected, visited, queue, v, w,
          deg1_count, deg2_count, parts, bipartite, color, ok, i;

    n_edges := CountEdges(qadj);
    max_val := Maximum(valencies);
    min_val := Minimum(valencies);

    # Single vertex
    if n_orbs = 1 then return "K1"; fi;

    # Single edge
    if n_orbs = 2 and n_edges = 1 then return "K2"; fi;

    # Check connectivity via BFS
    visited := BlistList([1..n_orbs], [1]);
    queue := [1];
    while Length(queue) > 0 do
        v := Remove(queue, 1);
        for w in qadj[v] do
            if not visited[w] then
                visited[w] := true;
                Add(queue, w);
            fi;
        od;
    od;
    is_connected := SizeBlist(visited) = n_orbs;

    # Complete graph: every vertex has degree n-1
    if min_val = n_orbs - 1 then
        return Concatenation("K", String(n_orbs));
    fi;

    # All degree 2 and connected → cycle
    if min_val = 2 and max_val = 2 and is_connected then
        return Concatenation("C", String(n_orbs));
    fi;

    # Path: connected, all deg ≤ 2, exactly 2 vertices of degree 1
    deg1_count := Length(Filtered(valencies, v -> v = 1));
    deg2_count := Length(Filtered(valencies, v -> v = 2));
    if is_connected and deg1_count = 2 and deg2_count = n_orbs - 2 and max_val <= 2 then
        return Concatenation("P", String(n_orbs));
    fi;

    # Check bipartite (for complete bipartite identification)
    color := List([1..n_orbs], i -> -1);
    color[1] := 0;
    queue := [1];
    bipartite := true;
    while Length(queue) > 0 and bipartite do
        v := Remove(queue, 1);
        for w in qadj[v] do
            if color[w] = -1 then
                color[w] := 1 - color[v];
                Add(queue, w);
            elif color[w] = color[v] then
                bipartite := false;
            fi;
        od;
    od;

    if bipartite and is_connected then
        parts := [Length(Filtered(color, c -> c = 0)),
                  Length(Filtered(color, c -> c = 1))];
        Sort(parts);
        # Complete bipartite: every vertex in one part adjacent to all in the other
        ok := true;
        for i in [1..n_orbs] do
            if color[i] = 0 and valencies[i] <> parts[2] then ok := false; fi;
            if color[i] = 1 and valencies[i] <> parts[1] then ok := false; fi;
        od;
        if ok then
            if parts[1] = 1 then
                return Concatenation("K1,", String(parts[2]), " (star)");
            else
                return Concatenation("K", String(parts[1]), ",", String(parts[2]));
            fi;
        fi;
    fi;

    # Petersen graph: 10 vertices, 15 edges, 3-regular
    if n_orbs = 10 and n_edges = 15 and min_val = 3 and max_val = 3 then
        return "Petersen?";
    fi;

    # Generic description
    if is_connected then
        return Concatenation(String(n_orbs), "-v ", String(n_edges), "-e");
    else
        return Concatenation(String(n_orbs), "-v ", String(n_edges),
                            "-e (disconnected)");
    fi;
end;

# ============================================================
# Exhaustive quotient enumeration
# ============================================================

EnumerateAllQuotients := function(data)
    local graph, aut, adj, n, conj_classes, results, seen,
          cc, H, orbs, orb_key, quot, orbit_sizes, qtype,
          n_edges, i;

    n := data.n;
    adj := data.adj;
    graph := BuildGraph(data);
    aut := AutGroupGraph(graph);
    data.aut_order := Size(aut);
    data.aut_structure := StructureDescription(aut);

    Print("================================================================\n");
    Print("  ", data.name, "  (n=", n, ",  |Aut| = ", Size(aut),
          ",  ", StructureDescription(aut), ")\n");
    Print("================================================================\n\n");

    # Get ALL conjugacy classes of subgroups
    conj_classes := ConjugacyClassesSubgroups(aut);
    Print("  Total conjugacy classes of subgroups: ", Length(conj_classes), "\n\n");

    results := [];
    seen := [];

    for cc in conj_classes do
        H := Representative(cc);

        # Skip trivial (gives the graph itself) and full group (gives K1)
        if Size(H) = 1 or Size(H) = Size(aut) then continue; fi;

        orbs := Orbits(H, [1..n]);
        orbs := List(orbs, Set);
        Sort(orbs, function(a,b) return Minimum(a) < Minimum(b); end);

        # Deduplicate by orbit partition
        orb_key := orbs;
        if orb_key in seen then continue; fi;
        Add(seen, orb_key);

        orbit_sizes := List(orbs, Length);
        quot := ComputeQuotient(orbs, adj, n);
        n_edges := CountEdges(quot.adj);
        qtype := IdentifyQuotient(quot.n_orbits, quot.adj, quot.valencies);

        Add(results, rec(
            subgroup_order := Size(H),
            subgroup_structure := StructureDescription(H),
            n_class := Size(cc),
            n_orbits := quot.n_orbits,
            orbit_sizes := orbit_sizes,
            orbits := orbs,
            quotient_adj := quot.adj,
            quotient_valencies := quot.valencies,
            max_valency := quot.max_valency,
            n_edges := n_edges,
            quotient_type := qtype,
            edge_counts := quot.edge_counts,
            internal_counts := quot.internal_counts
        ));
    od;

    # Sort by number of orbits (descending = most orbits first), then by |H|
    Sort(results, function(a, b)
        if a.n_orbits <> b.n_orbits then return a.n_orbits > b.n_orbits; fi;
        return a.subgroup_order < b.subgroup_order;
    end);

    # Print results
    Print("  ", Length(results), " distinct quotients (excluding trivial and collapsed)\n\n");
    Print("  #  | |H|   | H structure           | orbits | orbit sizes         | Q type      | deg seq\n");
    Print("  ---|-------|------------------------|--------|---------------------|-------------|--------\n");

    for i in [1..Length(results)] do
        Print("  ", String(i, 2), " | ");
        Print(String(results[i].subgroup_order, 5), " | ");
        # Pad structure description
        Print(String(results[i].subgroup_structure, -22), " | ");
        Print(String(results[i].n_orbits, 6), " | ");
        # Orbit sizes — compact
        Print(String(results[i].orbit_sizes, -19), " | ");
        Print(String(results[i].quotient_type, -11), " | ");
        Print(results[i].quotient_valencies, "\n");
    od;
    Print("\n");

    return results;
end;

# ============================================================
# Quotient enumeration via vertex stabilizer overgroups
# ============================================================
#
# The correct graph-theoretic quotient: Γ = (G, H) where G = Aut(Γ),
# H = vertex stabilizer. A quotient is (G, H') for H < H' < G.
# This guarantees equal-sized blocks and a vertex-transitive quotient.

EnumerateQuotientsViaStabilizer := function(data)
    local graph, aut, adj, n, stab, cc, results, seen,
          c, K, g, Kg, block, block_system,
          orbs, orbit_sizes, quot, qtype, n_edges, trans;

    n := data.n;
    adj := data.adj;
    graph := BuildGraph(data);
    aut := AutGroupGraph(graph);

    data.aut_order := Size(aut);
    data.aut_structure := StructureDescription(aut);

    stab := Stabilizer(aut, 1);

    Print("================================================================\n");
    Print("  ", data.name, "  (n=", n, ",  |Aut| = ", Size(aut),
          ",  ", StructureDescription(aut), ")\n");
    Print("  Vertex stabilizer: ", StructureDescription(stab),
          " (order ", Size(stab), ")\n");
    Print("================================================================\n\n");

    cc := ConjugacyClassesSubgroups(aut);

    results := [];
    seen := [];

    for c in cc do
        K := Representative(c);

        # Skip if too small or too large
        if Size(K) <= Size(stab) then continue; fi;
        if Size(K) >= Size(aut) then continue; fi;

        # Find a conjugate of K that contains stab
        Kg := fail;
        trans := RightTransversal(aut, Normalizer(aut, K));
        for g in trans do
            if IsSubgroup(K^g, stab) then
                Kg := K^g;
                break;
            fi;
        od;

        if Kg = fail then continue; fi;

        # The candidate block = orbit of vertex 1 under Kg
        block := Set(Orbit(Kg, 1));

        # Verify it's a genuine block: G-translates must partition V
        block_system := Orbit(aut, block, OnSets);
        if Length(block_system) * Length(block) <> n then
            continue;  # Not a block system — skip
        fi;

        # Use the verified block system as our partition
        orbs := List(block_system, Set);
        Sort(orbs, function(a,b) return Minimum(a) < Minimum(b); end);

        # Deduplicate
        if orbs in seen then continue; fi;
        Add(seen, orbs);

        orbit_sizes := List(orbs, Length);
        quot := ComputeQuotient(orbs, adj, n);
        n_edges := CountEdges(quot.adj);
        qtype := IdentifyQuotient(quot.n_orbits, quot.adj, quot.valencies);

        Add(results, rec(
            subgroup_order := Size(Kg),
            subgroup_structure := StructureDescription(Kg),
            n_orbits := quot.n_orbits,
            orbit_sizes := orbit_sizes,
            orbits := orbs,
            quotient_adj := quot.adj,
            quotient_valencies := quot.valencies,
            max_valency := quot.max_valency,
            n_edges := n_edges,
            quotient_type := qtype,
            edge_counts := quot.edge_counts,
            internal_counts := quot.internal_counts
        ));

        Print("  H' = ", StructureDescription(Kg),
              " (order ", Size(Kg), "): ");
        Print(quot.n_orbits, " blocks of size ", orbit_sizes[1]);
        Print(", quotient: ", qtype, "\n");
    od;

    # Sort by number of orbits (descending), then |H'|
    Sort(results, function(a, b)
        if a.n_orbits <> b.n_orbits then return a.n_orbits > b.n_orbits; fi;
        return a.subgroup_order < b.subgroup_order;
    end);

    Print("\n  ", Length(results), " proper quotients\n\n");

    return results;
end;

# ============================================================
# JSON output for Python consumption
# ============================================================

OutputAllQuotientsJSON := function(data, results)
    local r, i;

    Print("=== JSON OUTPUT ===\n");
    Print("{\n");
    Print("  \"graph\": \"", data.name, "\",\n");
    Print("  \"n\": ", data.n, ",\n");
    Print("  \"aut_order\": ", data.aut_order, ",\n");
    Print("  \"aut_structure\": \"", data.aut_structure, "\",\n");
    Print("  \"adj\": ", data.adj, ",\n");
    Print("  \"results\": [\n");

    for i in [1..Length(results)] do
        r := results[i];
        Print("    {\n");
        Print("      \"subgroup_order\": ", r.subgroup_order, ",\n");
        Print("      \"subgroup_structure\": \"", r.subgroup_structure, "\",\n");
        Print("      \"n_orbits\": ", r.n_orbits, ",\n");
        Print("      \"orbit_sizes\": ", r.orbit_sizes, ",\n");
        Print("      \"orbits\": ", r.orbits, ",\n");
        Print("      \"quotient_adj\": ", r.quotient_adj, ",\n");
        Print("      \"quotient_valencies\": ", r.quotient_valencies, ",\n");
        Print("      \"edge_counts\": ", r.edge_counts, ",\n");
        Print("      \"internal_counts\": ", r.internal_counts, ",\n");
        Print("      \"quotient_type\": \"", r.quotient_type, "\"\n");
        if i < Length(results) then
            Print("    },\n");
        else
            Print("    }\n");
        fi;
    od;

    Print("  ]\n");
    Print("}\n");
    Print("=== END JSON ===\n");
end;

# ============================================================
# Main: enumerate all graphs or a single one
# ============================================================

# Ordered by size for clean output
all_graph_names := ["k4", "k33", "cube", "petersen", "heawood",
                    "moebiuskantor", "pappus", "dodecahedron", "desargues",
                    "nauru", "f26a", "coxeter", "tuttecoxeter"];

# Override by defining target_graphs before Read()ing this file
# e.g.: target_graphs := ["petersen"]; Read("enumerate-quotients.g");
if not IsBound(target_graphs) then
    target_graphs := all_graph_names;
fi;

Print("\n");
Print("╔══════════════════════════════════════════════════════════════════╗\n");
Print("║  EXHAUSTIVE QUOTIENT ENUMERATION                               ║\n");
Print("║  Cubic arc-transitive graphs, n ≤ 30 (Foster census)           ║\n");
Print("╚══════════════════════════════════════════════════════════════════╝\n\n");

all_data := rec();
all_results := rec();

for gname in target_graphs do
    if IsBound(GraphLibrary.(gname)) then
        data := GraphLibrary.(gname)();
        if IsBound(json_mode) and json_mode = true then
            all_results.(gname) := EnumerateQuotientsViaStabilizer(data);
        else
            all_results.(gname) := EnumerateAllQuotients(data);
        fi;
        all_data.(gname) := data;
    else
        Print("Unknown graph: ", gname, "\n");
    fi;
od;

# JSON output mode
if IsBound(json_mode) and json_mode = true then
    for gname in target_graphs do
        if IsBound(all_results.(gname)) then
            OutputAllQuotientsJSON(all_data.(gname), all_results.(gname));
        fi;
    od;
fi;

# Final summary
Print("\n");
Print("╔══════════════════════════════════════════════════════════════════╗\n");
Print("║  SUMMARY                                                       ║\n");
Print("╚══════════════════════════════════════════════════════════════════╝\n\n");
Print("  Graph            |  n  | |Aut|  | # subgp classes | # quotients\n");
Print("  -----------------|-----|--------|-----------------|------------\n");

for gname in target_graphs do
    if IsBound(GraphLibrary.(gname)) then
        data := GraphLibrary.(gname)();
        graph := BuildGraph(data);
        aut := AutGroupGraph(graph);
        conj_classes := ConjugacyClassesSubgroups(aut);
        Print("  ", String(data.name, -17), "| ");
        Print(String(data.n, 3), " | ");
        Print(String(Size(aut), 6), " | ");
        Print(String(Length(conj_classes), 15), " | ");
        if IsBound(all_results.(gname)) then
            Print(String(Length(all_results.(gname)), 10));
        fi;
        Print("\n");
    fi;
od;

Print("\n");
