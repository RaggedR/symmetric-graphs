# symmetric-layout.g
# General tool: find subgroups of Aut(Γ) whose orbit quotients have low degree
# Output JSON for the Python layout engine
#
# Usage: gap -q symmetric-layout.g -- graphname [max_degree]
# e.g.:  gap -q symmetric-layout.g -- dodecahedron 3

LoadPackage("grape");
LoadPackage("json", false);  # may not be available, we'll format manually

# ============================================================
# Graph library
# ============================================================

GraphLibrary := rec();

GraphLibrary.petersen := function()
    return rec(
        name := "Petersen",
        n := 10,
        adj := [
            [2,5,6], [1,3,7], [2,4,8], [3,5,9], [1,4,10],
            [1,8,9], [2,9,10], [3,6,10], [4,6,7], [5,7,8]
        ]
    );
end;

GraphLibrary.dodecahedron := function()
    return rec(
        name := "Dodecahedron",
        n := 20,
        adj := [
            [2,5,6], [1,3,8], [2,4,10], [3,5,12], [1,4,14],
            [1,7,15], [6,8,16], [2,7,9], [8,10,17], [3,9,11],
            [10,12,18], [4,11,13], [12,14,19], [5,13,15], [6,14,20],
            [7,17,20], [9,16,18], [11,17,19], [13,18,20], [15,16,19]
        ]
    );
end;

GraphLibrary.cube := function()
    return rec(
        name := "Cube",
        n := 8,
        adj := [
            [2,4,5], [1,3,6], [2,4,7], [1,3,8],
            [1,6,8], [2,5,7], [3,6,8], [4,5,7]
        ]
    );
end;

GraphLibrary.icosahedron := function()
    return rec(
        name := "Icosahedron",
        n := 12,
        adj := [
            [2,3,4,5,6], [1,3,6,7,8], [1,2,4,8,9], [1,3,5,9,10],
            [1,4,6,10,11], [1,2,5,7,11], [2,6,8,11,12], [2,3,7,9,12],
            [3,4,8,10,12], [4,5,9,11,12], [5,6,7,10,12], [7,8,9,10,11]
        ]
    );
end;

GraphLibrary.heawood := function()
    # Heawood graph: incidence graph of Fano plane, 14 vertices, cubic, bipartite
    return rec(
        name := "Heawood",
        n := 14,
        adj := [
            [2,6,14], [1,3,11], [2,4,8], [3,5,13], [4,6,10],
            [1,5,7], [6,8,12], [3,7,9], [8,10,14], [5,9,11],
            [2,10,12], [7,11,13], [4,12,14], [1,9,13]
        ]
    );
end;

GraphLibrary.pappus := function()
    # Pappus graph: 18 vertices, cubic, bipartite, 3-arc-transitive
    # LCF notation: [5,7,-7,7,-7,-5]^3, |Aut| = 216
    return rec(
        name := "Pappus",
        n := 18,
        adj := [
            [2,6,18], [1,3,9], [2,4,14], [3,5,11], [4,6,16],
            [1,5,7], [6,8,12], [7,9,15], [2,8,10], [9,11,17],
            [4,10,12], [7,11,13], [12,14,18], [3,13,15], [8,14,16],
            [5,15,17], [10,16,18], [1,13,17]
        ]
    );
end;

GraphLibrary.desargues := function()
    # Desargues graph: 20 vertices, cubic, bipartite
    # It is the generalised Petersen graph GP(10,3)
    return rec(
        name := "Desargues",
        n := 20,
        adj := [
            [2,10,11], [1,3,12], [2,4,13], [3,5,14], [4,6,15],
            [5,7,16], [6,8,17], [7,9,18], [8,10,19], [1,9,20],
            [1,14,18], [2,15,19], [3,16,20], [4,11,17], [5,12,18],
            [6,13,19], [7,14,20], [8,11,15], [9,12,16], [10,13,17]
        ]
    );
end;

GraphLibrary.moebiuskantor := function()
    # Moebius-Kantor graph: generalised Petersen GP(8,3), 16 vertices, cubic
    return rec(
        name := "Moebius-Kantor",
        n := 16,
        adj := [
            [2,8,9], [1,3,10], [2,4,11], [3,5,12], [4,6,13],
            [5,7,14], [6,8,15], [1,7,16], [1,12,14], [2,13,15],
            [3,14,16], [4,9,15], [5,10,16], [6,9,11], [7,10,12],
            [8,11,13]
        ]
    );
end;

GraphLibrary.k33 := function()
    # Complete bipartite K_{3,3}
    return rec(
        name := "K33",
        n := 6,
        adj := [
            [4,5,6], [4,5,6], [4,5,6],
            [1,2,3], [1,2,3], [1,2,3]
        ]
    );
end;

# ============================================================
# Core functions
# ============================================================

# Build GRAPE graph from adjacency list
BuildGraph := function(data)
    return Graph(Group(()), [1..data.n], function(x,g) return x; end,
        function(x,y) return y in data.adj[x]; end, true);
end;

# Compute quotient from a list of orbits
ComputeQuotient := function(orbits, adj, n)
    local n_orbs, quot_adj, i, j, found, u, v, valencies, max_val,
          edge_counts, internal_counts;

    n_orbs := Length(orbits);
    quot_adj := List([1..n_orbs], i -> []);
    edge_counts := List([1..n_orbs], i -> List([1..n_orbs], j -> 0));
    internal_counts := List([1..n_orbs], i -> 0);

    for i in [1..n_orbs] do
        # Internal edges
        for u in orbits[i] do
            for v in orbits[i] do
                if u < v and v in adj[u] then
                    internal_counts[i] := internal_counts[i] + 1;
                fi;
            od;
        od;
        # Cross edges
        for j in [i+1..n_orbs] do
            found := 0;
            for u in orbits[i] do
                for v in orbits[j] do
                    if v in adj[u] then
                        found := found + 1;
                    fi;
                od;
            od;
            edge_counts[i][j] := found;
            edge_counts[j][i] := found;
            if found > 0 then
                Add(quot_adj[i], j);
                Add(quot_adj[j], i);
            fi;
        od;
    od;

    valencies := List(quot_adj, Length);
    max_val := Maximum(valencies);

    return rec(
        n_orbits := n_orbs,
        adj := quot_adj,
        valencies := valencies,
        max_valency := max_val,
        edge_counts := edge_counts,
        internal_counts := internal_counts
    );
end;

# Analyse angular offsets between orbits
# For each pair of adjacent orbits, determine how cross-edges connect positions
AnalyseOffsets := function(orbits, adj, generator, n)
    local n_orbs, offsets, i, j, orbit_order, pos_of, k, v, seq, start,
          u, w, pos_u, pos_w, diff, diffs, mode_diff;

    n_orbs := Length(orbits);

    # For each orbit, determine the ordering under the generator
    orbit_order := [];
    pos_of := List([1..n], x -> rec(orbit := 0, pos := 0));

    for i in [1..n_orbs] do
        start := Minimum(orbits[i]);
        seq := [start];
        v := start ^ generator;
        while v <> start do
            Add(seq, v);
            v := v ^ generator;
        od;
        Add(orbit_order, seq);
        for k in [1..Length(seq)] do
            pos_of[seq[k]] := rec(orbit := i, pos := k - 1);
        od;
    od;

    # For each pair of adjacent orbits, find the position offsets
    offsets := List([1..n_orbs], i -> List([1..n_orbs], j -> 0));
    for i in [1..n_orbs] do
        for j in [i+1..n_orbs] do
            diffs := [];
            for u in orbits[i] do
                for w in orbits[j] do
                    if w in adj[u] then
                        pos_u := pos_of[u].pos;
                        pos_w := pos_of[w].pos;
                        diff := (pos_w - pos_u) mod Length(orbits[i]);
                        Add(diffs, diff);
                    fi;
                od;
            od;
            if Length(diffs) > 0 then
                # Find most common offset (mode)
                Sort(diffs);
                mode_diff := diffs[1];  # simple: take smallest
                offsets[i][j] := mode_diff;
                offsets[j][i] := (-mode_diff) mod Length(orbits[i]);
            fi;
        od;
    od;

    return rec(
        orbit_order := orbit_order,
        pos_of := pos_of,
        offsets := offsets
    );
end;

# ============================================================
# Search for good quotients
# ============================================================

SearchQuotients := function(data, max_degree)
    local graph, aut, adj, n, results, seen, ProcessSubgroup,
          norms, conj_classes, cc, rep, H, g, ord, orders_seen, cyc;

    n := data.n;
    adj := data.adj;
    graph := BuildGraph(data);
    aut := AutGroupGraph(graph);

    Print("Graph: ", data.name, " (", n, " vertices)\n");
    Print("|Aut| = ", Size(aut), "  ", StructureDescription(aut), "\n\n");

    results := [];
    seen := [];  # track orbit sets we've already processed

    ProcessSubgroup := function(H, label)
        local orbs, orb_key, quot, offset_data, orbit_sizes, gen,
              has_pentagon, i, seq, start, v, is_cycle;

        if Size(H) = 1 or Size(H) = Size(aut) then return; fi;

        orbs := Orbits(H, [1..n]);
        orbs := List(orbs, Set);
        Sort(orbs, function(a,b) return Minimum(a) < Minimum(b); end);

        orb_key := orbs;
        if orb_key in seen then return; fi;
        Add(seen, orb_key);

        orbit_sizes := List(orbs, Length);

        quot := ComputeQuotient(orbs, adj, n);

        if quot.max_valency > max_degree then return; fi;
        # Allow 2 orbits if they're equal-sized (generalized Petersen pattern)
        if quot.n_orbits < 2 then return; fi;
        if quot.n_orbits = 2 and Length(Set(orbit_sizes)) > 1 then return; fi;

        # Check if quotient is a cycle
        is_cycle := ForAll(quot.valencies, v -> v = 2);

        # Find the highest-order element of H that preserves each orbit
        # (needed for cycle detection and offset analysis)
        gen := fail;
        for g in H do
            if Order(g) > 1 and ForAll(orbs, o -> Set(List(o, x -> x^g)) = o) then
                if gen = fail or Order(g) > Order(gen) then
                    gen := g;
                fi;
            fi;
        od;

        # Check which orbits form cycles under any generator of H
        has_pentagon := List([1..Length(orbs)], i -> false);
        if gen <> fail then
            for i in [1..Length(orbs)] do
                start := Minimum(orbs[i]);
                seq := [start];
                v := start ^ gen;
                while v <> start do
                    Add(seq, v);
                    v := v ^ gen;
                od;
                # Check if consecutive elements are adjacent
                if Length(seq) = Length(orbs[i]) then
                    has_pentagon[i] := ForAll([1..Length(seq)], function(k)
                        local next;
                        next := seq[(k mod Length(seq)) + 1];
                        return next in adj[seq[k]];
                    end);
                fi;
            od;
        fi;

        Add(results, rec(
            label := label,
            subgroup_order := Size(H),
            subgroup_structure := StructureDescription(H),
            n_orbits := quot.n_orbits,
            orbit_sizes := orbit_sizes,
            orbits := orbs,
            quotient_adj := quot.adj,
            quotient_valencies := quot.valencies,
            max_valency := quot.max_valency,
            edge_counts := quot.edge_counts,
            internal_counts := quot.internal_counts,
            is_cycle := is_cycle,
            has_internal_cycle := has_pentagon,
            generator := gen
        ));

        Print("  ", label, ": |H|=", Size(H), " (", StructureDescription(H), ")");
        Print(", ", quot.n_orbits, " orbits ", orbit_sizes);
        Print(", max_deg=", quot.max_valency);
        if is_cycle then Print(" [CYCLE]"); fi;
        Print("\n");
    end;

    # 1. Normal subgroups
    Print("--- Normal subgroups ---\n");
    norms := NormalSubgroups(aut);
    for H in norms do
        ProcessSubgroup(H, Concatenation("Normal |", String(Size(H)), "|"));
    od;

    # 2. Conjugacy classes of subgroups (for small groups)
    if Size(aut) <= 1000 then
        Print("\n--- Conjugacy classes of subgroups ---\n");
        conj_classes := ConjugacyClassesSubgroups(aut);
        for cc in conj_classes do
            rep := Representative(cc);
            ProcessSubgroup(rep, Concatenation("Conj |", String(Size(rep)), "| ",
                StructureDescription(rep)));
        od;
    else
        # For large groups, just try cyclic subgroups
        Print("\n--- Cyclic subgroups (group too large for full search) ---\n");
        orders_seen := [];
        for g in aut do
            ord := Order(g);
            if ord >= 2 and not ord in orders_seen then
                Add(orders_seen, ord);
                cyc := Group(g);
                ProcessSubgroup(cyc, Concatenation("Cyclic Z", String(ord)));
            fi;
        od;
    fi;

    return results;
end;

# ============================================================
# Output JSON
# ============================================================

OutputJSON := function(data, results)
    local r, i, j, offset_data;

    Print("\n\n=== JSON OUTPUT ===\n");
    Print("{\n");
    Print("  \"graph\": \"", data.name, "\",\n");
    Print("  \"n\": ", data.n, ",\n");
    Print("  \"adj\": ", data.adj, ",\n");
    Print("  \"results\": [\n");

    for i in [1..Length(results)] do
        r := results[i];
        Print("    {\n");
        Print("      \"label\": \"", r.label, "\",\n");
        Print("      \"subgroup_order\": ", r.subgroup_order, ",\n");
        Print("      \"subgroup_structure\": \"", r.subgroup_structure, "\",\n");
        Print("      \"n_orbits\": ", r.n_orbits, ",\n");
        Print("      \"orbit_sizes\": ", r.orbit_sizes, ",\n");
        Print("      \"orbits\": ", r.orbits, ",\n");
        Print("      \"quotient_adj\": ", r.quotient_adj, ",\n");
        Print("      \"quotient_valencies\": ", r.quotient_valencies, ",\n");
        Print("      \"max_valency\": ", r.max_valency, ",\n");
        Print("      \"edge_counts\": ", r.edge_counts, ",\n");
        Print("      \"internal_counts\": ", r.internal_counts, ",\n");
        Print("      \"is_cycle\": ", r.is_cycle, ",\n");
        Print("      \"has_internal_cycle\": ", r.has_internal_cycle, "\n");

        # Compute offsets if we have a generator
        if r.generator <> fail then
            offset_data := AnalyseOffsets(r.orbits, data.adj, r.generator, data.n);
            Print("      ,\"orbit_order\": ", offset_data.orbit_order, "\n");
            Print("      ,\"offsets\": ", offset_data.offsets, "\n");
        fi;

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
# Main
# ============================================================

# Parse arguments — set these before running
# Override by defining graph_name and max_deg before Read()ing this file
if not IsBound(graph_name) then
    graph_name := "dodecahedron";
fi;
if not IsBound(max_deg) then
    max_deg := 4;
fi;

if not IsBound(GraphLibrary.(graph_name)) then
    Print("Unknown graph: ", graph_name, "\n");
    Print("Available: ");
    for graph_name in RecNames(GraphLibrary) do
        Print(graph_name, " ");
    od;
    Print("\n");
else

data := GraphLibrary.(graph_name)();
Print("Searching for quotients of ", data.name, " with max degree <= ", max_deg, "\n\n");
results := SearchQuotients(data, max_deg);

Print("\n=== SUMMARY: ", Length(results), " quotients found with max degree <= ", max_deg, " ===\n\n");
for i in [1..Length(results)] do
    r := results[i];
    Print(i, ". ", r.label, "\n");
    Print("   ", r.n_orbits, " orbits of sizes ", r.orbit_sizes, "\n");
    Print("   Quotient degree sequence: ", r.quotient_valencies, "\n");
    Print("   Internal cycles: ", r.has_internal_cycle, "\n");
    if r.is_cycle then Print("   ** Quotient is a cycle! **\n"); fi;
    Print("\n");
od;

OutputJSON(data, results);

fi; # end if graph found
