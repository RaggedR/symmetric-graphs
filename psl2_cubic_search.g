# psl2_cubic_search.g
# Systematic search for cubic arc-transitive graphs from PSL(2,q)
# Try all conjugacy classes of subgroups H, all involutions a,
# require: valency 3, connected, a^2 = 1, a not in H.

LoadPackage("grape");

SabidussiGraphFull := function(G, H, a)
    local cosets, n, adj, i, j, dcoset, x, y;
    cosets := RightCosets(G, H);
    n := Length(cosets);
    adj := List([1..n], i -> []);
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
    return rec(n := n, adj := adj);
end;

# Check if an adjacency list gives a connected graph via BFS
IsConnectedAdj := function(adj, n)
    local visited, queue, v, w, count;
    if n = 0 then return true; fi;
    visited := BlistList([1..n], [1]);
    queue := [1];
    count := 1;
    while Length(queue) > 0 do
        v := Remove(queue, 1);
        for w in adj[v] do
            if not visited[w] then
                visited[w] := true;
                Add(queue, w);
                count := count + 1;
            fi;
        od;
    od;
    return count = n;
end;

# Compute girth of graph via BFS from each vertex
ComputeGirth := function(adj, n)
    local girth, v, dist, queue, u, w, d;
    girth := infinity;
    for v in [1..Minimum(n, 50)] do  # sample first 50 vertices for speed
        dist := List([1..n], i -> -1);
        dist[v] := 0;
        queue := [v];
        while Length(queue) > 0 do
            u := Remove(queue, 1);
            for w in adj[u] do
                if dist[w] = -1 then
                    dist[w] := dist[u] + 1;
                    Add(queue, w);
                elif w <> v and dist[w] >= dist[u] then
                    # Found a cycle of length dist[u] + dist[w] + 1
                    d := dist[u] + dist[w] + 1;
                    if d < girth then girth := d; fi;
                fi;
            od;
            if dist[u] >= girth/2 then break; fi;
        od;
    od;
    return girth;
end;


# ============================================================
# Main search
# ============================================================

Print("╔══════════════════════════════════════════════════════════╗\n");
Print("║  Cubic arc-transitive graphs from PSL(2,q)             ║\n");
Print("║  Searching all subgroups H and involutions a            ║\n");
Print("╚══════════════════════════════════════════════════════════╝\n\n");

all_results := [];

for q in [5, 7, 11, 13, 17, 19, 23] do
    Print("════════════════════════════════════════\n");
    Print("  q = ", q, "   PSL(2,", q, ")   |G| = ", Size(PSL(2,q)), "\n");
    Print("════════════════════════════════════════\n");

    G := PSL(2, q);

    # Get conjugacy classes of subgroups
    cc_subs := ConjugacyClassesSubgroups(G);

    # Get involutions (all conjugate for q odd)
    involutions := [];
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 2 then
            involutions := AsList(c);
            break;
        fi;
    od;
    Print("  Involutions: ", Length(involutions), "\n");

    found_any := false;

    for cc in cc_subs do
        H := Representative(cc);
        h_order := Size(H);

        # For cubic: need |H|/|intersection| = 3, so 3 divides |H|
        if h_order mod 3 <> 0 then continue; fi;
        # Skip trivial and full group
        if h_order <= 2 or h_order = Size(G) then continue; fi;
        # Skip if too many vertices (> 500)
        if Size(G) / h_order > 500 then continue; fi;

        # Try involutions outside H
        tested := 0;
        for a in involutions do
            if a in H then continue; fi;

            dcoset := DoubleCoset(H, a, H);
            val := Size(dcoset) / h_order;

            if val = 3 then
                # Build the graph and check connectivity
                result := SabidussiGraphFull(G, H, a);

                if IsConnectedAdj(result.adj, result.n) then
                    # Connected cubic graph! Compute properties.
                    graph := Graph(Group(()), [1..result.n],
                        function(x,g) return x; end,
                        function(x,y) return y in result.adj[x]; end, true);
                    aut := AutGroupGraph(graph);
                    girth := ComputeGirth(result.adj, result.n);

                    Print("  ★ FOUND: H = ", StructureDescription(H),
                          " (order ", h_order, ")\n");
                    Print("    |V| = ", result.n, ", valency = 3",
                          ", girth = ", girth, "\n");
                    Print("    |Aut| = ", Size(aut), "  ",
                          StructureDescription(aut), "\n");
                    Print("    |G|/|V| = ", Size(G)/result.n,
                          " (arc-stabilizer index)\n\n");

                    Add(all_results, rec(
                        q := q,
                        n := result.n,
                        valency := 3,
                        girth := girth,
                        h_structure := StructureDescription(H),
                        h_order := h_order,
                        aut_order := Size(aut),
                        aut_structure := StructureDescription(aut),
                        adj := result.adj
                    ));

                    found_any := true;
                    break;  # found one for this H, move to next subgroup class
                else
                    Print("  · H = ", StructureDescription(H),
                          " (order ", h_order, "): cubic but DISCONNECTED (",
                          result.n, " vertices)\n");
                    break;  # this H doesn't work
                fi;
            fi;

            tested := tested + 1;
            if tested >= 10 then break; fi;  # sample at most 10 involutions per H
        od;
    od;

    if not found_any then
        Print("  No connected cubic graph found for q = ", q, "\n\n");
    fi;
od;

# ============================================================
# Summary
# ============================================================

Print("\n\n");
Print("╔══════╦══════╦═════════╦═══════╦══════════════════╦═══════════════════════════╗\n");
Print("║  q   ║ |V|  ║ girth   ║ |Aut| ║ H                ║ Aut(Γ)                    ║\n");
Print("╠══════╬══════╬═════════╬═══════╬══════════════════╬═══════════════════════════╣\n");
for r in all_results do
    Print("║ ", String(r.q, 4), " ║ ",
          String(r.n, 4), " ║ ",
          String(r.girth, 7), " ║ ",
          String(r.aut_order, 5), " ║ ",
          String(r.h_structure, -16), " ║ ",
          r.aut_structure, "\n");
od;
Print("╚══════╩══════╩═════════╩═══════╩══════════════════╩═══════════════════════════╝\n");

# Output adjacency lists for drawing
for r in all_results do
    Print("\n=== ADJ PSL2_", r.q, " ===\n");
    Print("n := ", r.n, ";\n");
    Print("adj := ", r.adj, ";\n");
    Print("=== END ADJ ===\n");
od;
