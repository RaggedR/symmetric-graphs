# macbeath_graph.g
# Construct the Macbeath graph: the cubic arc-transitive graph on 168 vertices
# arising from the Macbeath surface (genus 7, Hurwitz surface)
#
# This is the 1-skeleton of the {7,3} tessellation of the Macbeath surface,
# constructed as the coset graph PSL(2,8)/C_3.
#
# Compare with the Klein graph (F056A): PSL(2,7)/C_3, 56 vertices.

LoadPackage("grape");

SabidussiGraph := function(G, H, a)
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

ComputeGirth := function(adj, n)
    local girth, v, dist, queue, u, w, d;
    girth := infinity;
    for v in [1..Minimum(n, 50)] do
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
# Build the Macbeath graph from PSL(2,8)
# ============================================================

Print("════════════════════════════════════════════════════════\n");
Print("  Macbeath graph: PSL(2,8)/C_3\n");
Print("  The cubic graph of the Macbeath surface (genus 7)\n");
Print("════════════════════════════════════════════════════════\n\n");

G := PSL(2, 8);
Print("G = PSL(2,8),  |G| = ", Size(G), "\n\n");

# Find all subgroups of order 3
cc_subs := ConjugacyClassesSubgroups(G);
Print("Conjugacy classes of subgroups: ", Length(cc_subs), "\n\n");

# Find involutions
involutions := [];
for c in ConjugacyClasses(G) do
    if Order(Representative(c)) = 2 then
        involutions := AsList(c);
        break;
    fi;
od;
Print("Involutions: ", Length(involutions), "\n\n");

# Search for cubic coset graphs
found := false;
for cc in cc_subs do
    H := Representative(cc);
    h_order := Size(H);

    if h_order <= 2 or h_order = Size(G) then continue; fi;

    for a in involutions do
        if a in H then continue; fi;

        dcoset := DoubleCoset(H, a, H);
        val := Size(dcoset) / h_order;

        if val = 3 then
            result := SabidussiGraph(G, H, a);

            if IsConnectedAdj(result.adj, result.n) then
                graph := Graph(Group(()), [1..result.n],
                    function(x,g) return x; end,
                    function(x,y) return y in result.adj[x]; end, true);
                aut := AutGroupGraph(graph);
                girth := ComputeGirth(result.adj, result.n);

                Print("★ FOUND: H = ", StructureDescription(H),
                      " (order ", h_order, ")\n");
                Print("  |V| = ", result.n, ", valency = 3",
                      ", girth = ", girth, "\n");
                Print("  |Aut| = ", Size(aut), "  ",
                      StructureDescription(aut), "\n");
                Print("  Connected: ", IsConnectedAdj(result.adj, result.n), "\n\n");

                if result.n = 168 then
                    Print("  *** This is the Macbeath graph! ***\n\n");

                    # Also build the dual: 7-regular graph on 72 vertices = PSL(2,8)/C_7
                    Print("Building dual graph: PSL(2,8)/C_7 (72 vertices, 7-regular)...\n");
                    found := true;
                fi;

                # Output adjacency list
                Print("\n=== ADJ MACBEATH ===\n");
                Print("n := ", result.n, ";\n");
                Print("adj := ", result.adj, ";\n");
                Print("=== END ADJ ===\n");

                break;
            else
                Print("  H = ", StructureDescription(H),
                      " (order ", h_order, "): cubic but DISCONNECTED\n");
                break;
            fi;
        fi;
    od;

    if found then break; fi;
od;

if not found then
    Print("No Macbeath graph found — trying all subgroup classes...\n");
fi;


# ============================================================
# Also build the dual: 7-regular graph on 72 vertices
# ============================================================

Print("\n════════════════════════════════════════════════════════\n");
Print("  Dual graph: PSL(2,8)/C_7 (72 vertices, 7-regular)\n");
Print("════════════════════════════════════════════════════════\n\n");

for cc in cc_subs do
    H := Representative(cc);
    if Size(H) = 7 then
        Print("Found C_7 subgroup: ", StructureDescription(H), "\n");
        n_verts := Size(G) / Size(H);
        Print("  Vertices: ", n_verts, "\n");

        # Build coset graph with all involutions to find 7-regular one
        for a in involutions do
            if a in H then continue; fi;
            dcoset := DoubleCoset(H, a, H);
            val := Size(dcoset) / Size(H);
            Print("  Involution gives valency: ", val, "\n");

            if val = 7 then
                result := SabidussiGraph(G, H, a);
                if IsConnectedAdj(result.adj, result.n) then
                    girth := ComputeGirth(result.adj, result.n);
                    Print("★ Dual graph: |V| = ", result.n,
                          ", valency = 7, girth = ", girth, "\n");
                    Print("\n=== ADJ MACBEATH_DUAL ===\n");
                    Print("n := ", result.n, ";\n");
                    Print("adj := ", result.adj, ";\n");
                    Print("=== END ADJ ===\n");
                fi;
            fi;
            break;
        od;
        break;
    fi;
od;


# ============================================================
# Summary comparison
# ============================================================

Print("\n\n");
Print("════════════════════════════════════════════════════════\n");
Print("  Comparison: Klein quartic vs Macbeath surface\n");
Print("════════════════════════════════════════════════════════\n");
Print("\n");
Print("                    Klein quartic     Macbeath surface\n");
Print("  Genus              3                 7\n");
Print("  Group              PSL(2,7)          PSL(2,8)\n");
Print("  |Aut|              168               504\n");
Print("  Cubic graph        56 vertices       168 vertices\n");
Print("  r-regular graph    24 verts, 7-reg   72 verts, 7-reg\n");
Print("  Triangle group     (2,3,7)           (2,3,7)\n");
Print("  Both Hurwitz       yes               yes\n");
Print("\n");
