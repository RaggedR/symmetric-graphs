# bolza_graphs.g
# Construct all Sabidussi graphs of GL(2,3) — the automorphism group
# of the Bolza surface (genus 2, from Delta(2,3,8)).
#
# For each subgroup H of GL(2,3), build the coset graph GL(2,3)/H
# and determine: vertices, valency, girth, connectivity, automorphism group.
#
# The key graphs for the Bolza surface tessellation:
#   GL(2,3)/C_3 = 16 vertices (the Moebius-Kantor graph F016)
#   GL(2,3)/C_8 =  6 vertices (the {3,8} dual — expected multigraph)

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

# Multigraph version: count multiplicities
SabidussiMultigraph := function(G, H, involutions)
    local cosets, n, mult, i, j, x, y, dcoset, a, val, total_edges;
    cosets := RightCosets(G, H);
    n := Length(cosets);
    # mult[i][j] counts edges between vertex i and vertex j
    mult := List([1..n], i -> List([1..n], j -> 0));

    # Use ALL involutions to build the full connection
    for a in involutions do
        if a in H then continue; fi;
        dcoset := DoubleCoset(H, a, H);
        for i in [1..n] do
            x := Representative(cosets[i]);
            for j in [i+1..n] do
                y := Representative(cosets[j]);
                if x * y^(-1) in dcoset then
                    mult[i][j] := mult[i][j] + 1;
                    mult[j][i] := mult[j][i] + 1;
                fi;
            od;
        od;
    od;
    return rec(n := n, mult := mult);
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
    for v in [1..n] do
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
# Build GL(2,3) and find its structure
# ============================================================

Print("════════════════════════════════════════════════════════\n");
Print("  Sabidussi graphs of GL(2,3)\n");
Print("  Automorphism group of the Bolza surface (genus 2)\n");
Print("════════════════════════════════════════════════════════\n\n");

G := GL(2, 3);
Print("G = GL(2,3),  |G| = ", Size(G), "\n");
Print("Structure: ", StructureDescription(G), "\n\n");

# Conjugacy classes
Print("Conjugacy classes:\n");
for c in ConjugacyClasses(G) do
    Print("  Order ", Order(Representative(c)),
          ", size ", Size(c), "\n");
od;
Print("\n");

# Find involutions
involutions := [];
for c in ConjugacyClasses(G) do
    if Order(Representative(c)) = 2 then
        Append(involutions, AsList(c));
    fi;
od;
Print("Total involutions: ", Length(involutions), "\n\n");

# Subgroup classes
cc_subs := ConjugacyClassesSubgroups(G);
Print("Conjugacy classes of subgroups: ", Length(cc_subs), "\n\n");


# ============================================================
# Build all Sabidussi coset graphs
# ============================================================

Print("════════════════════════════════════════════════════════\n");
Print("  All coset graphs GL(2,3)/H\n");
Print("════════════════════════════════════════════════════════\n\n");

for cc in cc_subs do
    H := Representative(cc);
    h_order := Size(H);
    n_verts := Size(G) / h_order;

    # Skip trivial cases
    if h_order = Size(G) then continue; fi;
    if h_order = 1 then continue; fi;

    Print("────────────────────────────────────────\n");
    Print("H = ", StructureDescription(H),
          " (order ", h_order, ")\n");
    Print("  Vertices: ", n_verts,
          "  (conjugacy classes: ", Size(cc), ")\n");

    # Try each involution
    found := false;
    for a in involutions do
        if a in H then continue; fi;

        dcoset := DoubleCoset(H, a, H);
        val := Size(dcoset) / h_order;

        result := SabidussiGraph(G, H, a);

        if IsConnectedAdj(result.adj, result.n) then
            girth := ComputeGirth(result.adj, result.n);
            deg_seq := List([1..result.n], i -> Length(result.adj[i]));

            Print("  Valency: ", deg_seq[1],
                  ", girth: ", girth,
                  ", connected: yes\n");

            # Check if regular
            if ForAll(deg_seq, d -> d = deg_seq[1]) then
                Print("  Regular: yes (", deg_seq[1], "-regular)\n");
            else
                Print("  Regular: no, deg seq = ", Set(deg_seq), "\n");
            fi;

            # Compute automorphism group for small graphs
            if result.n <= 200 then
                graph := Graph(Group(()), [1..result.n],
                    function(x,g) return x; end,
                    function(x,y) return y in result.adj[x]; end, true);
                aut := AutGroupGraph(graph);
                Print("  |Aut| = ", Size(aut), "  ",
                      StructureDescription(aut), "\n");
            fi;

            Print("  adj = ", result.adj, "\n");

            found := true;
            break;
        fi;
    od;

    if not found then
        Print("  No connected coset graph found with any involution.\n");

        # Check if it works as multigraph (for the 6-vertex case)
        if n_verts <= 10 then
            Print("  Checking multigraph structure...\n");
            mg := SabidussiMultigraph(G, H, involutions);
            Print("  Multiplicity matrix:\n");
            for i in [1..mg.n] do
                Print("    ", mg.mult[i], "\n");
            od;
            # Underlying simple graph
            Print("  Underlying simple graph adjacency:\n");
            for i in [1..mg.n] do
                nbrs := [];
                for j in [1..mg.n] do
                    if i <> j and mg.mult[i][j] > 0 then
                        Add(nbrs, j);
                    fi;
                od;
                Print("    ", i, " -> ", nbrs,
                      " (degree ", Length(nbrs), ")\n");
            od;
        fi;
    fi;
    Print("\n");
od;
