# psl2_cubic_census.g
# Census of cubic arc-transitive graphs from PSL(2,q) for all primes up to 200.
# For each q, try ALL conjugacy classes of subgroups H and involutions a.
# Report: q, |V|, girth, H, |Aut|.

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
    for v in [1..Minimum(n, 30)] do
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
# Main census
# ============================================================

Print("q,n_vertices,girth,H_structure,H_order,aut_order,aut_structure\n");

for q in Filtered([5..200], IsPrimeInt) do
    G := PSL(2, q);
    g_order := Size(G);

    # Get involutions
    involutions := [];
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 2 then
            involutions := AsList(c);
            break;
        fi;
    od;

    if Length(involutions) = 0 then continue; fi;

    # Track which vertex counts we've already found (avoid duplicates)
    found_sizes := [];

    cc_subs := ConjugacyClassesSubgroups(G);

    for cc in cc_subs do
        H := Representative(cc);
        h_order := Size(H);

        # Need 3 | |H| for cubic
        if h_order mod 3 <> 0 then continue; fi;
        if h_order <= 2 or h_order = g_order then continue; fi;
        # Skip if too many vertices
        if g_order / h_order > 1000 then continue; fi;

        for a in involutions do
            if a in H then continue; fi;

            dcoset := DoubleCoset(H, a, H);
            val := Size(dcoset) / h_order;

            if val = 3 then
                result := SabidussiGraph(G, H, a);
                if IsConnectedAdj(result.adj, result.n) then
                    # Avoid reporting same graph twice
                    if not result.n in found_sizes then
                        Add(found_sizes, result.n);
                        girth := ComputeGirth(result.adj, result.n);

                        # Only compute Aut for manageable sizes
                        if result.n <= 500 then
                            graph := Graph(Group(()), [1..result.n],
                                function(x,g) return x; end,
                                function(x,y) return y in result.adj[x]; end, true);
                            aut := AutGroupGraph(graph);
                            Print(q, ",", result.n, ",", girth, ",",
                                  StructureDescription(H), ",", h_order, ",",
                                  Size(aut), ",", StructureDescription(aut), "\n");
                        else
                            Print(q, ",", result.n, ",", girth, ",",
                                  StructureDescription(H), ",", h_order, ",,\n");
                        fi;
                    fi;
                fi;
                break;  # one involution per H is enough
            fi;
        od;
    od;
od;
