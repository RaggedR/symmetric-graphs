# klein_quartic_codes.g
# Construct the Klein quartic tilings: {7,3} and {3,7} on genus-3 surface
# Triangle group Δ(2,3,7), quotient PSL(2,7), order 168
#
# {7,3}: 56 vertices (degree 3), 84 edges, 24 heptagonal faces
# {3,7}: 24 vertices (degree 7), 84 edges, 56 triangular faces
#
# Usage: gap klein_quartic_codes.g

LoadPackage("grape");

G := PSL(2, 7);
Print("G = PSL(2,7)\n");
Print("|G| = ", Size(G), "\n\n");

# ============================================================
# Helper
# ============================================================

SabidussiGraph := function(G, H, a)
    local cosets, n, adj, i, j, dc, x, y;
    cosets := RightCosets(G, H);
    n := Length(cosets);
    adj := List([1..n], i -> []);
    dc := DoubleCoset(H, a, H);
    for i in [1..n] do
        x := Representative(cosets[i]);
        for j in [i+1..n] do
            y := Representative(cosets[j]);
            if x * y^(-1) in dc then
                Add(adj[i], j);
                Add(adj[j], i);
            fi;
        od;
    od;
    for i in [1..n] do adj[i] := Set(adj[i]); od;
    return rec(n := n, adj := adj);
end;

ComputeGirth := function(adj, n)
    local girth, v, dist, queue, u, w, d, parent;
    girth := infinity;
    for v in [1..n] do
        dist := List([1..n], i -> -1);
        parent := List([1..n], i -> -1);
        dist[v] := 0;
        queue := [v];
        while Length(queue) > 0 do
            u := Remove(queue, 1);
            for w in adj[u] do
                if dist[w] = -1 then
                    dist[w] := dist[u] + 1;
                    parent[w] := u;
                    Add(queue, w);
                elif parent[u] <> w and parent[w] <> u then
                    d := dist[u] + dist[w] + 1;
                    if d < girth then girth := d; fi;
                fi;
            od;
            if dist[u] >= girth/2 then break; fi;
        od;
    od;
    return girth;
end;

CheckConnected := function(adj, n)
    local visited, queue, v, w, count;
    visited := BlistList([1..n], [1]);
    queue := [1]; count := 1;
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

# ============================================================
# Find {7,3} tiling: Sab(PSL(2,7), C_3, a) — 56 vertices, cubic
# ============================================================

Print("=== {7,3} Klein quartic tiling ===\n");
Print("Target: 56 vertices, degree 3\n\n");

cc := ConjugacyClasses(G);
inv_reps := List(Filtered(cc, c -> Order(Representative(c)) = 2), Representative);
Print("Involution classes: ", Length(inv_reps), "\n");

el3 := Filtered(Elements(G), g -> Order(g) = 3);
found73 := false;
for g3 in el3 do
    H3 := Group(g3);
    if Size(H3) <> 3 then continue; fi;
    for a in inv_reps do
        if a in H3 then continue; fi;
        dc := DoubleCoset(H3, a, H3);
        val := Size(dc) / Size(H3);
        if val = 3 then
            graph73 := SabidussiGraph(G, H3, a);
            if CheckConnected(graph73.adj, graph73.n) then
                degs := Set(List(graph73.adj, Length));
                nedges := Sum(List(graph73.adj, Length)) / 2;
                girth := ComputeGirth(graph73.adj, graph73.n);
                Print("Found {7,3}: ", graph73.n, " vertices, degree ", degs, "\n");
                Print("  Edges: ", nedges, "\n");
                Print("  Girth: ", girth, "\n");
                found73 := true;
                break;
            fi;
        fi;
    od;
    if found73 then break; fi;
od;

# ============================================================
# Find {3,7} tiling: Sab(PSL(2,7), C_7, a) — 24 vertices, degree 7
# ============================================================

Print("\n=== {3,7} Klein quartic dual tiling ===\n");
Print("Target: 24 vertices, degree 7\n\n");

el7 := Filtered(Elements(G), g -> Order(g) = 7);
found37 := false;
for g7 in el7 do
    H7 := Group(g7);
    if Size(H7) <> 7 then continue; fi;
    for a in inv_reps do
        if a in H7 then continue; fi;
        dc := DoubleCoset(H7, a, H7);
        val := Size(dc) / Size(H7);
        if val = 7 then
            graph37 := SabidussiGraph(G, H7, a);
            if CheckConnected(graph37.adj, graph37.n) then
                degs := Set(List(graph37.adj, Length));
                nedges := Sum(List(graph37.adj, Length)) / 2;
                girth := ComputeGirth(graph37.adj, graph37.n);
                Print("Found {3,7}: ", graph37.n, " vertices, degree ", degs, "\n");
                Print("  Edges: ", nedges, "\n");
                Print("  Girth: ", girth, "\n");
                found37 := true;
                break;
            fi;
        fi;
    od;
    if found37 then break; fi;
od;

# ============================================================
# Output JSON
# ============================================================

Print("\n\n=== JSON OUTPUT ===\n");
Print("{\n");

if found73 then
    Print("  \"klein73\": {\n");
    Print("    \"name\": \"Klein-73\",\n");
    Print("    \"description\": \"{7,3} tiling of Klein quartic (genus 3)\",\n");
    Print("    \"n\": ", graph73.n, ",\n");
    Print("    \"adj\": [\n");
    for i in [1..graph73.n] do
        Print("      ", graph73.adj[i]);
        if i < graph73.n then Print(","); fi;
        Print("\n");
    od;
    Print("    ]\n");
    Print("  }");
    if found37 then Print(","); fi;
    Print("\n");
fi;

if found37 then
    Print("  \"klein37\": {\n");
    Print("    \"name\": \"Klein-37\",\n");
    Print("    \"description\": \"{3,7} tiling of Klein quartic (genus 3)\",\n");
    Print("    \"n\": ", graph37.n, ",\n");
    Print("    \"adj\": [\n");
    for i in [1..graph37.n] do
        Print("      ", graph37.adj[i]);
        if i < graph37.n then Print(","); fi;
        Print("\n");
    od;
    Print("    ]\n");
    Print("  }\n");
fi;

Print("}\n");
Print("=== END JSON ===\n");
