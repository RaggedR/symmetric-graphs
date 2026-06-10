# clayworth_build.g
# Construct the Clayworth graphs: {12,3} and {3,12} tilings of the genus-505 surface
#
# {12,3}: Sab(G_2(2), C_3) — 4032 vertices, cubic, 6048 edges, 1008 dodecagonal faces
# {3,12}: Sab(G_2(2), C_12) — 1008 vertices, degree 12, 6048 edges, 4032 triangular faces
#
# Both live on a genus-505 surface (V - E + F = -1008 = 2 - 2*505)
# Together they define a CSS quantum code.
#
# Usage: gap clayworth_build.g

LoadPackage("grape");

# Build G_2(2) as Aut of the Tutte 12-cage
# We need to source the graph library
skip_main := true;
Read("symmetric-layout.g");

Print("Building G_2(2) from Tutte 12-cage automorphisms...\n");
data := GraphLibrary.tuttecage12();
gamma := Graph(Group(()), [1..data.n], function(x,g) return x; end,
    function(x,y) return y in data.adj[x]; end, true);
G := AutGroupGraph(gamma);

Print("|G| = ", Size(G), "\n");
Print("G = ", StructureDescription(G), "\n\n");

# ============================================================
# Helper: build Sabidussi graph Sab(G, H, a)
# ============================================================

SabidussiGraph := function(G, H, a)
    local cosets, n, adj, i, j, dc, x, y;
    cosets := RightCosets(G, H);
    n := Length(cosets);
    Print("  Building coset graph: ", n, " vertices\n");
    adj := List([1..n], i -> []);
    dc := DoubleCoset(H, a, H);
    Print("  |HaH|/|H| = ", Size(dc)/Size(H), " (valency)\n");
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

# ============================================================
# Find the right C_3 and C_12 subgroups and involutions
# ============================================================

# Find all involutions (up to conjugacy)
cc := ConjugacyClasses(G);
inv_reps := List(Filtered(cc, c -> Order(Representative(c)) = 2), Representative);
Print("Involution classes: ", Length(inv_reps), "\n\n");

# Find C_3 subgroups
Print("=== Clayworth-3 graph: Sab(G_2(2), C_3, a) ===\n");
Print("Target: 4032 vertices, degree 3\n\n");

# We need a C_3 subgroup and an involution giving valency 3
el3 := Filtered(Elements(G), g -> Order(g) = 3);
Print("Elements of order 3: ", Length(el3), "\n");

# Try first element of order 3
found_c3 := false;
for g3 in el3 do
    H3 := Group(g3);
    if Size(H3) <> 3 then continue; fi;
    for a in inv_reps do
        if a in H3 then continue; fi;
        dc := DoubleCoset(H3, a, H3);
        val := Size(dc) / Size(H3);
        if val = 3 then
            Print("Found! H = <g3>, |HaH|/|H| = ", val, "\n");
            Print("Building graph...\n");
            cw3 := SabidussiGraph(G, H3, a);

            # Check connectivity
            visited := BlistList([1..cw3.n], [1]);
            queue := [1]; count := 1;
            while Length(queue) > 0 do
                v := Remove(queue, 1);
                for w in cw3.adj[v] do
                    if not visited[w] then
                        visited[w] := true;
                        Add(queue, w);
                        count := count + 1;
                    fi;
                od;
            od;

            if count = cw3.n then
                Print("Connected: YES\n");

                # Verify degrees
                degs := Set(List(cw3.adj, Length));
                Print("Degree set: ", degs, "\n");

                # Edge count
                nedges := Sum(List(cw3.adj, Length)) / 2;
                Print("Edges: ", nedges, "\n");

                # Euler characteristic check
                # For {12,3}: V=4032, E=6048, F=1008 → chi = V-E+F = -1008 = 2-2*505
                Print("V - E = ", cw3.n - nedges, "\n");
                Print("Expected V=4032, E=6048\n");

                # Compute girth
                girth := infinity;
                for start in [1..Minimum(cw3.n, 50)] do
                    dist := List([1..cw3.n], i -> -1);
                    dist[start] := 0;
                    bfs := [start];
                    while Length(bfs) > 0 do
                        u := Remove(bfs, 1);
                        for w in cw3.adj[u] do
                            if dist[w] = -1 then
                                dist[w] := dist[u] + 1;
                                Add(bfs, w);
                            elif dist[w] >= dist[u] then
                                d := dist[u] + dist[w] + 1;
                                if d < girth then girth := d; fi;
                            fi;
                        od;
                        if dist[u] >= girth/2 then break; fi;
                    od;
                od;
                Print("Girth: ", girth, "\n");

                found_c3 := true;
                break;
            else
                Print("Connected: NO (", count, "/", cw3.n, ")\n");
            fi;
        fi;
    od;
    if found_c3 then break; fi;
od;

if not found_c3 then
    Print("WARNING: No cubic Clayworth graph found with C_3!\n");
fi;

Print("\n");

# ============================================================
# Clayworth-12: Sab(G_2(2), C_12, a) — the dual tiling
# ============================================================

Print("=== Clayworth-12 graph: Sab(G_2(2), C_12, a) ===\n");
Print("Target: 1008 vertices, degree 12\n\n");

# Find C_12 subgroups
el12 := Filtered(Elements(G), g -> Order(g) = 12);
Print("Elements of order 12: ", Length(el12), "\n");

found_c12 := false;
for g12 in el12 do
    H12 := Group(g12);
    if Size(H12) <> 12 then continue; fi;
    for a in inv_reps do
        if a in H12 then continue; fi;
        dc := DoubleCoset(H12, a, H12);
        val := Size(dc) / Size(H12);
        if val = 12 then
            Print("Found! H = <g12>, |HaH|/|H| = ", val, "\n");
            Print("Building graph...\n");
            cw12 := SabidussiGraph(G, H12, a);

            # Check connectivity
            visited := BlistList([1..cw12.n], [1]);
            queue := [1]; count := 1;
            while Length(queue) > 0 do
                v := Remove(queue, 1);
                for w in cw12.adj[v] do
                    if not visited[w] then
                        visited[w] := true;
                        Add(queue, w);
                        count := count + 1;
                    fi;
                od;
            od;

            if count = cw12.n then
                Print("Connected: YES\n");
                degs := Set(List(cw12.adj, Length));
                Print("Degree set: ", degs, "\n");
                nedges := Sum(List(cw12.adj, Length)) / 2;
                Print("Edges: ", nedges, "\n");
                Print("V - E = ", cw12.n - nedges, "\n");
                Print("Expected V=1008, E=6048\n");

                # Girth
                girth := infinity;
                for start in [1..Minimum(cw12.n, 50)] do
                    dist := List([1..cw12.n], i -> -1);
                    dist[start] := 0;
                    bfs := [start];
                    while Length(bfs) > 0 do
                        u := Remove(bfs, 1);
                        for w in cw12.adj[u] do
                            if dist[w] = -1 then
                                dist[w] := dist[u] + 1;
                                Add(bfs, w);
                            elif dist[w] >= dist[u] then
                                d := dist[u] + dist[w] + 1;
                                if d < girth then girth := d; fi;
                            fi;
                        od;
                        if dist[u] >= girth/2 then break; fi;
                    od;
                od;
                Print("Girth: ", girth, "\n");

                found_c12 := true;
                break;
            else
                Print("Connected: NO (", count, "/", cw12.n, ")\n");
            fi;
        fi;
    od;
    if found_c12 then break; fi;
od;

if not found_c12 then
    Print("WARNING: No degree-12 Clayworth graph found with C_12!\n");
fi;

# ============================================================
# Output JSON for Python code analysis
# ============================================================

Print("\n\n=== JSON OUTPUT ===\n");
Print("{\n");

if found_c3 then
    Print("  \"clayworth3\": {\n");
    Print("    \"name\": \"Clayworth-3\",\n");
    Print("    \"description\": \"{12,3} tiling of genus-505 surface\",\n");
    Print("    \"n\": ", cw3.n, ",\n");
    Print("    \"adj\": [\n");
    for i in [1..cw3.n] do
        Print("      ", cw3.adj[i]);
        if i < cw3.n then Print(","); fi;
        Print("\n");
    od;
    Print("    ]\n");
    Print("  }");
    if found_c12 then Print(","); fi;
    Print("\n");
fi;

if found_c12 then
    Print("  \"clayworth12\": {\n");
    Print("    \"name\": \"Clayworth-12\",\n");
    Print("    \"description\": \"{3,12} tiling of genus-505 surface\",\n");
    Print("    \"n\": ", cw12.n, ",\n");
    Print("    \"adj\": [\n");
    for i in [1..cw12.n] do
        Print("      ", cw12.adj[i]);
        if i < cw12.n then Print(","); fi;
        Print("\n");
    od;
    Print("    ]\n");
    Print("  }\n");
fi;

Print("}\n");
Print("=== END JSON ===\n");
