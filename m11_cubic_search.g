# m11_cubic_search.g
# Search M11 for cubic Sabidussi graphs Sab(M11, H, HaH)
# A cubic graph requires |HaH|/|H| = 3 (valency 3)
#
# Usage: gap m11_cubic_search.g

LoadPackage("grape");

G := MathieuGroup(11);
Print("M11: order ", Size(G), "\n\n");

# Find all conjugacy classes of involutions
cc := ConjugacyClasses(G);
involution_classes := Filtered(cc, c -> Order(Representative(c)) = 2);
Print("Number of involution classes: ", Length(involution_classes), "\n");
for c in involution_classes do
    Print("  Class size: ", Size(c), ", representative: ", Representative(c), "\n");
od;
Print("\n");

# Find all conjugacy classes of subgroups
subs := ConjugacyClassesSubgroups(G);
Print("Number of conjugacy classes of subgroups: ", Length(subs), "\n\n");

# For each subgroup H and each involution a (not in H),
# check if |HaH|/|H| = 3
Print("Searching for cubic Sabidussi triples (H, a) with |HaH|/|H| = 3...\n");
Print("================================================================\n\n");

found := 0;
for sc in subs do
    H := Representative(sc);
    ordH := Size(H);
    if ordH < 2 then continue; fi;  # skip trivial

    n_vertices := Size(G) / ordH;

    for ic in involution_classes do
        a := Representative(ic);

        # Check if a is NOT in H (need coset graph to be non-trivial)
        if a in H then continue; fi;

        # Compute double coset HaH
        dc := DoubleCoset(H, a, H);
        valency := Size(dc) / ordH;

        if valency = 3 then
            found := found + 1;
            Print("FOUND cubic triple #", found, ":\n");
            Print("  H = ", StructureDescription(H), ", |H| = ", ordH, "\n");
            Print("  Involution class size: ", Size(ic), "\n");
            Print("  Vertices: ", n_vertices, "\n");
            Print("  Valency: ", valency, "\n");

            # Check if a^2 is in H (needed for undirected graph)
            if a^2 in H then
                Print("  a^2 in H: YES (graph is undirected)\n");
            else
                Print("  a^2 in H: NO (directed graph)\n");
            fi;

            # Build the actual graph to check connectivity
            cosets := RightCosets(G, H);
            n := Length(cosets);
            adj := List([1..n], i -> []);
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

            # Check connectivity via BFS
            visited := BlistList([1..n], [1]);
            queue := [1];
            while Length(queue) > 0 do
                v := Remove(queue, 1);
                for w in adj[v] do
                    if not visited[w] then
                        visited[w] := true;
                        Add(queue, w);
                    fi;
                od;
            od;
            n_visited := SizeBlist(visited);

            if n_visited = n then
                Print("  Connected: YES\n");

                # Compute girth
                girth := infinity;
                for v in [1..n] do
                    dist := ListWithIdenticalEntries(n, -1);
                    dist[v] := 0;
                    q := [v];
                    found_cycle := false;
                    while Length(q) > 0 and not found_cycle do
                        u := Remove(q, 1);
                        for w in adj[u] do
                            if dist[w] = -1 then
                                dist[w] := dist[u] + 1;
                                Add(q, w);
                            elif dist[w] >= dist[u] then
                                cycle_len := dist[u] + dist[w] + 1;
                                if cycle_len < girth then
                                    girth := cycle_len;
                                fi;
                                found_cycle := true;
                            fi;
                        od;
                    od;
                od;
                Print("  Girth: ", girth, "\n");

                # Compute automorphism group using GRAPE
                gamma := Graph(Group(()), [1..n], OnPoints,
                    function(x, y) return y in adj[x]; end, true);
                aut := AutGroupGraph(gamma);
                Print("  |Aut(G)| = ", Size(aut), "\n");
                Print("  Aut(G) = ", StructureDescription(aut), "\n");

            else
                Print("  Connected: NO (largest component: ", n_visited, "/", n, ")\n");
            fi;
            Print("\n");
        fi;
    od;
od;

if found = 0 then
    Print("No cubic Sabidussi graphs found for M11.\n");
fi;

Print("\nDone. Found ", found, " cubic triples total.\n");
