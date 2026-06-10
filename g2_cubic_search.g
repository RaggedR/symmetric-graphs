# g2_cubic_search.g
# Search G_2(2) = Aut(GH(2,2)) for cubic and 6-regular Sabidussi graphs
# Sab(G, H, HaH) where |HaH|/|H| = 3 or 6
#
# G_2(2) has order 12096 = 2^6 * 3^3 * 7
# PSU(3,3) = G_2(2)' has order 6048
#
# Usage: gap g2_cubic_search.g

LoadPackage("grape");

# Build G_2(2) as Aut of the GH(2,2) incidence graph (Tutte 12-cage)

Read("symmetric-layout.g");
skip_main := true;

data := GraphLibrary.tuttecage12();
gamma := Graph(Group(()), [1..data.n], function(x,g) return x; end,
    function(x,y) return y in data.adj[x]; end, true);
G := AutGroupGraph(gamma);

Print("G = Aut(Tutte-12-cage)\n");
Print("|G| = ", Size(G), "\n");
Print("Structure: ", StructureDescription(G), "\n\n");

# Derived subgroup = PSU(3,3)
Gprime := DerivedSubgroup(G);
Print("|G'| = ", Size(Gprime), "  ", StructureDescription(Gprime), "\n\n");

# Find involution classes
cc := ConjugacyClasses(G);
inv_classes := Filtered(cc, c -> Order(Representative(c)) = 2);
Print("Involution classes in G: ", Length(inv_classes), "\n");
for c in inv_classes do
    Print("  size ", Size(c), ", rep in G': ", Representative(c) in Gprime, "\n");
od;
Print("\n");

# Get conjugacy classes of subgroups
# G_2(2) has order 12096 — small enough for full enumeration
Print("Computing conjugacy classes of subgroups...\n");
subs := ConjugacyClassesSubgroups(G);
Print("Number of conjugacy classes: ", Length(subs), "\n\n");

# Search for valency 3 and 6
Print("Searching for Sabidussi triples with valency 3 or 6...\n");
Print("================================================================\n\n");

found := 0;
for sc in subs do
    H := Representative(sc);
    ordH := Size(H);
    if ordH < 2 then continue; fi;

    n_vertices := Size(G) / ordH;
    if n_vertices < 10 then continue; fi;  # skip tiny graphs

    for ic in inv_classes do
        a := Representative(ic);
        if a in H then continue; fi;

        dc := DoubleCoset(H, a, H);
        valency := Size(dc) / ordH;

        if valency in [3, 6] then
            found := found + 1;
            Print("FOUND #", found, ": valency ", valency, "\n");
            Print("  H = ", StructureDescription(H), ", |H| = ", ordH, "\n");
            Print("  |V| = ", n_vertices, "\n");

            # Build graph
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

            # Check degrees
            degs := List(adj, Length);
            deg_set := Set(degs);
            Print("  Degree set: ", deg_set, "\n");

            # Check connectivity
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

            if count = n then
                Print("  Connected: YES\n");

                # Compute girth (BFS from first few vertices)
                girth := infinity;
                for start in [1..Minimum(n, 30)] do
                    dist := List([1..n], i -> -1);
                    dist[start] := 0;
                    bfs := [start];
                    while Length(bfs) > 0 do
                        u := Remove(bfs, 1);
                        for w in adj[u] do
                            if dist[w] = -1 then
                                dist[w] := dist[u] + 1;
                                Add(bfs, w);
                            elif dist[w] >= dist[u] and w <> start then
                                d := dist[u] + dist[w] + 1;
                                if d < girth then girth := d; fi;
                            fi;
                        od;
                        if dist[u] >= girth/2 then break; fi;
                    od;
                od;
                Print("  Girth: ", girth, "\n");

                # Output adjacency list for Python
                Print("  ADJ_START\n");
                for i in [1..n] do
                    Print("  ", adj[i], "\n");
                od;
                Print("  ADJ_END\n");
            else
                Print("  Connected: NO (largest component: ", count, "/", n, ")\n");
            fi;
            Print("\n");
        fi;
    od;
od;

Print("\nTotal found: ", found, "\n");
