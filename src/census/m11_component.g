# m11_component.g
# Investigate the disconnected Sabidussi graph from M11/C3
# What is the connected component?

LoadPackage("grape");

G := MathieuGroup(11);
H := Group([(1,2,3)]);  # some C3 — we need the RIGHT one

# Find a C3 subgroup and the involution class
cc := ConjugacyClasses(G);
inv_class := First(cc, c -> Order(Representative(c)) = 2);
a := Representative(inv_class);

# Find all conjugacy classes of subgroups of order 3
subs := ConjugacyClassesSubgroups(G);
c3_classes := Filtered(subs, sc -> Size(Representative(sc)) = 3);
Print("Number of conjugacy classes of C3: ", Length(c3_classes), "\n");

for sc in c3_classes do
    H := Representative(sc);
    dc := DoubleCoset(H, a, H);
    val := Size(dc) / Size(H);
    Print("C3 class, valency = ", val, "\n");

    if val = 3 then
        # Build the graph
        cosets := RightCosets(G, H);
        n := Length(cosets);
        adj_list := List([1..n], i -> []);
        for i in [1..n] do
            x := Representative(cosets[i]);
            for j in [i+1..n] do
                y := Representative(cosets[j]);
                if x * y^(-1) in dc then
                    Add(adj_list[i], j);
                    Add(adj_list[j], i);
                fi;
            od;
        od;
        for i in [1..n] do adj_list[i] := Set(adj_list[i]); od;

        # Find connected component containing vertex 1
        visited := BlistList([1..n], [1]);
        queue := [1];
        component := [1];
        while Length(queue) > 0 do
            v := Remove(queue, 1);
            for w in adj_list[v] do
                if not visited[w] then
                    visited[w] := true;
                    Add(queue, w);
                    Add(component, w);
                fi;
            od;
        od;

        Print("Component size: ", Length(component), "\n");
        Print("Number of components: ", n / Length(component), "\n\n");

        # Extract the component as a graph and compute its automorphism group
        comp := Set(component);
        m := Length(comp);
        # Renumber vertices 1..m
        comp_adj := List([1..m], i -> []);
        for i in [1..m] do
            for w in adj_list[comp[i]] do
                j := Position(comp, w);
                if j <> fail then
                    AddSet(comp_adj[i], j);
                fi;
            od;
        od;

        # Compute girth of component
        girth := infinity;
        for v in [1..m] do
            dist := ListWithIdenticalEntries(m, -1);
            dist[v] := 0;
            q := [v];
            found_cycle := false;
            while Length(q) > 0 and not found_cycle do
                u := Remove(q, 1);
                for w in comp_adj[u] do
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
        Print("Component girth: ", girth, "\n");

        # Use GRAPE for automorphism group of component
        gamma := Graph(Group(()), [1..m], OnPoints,
            function(x, y) return y in comp_adj[x]; end, true);
        aut := AutGroupGraph(gamma);
        Print("Component |Aut| = ", Size(aut), "\n");
        Print("Component Aut = ", StructureDescription(aut), "\n");

        # Check: is this the dodecahedron? (20 vertices, cubic, girth 5)
        Print("\nFor reference: dodecahedron has 20 vertices, cubic, girth 5, |Aut| = 120 = A5 x C2\n");
    fi;
od;
