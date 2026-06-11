# sp_cubic_census.g — v2
# Census of arc-transitive graphs from PSp(2n, q).
#
# Key insight: must use PROJECTIVE symplectic group (quotient by center),
# analogous to using PSL not SL. PSp(4,q) is simple for q >= 3.
# PSp(4,2)' = A_6 is the simple group in the q=2 case.
#
# Wan's "Geometry of Classical Groups" gives orbit structure.
# Searching valencies 3–7 since the symplectic geometry is rich.
#
# Usage: gap -q sp_cubic_census.g

LoadPackage("grape");

# ─────────────────────────────────────────────
# Utility functions
# ─────────────────────────────────────────────

SabidussiGraph := function(G, H, a)
    local cosets, n, adj, i, j, dcoset, x, y;
    cosets := RightCosets(G, H);
    n := Length(cosets);
    if n > 3000 then return fail; fi;
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
    local girth, start, dist, bfs, u, w, d;
    girth := infinity;
    for start in [1..Minimum(n, 50)] do
        dist := List([1..n], i -> -1);
        dist[start] := 0;
        bfs := [start];
        while Length(bfs) > 0 do
            u := Remove(bfs, 1);
            for w in adj[u] do
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
    return girth;
end;

# ─────────────────────────────────────────────
# Search a group G for Sabidussi graphs
# with valency in target_vals
# ─────────────────────────────────────────────

SearchGraphs := function(G, group_name, target_vals, max_vertices)
    local cc, subs, sc, H, ordH, n_vertices, ic, inv_classes,
          a, dc, valency, result, girth, found, gamma, aut_graph,
          aut_size, s_arc;

    Print("═══════════════════════════════════════════════════\n");
    Print("Group: ", group_name, "\n");
    Print("|G| = ", Size(G), "\n");
    Print("Structure: ", StructureDescription(G), "\n");
    Print("Simple: ", IsSimpleGroup(G), "\n");
    Print("Target valencies: ", target_vals, "\n\n");

    cc := ConjugacyClasses(G);
    inv_classes := Filtered(cc, c -> Order(Representative(c)) = 2);
    Print("Involution classes: ", Length(inv_classes));
    for ic in inv_classes do
        Print(" (size ", Size(ic), ")");
    od;
    Print("\n");

    Print("Computing conjugacy classes of subgroups...\n");
    subs := ConjugacyClassesSubgroups(G);
    Print("Subgroup classes: ", Length(subs), "\n\n");

    found := 0;

    for sc in subs do
        H := Representative(sc);
        ordH := Size(H);
        if ordH < 2 then continue; fi;

        n_vertices := Size(G) / ordH;
        if n_vertices < 6 or n_vertices > max_vertices then continue; fi;

        for ic in inv_classes do
            a := Representative(ic);
            if a in H then continue; fi;

            dc := DoubleCoset(H, a, H);
            valency := Size(dc) / ordH;

            if not valency in target_vals then continue; fi;

            result := SabidussiGraph(G, H, a);
            if result = fail then continue; fi;

            if not IsConnectedAdj(result.adj, result.n) then
                # only report disconnected cubic for debugging
                if valency = 3 then
                    Print("  [disc] H=", StructureDescription(H),
                          " |V|=", result.n, " val=", valency, "\n");
                fi;
                continue;
            fi;

            found := found + 1;
            girth := ComputeGirth(result.adj, result.n);

            # Full automorphism group
            gamma := Graph(Group(()), [1..result.n],
                function(x,g) return x; end,
                function(x,y) return y in result.adj[x]; end, true);
            aut_graph := AutGroupGraph(gamma);
            aut_size := Size(aut_graph);

            Print("────────────────────────────────────\n");
            Print("GRAPH #", found, "\n");
            Print("  Group: ", group_name, "\n");
            Print("  H = ", StructureDescription(H), ", |H| = ", ordH, "\n");
            Print("  Vertices: ", result.n, "\n");
            Print("  Valency: ", valency, "\n");
            Print("  Girth: ", girth, "\n");
            Print("  |Aut| = ", aut_size, "\n");
            Print("  Aut = ", StructureDescription(aut_graph), "\n");

            # Arc-transitivity check
            # For k-regular: vertex-trans needs |Aut| divisible by n
            # arc-trans needs |Aut| divisible by k*n
            if aut_size mod (valency * result.n) = 0 then
                s_arc := aut_size / (valency * result.n);
                Print("  Arc-transitive: YES\n");
                if s_arc > 1 then
                    Print("  |Aut|/(k*n) = ", s_arc, "\n");
                fi;
            elif aut_size mod result.n = 0 then
                Print("  Vertex-transitive: YES (not arc-transitive)\n");
            else
                Print("  Not vertex-transitive (interesting!)\n");
            fi;

            # Is it bipartite? (check: girth even and girth = diameter * 2?)
            # Simple check: odd girth => not bipartite
            if girth mod 2 = 0 then
                Print("  Bipartite: possibly\n");
            fi;

            Print("\n");
        od;
    od;

    Print("Total connected graphs from ", group_name, ": ", found, "\n\n");
    return found;
end;


# ─────────────────────────────────────────────
# Main census
# ─────────────────────────────────────────────

Print("╔═══════════════════════════════════════════════════╗\n");
Print("║  Graphs from Projective Symplectic Groups         ║\n");
Print("║  Wan's Anzahl theorems → Hecke algebra H(G,H)     ║\n");
Print("╚═══════════════════════════════════════════════════╝\n\n");

total := 0;
vals := [3, 4, 5, 6, 7];

# ── PSp(4,2)' = A_6, order 360 (simple) ──
# Exceptional: PSp(4,2) ≅ S_6, derived subgroup is A_6 ≅ PSL(2,9)
Print(">>> PSp(4,2)' = A_6 ≅ PSL(2,9) <<<\n\n");
G := DerivedSubgroup(Sp(4, 2));
total := total + SearchGraphs(G, "PSp(4,2)' = A_6", vals, 3000);

# ── PSp(4,2) = S_6, order 720 ──
Print(">>> PSp(4,2) = S_6 <<<\n\n");
G := Sp(4, 2);
total := total + SearchGraphs(G, "Sp(4,2) = S_6", vals, 3000);

# ── PSp(4,3), order 25920 ──
# Sp(4,3) has center {±I} of order 2. PSp(4,3) = Sp(4,3)/{±I}.
# PSp(4,3) is simple, isomorphic to O(5,3) (Omega, not full orthogonal).
Print(">>> PSp(4,3) = Omega(5,3), order 25920 <<<\n\n");
G := Sp(4, 3);
ctr := Center(G);
if Size(ctr) > 1 then
    G := G / ctr;
fi;
total := total + SearchGraphs(G, "PSp(4,3)", vals, 3000);

# ── PSp(4,4), order 979200 ──
# Too large for full subgroup search. Targeted search only.
Print(">>> PSp(4,4) — targeted search <<<\n\n");
G := Sp(4, 4);
# Center is trivial for q even, so Sp(4,4) = PSp(4,4)
Print("═══════════════════════════════════════════════════\n");
Print("Group: Sp(4,4) = PSp(4,4)\n");
Print("|G| = ", Size(G), "\n");
Print("Structure: ", StructureDescription(G), "\n\n");

cc := ConjugacyClasses(G);
inv_classes := Filtered(cc, c -> Order(Representative(c)) = 2);
Print("Involution classes: ", Length(inv_classes), "\n");

found_sp44 := 0;
# Try C_3 subgroups
for c in cc do
    if Order(Representative(c)) = 3 then
        H := Group(Representative(c));
        n_vert := Size(G) / Size(H);
        Print("  C_3: |V| = ", n_vert);
        if n_vert > 3000 then
            Print(" — skipping\n");
            continue;
        fi;
        for ic in inv_classes do
            a := Representative(ic);
            if a in H then continue; fi;
            dc := DoubleCoset(H, a, H);
            valency := Size(dc) / Size(H);
            if valency in vals then
                Print(" val=", valency);
                result := SabidussiGraph(G, H, a);
                if result <> fail and IsConnectedAdj(result.adj, result.n) then
                    found_sp44 := found_sp44 + 1;
                    girth := ComputeGirth(result.adj, result.n);
                    gamma := Graph(Group(()), [1..result.n],
                        function(x,g) return x; end,
                        function(x,y) return y in result.adj[x]; end, true);
                    aut_graph := AutGroupGraph(gamma);
                    Print("\n────────────────────────────────────\n");
                    Print("GRAPH from Sp(4,4)\n");
                    Print("  H = C_3, |V| = ", result.n, ", val = ", valency, "\n");
                    Print("  Girth: ", girth, "\n");
                    Print("  |Aut| = ", Size(aut_graph), "\n");
                    Print("  Aut = ", StructureDescription(aut_graph), "\n\n");
                fi;
            fi;
        od;
        Print("\n");
        break;  # just try one class of C_3
    fi;
od;

# Try S_3 and S_4 in Sp(4,4) via normalizers
for c in cc do
    if Order(Representative(c)) = 3 then
        H := Normalizer(G, Group(Representative(c)));
        n_vert := Size(G) / Size(H);
        Print("  N(C_3): |H| = ", Size(H), ", H = ", StructureDescription(H),
              ", |V| = ", n_vert);
        if n_vert > 3000 then
            Print(" — skipping\n");
        else
            Print("\n");
            for ic in inv_classes do
                a := Representative(ic);
                if a in H then continue; fi;
                dc := DoubleCoset(H, a, H);
                valency := Size(dc) / Size(H);
                if valency in vals then
                    result := SabidussiGraph(G, H, a);
                    if result <> fail and IsConnectedAdj(result.adj, result.n) then
                        found_sp44 := found_sp44 + 1;
                        girth := ComputeGirth(result.adj, result.n);
                        gamma := Graph(Group(()), [1..result.n],
                            function(x,g) return x; end,
                            function(x,y) return y in result.adj[x]; end, true);
                        aut_graph := AutGroupGraph(gamma);
                        Print("────────────────────────────────────\n");
                        Print("GRAPH from Sp(4,4)\n");
                        Print("  H = ", StructureDescription(H),
                              ", |V| = ", result.n, ", val = ", valency, "\n");
                        Print("  Girth: ", girth, "\n");
                        Print("  |Aut| = ", Size(aut_graph), "\n");
                        Print("  Aut = ", StructureDescription(aut_graph), "\n\n");
                    fi;
                fi;
            od;
        fi;
        break;
    fi;
od;

total := total + found_sp44;

Print("════════════════════════════════════════════════════\n");
Print("TOTAL connected graphs found: ", total, "\n");
