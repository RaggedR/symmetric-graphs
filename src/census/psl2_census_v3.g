# psl2_census_v3.g
# Census of cubic graphs from PSL(2,q).
# Strategy: for each q, construct C3, S3, and (when applicable) S4 directly.
# Use ALL involutions for the double-coset valency check (cheap).
# Only build the Sabidussi graph if we find valency 3.

LoadPackage("grape");

SabidussiGraph := function(G, H, a)
    local cosets, n, adj, i, j, dcoset, x, y;
    cosets := RightCosets(G, H);
    n := Length(cosets);
    if n > 2000 then return fail; fi;
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

# Try to build a cubic graph. Uses ALL involutions for the valency check.
TryCubic := function(G, H, q, label)
    local a, dcoset, val, result, c, inv_class;
    if Size(G) / Size(H) > 2000 then return fail; fi;

    # Get ALL involutions
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 2 then
            inv_class := c;
            break;
        fi;
    od;

    for a in inv_class do
        if a in H then continue; fi;
        dcoset := DoubleCoset(H, a, H);
        val := Size(dcoset) / Size(H);
        if val = 3 then
            result := SabidussiGraph(G, H, a);
            if result = fail then return fail; fi;
            if IsConnectedAdj(result.adj, result.n) then
                Print(q, ",", result.n, ",", label, ",", Size(H), "\n");
                return result.n;
            else
                # cubic but disconnected — no point trying more involutions
                return fail;
            fi;
        fi;
    od;
    return fail;
end;

Print("q,n_vertices,H_type,H_order\n");

for q in Filtered([5..100], IsPrimeInt) do
    G := PSL(2, q);

    # --- C_3 ---
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 3 then
            TryCubic(G, Group(Representative(c)), q, "C3");

            # --- S_3 = Normalizer(G, C_3) ---
            TryCubic(G, Normalizer(G, Group(Representative(c))), q, "S3");
            break;
        fi;
    od;

    # --- S_4: exists iff q ≡ ±1 mod 8 ---
    if q mod 8 in [1, 7] then
        # S_4 = Normalizer(G, V4) where V4 is a Klein 4-group.
        # Find V4: two commuting involutions generating a group of order 4.
        for c in ConjugacyClasses(G) do
            if Order(Representative(c)) = 2 then
                for a2 in c do
                    if a2 <> Representative(c) and
                       a2 * Representative(c) = Representative(c) * a2 then
                        H := Normalizer(G, Group(Representative(c), a2));
                        if Size(H) = 24 then
                            TryCubic(G, H, q, "S4");
                        fi;
                        break;
                    fi;
                od;
                break;
            fi;
        od;
    fi;

    # --- A_5: exists iff q ≡ ±1 mod 5 or q = 5 ---
    if q = 5 or q mod 5 in [1, 4] then
        # A_5 contains elements of orders 2, 3, 5.
        # Find it as a subgroup: take elements of order 3 and 5 that generate
        # a group of order 60.
        for c3 in ConjugacyClasses(G) do
            if Order(Representative(c3)) = 3 then
                for c5 in ConjugacyClasses(G) do
                    if Order(Representative(c5)) = 5 then
                        H := Group(Representative(c3), Representative(c5));
                        if Size(H) = 60 then
                            TryCubic(G, H, q, "A5");
                        fi;
                        break;
                    fi;
                od;
                break;
            fi;
        od;
    fi;

    # --- D_12 (dihedral of order 12) ---
    # Normalizer of C_6 when 6 | (q±1)/2
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 6 then
            H := Normalizer(G, Group(Representative(c)));
            if Size(H) = 12 then
                TryCubic(G, H, q, "D12");
            fi;
            break;
        fi;
    od;

od;
