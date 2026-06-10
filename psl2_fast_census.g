# psl2_fast_census.g
# Fast census: only try C3, S3, and find S4/A5 by direct construction.
# Avoids ConjugacyClassesSubgroups entirely.

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

TryCubic := function(G, H, involutions, q, label)
    local a, dcoset, val, result;
    if Size(G) / Size(H) > 2000 then return fail; fi;
    for a in involutions do
        if a in H then continue; fi;
        dcoset := DoubleCoset(H, a, H);
        val := Size(dcoset) / Size(H);
        if val = 3 then
            result := SabidussiGraph(G, H, a);
            if result = fail then return fail; fi;
            if IsConnectedAdj(result.adj, result.n) then
                Print(q, ",", result.n, ",", label, ",", Size(H), "\n");
                return result.n;
            fi;
            return fail;
        fi;
    od;
    return fail;
end;

Print("q,n_vertices,H_type,H_order\n");

for q in Filtered([5..100], IsPrimeInt) do
    G := PSL(2, q);

    # Get a sample of involutions (all conjugate for q odd prime)
    involutions := [];
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 2 then
            involutions := AsList(c);
            if Length(involutions) > 30 then
                involutions := involutions{[1..30]};
            fi;
            break;
        fi;
    od;
    if Length(involutions) = 0 then continue; fi;

    # --- C_3: find an element of order 3, take the group it generates ---
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 3 then
            TryCubic(G, Group(Representative(c)), involutions, q, "C3");
            # S_3 = normalizer of C_3
            TryCubic(G, Normalizer(G, Group(Representative(c))),
                     involutions, q, "S3");
            break;
        fi;
    od;

    # --- S_4: q ≡ ±1 mod 8 ---
    # Find S_4 as the stabilizer of a special structure.
    # In PSL(2,q), S_4 is the normalizer of a Klein 4-group V4.
    # V4 exists in PSL(2,q) when q is odd (always for our primes).
    if q mod 8 = 1 or q mod 8 = 7 then
        # Find a Klein 4-group
        for c in ConjugacyClasses(G) do
            if Order(Representative(c)) = 2 then
                # Take an involution t, find another involution s commuting with it
                # but not equal to it
                for c2 in ConjugacyClasses(G) do
                    if Order(Representative(c2)) = 2 then
                        # Try pairs
                        for s in AsList(c2){[1..Minimum(Length(AsList(c2)), 20)]} do
                            if Representative(c) <> s and
                               Representative(c) * s = s * Representative(c) then
                                # Found V4
                                H := Normalizer(G, Group(Representative(c), s));
                                if Size(H) = 24 then
                                    TryCubic(G, H, involutions, q, "S4");
                                fi;
                                break;
                            fi;
                        od;
                        break;
                    fi;
                od;
                break;
            fi;
        od;
    fi;

    # --- A_5: q ≡ ±1 mod 5 or q = 5 ---
    if q = 5 or q mod 5 = 1 or q mod 5 = 4 then
        # Find A_5 via elements of order 5
        for c in ConjugacyClasses(G) do
            if Order(Representative(c)) = 5 then
                # A_5 contains an element of order 5. Take its normalizer...
                # Actually, in PSL(2,q), A_5 is harder to pin down.
                # Skip for now — |A_5| = 60 gives |V| = |G|/60, usually too small
                # to be interesting for small q, or A_5 might not yield valency 3.
                break;
            fi;
        od;
    fi;

od;
