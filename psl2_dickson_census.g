# psl2_dickson_census.g
#
# Targeted search for cubic graphs from PSL(2,q) using Dickson's theorem.
# Instead of enumerating ALL subgroups (expensive), only try the subgroups
# that Dickson's theorem says can exist:
#   - C_n for n | (q±1)/2
#   - D_n for n | (q±1)/2
#   - A_4 (when p > 2 and p ≢ ±3 mod 8, or p=3)
#   - S_4 (when q ≡ ±1 mod 8)
#   - A_5 (when q ≡ ±1 mod 5 or q = 5)
#   - Borel subgroup (upper triangular)
#   - PSL(2, q_0) / PGL(2, q_0) for subfields
#
# For cubic: need |H| divisible by 3.

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

# Try a specific subgroup H: find involution giving valency 3, check connectivity
TryCubicGraph := function(G, H, involutions, q, h_label)
    local a, dcoset, val, result, n_vertices;
    n_vertices := Size(G) / Size(H);
    if n_vertices > 2000 then return fail; fi;

    for a in involutions do
        if a in H then continue; fi;
        dcoset := DoubleCoset(H, a, H);
        val := Size(dcoset) / Size(H);
        if val = 3 then
            result := SabidussiGraph(G, H, a);
            if IsConnectedAdj(result.adj, result.n) then
                Print(q, ",", result.n, ",", h_label, ",", Size(H), "\n");
                return result.n;
            fi;
            return fail;  # cubic but disconnected
        fi;
    od;
    return fail;  # no cubic involution
end;

# ============================================================
# Main census using Dickson's theorem
# ============================================================

Print("q,n_vertices,H_type,H_order\n");

for q in Filtered([5..500], IsPrimeInt) do
    G := PSL(2, q);

    # Get involutions (all conjugate for q odd prime)
    involutions := [];
    for c in ConjugacyClasses(G) do
        if Order(Representative(c)) = 2 then
            # Take a sample — don't need all of them
            involutions := AsList(c);
            if Length(involutions) > 50 then
                involutions := involutions{[1..50]};
            fi;
            break;
        fi;
    od;

    if Length(involutions) = 0 then continue; fi;

    found_sizes := [];

    # --- Try C_3 ---
    # C_3 exists in PSL(2,q) when 3 | (q-1)/2 or 3 | (q+1)/2,
    # i.e., when q ≢ 0 mod 3 (which is always true for q prime > 3)
    if q > 3 then
        for c_class in ConjugacyClasses(G) do
            if Order(Representative(c_class)) = 3 then
                H := Group(Representative(c_class));
                n := TryCubicGraph(G, H, involutions, q, "C3");
                if n <> fail then Add(found_sizes, n); fi;
                break;
            fi;
        od;
    fi;

    # --- Try S_3 = D_3 ---
    # S_3 exists when C_3 exists (normalizer of C_3 in PSL(2,q))
    if q > 3 then
        for c_class in ConjugacyClasses(G) do
            if Order(Representative(c_class)) = 3 then
                H := Normalizer(G, Group(Representative(c_class)));
                if Size(H) = 6 then
                    n := TryCubicGraph(G, H, involutions, q, "S3");
                    if n <> fail and not n in found_sizes then
                        Add(found_sizes, n);
                    fi;
                fi;
                break;
            fi;
        od;
    fi;

    # --- Try A_4 ---
    # A_4 <= PSL(2,q) when q ≡ ±1 mod 8 or q = 3 or q ≡ ±3 mod 8 with some conditions
    # Actually A_4 is always in PSL(2,q) for q odd prime >= 5 (order 12, divisible by 3)
    # But we just try to find it
    if q >= 5 then
        for cc in ConjugacyClassesSubgroups(G) do
            H := Representative(cc);
            if Size(H) = 12 and StructureDescription(H) = "A4" then
                n := TryCubicGraph(G, H, involutions, q, "A4");
                if n <> fail and not n in found_sizes then
                    Add(found_sizes, n);
                fi;
                break;
            fi;
        od;
    fi;

    # --- Try S_4 ---
    # S_4 <= PSL(2,q) iff q ≡ ±1 mod 8
    if q mod 8 = 1 or q mod 8 = 7 then
        for cc in ConjugacyClassesSubgroups(G) do
            H := Representative(cc);
            if Size(H) = 24 and StructureDescription(H) = "S4" then
                n := TryCubicGraph(G, H, involutions, q, "S4");
                if n <> fail and not n in found_sizes then
                    Add(found_sizes, n);
                fi;
                break;
            fi;
        od;
    fi;

    # --- Try A_5 ---
    # A_5 <= PSL(2,q) iff q ≡ ±1 mod 5 or q = 5
    if q = 5 or q mod 5 = 1 or q mod 5 = 4 then
        for cc in ConjugacyClassesSubgroups(G) do
            H := Representative(cc);
            if Size(H) = 60 and StructureDescription(H) = "A5" then
                n := TryCubicGraph(G, H, involutions, q, "A5");
                if n <> fail and not n in found_sizes then
                    Add(found_sizes, n);
                fi;
                break;
            fi;
        od;
    fi;

    # --- Try D_6 (dihedral of order 12, contains C_3) ---
    if q > 5 then
        for cc in ConjugacyClassesSubgroups(G) do
            H := Representative(cc);
            if Size(H) = 12 and (StructureDescription(H) = "D12" or
                                  StructureDescription(H) = "C6 : C2") then
                n := TryCubicGraph(G, H, involutions, q, "D12");
                if n <> fail and not n in found_sizes then
                    Add(found_sizes, n);
                fi;
                break;
            fi;
        od;
    fi;

od;
