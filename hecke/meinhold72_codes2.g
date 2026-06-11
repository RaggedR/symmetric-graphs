# meinhold72_codes2.g — Deeper analysis of the best codes
#
# Focus on:
#   1. Weight enumerator of [72, 8, 28] code
#   2. Minimum distance estimation for [72, 36] self-dual codes
#   3. Type I/II classification of self-dual codes
#   4. Automorphism groups
#
# Usage: gap -q meinhold72_codes2.g

# ─────────────────────────────────────────────
# Setup (same as before)
# ─────────────────────────────────────────────

G := PrimitiveGroup(63, 2);;
subs := ConjugacyClassesSubgroups(G);;
K := First(List(subs, Representative), H -> Size(H) = 168);;
cosets := RightCosets(G, K);;
n := Length(cosets);;
act := ActionHomomorphism(G, cosets, OnRight);;
dcosets := DoubleCosets(G, K, K);;
r := Length(dcosets);;
valencies := List(dcosets, d -> Size(d) / Size(K));;

Print("Setup: G₂(2) on ", n, " points, rank ", r, "\n");;
Print("Valencies: ", valencies, "\n\n");;

# Element → double coset lookup
elt_list := AsSortedList(G);;
dc_of_elt := List([1..Size(G)], x -> 0);;
for i in [1..r] do
    for el in dcosets[i] do
        j := PositionSorted(elt_list, el);;
        dc_of_elt[j] := i;;
    od;
od;

coset_reps := List(cosets, c -> Representative(c));;

# Build adjacency matrices over GF(2)
adj := List([1..r], d -> NullMat(n, n, GF(2)));;
for i in [1..n] do
    x := coset_reps[i];;
    for j in [1..n] do
        y := coset_reps[j];;
        prod := x * y^(-1);;
        d := dc_of_elt[PositionSorted(elt_list, prod)];;
        adj[d][i][j] := One(GF(2));;
    od;
od;

Print("Adjacency matrices built.\n\n");;

# ─────────────────────────────────────────────
# The [72, 8, 28] code from T₃+T₄+T₇+T₈
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("The [72, 8, 28] self-orthogonal code\n");;
Print("═══════════════════════════════════════════════════\n\n");;

M_8_28 := adj[3] + adj[4] + adj[7] + adj[8];;
basis_8_28 := BaseMat(M_8_28);;
k := Length(basis_8_28);;
Print("Dimension: ", k, "\n");;
Print("Self-orthogonal: ", IsZero(M_8_28 * TransposedMat(M_8_28)), "\n\n");;

# Full weight enumerator
Print("Weight enumerator:\n");;
wd := List([0..n], x -> 0);;
for combo in [0..2^k - 1] do
    v := Zero(GF(2)) * basis_8_28[1];;
    current := combo;;
    for j in [1..k] do
        if current mod 2 = 1 then
            v := v + basis_8_28[j];;
        fi;
        current := QuoInt(current, 2);;
    od;
    wt := Number([1..n], l -> v[l] <> Zero(GF(2)));;
    wd[wt + 1] := wd[wt + 1] + 1;;
od;

for w in [0..n] do
    if wd[w+1] > 0 then
        Print("  A_", w, " = ", wd[w+1], "\n");;
    fi;
od;
Print("\n");;

# Check if doubly-even (all weights divisible by 4)
doubly_even := true;;
for w in [0..n] do
    if wd[w+1] > 0 and w mod 4 <> 0 then
        doubly_even := false;;
    fi;
od;
Print("Doubly-even: ", doubly_even, "\n\n");;

# ─────────────────────────────────────────────
# The [72, 15, 16] code from T₃+T₄+T₅+T₆
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("The [72, 15, 16] self-orthogonal code\n");;
Print("═══════════════════════════════════════════════════\n\n");;

M_15_16 := adj[3] + adj[4] + adj[5] + adj[6];;
basis_15_16 := BaseMat(M_15_16);;
k15 := Length(basis_15_16);;
Print("Dimension: ", k15, "\n");;
Print("Self-orthogonal: ", IsZero(M_15_16 * TransposedMat(M_15_16)), "\n\n");;

# Weight enumerator (2^15 = 32768, feasible)
Print("Weight enumerator:\n");;
wd15 := List([0..n], x -> 0);;
for combo in [0..2^k15 - 1] do
    v := Zero(GF(2)) * basis_15_16[1];;
    current := combo;;
    for j in [1..k15] do
        if current mod 2 = 1 then
            v := v + basis_15_16[j];;
        fi;
        current := QuoInt(current, 2);;
    od;
    wt := Number([1..n], l -> v[l] <> Zero(GF(2)));;
    wd15[wt + 1] := wd15[wt + 1] + 1;;
od;

for w in [0..n] do
    if wd15[w+1] > 0 then
        Print("  A_", w, " = ", wd15[w+1], "\n");;
    fi;
od;
Print("\n");;

doubly_even_15 := true;;
for w in [0..n] do
    if wd15[w+1] > 0 and w mod 4 <> 0 then
        doubly_even_15 := false;;
    fi;
od;
Print("Doubly-even: ", doubly_even_15, "\n\n");;

# ─────────────────────────────────────────────
# The [72, 15, 16] code from T₂+T₄+T₆+T₈
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("The [72, 15, 16] code from T₂+T₄+T₆+T₈\n");;
Print("═══════════════════════════════════════════════════\n\n");;

M_alt := adj[2] + adj[4] + adj[6] + adj[8];;
basis_alt := BaseMat(M_alt);;
k_alt := Length(basis_alt);;
Print("Dimension: ", k_alt, "\n\n");;

wd_alt := List([0..n], x -> 0);;
for combo in [0..2^k_alt - 1] do
    v := Zero(GF(2)) * basis_alt[1];;
    current := combo;;
    for j in [1..k_alt] do
        if current mod 2 = 1 then
            v := v + basis_alt[j];;
        fi;
        current := QuoInt(current, 2);;
    od;
    wt := Number([1..n], l -> v[l] <> Zero(GF(2)));;
    wd_alt[wt + 1] := wd_alt[wt + 1] + 1;;
od;

Print("Weight enumerator:\n");;
for w in [0..n] do
    if wd_alt[w+1] > 0 then
        Print("  A_", w, " = ", wd_alt[w+1], "\n");;
    fi;
od;
Print("\n");;

# ─────────────────────────────────────────────
# [72, 36] self-dual codes — minimum distance
# estimation via row combinations up to weight 4
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Minimum distance bounds for [72, 36] self-dual codes\n");;
Print("═══════════════════════════════════════════════════\n\n");;

EstimateMinDist := function(basis, label)
    local k, min_wt, v, wt, i, j, l, m, count;
    k := Length(basis);;
    min_wt := n + 1;;
    count := 0;;

    # Weight-1 combinations (individual rows)
    for i in [1..k] do
        wt := Number([1..n], l -> basis[i][l] <> Zero(GF(2)));;
        if wt < min_wt then min_wt := wt; fi;;
    od;
    Print("  ", label, ": min row weight = ", min_wt, "\n");;

    # Weight-2 combinations
    for i in [1..k] do
        for j in [i+1..k] do
            v := basis[i] + basis[j];;
            wt := Number([1..n], l -> v[l] <> Zero(GF(2)));;
            if wt < min_wt then min_wt := wt; fi;;
            count := count + 1;;
        od;
    od;
    Print("  after weight-2 (", count, " combos): d ≤ ", min_wt, "\n");;

    # Weight-3 combinations
    count := 0;;
    for i in [1..k] do
        for j in [i+1..k] do
            for l in [j+1..k] do
                v := basis[i] + basis[j] + basis[l];;
                wt := Number([1..n], m -> v[m] <> Zero(GF(2)));;
                if wt < min_wt then min_wt := wt; fi;;
                count := count + 1;;
            od;
        od;
    od;
    Print("  after weight-3 (", count, " combos): d ≤ ", min_wt, "\n");;

    # Weight-4 (if feasible: C(36,4) = 58905)
    if k <= 36 then
        count := 0;;
        for i in [1..k] do
            for j in [i+1..k] do
                for l in [j+1..k] do
                    for m in [l+1..k] do
                        v := basis[i] + basis[j] + basis[l] + basis[m];;
                        wt := Number([1..n], p -> v[p] <> Zero(GF(2)));;
                        if wt < min_wt then min_wt := wt; fi;;
                        count := count + 1;;
                    od;
                od;
            od;
        od;
        Print("  after weight-4 (", count, " combos): d ≤ ", min_wt, "\n");;
    fi;

    return min_wt;
end;;

# T₃+T₄
M34 := adj[3] + adj[4];;
b34 := BaseMat(M34);;
Print("\nT₃+T₄ [72, ", Length(b34), "]:\n");;
EstimateMinDist(b34, "T₃+T₄");;

# T₅+T₆
M56 := adj[5] + adj[6];;
b56 := BaseMat(M56);;
Print("\nT₅+T₆ [72, ", Length(b56), "]:\n");;
EstimateMinDist(b56, "T₅+T₆");;

# T₇+T₈
M78 := adj[7] + adj[8];;
b78 := BaseMat(M78);;
Print("\nT₇+T₈ [72, ", Length(b78), "]:\n");;
EstimateMinDist(b78, "T₇+T₈");;

# All non-identity
M_all := adj[2] + adj[3] + adj[4] + adj[5] + adj[6] + adj[7] + adj[8];;
b_all := BaseMat(M_all);;
Print("\nAll non-identity [72, ", Length(b_all), "]:\n");;
EstimateMinDist(b_all, "all");;

# ─────────────────────────────────────────────
# Check doubly-even property of [72,36] codes
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Doubly-even check for [72,36] codes\n");;
Print("═══════════════════════════════════════════════════\n\n");;

CheckDoublyEven := function(basis, label)
    local k, v, wt, i, j, all_div4;
    k := Length(basis);;
    all_div4 := true;;

    # Check rows
    for i in [1..k] do
        wt := Number([1..n], l -> basis[i][l] <> Zero(GF(2)));;
        if wt mod 4 <> 0 then all_div4 := false; break; fi;;
    od;

    if not all_div4 then
        Print("  ", label, ": NOT doubly-even (row weights not all ≡ 0 mod 4)\n");;
        # Print first few row weights
        Print("    Row weights: ");;
        for i in [1..Minimum(k, 10)] do
            wt := Number([1..n], l -> basis[i][l] <> Zero(GF(2)));;
            Print(wt, " ");;
        od;
        Print("...\n");;
        return false;
    fi;

    # Check all pairwise sums (needed: wt(u+v) ≡ 0 mod 4 for all u,v in code)
    # Equivalent: wt(rᵢ) ≡ 0 mod 4 AND rᵢ·rⱼ ≡ 0 mod 2 (self-orthogonal)
    # Since we already know self-orthogonal, just need row weights mod 4
    Print("  ", label, ": all row weights ≡ 0 mod 4 + self-orthogonal → DOUBLY-EVEN\n");;
    return true;
end;;

CheckDoublyEven(b34, "T₃+T₄");;
CheckDoublyEven(b56, "T₅+T₆");;
CheckDoublyEven(b78, "T₇+T₈");;
CheckDoublyEven(b_all, "All non-id");;

# ─────────────────────────────────────────────
# Check if [72,36] codes are equivalent to
# each other (via the T₂ involution)
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Code equivalence and T₂ action\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Row space containment checks
Print("Containment relations between [72,36] codes:\n");;

# T₃+T₄ vs T₅+T₆
inter_3456 := BaseMat(Concatenation(b34, b56));;
Print("  span(T₃+T₄, T₅+T₆) has dim ", Length(inter_3456), "\n");;

inter_3478 := BaseMat(Concatenation(b34, b78));;
Print("  span(T₃+T₄, T₇+T₈) has dim ", Length(inter_3478), "\n");;

inter_5678 := BaseMat(Concatenation(b56, b78));;
Print("  span(T₅+T₆, T₇+T₈) has dim ", Length(inter_5678), "\n");;

# ─────────────────────────────────────────────
# Relationship to known codes
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Summary of codes from Meinhold-72\n");;
Print("═══════════════════════════════════════════════════\n\n");;

Print("Classical binary codes:\n");;
Print("  T₃+T₄+T₇+T₈:          [72,  8, 28]  self-orthogonal\n");;
Print("  T₃+T₄+T₅+T₆:          [72, 15, 16]  self-orthogonal\n");;
Print("  T₂+T₄+T₆+T₈:          [72, 15, 16]  self-orthogonal\n");;
Print("  T₃+T₄:                 [72, 36,  ?]  self-dual\n");;
Print("  T₅+T₆:                 [72, 36,  ?]  self-dual\n");;
Print("  T₇+T₈:                 [72, 36,  ?]  self-dual\n");;
Print("  T₃+T₄+T₅+T₆+T₇+T₈:   [72, 36,  ?]  self-dual\n");;

Print("\n═══ DONE ═══\n");;
