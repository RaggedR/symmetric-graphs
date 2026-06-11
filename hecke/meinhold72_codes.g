# meinhold72_codes.g — Binary codes from the Meinhold-72 association scheme
#
# Constructs orbital adjacency matrices over GF(2), computes:
#   - [n, k, d] parameters for each orbital code
#   - Codes from sums of orbitals
#   - CSS orthogonality conditions
#   - Self-orthogonal codes
#
# Usage: gap -q meinhold72_codes.g

# ─────────────────────────────────────────────
# Construct G₂(2) and the coset action
# ─────────────────────────────────────────────

G := PrimitiveGroup(63, 2);;   # G₂(2) = PSU(3,3):C₂, order 12096
subs := ConjugacyClassesSubgroups(G);;
K := First(List(subs, Representative), H -> Size(H) = 168);;

Print("G = ", StructureDescription(G), ", |G| = ", Size(G), "\n");;
Print("K = ", StructureDescription(K), ", |K| = ", Size(K), "\n");;

cosets := RightCosets(G, K);;
n := Length(cosets);;
Print("n = |G/K| = ", n, "\n\n");;

act := ActionHomomorphism(G, cosets, OnRight);;

# ─────────────────────────────────────────────
# Double cosets and orbital adjacency matrices
# ─────────────────────────────────────────────

dcosets := DoubleCosets(G, K, K);;
r := Length(dcosets);;
valencies := List(dcosets, d -> Size(d) / Size(K));;
Print("rank = ", r, "\n");;
Print("valencies = ", valencies, "\n\n");;

# Element → double coset lookup
elt_list := AsSortedList(G);;
dc_of_elt := List([1..Size(G)], x -> 0);;
for i in [1..r] do
    for el in dcosets[i] do
        j := PositionSorted(elt_list, el);;
        dc_of_elt[j] := i;;
    od;
od;

# Coset representatives
coset_reps := List(cosets, c -> Representative(c));;

# Build n×n adjacency matrices over GF(2)
Print("Building adjacency matrices over GF(2)...\n");;
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
Print("Done.\n\n");;

# ─────────────────────────────────────────────
# Helper: compute minimum distance of a code
# given by a generator matrix (list of basis vectors)
# ─────────────────────────────────────────────

MinDistCode := function(basis)
    local k, min_wt, combo, v, wt, i, j, l, current;
    k := Length(basis);;
    if k = 0 then return 0; fi;;
    if k > 20 then return -1; fi;;  # too large to enumerate

    min_wt := Length(basis[1]) + 1;;

    # Enumerate all non-zero linear combinations
    for combo in [1..2^k - 1] do
        v := Zero(GF(2)) * basis[1];;
        current := combo;;
        for j in [1..k] do
            if current mod 2 = 1 then
                v := v + basis[j];;
            fi;
            current := QuoInt(current, 2);;
        od;
        wt := Number([1..Length(v)], l -> v[l] <> Zero(GF(2)));;
        if wt < min_wt then
            min_wt := wt;;
        fi;
    od;
    return min_wt;
end;;

# ─────────────────────────────────────────────
# Analyse each orbital code
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Individual orbital codes over GF(2)\n");;
Print("═══════════════════════════════════════════════════\n\n");;

for d in [1..r] do
    basis_d := BaseMat(adj[d]);;
    k := Length(basis_d);;
    dist := MinDistCode(basis_d);;

    # Check self-orthogonality: C ⊆ C⊥ iff A·Aᵀ = 0 over GF(2)
    prod := adj[d] * TransposedMat(adj[d]);;
    self_orth := IsZero(prod);;

    # Check A = Aᵀ (symmetric)
    is_sym := adj[d] = TransposedMat(adj[d]);;

    Print("T_", d, " [valency ", valencies[d], "]:  ");;
    Print("[", n, ", ", k, ", ", dist, "]");;
    if self_orth then Print("  SELF-ORTHOGONAL"); fi;;
    if is_sym then Print("  symmetric"); fi;;
    Print("\n");;
od;

# ─────────────────────────────────────────────
# Codes from sums of orbitals
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Codes from sums of orbital adjacency matrices\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Try all pairs and interesting combinations
# Focus on: paired orbitals (T₃+T₄, T₅+T₆, T₇+T₈),
# same-valency sums, and "distance shells"

interesting := [
    [ [3,4], "T₃+T₄ (paired val-7)" ],
    [ [5,6], "T₅+T₆ (paired val-7)" ],
    [ [7,8], "T₇+T₈ (paired val-21)" ],
    [ [3,5], "T₃+T₅ (two val-7)" ],
    [ [3,6], "T₃+T₆ (val-7 cross)" ],
    [ [3,4,5,6], "T₃+T₄+T₅+T₆ (all val-7)" ],
    [ [3,7], "T₃+T₇ (val-7 + val-21)" ],
    [ [3,4,7,8], "T₃+T₄+T₇+T₈ (paired sums)" ],
    [ [3,4,5,6,7,8], "all non-identity" ],
    [ [2,4,6,8], "T₂+T₄+T₆+T₈ (one from each pair)" ],
    [ [2,3,5,7], "T₂+T₃+T₅+T₇ (other from each pair)" ]
];;

for info in interesting do
    indices := info[1];;
    label := info[2];;
    M := NullMat(n, n, GF(2));;
    for d in indices do
        M := M + adj[d];;
    od;
    basis_M := BaseMat(M);;
    k := Length(basis_M);;
    dist := MinDistCode(basis_M);;

    # Self-orthogonality
    prod := M * TransposedMat(M);;
    self_orth := IsZero(prod);;

    Print(label, ":  [", n, ", ", k, ", ", dist, "]");;
    if self_orth then Print("  SELF-ORTHOGONAL"); fi;;
    Print("\n");;
od;

# ─────────────────────────────────────────────
# CSS conditions: find pairs (C_X, C_Z) with
# C_Z⊥ ⊆ C_X, equivalently C_X · C_Zᵀ = 0
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("CSS quantum code search\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# For each self-orthogonal code C (A·Aᵀ = 0),
# we get a CSS code [[n, k_X + k_Z - n, min(d_X, d_Z)]]
# where C_X = C_Z = C.

Print("Self-dual CSS codes (C_X = C_Z = C, needs A·Aᵀ = 0 mod 2):\n\n");;

# Check all individual orbitals
for d in [1..r] do
    prod := adj[d] * TransposedMat(adj[d]);;
    if IsZero(prod) then
        basis_d := BaseMat(adj[d]);;
        k := Length(basis_d);;
        dist := MinDistCode(basis_d);;
        # CSS parameters: [[n, 2k - n, d]]
        css_k := 2 * k - n;;
        if css_k > 0 then
            Print("  T_", d, ": [[", n, ", ", css_k, ", ", dist, "]]\n");;
        else
            Print("  T_", d, ": self-orth but k=", k, " < n/2=", n/2, " (degenerate)\n");;
        fi;
    fi;
od;

# Check sums
for info in interesting do
    indices := info[1];;
    label := info[2];;
    M := NullMat(n, n, GF(2));;
    for d in indices do
        M := M + adj[d];;
    od;
    prod := M * TransposedMat(M);;
    if IsZero(prod) then
        basis_M := BaseMat(M);;
        k := Length(basis_M);;
        dist := MinDistCode(basis_M);;
        css_k := 2 * k - n;;
        if css_k > 0 then
            Print("  ", label, ": [[", n, ", ", css_k, ", ", dist, "]]\n");;
        else
            Print("  ", label, ": self-orth but k=", k, " (degenerate)\n");;
        fi;
    fi;
od;

# ─────────────────────────────────────────────
# Cross-pair CSS: C_X from one set, C_Z from another
# Condition: M_X · M_Z^T = 0 over GF(2)
# ─────────────────────────────────────────────

Print("\nCross-pair CSS codes (C_X ≠ C_Z, need M_X · M_Zᵀ = 0):\n\n");;

# Build a list of candidate matrices
candidates := [];;
for d in [1..r] do
    Add(candidates, rec(mat := adj[d], label := Concatenation("T_", String(d)),
                        val := valencies[d]));;
od;
for info in interesting do
    M := NullMat(n, n, GF(2));;
    for d in info[1] do M := M + adj[d]; od;;
    Add(candidates, rec(mat := M, label := info[2]));;
od;

# Check all pairs
css_found := 0;;
for i in [1..Length(candidates)] do
    for j in [i..Length(candidates)] do
        prod := candidates[i].mat * TransposedMat(candidates[j].mat);;
        if IsZero(prod) then
            bX := BaseMat(candidates[i].mat);;
            bZ := BaseMat(candidates[j].mat);;
            kX := Length(bX);;
            kZ := Length(bZ);;
            css_k := kX + kZ - n;;
            if css_k > 0 then
                dX := MinDistCode(bX);;
                dZ := MinDistCode(bZ);;
                css_d := Minimum(dX, dZ);;
                Print("  C_X = ", candidates[i].label,
                      " [", n, ",", kX, ",", dX, "]");;
                Print("  ×  C_Z = ", candidates[j].label,
                      " [", n, ",", kZ, ",", dZ, "]");;
                Print("  →  [[", n, ", ", css_k, ", ", css_d, "]]\n");;
                css_found := css_found + 1;;
            fi;
        fi;
    od;
od;

if css_found = 0 then
    Print("  (no non-degenerate cross-pair CSS codes found)\n");;
fi;

# ─────────────────────────────────────────────
# Dual codes: for each orbital code C, compute C⊥
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Dual code parameters\n");;
Print("═══════════════════════════════════════════════════\n\n");;

for d in [1..r] do
    basis_d := BaseMat(adj[d]);;
    k := Length(basis_d);;
    # Dual: null space of adj[d]
    null_d := NullspaceMat(TransposedMat(adj[d]));;
    k_perp := Length(null_d);;
    dist_perp := MinDistCode(null_d);;
    Print("T_", d, "⊥ [valency ", valencies[d], "]:  ");;
    Print("[", n, ", ", k_perp, ", ", dist_perp, "]");;
    Print("  (k + k⊥ = ", k, " + ", k_perp, " = ", k + k_perp, ")\n");;
od;

# ─────────────────────────────────────────────
# Weight enumerators for the most interesting codes
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Weight distributions of small codes\n");;
Print("═══════════════════════════════════════════════════\n\n");;

WeightDistribution := function(basis)
    local k, wd, combo, v, wt, j, current;
    k := Length(basis);;
    if k = 0 or k > 20 then return []; fi;;
    wd := List([0..Length(basis[1])], x -> 0);;
    for combo in [0..2^k - 1] do
        v := Zero(GF(2)) * basis[1];;
        current := combo;;
        for j in [1..k] do
            if current mod 2 = 1 then
                v := v + basis[j];;
            fi;
            current := QuoInt(current, 2);;
        od;
        wt := Number([1..Length(v)], l -> v[l] <> Zero(GF(2)));;
        wd[wt + 1] := wd[wt + 1] + 1;;
    od;
    return wd;
end;;

for d in [1..r] do
    basis_d := BaseMat(adj[d]);;
    k := Length(basis_d);;
    if k <= 15 then
        wd := WeightDistribution(basis_d);;
        Print("T_", d, " [", n, ", ", k, "]:  ");;
        for w in [0..n] do
            if wd[w+1] > 0 then
                Print("A_", w, "=", wd[w+1], " ");;
            fi;
        od;
        Print("\n");;
    fi;
od;

Print("\n═══ DONE ═══\n");;
