# meinhold72_irreps.g — Irreducible decomposition of C[G₂(2)/PSL(3,2)]
# and Wedderburn structure of H(G₂(2), PSL(3,2))
#
# Usage: gap -q meinhold72_irreps.g

# ─────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────

G := PrimitiveGroup(63, 2);;   # G₂(2) = PSU(3,3):C₂, order 12096
subs := ConjugacyClassesSubgroups(G);;
K := First(List(subs, Representative), H -> Size(H) = 168);;

Print("G = ", StructureDescription(G), ", |G| = ", Size(G), "\n");;
Print("K = ", StructureDescription(K), ", |K| = ", Size(K), "\n");;
Print("|G/K| = ", Size(G)/Size(K), "\n\n");;

# ─────────────────────────────────────────────
# Character table
# ─────────────────────────────────────────────

Print("Computing character table of G₂(2)...\n");;
ct := CharacterTable(G);;
irrs := Irr(ct);;
Print("Number of irreducible characters: ", Length(irrs), "\n");;
Print("Degrees: ", List(irrs, chi -> chi[1]), "\n");;
Print("Class sizes: ", SizesConjugacyClasses(ct), "\n\n");;

# ─────────────────────────────────────────────
# Permutation character via the action on cosets
# ─────────────────────────────────────────────

Print("Computing permutation representation on G/K...\n");;
cosets := RightCosets(G, K);;
n := Length(cosets);;  # should be 72

# Build the permutation action
act := ActionHomomorphism(G, cosets, OnRight);;
img := Image(act);;
Print("Permutation group on ", NrMovedPoints(img), " points\n");;

# Permutation character: for each g, count fixed points
perm_char_values := List(ConjugacyClasses(G), c ->
    NrMovedPoints(img) - NrMovedPoints(Image(act, Representative(c)))
    + Number([1..NrMovedPoints(img)], p -> p^Image(act, Representative(c)) = p)
);;

perm_char_values := [];;
for c in ConjugacyClasses(G) do
    g := Image(act, Representative(c));;
    fix := Number([1..n], p -> p^g = p);;
    Add(perm_char_values, fix);;
od;

Print("Permutation character values: ", perm_char_values, "\n");;
Print("Check: perm_char(e) = ", perm_char_values[1], " (should be 72)\n\n");;

# Create as a proper class function
perm_char := ClassFunction(ct, perm_char_values);;

# ─────────────────────────────────────────────
# Decompose into irreducibles
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Decomposition of 1_K^G into irreducibles of G₂(2)\n");;
Print("═══════════════════════════════════════════════════\n\n");;

multiplicities := [];;
dim_check := 0;;
for i in [1..Length(irrs)] do
    mi := ScalarProduct(ct, perm_char, irrs[i]);;
    Add(multiplicities, mi);;
    if mi > 0 then
        Print("  chi_", i, ": degree ", irrs[i][1],
              ", multiplicity ", mi, "\n");;
        dim_check := dim_check + mi * irrs[i][1];;
    fi;
od;

Print("\nDimension check: Σ mᵢ·dᵢ = ", dim_check, " (should be 72)\n");;

# ─────────────────────────────────────────────
# Wedderburn decomposition
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Wedderburn decomposition: H ≅ ⊕ M_{mᵢ}(C)\n");;
Print("═══════════════════════════════════════════════════\n\n");;

dim_sq_check := 0;;
active_irreps := [];;
for i in [1..Length(irrs)] do
    if multiplicities[i] > 0 then
        Add(active_irreps, i);;
        Print("  M_", multiplicities[i], "(C)  from chi_", i,
              " (degree ", irrs[i][1], ")\n");;
        dim_sq_check := dim_sq_check + multiplicities[i]^2;;
    fi;
od;

Print("\nDimension check: Σ mᵢ² = ", dim_sq_check, " (should be 8)\n");;
Print("Number of simple summands: ", Length(active_irreps), "\n\n");;

# ─────────────────────────────────────────────
# Build adjacency matrices for the 8 orbitals
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Building orbital adjacency matrices\n");;
Print("═══════════════════════════════════════════════════\n\n");;

dcosets := DoubleCosets(G, K, K);;
r := Length(dcosets);;
valencies := List(dcosets, d -> Size(d) / Size(K));;
Print("Double cosets: ", r, ", valencies: ", valencies, "\n\n");;

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
coset_reps := List(cosets, Representative);;

# Build n×n adjacency matrices
adj := List([1..r], d -> NullMat(n, n));;
for i in [1..n] do
    x := coset_reps[i];;
    for j in [1..n] do
        y := coset_reps[j];;
        prod := x * y^(-1);;
        d := dc_of_elt[PositionSorted(elt_list, prod)];;
        adj[d][i][j] := 1;;
    od;
od;

# Verify row sums
for d in [1..r] do
    rs := Set(List([1..n], i -> Sum(adj[d][i])));;
    if rs <> [valencies[d]] then
        Print("ERROR: T_", d, " row sums = ", rs, "\n");;
    fi;
od;
Print("All adjacency matrices verified (row sums match valencies).\n\n");;

# ─────────────────────────────────────────────
# Eigenvalues of each adjacency matrix
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Eigenvalues of orbital adjacency matrices\n");;
Print("═══════════════════════════════════════════════════\n\n");;

for d in [1..r] do
    Print("T_", d, " (valency ", valencies[d], "):\n");;
    cp := CharacteristicPolynomial(Rationals, Rationals, adj[d], 1);;
    facs := Collected(Factors(cp));;
    Print("  Eigenvalues: ");;
    for f in facs do
        deg := DegreeOfLaurentPolynomial(f[1]);;
        mult := f[2];;
        if deg = 1 then
            coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
            ev := -coeffs[1];;
            Print(ev, " (×", mult, ")  ");;
        else
            Print("[irred deg ", deg, "] (×", mult, ")  ");;
        fi;
    od;
    Print("\n");;
od;

# ─────────────────────────────────────────────
# K-averaging projector and isotypic projectors
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("K-fixed vectors and Hecke algebra representation matrices\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Build K-averaging projector Q = (1/|K|) Σ_{k∈K} P_k
Print("Building K-averaging projector...\n");;

# Precompute: for each k in K, its permutation on cosets
# k sends coset i to coset sigma_k(i) where k*x_i*K = x_{sigma(i)}*K

# Use the action homomorphism we already have
Q := NullMat(n, n);;
for k in K do
    perm := Image(act, k);;
    for i in [1..n] do
        Q[i][i^perm] := Q[i][i^perm] + 1;;
    od;
od;
Q := Q / Size(K);;
Print("Rank(Q) = ", RankMat(Q), " (should be ", r, " = rank of H)\n\n");;

# Build isotypic projectors P_i = (dᵢ/|G|) Σ_{g∈G} conj(χᵢ(g)) P_g
# Optimization: sum over conjugacy classes instead of all elements

Print("Building isotypic projectors and extracting Hecke matrices...\n\n");;

cc_list := ConjugacyClasses(G);;

for idx in [1..Length(active_irreps)] do
    i := active_irreps[idx];;
    chi := irrs[i];;
    di := chi[1];;
    mi := multiplicities[i];;

    Print("────────────────────────────────────────\n");;
    Print("chi_", i, ": degree ", di, ", multiplicity ", mi, "\n");;
    Print("────────────────────────────────────────\n");;

    # P_i = (dᵢ/|G|) Σ_g conj(χᵢ(g)) P_g
    Pi := NullMat(n, n);;
    for j in [1..Length(cc_list)] do
        chi_val := ComplexConjugate(
            ValuesOfClassFunction(chi)[j]
        );;
        if IsZero(chi_val) then continue; fi;;
        for g in cc_list[j] do
            perm := Image(act, g);;
            for ii in [1..n] do
                Pi[ii][ii^perm] := Pi[ii][ii^perm] + chi_val;;
            od;
        od;
    od;
    Pi := Pi * di / Size(G);;

    # K-fixed vectors in isotypic component = image of P_i * Q
    PiQ := Pi * Q;;
    basis_raw := SemiEchelonMat(PiQ);;
    basis_vecs := basis_raw.vectors;;
    Print("  dim(V_i^K) = ", Length(basis_vecs), " (should be ", mi, ")\n\n");;

    if mi = 1 then
        # Multiplicity 1: each T_d acts as a scalar
        Print("  Eigenvalues (T_d acts as scalar on V_i):\n");;
        for d in [1..r] do
            # T_d * v = lambda * v, so lambda = (v * adj[d] * v^T) / (v * v^T)
            v := basis_vecs[1];;
            img_v := v * adj[d];;
            lambda := (img_v * v) / (v * v);;
            Print("    T_", d, " (val ", valencies[d], "): λ = ", lambda, "\n");;
        od;
    else
        # Multiplicity m: T_d acts as an m×m matrix
        Print("  Hecke representation matrices (T_d on V_i^K):\n");;
        gram := basis_vecs * TransposedMat(basis_vecs);;
        gram_inv := gram^(-1);;

        for d in [1..r] do
            images := basis_vecs * adj[d];;
            M_d := images * TransposedMat(basis_vecs) * gram_inv;;
            Print("    T_", d, " (val ", valencies[d], "): ");;
            if mi = 2 then
                Print("[[", M_d[1][1], ", ", M_d[1][2], "], [",
                      M_d[2][1], ", ", M_d[2][2], "]]");;
            else
                Print(M_d);;
            fi;
            Print("\n");;
        od;

        # Eigenvalues of each block
        Print("\n  Block eigenvalues:\n");;
        for d in [1..r] do
            images := basis_vecs * adj[d];;
            M_d := images * TransposedMat(basis_vecs) * gram_inv;;
            cp := CharacteristicPolynomial(Rationals, Rationals, M_d, 1);;
            facs := Collected(Factors(cp));;
            Print("    T_", d, ": ");;
            for f in facs do
                deg := DegreeOfLaurentPolynomial(f[1]);;
                if deg = 1 then
                    coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
                    Print(-coeffs[1], " ");;
                else
                    Print("[irred deg ", deg, "] ");;
                fi;
            od;
            Print("\n");;
        od;
    fi;
    Print("\n");;
od;

# ─────────────────────────────────────────────
# Summary: eigenvalue table
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Eigenvalue table (rows = irreps, cols = orbitals)\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Header
Print("              ");;
for d in [1..r] do
    Print("T_", d, "(", valencies[d], ")  ");;
od;
Print("\n");;

for idx in [1..Length(active_irreps)] do
    i := active_irreps[idx];;
    chi := irrs[i];;
    di := chi[1];;
    mi := multiplicities[i];;

    Pi := NullMat(n, n);;
    for j in [1..Length(cc_list)] do
        chi_val := ComplexConjugate(ValuesOfClassFunction(chi)[j]);;
        if IsZero(chi_val) then continue; fi;;
        for g in cc_list[j] do
            perm := Image(act, g);;
            for ii in [1..n] do
                Pi[ii][ii^perm] := Pi[ii][ii^perm] + chi_val;;
            od;
        od;
    od;
    Pi := Pi * di / Size(G);;

    PiQ := Pi * Q;;
    basis_vecs := SemiEchelonMat(PiQ).vectors;;

    Print("d=", di, " m=", mi, ":  ");;

    if mi = 1 then
        v := basis_vecs[1];;
        for d in [1..r] do
            img_v := v * adj[d];;
            lambda := (img_v * v) / (v * v);;
            Print(lambda, "  ");;
        od;
    else
        gram := basis_vecs * TransposedMat(basis_vecs);;
        gram_inv := gram^(-1);;
        for d in [1..r] do
            images := basis_vecs * adj[d];;
            M_d := images * TransposedMat(basis_vecs) * gram_inv;;
            Print("[", TraceMat(M_d), "/", mi, "]  ");;
        od;
    fi;
    Print("\n");;
od;

Print("\n═══ DONE ═══\n");;
