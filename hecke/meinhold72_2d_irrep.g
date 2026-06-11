# meinhold72_2d_irrep.g — Extract the 2D irrep of H(G₂(2), PSL(3,2))
# from the left regular representation.
#
# Strategy: L(T₇) has a 3D eigenspace at λ=3. The 2D irrep lives inside
# this 3D space (the third dimension is from a 1D irrep where T₇→3).
# Find the 1D sub (common eigenvector of all L(Tᵢ)) and take its complement.
#
# Usage: gap -q meinhold72_2d_irrep.g

Read("/tmp/sc72.g");;
sc := struct_const_72;;
r := 8;;
vals := valencies_72;;

Print("Valencies: ", vals, "\n\n");;

# ─────────────────────────────────────────────
# Left regular representation
# L(Tᵢ)[j][k] = c_{ij}^k
# ─────────────────────────────────────────────

L := List([1..r], i -> List([1..r], j -> List([1..r], k -> sc[i][j][k])));;

# ─────────────────────────────────────────────
# Find the 3D eigenspace of L(T₇) at λ = 3
# ─────────────────────────────────────────────

E3 := NullspaceMat(L[7] - 3 * IdentityMat(r));;
Print("dim Eigenspace(L(T₇), 3) = ", Length(E3), "\n");;
Print("Basis:\n");;
for v in E3 do Print("  ", v, "\n"); od;
Print("\n");;

# ─────────────────────────────────────────────
# Within this 3D space, find the 1D irrep(s).
# A common eigenvector of all L(Tᵢ) restricted
# to this eigenspace is a 1D irrep.
# ─────────────────────────────────────────────

# Restrict each L(Tᵢ) to the 3D eigenspace
gram := E3 * TransposedMat(E3);;
gram_inv := gram^(-1);;

R := List([1..r], d ->
    (E3 * L[d]) * TransposedMat(E3) * gram_inv
);;

Print("Restricted matrices (3×3):\n");;
for d in [1..r] do
    Print("  R(T_", d, ") = ", R[d], "\n");;
od;
Print("\n");;

# Find common eigenvectors of all R[d].
# R[2] is an involution (R[2]^2 should be I₃).
# Use R[2] eigenspaces to decompose.

Print("R(T₂)^2 = I: ", R[2]^2 = IdentityMat(3), "\n");;

ev2_plus := NullspaceMat(R[2] - IdentityMat(3));;
ev2_minus := NullspaceMat(R[2] + IdentityMat(3));;
Print("R(T₂) eigenspace +1: dim ", Length(ev2_plus), "\n");;
Print("R(T₂) eigenspace -1: dim ", Length(ev2_minus), "\n\n");;

# The 1D irreps at T₇=3 are χ₁₀ (where T₂→1) and χ₁₁ (where T₂→-1).
# The 2D irrep should be the remaining space.

# In the +1 eigenspace, check if it's a common eigenspace for all R[d]
Print("═══ +1 eigenspace of R(T₂) ═══\n");;
for v in ev2_plus do
    Print("  v = ", v, "\n  eigenvalues: ");;
    is_eigenvec := true;;
    for d in [1..r] do
        w := v * R[d];;
        j := First([1..3], i -> v[i] <> 0);;
        ratio := w[j] / v[j];;
        if w = ratio * v then
            Print(ratio, " ");;
        else
            Print("* ");;
            is_eigenvec := false;;
        fi;
    od;
    Print("\n  common eigenvector: ", is_eigenvec, "\n\n");;
od;

Print("═══ -1 eigenspace of R(T₂) ═══\n");;
for v in ev2_minus do
    Print("  v = ", v, "\n  eigenvalues: ");;
    is_eigenvec := true;;
    for d in [1..r] do
        w := v * R[d];;
        j := First([1..3], i -> v[i] <> 0);;
        ratio := w[j] / v[j];;
        if w = ratio * v then
            Print(ratio, " ");;
        else
            Print("* ");;
            is_eigenvec := false;;
        fi;
    od;
    Print("\n  common eigenvector: ", is_eigenvec, "\n\n");;
od;

# ─────────────────────────────────────────────
# Also look at the λ = -3 eigenspace of L(T₇)
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Eigenspace of L(T₇) at λ = -3\n");;
Print("═══════════════════════════════════════════════════\n\n");;

Em3 := NullspaceMat(L[7] + 3 * IdentityMat(r));;
Print("dim = ", Length(Em3), "\n");;

gram_m := Em3 * TransposedMat(Em3);;
gram_m_inv := gram_m^(-1);;
Rm := List([1..r], d ->
    (Em3 * L[d]) * TransposedMat(Em3) * gram_m_inv
);;

Print("Restricted matrices (3×3):\n");;
for d in [1..r] do
    Print("  R(T_", d, ") = ", Rm[d], "\n");;
od;
Print("\n");;

# ─────────────────────────────────────────────
# Approach 2: Use L(T₈) which has eigenvalue 21
# with multiplicity 2 — only the two 1D irreps
# where T₈ = 21 (χ₁: all vals, χ₁₀: T₈=3??)
# Actually, let's use T₇ AND T₈ together
# to separate the 1D and 2D components.
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Using multiple operators to isolate the 2D irrep\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# From the 1D irreps (common eigenvectors of all L[d]):
# They must be in the λ=21 eigenspace of L[7] or λ=-21 eigenspace.
# We saw: λ=21 gives [1,1,1,1,1,1,1,1] with eigenvalues (1,1,7,7,7,7,21,21)
#          λ=-21 gives [1,-1,-1,1,-1,1,-1,1] with eigenvalues (1,-1,-7,7,-7,7,-21,21)
# And [0,0,-1,-1,1,1,0,0] at λ=3 with eigenvalues (1,1,*,...) is not a common eigenvec.

# So χ₁ and χ₂ have T₇ eigenvalues 21 and -21.
# χ₁₀ and χ₁₁ have T₇ eigenvalue... let me check.
# From the output: λ=3 has dim 3 (containing 1D irreps χ₁₀, χ₁₁ + part of 2D)
# λ=-3 has dim 3 (containing part of 2D + maybe something else)

# Let me find ALL common eigenvectors across all 8 operators simultaneously.
# These are the 1D irreps of H.

Print("Finding all 1D irreps (common eigenvectors of all L[d])...\n\n");;

# Method: intersect eigenspaces progressively.
# Start with L[2] eigenspace +1 (dim 4), then intersect with eigenspaces of L[3], etc.

# L[2] eigenspace +1
space_plus := NullspaceMat(L[2] - IdentityMat(r));;
Print("L[2]=+1: dim ", Length(space_plus), "\n");;

# Restrict L[8] to this 4D space and find eigenspaces
gr := space_plus * TransposedMat(space_plus);;
gri := gr^(-1);;
L8_rest := (space_plus * L[8]) * TransposedMat(space_plus) * gri;;
Print("L[8] restricted to L[2]=+1 eigenspace:\n");;
Print("  ", L8_rest, "\n");;
cp8 := CharacteristicPolynomial(Rationals, Rationals, L8_rest, 1);;
Print("  char poly: ", cp8, "\n");;
facs8 := Collected(Factors(cp8));;
for f in facs8 do
    deg := DegreeOfLaurentPolynomial(f[1]);;
    mult := f[2];;
    if deg = 1 then
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        ev := -coeffs[1];;
        Print("  eigenvalue ", ev, " (×", mult, ")\n");;
        sub := NullspaceMat(L8_rest - ev * IdentityMat(4));;
        # Map back to H
        for v in sub do
            w := v * space_plus;;
            Print("    ", w, " → eigenvalues: ");;
            for d in [1..r] do
                u := w * L[d];;
                j := First([1..r], i -> w[i] <> 0);;
                ratio := u[j] / w[j];;
                if u = ratio * w then Print(ratio, " ");
                else Print("* "); fi;;
            od;
            Print("\n");;
        od;
    else
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        Print("  irred factor deg ", deg, ": ");;
        Print(f[1], " (×", mult, ")\n");;
        sub := NullspaceMat(Value(f[1], L8_rest));;
        Print("  kernel dim = ", Length(sub), "\n");;
        # This is the 2D irrep restricted to the +1 eigenspace of T₂
        # Map back to H
        basis_2d := sub * space_plus;;
        Print("  Basis vectors in H:\n");;
        for v in basis_2d do
            Print("    ", v, "\n");;
        od;
    fi;
od;

Print("\n");;

# Similarly for L[2] eigenspace -1
space_minus := NullspaceMat(L[2] + IdentityMat(r));;
Print("L[2]=-1: dim ", Length(space_minus), "\n");;

grm := space_minus * TransposedMat(space_minus);;
grmi := grm^(-1);;
L8_rest_m := (space_minus * L[8]) * TransposedMat(space_minus) * grmi;;
cp8m := CharacteristicPolynomial(Rationals, Rationals, L8_rest_m, 1);;
facs8m := Collected(Factors(cp8m));;
for f in facs8m do
    deg := DegreeOfLaurentPolynomial(f[1]);;
    mult := f[2];;
    if deg = 1 then
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        ev := -coeffs[1];;
        Print("  eigenvalue ", ev, " (×", mult, ")\n");;
        sub := NullspaceMat(L8_rest_m - ev * IdentityMat(4));;
        for v in sub do
            w := v * space_minus;;
            Print("    ", w, " → eigenvalues: ");;
            for d in [1..r] do
                u := w * L[d];;
                j := First([1..r], i -> w[i] <> 0);;
                ratio := u[j] / w[j];;
                if u = ratio * w then Print(ratio, " ");
                else Print("* "); fi;;
            od;
            Print("\n");;
        od;
    else
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        Print("  irred factor deg ", deg, ": ");;
        Print(f[1], " (×", mult, ")\n");;
        sub := NullspaceMat(Value(f[1], L8_rest_m));;
        Print("  kernel dim = ", Length(sub), "\n");;
        basis_2d_m := sub * space_minus;;
        Print("  Basis vectors in H:\n");;
        for v in basis_2d_m do
            Print("    ", v, "\n");;
        od;
    fi;
od;

# ─────────────────────────────────────────────
# Extract the 2D irrep from the 4D invariant
# subspace (combine the two 2D pieces from
# the ±1 eigenspaces of T₂)
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("The 2D irrep over Q\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# The 2D irrep V of H is a simple left H-module.
# In the left regular representation (H acting on itself),
# V appears with multiplicity 2 (Artin-Wedderburn: simple
# module of M₂(C) has dim 2, appears dim=2 times in the
# regular rep of M₂(C)).
#
# Over Q, the two copies of V fuse into a 4D Q-irreducible.
# We found 2D pieces in each T₂-eigenspace.
# Take the 2D piece from the +1 eigenspace.

if IsBound(basis_2d) then
    Print("2D invariant subspace (from T₂ = +1 block):\n");;
    for v in basis_2d do Print("  ", v, "\n"); od;;
    Print("\n");;

    gr2 := basis_2d * TransposedMat(basis_2d);;
    gr2i := gr2^(-1);;

    Print("Representation matrices (2×2 over Q):\n\n");;
    for d in [1..r] do
        images := basis_2d * L[d];;
        M := images * TransposedMat(basis_2d) * gr2i;;
        Print("  ρ(T_", d, ") = ", M, "  [val ", vals[d], "]\n");;
    od;

    Print("\nVerification:\n");;
    # Check ρ(T₂)
    im2 := basis_2d * L[2];;
    M2 := im2 * TransposedMat(basis_2d) * gr2i;;
    Print("  ρ(T₂)² = ", M2^2, " (should be I)\n");;

    # Check some products
    for i in [2..4] do
        for j in [i..5] do
            imi := basis_2d * L[i];;
            Mi := imi * TransposedMat(basis_2d) * gr2i;;
            imj := basis_2d * L[j];;
            Mj := imj * TransposedMat(basis_2d) * gr2i;;
            prod_ij := Mi * Mj;;
            # What should TᵢTⱼ be?
            sum_mat := NullMat(2,2);;
            for k in [1..r] do
                imk := basis_2d * L[k];;
                Mk := imk * TransposedMat(basis_2d) * gr2i;;
                sum_mat := sum_mat + sc[i][j][k] * Mk;;
            od;
            if prod_ij = sum_mat then
                Print("  ρ(T_", i, ")·ρ(T_", j, ") = ρ(T_", i, "·T_", j, ") ✓\n");;
            else
                Print("  ρ(T_", i, ")·ρ(T_", j, ") ≠ ρ(T_", i, "·T_", j, ") ✗\n");;
                Print("    product: ", prod_ij, "\n");;
                Print("    expected: ", sum_mat, "\n");;
            fi;
        od;
    od;

    # Characteristic polynomials of each matrix
    Print("\nMinimal/characteristic polynomials of ρ(Tᵢ):\n");;
    for d in [1..r] do
        images := basis_2d * L[d];;
        M := images * TransposedMat(basis_2d) * gr2i;;
        cp := CharacteristicPolynomial(Rationals, Rationals, M, 1);;
        Print("  ρ(T_", d, "): ", cp, "  → eigenvalues: ");;
        facs := Factors(cp);;
        for f in facs do
            deg := DegreeOfLaurentPolynomial(f);;
            if deg = 1 then
                coeffs := CoefficientsOfUnivariatePolynomial(f);;
                Print(-coeffs[1], " ");;
            else
                disc := CoefficientsOfUnivariatePolynomial(f);;
                a := disc[1]; b := disc[2];; # f = a + bx + x²
                Print("(", -b/2, " ± √(", b^2/4 - a, ")) ");;
            fi;
        od;
        Print("\n");;
    od;

    # Commutator [ρ(T₃), ρ(T₄)]
    Print("\nCommutators in the 2×2 rep:\n");;
    for i in [2..6] do
        for j in [i+1..r] do
            imi := basis_2d * L[i];;
            Mi := imi * TransposedMat(basis_2d) * gr2i;;
            imj := basis_2d * L[j];;
            Mj := imj * TransposedMat(basis_2d) * gr2i;;
            comm := Mi * Mj - Mj * Mi;;
            if not IsZero(comm) then
                Print("  [ρ(T_", i, "), ρ(T_", j, ")] = ", comm, "\n");;
            fi;
        od;
    od;
fi;

Print("\n═══ DONE ═══\n");;
