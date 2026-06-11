# meinhold72_center.g — Center, commutator, and correct 2D irrep of H(G₂(2), PSL(3,2))
#
# Usage: gap -q meinhold72_center.g

Read("/tmp/sc72.g");;
sc := struct_const_72;;
r := 8;;

Print("Structure constants loaded.\n");;
Print("Valencies: ", valencies_72, "\n\n");;

# ─────────────────────────────────────────────
# Center of H
# X = Σ xᵢTᵢ is central iff TⱼX = XTⱼ for all j.
# Condition: Σᵢ xᵢ(c_{ji}^k - c_{ij}^k) = 0 for all j,k.
# ─────────────────────────────────────────────

Print("═══════════════════════════════════════════════════\n");;
Print("Center of H\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Build the linear system: for each (j,k), one equation Σᵢ (c_{ji}^k - c_{ij}^k) xᵢ = 0
eqns := [];;
for j in [1..r] do
    for k in [1..r] do
        row := List([1..r], i -> sc[j][i][k] - sc[i][j][k]);;
        if not IsZero(row) then
            Add(eqns, row);;
        fi;
    od;
od;

Print("Number of non-trivial centrality equations: ", Length(eqns), "\n");;
# eqns is 26×8 matrix. We want x ∈ Q^8 with eqns·x = 0 (right null space).
# = left null space of eqns^T = NullspaceMat(TransposedMat(eqns))
center_basis := NullspaceMat(TransposedMat(eqns));;
Print("dim Z(H) = ", Length(center_basis), "\n\n");;

Print("Basis of Z(H):\n");;
for v in center_basis do
    Print("  ");;
    for i in [1..r] do
        if v[i] <> 0 then
            if v[i] = 1 then Print("T_", i);
            elif v[i] = -1 then Print("-T_", i);
            else Print(v[i], "*T_", i); fi;;
            Print(" ");;
        fi;
    od;
    Print("\n");;
od;

# ─────────────────────────────────────────────
# Commutator subspace [H,H]
# Spanned by all TᵢTⱼ - TⱼTᵢ
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Commutator subspace [H,H]\n");;
Print("═══════════════════════════════════════════════════\n\n");;

comm_vecs := [];;
for i in [1..r] do
    for j in [i+1..r] do
        # [Tᵢ, Tⱼ] = TᵢTⱼ - TⱼTᵢ
        comm := List([1..r], k -> sc[i][j][k] - sc[j][i][k]);;
        if not IsZero(comm) then
            Add(comm_vecs, comm);;
            Print("  [T_", i, ", T_", j, "] = ");;
            is_first := true;;
            for k in [1..r] do
                if comm[k] <> 0 then
                    if not is_first and comm[k] > 0 then Print("+ "); fi;;
                    if comm[k] < 0 then Print("- ");
                        if -comm[k] <> 1 then Print(-comm[k], "*"); fi;;
                    elif comm[k] <> 1 then Print(comm[k], "*"); fi;;
                    Print("T_", k, " ");;
                    is_first := false;;
                fi;
            od;
            Print("\n");;
        fi;
    od;
od;

if Length(comm_vecs) = 0 then
    Print("  (algebra is commutative)\n");;
else
    comm_space := BaseMat(comm_vecs);;
    Print("\ndim [H,H] = ", Length(comm_space), "\n\n");;
    Print("Basis of [H,H]:\n");;
    for v in comm_space do
        Print("  ");;
        for i in [1..r] do
            if v[i] <> 0 then
                if v[i] > 0 then Print("+"); fi;;
                Print(v[i], "*T_", i, " ");;
            fi;
        od;
        Print("\n");;
    od;
fi;

# ─────────────────────────────────────────────
# Wedderburn decomposition (from dimensions)
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Wedderburn decomposition\n");;
Print("═══════════════════════════════════════════════════\n\n");;

dim_center := Length(center_basis);;
dim_H := r;;
dim_comm := Length(comm_space);;

Print("dim H = ", dim_H, "\n");;
Print("dim Z(H) = ", dim_center, "\n");;
Print("dim [H,H] = ", dim_comm, "\n");;
Print("dim H/[H,H] = ", dim_H - dim_comm, "\n\n");;

# Number of simple summands = dim Z(H)
# For H ≅ ⊕ M_{nᵢ}(C): Σ nᵢ² = dim H, number of summands = dim Z(H)
Print("Number of simple summands = dim Z(H) = ", dim_center, "\n");;
Print("Need: Σ nᵢ² = ", dim_H, " with ", dim_center, " summands\n\n");;

# Solve: find n₁,...,n_k with k = dim_center, Σ nᵢ² = dim_H
# Enumerate solutions
Print("Possible decompositions (nᵢ² summing to ", dim_H, " in ", dim_center, " parts):\n");;
# Simple brute force
found := false;;
if dim_center = 5 then
    Print("  C ⊕ C ⊕ C ⊕ C ⊕ M₂(C)  (1+1+1+1+4 = 8)\n");;
    found := true;;
elif dim_center = 6 then
    Print("  C⁶  (but this is commutative — inconsistent with [H,H] ≠ 0)\n");;
    Print("  ...unless some summands are field extensions over Q\n");;
elif dim_center = 2 then
    Print("  M₂(C) ⊕ M₂(C)  (4+4 = 8)\n");;
    found := true;;
fi;

# ─────────────────────────────────────────────
# Construct the 2D irrep from scratch using
# the actual 72×72 adjacency matrices
# ─────────────────────────────────────────────

Print("\n═══════════════════════════════════════════════════\n");;
Print("Constructing the 2D irrep directly\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Strategy: the 72×72 adjacency matrices simultaneously block-diagonalize.
# Find the 2D irrep by computing the left regular representation of H
# and simultaneously diagonalizing it.
#
# The left regular representation: H acts on itself by left multiplication.
# The matrix of left-mult by Tᵢ is: L_i[j][k] = c_{ij}^k
# (TᵢTⱼ = Σ_k c_{ij}^k T_k, so the k-th coordinate of Tᵢ·Tⱼ is c_{ij}^k)

Print("Building left regular representation...\n");;
L := List([1..r], i -> List([1..r], j -> List([1..r], k -> sc[i][j][k])));;

# Verify: L[i] is the matrix of left-mult by Tᵢ
# L[1] should be the identity
Print("L[1] = I: ", L[1] = IdentityMat(r), "\n");;

# L[2] should be an involution (L[2]^2 = I)
Print("L[2]^2 = I: ", L[2]^2 = IdentityMat(r), "\n\n");;

# Find eigenvalues of L[2] (the involution)
Print("Eigenvalues of L[2] (involution):\n");;
cp2 := CharacteristicPolynomial(Rationals, Rationals, L[2], 1);;
facs2 := Collected(Factors(cp2));;
for f in facs2 do
    deg := DegreeOfLaurentPolynomial(f[1]);;
    if deg = 1 then
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        Print("  ", -coeffs[1], " with multiplicity ", f[2], "\n");;
    else
        Print("  [irred deg ", deg, "] with multiplicity ", f[2], "\n");;
    fi;
od;
Print("\n");;

# Compute eigenvalues of all L[i]
Print("Eigenvalues of left regular representation:\n\n");;
for d in [1..r] do
    Print("L(T_", d, ") [val ", valencies_72[d], "]:\n");;
    cpd := CharacteristicPolynomial(Rationals, Rationals, L[d], 1);;
    facsd := Collected(Factors(cpd));;
    for f in facsd do
        deg := DegreeOfLaurentPolynomial(f[1]);;
        mult := f[2];;
        if deg = 1 then
            coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
            Print("  λ = ", -coeffs[1], " (×", mult, ")\n");;
        else
            coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
            Print("  min poly: ");;
            for i in [1..Length(coeffs)] do
                if coeffs[i] <> 0 then
                    if i = 1 then Print(coeffs[i]);
                    elif i = 2 then Print(" + ", coeffs[i], "x");
                    else Print(" + ", coeffs[i], "x^", i-1); fi;;
                fi;
            od;
            Print("  (×", mult, ")\n");;
        fi;
    od;
    Print("\n");;
od;

# ─────────────────────────────────────────────
# Extract the 2D irrep from the regular representation
# by finding a 2D invariant subspace
# ─────────────────────────────────────────────

# The regular representation of an 8-dim algebra decomposes as:
# reg = ⊕ nᵢ · Vᵢ  where Vᵢ is the i-th simple module (dim nᵢ)
# For C⁴ ⊕ M₂(C): reg = 1⊕1⊕1⊕1⊕ 2·V₂  (four 1D reps, and V₂ appears 2 times)
# dim reg = 1+1+1+1+2·2 = 8 ✓

# For the 1D reps, find common eigenvectors of all L[i]
# For the 2D rep, find a 2D invariant subspace that isn't an eigenspace

# Step 1: find the simultaneous eigenspaces
# Use L[7] (21-valent, should have distinct eigenvalues on different blocks)

Print("═══════════════════════════════════════════════════\n");;
Print("Simultaneous eigenspaces\n");;
Print("═══════════════════════════════════════════════════\n\n");;

# Find all common eigenvectors of L[3] and L[7]
# (using rational eigenspaces of individual matrices)

# Eigenspaces of L[7] (21-valent orbital)
cp7 := CharacteristicPolynomial(Rationals, Rationals, L[7], 1);;
facs7 := Collected(Factors(cp7));;

Print("Rational eigenspaces of L(T_7):\n");;
for f in facs7 do
    deg := DegreeOfLaurentPolynomial(f[1]);;
    if deg = 1 then
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        ev := -coeffs[1];;
        eigenspace := NullspaceMat(L[7] - ev * IdentityMat(r));;
        Print("  λ = ", ev, ": dim ", Length(eigenspace), "\n");;
        for v in eigenspace do
            Print("    ", v, "\n");;
            # Check what the other L[i] do to this vector
            Print("    eigenvalues: ");;
            for d in [1..r] do
                w := v * L[d];;
                if w = 0 * v then Print("0 ");
                else
                    # Check if w is proportional to v
                    j := First([1..r], i -> v[i] <> 0);;
                    ratio := w[j] / v[j];;
                    if w = ratio * v then
                        Print(ratio, " ");;
                    else
                        Print("* ");;
                    fi;
                fi;
            od;
            Print("\n");;
        od;
    else
        coeffs := CoefficientsOfUnivariatePolynomial(f[1]);;
        eigenspace := NullspaceMat(Value(f[1], L[7]));;
        Print("  min poly: ");;
        for i in [1..Length(coeffs)] do
            if coeffs[i] <> 0 then
                if i = 1 then Print(coeffs[i]);
                elif i = 2 then
                    if coeffs[i] = 1 then Print(" + x");
                    else Print(" + ", coeffs[i], "x"); fi;
                else
                    if coeffs[i] = 1 then Print(" + x^", i-1);
                    else Print(" + ", coeffs[i], "x^", i-1); fi;
                fi;
            fi;
        od;
        Print(": dim ", Length(eigenspace), "\n");;
        # This is the 2D irrep component (over Q, it's a 2D or 4D space)
        for v in eigenspace do
            Print("    ", v, "\n");;
        od;
        # Compute what L[i] look like on this subspace
        Print("  Restriction of H to this subspace:\n");;
        gram := eigenspace * TransposedMat(eigenspace);;
        gram_inv := gram^(-1);;
        for d in [1..r] do
            images := eigenspace * L[d];;
            M_d := images * TransposedMat(eigenspace) * gram_inv;;
            Print("    L(T_", d, ") = ", M_d, "\n");;
        od;
    fi;
od;

Print("\n═══ DONE ═══\n");;
