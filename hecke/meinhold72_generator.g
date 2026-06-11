# meinhold72_generator.g — Extract the 8×72 generator matrix of the Meinhold code
#
# Usage: gap -q meinhold72_generator.g

G := PrimitiveGroup(63, 2);;
subs := ConjugacyClassesSubgroups(G);;
K := First(List(subs, Representative), H -> Size(H) = 168);;
cosets := RightCosets(G, K);;
n := Length(cosets);;
act := ActionHomomorphism(G, cosets, OnRight);;
dcosets := DoubleCosets(G, K, K);;
r := Length(dcosets);;
valencies := List(dcosets, d -> Size(d) / Size(K));;

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

# ─────────────────────────────────────────────
# Find the [72, 8, 32] code
# Try all 2^r subsets of orbitals, look for [72, 8, 32]
# ─────────────────────────────────────────────

Print("Searching for the [72, 8, 32] code...\n\n");;

for mask in [1..2^r - 1] do
    M := NullMat(n, n, GF(2));;
    current := mask;;
    label := "";;
    for d in [1..r] do
        if current mod 2 = 1 then
            M := M + adj[d];;
            Append(label, Concatenation("+T", String(d), "(", String(valencies[d]), ")"));;
        fi;
        current := QuoInt(current, 2);;
    od;
    basis := BaseMat(M);;
    k := Length(basis);;
    if k = 8 then
        # Check minimum distance
        min_wt := n + 1;;
        for combo in [1..2^k - 1] do
            v := Zero(GF(2)) * basis[1];;
            c2 := combo;;
            for j in [1..k] do
                if c2 mod 2 = 1 then v := v + basis[j]; fi;;
                c2 := QuoInt(c2, 2);;
            od;
            wt := Number([1..n], l -> v[l] <> Zero(GF(2)));;
            if wt < min_wt then min_wt := wt; fi;;
        od;
        if min_wt = 32 then
            Print("FOUND: mask=", mask, ", orbitals: ", label, "\n");;
            Print("Parameters: [72, 8, 32]\n\n");;

            # Print the generator matrix
            Print("Generator matrix G (8 × 72 over GF(2)):\n\n");;
            for i in [1..k] do
                Print("  ");;
                for j in [1..n] do
                    if basis[i][j] = One(GF(2)) then Print("1");
                    else Print("0"); fi;;
                od;
                Print("\n");;
            od;
            Print("\n");;

            # Row weights
            Print("Row weights: ");;
            for i in [1..k] do
                wt := Number([1..n], l -> basis[i][l] <> Zero(GF(2)));;
                Print(wt, " ");;
            od;
            Print("\n\n");;

            # Verify self-orthogonality
            Print("G · G^T = 0: ", IsZero(basis * TransposedMat(basis)), "\n\n");;

            # Print as a list of support sets (positions of 1s)
            Print("Support sets (1-indexed positions of 1s in each row):\n\n");;
            for i in [1..k] do
                Print("  g", i, " = { ");;
                for j in [1..n] do
                    if basis[i][j] = One(GF(2)) then Print(j, " "); fi;;
                od;
                Print("}\n");;
            od;
            Print("\n");;

            # Parity check matrix H (64 × 72)
            # H = null space of G (right kernel)
            null := NullspaceMat(TransposedMat(basis));;
            Print("Parity check matrix dimensions: ", Length(null), " × ", n, "\n");;
            Print("Verify H·G^T = 0: ", IsZero(null * TransposedMat(basis)), "\n\n");;

            break;
        fi;
    fi;
od;

Print("═══ DONE ═══\n");;
