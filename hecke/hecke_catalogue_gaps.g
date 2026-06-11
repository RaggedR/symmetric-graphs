# hecke_catalogue_gaps.g — Hecke data for the catalogue pages that were
# skipped in the first merge:
#   CR(19)  = orbital/coset graph of (C19:C6)/C3,  n = 38   (F26A sibling)
#   CR(31)  = orbital/coset graph of (C31:C6)/C3,  n = 62   (F26A sibling)
#   Klein {3,7} graph (24 vertices, valency 7)
#   Klein {7,3} graph (56 vertices, cubic)
#   GH(2,8) collinearity graph, 3D4(2)/P1, n = 819 (character table only)
#
# Usage: gap -q hecke_catalogue_gaps.g > hecke_catalogue_gaps.out 2>&1

LoadPackage("grape");;

# ───── machinery (from hecke_irreps3.g) ─────
suborbitdata := function(im, n)
    local st, orbs, mats, prs, o, A, L, p;
    st := Stabilizer(im, 1);;
    orbs := ShallowCopy(Orbits(st, [1..n]));;
    SortBy(orbs, Length);;
    mats := [];;  prs := [];;
    for o in orbs do
        A := NullMat(n, n);;
        L := Orbit(im, [1, o[1]], OnPairs);;
        for p in L do A[p[1]][p[2]] := 1; od;
        Add(mats, A);;  Add(prs, L);;
    od;
    return [mats, prs];
end;;

spin := function(v, mats)
    local V, i, M, w, changed;
    V := [v];;
    changed := true;;
    while changed do
        changed := false;;
        for i in [1..Length(V)] do
            for M in mats do
                w := V[i] * M;;
                if SolutionMat(V, w) = fail then
                    Add(V, w);; changed := true;;
                fi;
            od;
        od;
    od;
    return V;
end;;

rootsfan := function(f)
    local c, b0, b1, d;
    c := CoefficientsOfUnivariatePolynomial(f);;
    if Length(c) = 2 then return [ -c[1]/c[2] ]; fi;
    if Length(c) = 3 then
        b1 := c[2]/c[3];;  b0 := c[1]/c[3];;
        d := b1^2 - 4*b0;;
        return [ (-b1 + Sqrt(d))/2, (-b1 - Sqrt(d))/2 ];
    fi;
    return [];
end;;

analyze := function(name, im, n, mmax)
    local tbl, irr, ccl, csums, C, M, g, x, sd, mats, prs, r, lens,
          pivals, pi, idx, chi, m, P, W, d, B, vals, a, cp, f, lam, E,
          v, V, rho, rhos, evecs, mu, i, j, c2, S, w;
    Print("############ ", name, " ############\n");;
    tbl := CharacterTable(im);;
    irr := Irr(tbl);;
    ccl := ConjugacyClasses(tbl);;
    csums := [];;
    for C in ccl do
        M := NullMat(n, n);;
        for g in Elements(C) do
            for x in [1..n] do M[x][x^g] := M[x][x^g] + 1; od;
        od;
        Add(csums, M);;
    od;
    sd := suborbitdata(im, n);;
    mats := sd[1];;  prs := sd[2];;
    r := Length(mats);;
    lens := List(mats, M -> Sum(M[1]));;
    Print("rank ", r, ",  suborbit lengths ", Collected(lens), "\n\n");;
    pivals := List(ccl, C -> Number([1..n], x -> x^Representative(C) = x));;
    pi := ClassFunction(tbl, pivals);;
    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m = 0 then continue; fi;
        if m > mmax then
            Print("chi_", idx, " (deg ", chi[1], "): mult ", m,
                  " — skipped (> ", mmax, ")\n\n");;
            continue;;
        fi;
        P := (chi[1] / Size(im)) * Sum([1..Length(ccl)],
                 i -> ComplexConjugate(chi[i]) * csums[i]);;
        if m = 1 then
            vals := List(prs, L -> Sum(L, p -> P[p[2]][p[1]]) / chi[1]);;
            Print("1-dim irrep on chi_", idx, " (deg ", chi[1],
                  "):\n  T_i -> ", vals, "\n\n");;
            continue;;
        fi;
        Print(">>> ", m, "-dim irrep on chi_", idx, " (deg ", chi[1],
              ")\n");;
        W := BaseMat(P);;
        d := Length(W);;
        Print("  isotypic dim ", d, " (expect ", m * chi[1], ")\n");;
        B := List(mats, M -> List(W, w -> SolutionMat(W, w * M)));;
        V := fail;;
        for a in [2..r] do
            cp := CharacteristicPolynomial(B[a]);;
            for f in Set(Factors(cp)) do
                for lam in rootsfan(f) do
                    E := NullspaceMat(B[a] - lam * IdentityMat(d));;
                    if Length(E) <> chi[1] then continue; fi;
                    v := E[1];;
                    V := spin(v, B);;
                    if Length(V) = m then
                        Print("  module: eigenvector of T_", a,
                              " for eigenvalue ", lam, "\n");;
                        break;;
                    else
                        V := fail;;
                    fi;
                od;
                if V <> fail then break; fi;
            od;
            if V <> fail then break; fi;
        od;
        if V = fail then
            Print("  !! no irreducible module found by quadratic eigenvalues\n\n");;
            continue;;
        fi;
        rho := List(V, w -> SolutionMat(V, w * B[a]));;
        cp := CharacteristicPolynomial(rho);;
        evecs := [];;
        for f in Set(Factors(cp)) do
            for mu in rootsfan(f) do
                Append(evecs, NullspaceMat(rho - mu * IdentityMat(m)));;
            od;
        od;
        if Length(evecs) = m then
            V := List(evecs, e -> e * V);;
            Print("  basis: eigenbasis of rho(T_", a, ")\n");;
        else
            Print("  basis: spin basis (rho(T_", a,
                  ") not split over cyclotomics)\n");;
        fi;
        if m = 2 then
            for i in [2..r] do
                rho := List(V, w -> SolutionMat(V, w * B[i]));;
                if rho[1][2] <> 0 and rho[2][1] <> 0 then
                    c2 := rho[1][2] / rho[2][1];;
                    if IsRat(c2) and c2 > 0 then
                        V[2] := Sqrt(c2) * V[2];;
                        Print("  normalized: rho(T_", i, ") symmetric\n");;
                    fi;
                    break;;
                fi;
            od;
        fi;
        rhos := List(B, M -> List(V, w -> SolutionMat(V, w * M)));;
        Print("  matrices (suborbit lengths ", lens, "):\n");;
        for i in [1..r] do
            Print("  T_", i, " (val ", lens[i], ", charpoly ",
                  CharacteristicPolynomial(rhos[i]), "):\n");;
            PrintArray(rhos[i]);;
        od;
        S := [];;
        for i in [1..r] do
            w := Flat(rhos[i]);;
            if Length(S) = 0 or SolutionMat(S, w) = fail then Add(S, w); fi;
            for j in [1..r] do
                w := Flat(rhos[i] * rhos[j]);;
                if SolutionMat(S, w) = fail then
                    Add(S, w);;
                fi;
            od;
        od;
        Print("  check: dim of algebra generated = ", Length(S),
              "  (expect ", m^2, ")\n");;
        Print("  check: rho(sum T_i) = ");;
        PrintArray(Sum(rhos));;
        Print("\n");;
    od;
end;;

# ───── Frobenius pair builder: (C_p : C_6) / C_3 ─────
# mult = element of order 6 in (Z/p)^*; points are residues 0..p-1 -> 1..p
frobpair := function(p, mult)
    local trans, m, F, H;
    trans := PermList(List([0..p-1], x -> ((x+1) mod p) + 1));;
    m := PermList(List([0..p-1], x -> ((mult*x) mod p) + 1));;
    F := Group(trans, m);;
    H := Subgroup(F, [m^2]);;   # order 3
    Print("|F| = ", Size(F), ",  |H| = ", Size(H),
          ",  index = ", Index(F, H), "\n");;
    return Image(FactorCosetAction(F, H));
end;;

# ════════ CR(19): (C19:C6)/C3, n = 38 ════════
im := frobpair(19, 8);;     # ord(8 mod 19) = 6
analyze("CR(19) pair (C19:C6)/C3   (n = 38)", im, 38, 3);;

# ════════ CR(31): (C31:C6)/C3, n = 62 ════════
im := frobpair(31, 26);;    # ord(26 mod 31) = 6  (26 = -5)
analyze("CR(31) pair (C31:C6)/C3   (n = 62)", im, 62, 3);;

# ════════ Klein graphs ════════
Read("/Users/robin/git/hecke/klein_adj.g");;

kleingraph := function(name, adj)
    local n, gamma, A, orbs;
    n := Length(adj);;
    gamma := Graph(Group(()), [1..n], OnPoints,
                   function(x, y) return y in adj[x]; end, true);;
    A := AutGroupGraph(gamma);;
    Print("### ", name, ": n = ", n, ", |Aut| = ", Size(A),
          ", transitive = ", IsTransitive(A, [1..n]), "\n");;
    Print("girth = ", Girth(gamma), ", diameter = ", Diameter(gamma), "\n");;
    if IsTransitive(A, [1..n]) then
        analyze(name, A, n, 3);;
    else
        orbs := Orbits(A, [1..n]);;
        Print("NOT vertex-transitive: orbits ", List(orbs, Length), "\n\n");;
    fi;
end;;

kleingraph("Klein {3,7} graph (24 vertices)", klein37_adj);;
kleingraph("Klein {7,3} graph (56 vertices, cubic)", klein73_adj);;

# ════════ GH(2,8): 3D4(2)/P1, n = 819 — character table only ════════
Print("############ GH(2,8) = 3D4(2)/P1   (n = 819) ############\n");;
t := CharacterTable("3D4(2)");;
Print("maxes: ", Maxes(t), "\n");;
for nm in Maxes(t) do
    s := CharacterTable(nm);;
    if Size(t) / Size(s) = 819 then
        Print("subgroup ", nm, ", index 819\n");;
        pi := TrivialCharacter(s)^t;;
        mults := MatScalarProducts(t, Irr(t), [pi])[1];;
        Print("decomposition (deg, mult): ",
              Filtered(List([1..Length(mults)],
                       i -> [Irr(t)[i][1], mults[i]]),
                       p -> p[2] > 0), "\n");;
        Print("rank = ", Sum(mults, m -> m^2),
              ",  Gelfand = ", ForAll(mults, m -> m <= 1), "\n");;
    fi;
od;
# eigenvalues of the collinearity graph from intersection array {18,16,16;1,1,9}
B1 := [[0,18,0,0],[1,1,16,0],[0,1,1,16],[0,0,9,9]];;
Print("intersection matrix eigenvalue charpoly factors: ",
      Factors(CharacteristicPolynomial(B1)), "\n");;
Print("════════ DONE ════════\n");;
