# hecke_irreps.g — explicit irreps (matrix entries) of the non-Gelfand
# Hecke algebras found by hecke_census.g.
#
# Method: orbital matrices A_i on C[G/H]; isotypic projector
#   P_chi = (chi(1)/|G|) sum_g chi(g^-1) pi(g);
# the chi-isotypic W = (m-dim irrep) (x) C^chi(1).  An eigenvector of an
# orbital matrix, for an eigenvalue simple in the irrep, is a pure tensor,
# so its module closure under {A_i} is ONE m-dim copy.  Restrict, get
# exact m x m matrices over cyclotomics.  1-dim irreps: spherical values
# trace(A_i P)/chi(1).
#
# Usage: gap -q hecke_irreps.g > hecke_irreps.out 2>&1

LoadPackage("grape");;

lcfedges := function(code, reps)
    local full, m, edges, i, j;
    full := Concatenation(List([1..reps], r -> code));;
    m := Length(full);;
    edges := [];;
    for i in [1..m] do
        Add(edges, [i, (i mod m) + 1]);;
        j := ((i - 1 + full[i]) mod m) + 1;;
        if i < j then Add(edges, [i, j]); else Add(edges, [j, i]); fi;
    od;
    return Set(edges);;
end;;

suborbitmats := function(im, n)
    local st, orbs, mats, o, A, p;
    st := Stabilizer(im, 1);;
    orbs := ShallowCopy(Orbits(st, [1..n]));;
    SortBy(orbs, Length);;
    mats := [];;
    for o in orbs do
        A := NullMat(n, n);;
        for p in Orbit(im, [1, o[1]], OnPairs) do A[p[1]][p[2]] := 1; od;
        Add(mats, A);;
    od;
    return mats;
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

# roots of a rational polynomial factor, degree <= 2 (else [])
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
    local tbl, irr, ccl, csums, C, M, g, x, mats, r, lens, pivals, pi,
          idx, chi, m, P, W, vals, a, cp, facs, f, lam, E, I2, v, V,
          rho, rhos, evecs, mu, i, j, c2, S, w;
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
    mats := suborbitmats(im, n);;
    r := Length(mats);;
    lens := List(mats, M -> Sum(M[1]));;
    Print("rank ", r, ",  suborbit lengths ", lens, "\n\n");;
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
            vals := List(mats, M -> TraceMat(M * P) / chi[1]);;
            Print("1-dim irrep on chi_", idx, " (deg ", chi[1],
                  "):\n  T_i -> ", vals, "\n\n");;
            continue;;
        fi;
        Print(">>> ", m, "-dim irrep on chi_", idx, " (deg ", chi[1],
              ")\n");;
        W := BaseMat(P);;
        V := fail;;
        for a in [2..r] do
            cp := CharacteristicPolynomial(mats[a]);;
            for f in Set(Factors(cp)) do
                for lam in rootsfan(f) do
                    E := NullspaceMat(mats[a] - lam * IdentityMat(n));;
                    if Length(E) = 0 then continue; fi;
                    I2 := SumIntersectionMat(E, W)[2];;
                    if Length(I2) = chi[1] then
                        v := I2[1];;
                        V := spin(v, mats);;
                        if Length(V) = m then
                            Print("  module: eigenvector of T_", a,
                                  " for eigenvalue ", lam, "\n");;
                            break;;
                        else
                            V := fail;;
                        fi;
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
        # diagonalize rho(T_a) if its eigenvalues are (at worst) quadratic
        rho := List(V, w -> SolutionMat(V, w * mats[a]));;
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
        # symmetric normalization for m = 2
        if m = 2 then
            for i in [2..r] do
                rho := List(V, w -> SolutionMat(V, w * mats[i]));;
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
        rhos := List(mats, M -> List(V, w -> SolutionMat(V, w * M)));;
        Print("  matrices (suborbit lengths ", lens, "):\n");;
        for i in [1..r] do
            Print("  T_", i, " (val ", lens[i], ", charpoly ",
                  CharacteristicPolynomial(rhos[i]), "):\n");;
            PrintArray(rhos[i]);;
        od;
        # checks: sum of rhos = J-image? algebra dimension = m^2
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

# ════════ Case 1: PSL(2,5)/C3, n = 20 ════════
G := PSL(2,5);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
im := Image(FactorCosetAction(G, H));;
analyze("PSL(2,5) / C3   (n = 20)", im, Index(G, H), 3);;

# ════════ Case 2: PSL(2,7)/C3, n = 56 ════════
G := PSL(2,7);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
im := Image(FactorCosetAction(G, H));;
analyze("PSL(2,7) / C3   (n = 56)", im, Index(G, H), 3);;

# ════════ Case 3: F26A, n = 26 ════════
edges := lcfedges([-7, 7], 13);;
dir := Concatenation(edges, List(edges, e -> [e[2], e[1]]));;
gamma := EdgeOrbitsGraph(Group(()), dir, 26);;
A := AutGroupGraph(gamma);;
analyze("F26A   (Aut = C13:C6, n = 26)", A, 26, 3);;

# ════════ Case 4: G2(2)/PSL(3,2) = Meinhold-72 (cross-check vs why13) ════════
G63 := PrimitiveGroup(63, 2);;
K := Group(
  ( 1,33)( 2,61)( 3,57)( 4,27)( 5,47)( 6,36)( 7,60)( 8,32)( 9,59)(10,54)
    (12,37)(13,55)(14,29)(16,48)(17,20)(18,26)(19,31)(21,40)(23,56)(24,63)
    (25,39)(28,51)(30,41)(35,44)(38,58)(42,46)(43,62)(49,53),
  ( 1,12)( 2,56)( 3,55)( 4,51)( 5,59)( 6,61)( 7,39)( 8,34)( 9,57)(10,54)
    (11,22)(13,47)(14,38)(15,31)(16,36)(18,46)(19,52)(20,30)(23,48)(24,53)
    (25,27)(26,37)(28,60)(32,45)(33,42)(35,49)(43,63)(44,62),
  ( 2, 4)( 3,27)( 5,28)( 6, 9)( 7,36)( 8,42)(10,20)(11,14)(13,48)(15,44)
    (16,39)(17,41)(18,46)(19,24)(21,50)(22,38)(23,47)(25,55)(26,32)(30,54)
    (31,62)(33,34)(37,45)(43,63)(51,56)(52,53)(57,61)(59,60),
  ( 1,38)( 2,51)( 3,20)( 4,43)( 5,41)( 6,35)( 7,53)( 8,40)(10,31)(11,50)
    (12,37)(13,55)(14,26)(15,52)(16,25)(17,57)(18,29)(19,54)(21,32)(23,24)
    (27,62)(28,61)(30,47)(33,58)(36,44)(39,48)(49,60)(56,63) );;
im := Image(FactorCosetAction(G63, K));;
analyze("G2(2) / PSL(3,2)   (Meinhold-72, n = 72)", im, 72, 2);;

Print("════════ DONE ════════\n");;
