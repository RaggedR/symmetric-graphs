# hecke_psl2_11.g — Hecke algebra of PSL(2,11)/C3
#
# PSL(2,11) has order 660, C3 has index 220.
# Expected: rank ~20-30, non-Gelfand, 2D irreps with B^2 = N*I.
# The question: what is N for q = 11?
#
# Usage: gap -q hecke_psl2_11.g 2>&1 | tee hecke_psl2_11.out

LoadPackage("grape");;

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
    Print("rank ", r, ",  suborbit lengths ", Collected(lens), "\n\n");;
    pivals := List(ccl, C -> Number([1..n], x -> x^Representative(C) = x));;
    pi := ClassFunction(tbl, pivals);;

    # decomposition summary
    Print("Permutation character decomposition:\n");;
    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m > 0 then
            Print("  chi_", idx, " (deg ", chi[1], "): mult ", m, "\n");;
        fi;
    od;
    Print("\nWedderburn: ");;
    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m > 0 then
            if m = 1 then Print("C"); else Print("M_", m); fi;
            Print("(", chi[1], ") ");;
        fi;
    od;
    Print("\n\n");;

    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m = 0 then continue; fi;
        if m > mmax then
            Print("chi_", idx, " (deg ", chi[1], "): mult ", m,
                  " -- skipped (> ", mmax, ")\n\n");;
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
        # Find traceless generators and report B^2 = N*I
        Print("\n  === B^2 = N*I analysis ===\n");;
        for i in [2..r] do
            if TraceMat(rhos[i]) = 0 then
                rho := rhos[i] * rhos[i];;
                if rho[1][2] = 0 and rho[2][1] = 0 and rho[1][1] = rho[2][2] then
                    Print("  T_", i, " is traceless, T_", i, "^2 = ",
                          rho[1][1], " * I\n");;
                fi;
            fi;
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

# ════════ PSL(2,11) / C3, n = 220 ════════
Print("Computing PSL(2,11)/C3...\n\n");;
G := PSL(2,11);;
Print("|G| = ", Size(G), "\n");;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
Print("|H| = ", Size(H), ",  n = [G:H] = ", Index(G, H), "\n\n");;
im := Image(FactorCosetAction(G, H));;
analyze("PSL(2,11) / C3   (n = 220)", im, Index(G, H), 3);;

Print("════════ DONE ════════\n");;
