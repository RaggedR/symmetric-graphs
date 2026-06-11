# hecke_torsion.g — Torsion elements and irrationalities in H(PSL(2,q), C3)
#
# For q = 7, 11, 13: extract all non-commutative blocks (up to mult 5),
# catalog characteristic polynomials, identify finite-order elements and
# the algebraic number fields they generate.
#
# Usage: gap -q hecke_torsion.g 2>&1 | tee hecke_torsion.out

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

# Check if a matrix has finite order: M^k = I for some small k
finiteorder := function(M, maxord)
    local k, P, m;
    m := Length(M);;
    P := IdentityMat(m);;
    for k in [1..maxord] do
        P := P * M;;
        if P = IdentityMat(m) then return k; fi;
    od;
    return 0;  # not finite order up to maxord
end;;

analyzetorsion := function(name, im, n, mmax)
    local tbl, irr, ccl, csums, C, M, g, x, mats, r, lens, pivals, pi,
          idx, chi, m, P, W, vals, a, cp, facs, f, lam, E, I2, v, V,
          rho, rhos, evecs, mu, i, j, c2, S, w, cpset, ord, cpoly;
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
    pivals := List(ccl, C -> Number([1..n], x -> x^Representative(C) = x));;
    pi := ClassFunction(tbl, pivals);;

    Print("rank ", r, "\n\n");;

    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m = 0 or m = 1 then continue; fi;
        if m > mmax then
            Print("chi_", idx, " (deg ", chi[1], "): mult ", m,
                  " -- skipped (> ", mmax, ")\n\n");;
            continue;;
        fi;
        Print(">>> ", m, "-dim irrep on chi_", idx, " (deg ", chi[1],
              ")\n");;
        P := (chi[1] / Size(im)) * Sum([1..Length(ccl)],
                 i -> ComplexConjugate(chi[i]) * csums[i]);;
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
            Print("  !! no module found\n\n");; continue;;
        fi;
        # Diagonalize if possible
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
        fi;
        # Symmetric normalization for m = 2
        if m = 2 then
            for i in [2..r] do
                rho := List(V, w -> SolutionMat(V, w * mats[i]));;
                if rho[1][2] <> 0 and rho[2][1] <> 0 then
                    c2 := rho[1][2] / rho[2][1];;
                    if IsRat(c2) and c2 > 0 then
                        V[2] := Sqrt(c2) * V[2];;
                    fi;
                    break;;
                fi;
            od;
        fi;
        rhos := List(mats, M -> List(V, w -> SolutionMat(V, w * M)));;

        # Collect distinct charpolys and analyze
        Print("  Distinct characteristic polynomials:\n");;
        cpset := [];;
        for i in [1..r] do
            cpoly := CharacteristicPolynomial(rhos[i]);;
            if not cpoly in List(cpset, p -> p[1]) then
                Add(cpset, [cpoly, [i]]);;
            else
                First(cpset, p -> p[1] = cpoly)[2] :=
                    Concatenation(First(cpset, p -> p[1] = cpoly)[2], [i]);;
            fi;
        od;
        for p in cpset do
            Print("    ", p[1], "  (T_", p[2], ")\n");;
        od;

        # Check for finite-order elements
        Print("\n  Torsion elements (order <= 30):\n");;
        for i in [1..r] do
            ord := finiteorder(rhos[i], 30);;
            if ord > 0 and ord > 1 then
                Print("    T_", i, ": order ", ord,
                      "  (charpoly ", CharacteristicPolynomial(rhos[i]), ")\n");;
            fi;
        od;

        # B^2 = N*I for traceless elements (m=2 only)
        if m = 2 then
            Print("\n  B^2 = N*I for traceless generators:\n");;
            for i in [2..r] do
                if TraceMat(rhos[i]) = 0 then
                    rho := rhos[i] * rhos[i];;
                    if rho[1][2] = 0 and rho[2][1] = 0
                       and rho[1][1] = rho[2][2] then
                        Print("    T_", i, ": B^2 = ", rho[1][1],
                              " * I\n");;
                    fi;
                fi;
            od;
        fi;

        # Check
        S := [];;
        for i in [1..r] do
            w := Flat(rhos[i]);;
            if Length(S) = 0 or SolutionMat(S, w) = fail then Add(S, w); fi;
            for j in [1..r] do
                w := Flat(rhos[i] * rhos[j]);;
                if SolutionMat(S, w) = fail then Add(S, w);; fi;
            od;
        od;
        Print("  dim check: ", Length(S), " (expect ", m^2, ")\n\n");;
    od;
end;;

# ════════ PSL(2,7)/C3 ════════
Print("=== PSL(2,7)/C3 ===\n\n");;
G := PSL(2,7);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
im := Image(FactorCosetAction(G, H));;
analyzetorsion("PSL(2,7)/C3  (n=56)", im, Index(G,H), 3);;

# ════════ PSL(2,11)/C3 ════════
Print("=== PSL(2,11)/C3 ===\n\n");;
G := PSL(2,11);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
im := Image(FactorCosetAction(G, H));;
analyzetorsion("PSL(2,11)/C3  (n=220)", im, Index(G,H), 3);;

# ════════ PSL(2,13)/C3 ════════
Print("=== PSL(2,13)/C3 ===\n\n");;
G := PSL(2,13);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
im := Image(FactorCosetAction(G, H));;
analyzetorsion("PSL(2,13)/C3  (n=364)", im, Index(G,H), 5);;

Print("════════ DONE ════════\n");;
