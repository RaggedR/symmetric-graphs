# hecke_psl2_survey.g — Wedderburn decomposition of H(PSL(2,q), C3) for primes q
#
# Just the multiplicity data: which blocks appear as q varies?
# Usage: gap -q hecke_psl2_survey.g 2>&1 | tee hecke_psl2_survey.out

for q in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
    G := PSL(2, q);;
    cc := ConjugacyClasses(G);;
    c3rep := First(cc, c -> Order(Representative(c)) = 3);;
    if c3rep = fail then
        Print("q = ", q, ": no element of order 3\n\n");;
        continue;;
    fi;
    H := Group(Representative(c3rep));;
    n := Index(G, H);;
    im := Image(FactorCosetAction(G, H));;
    tbl := CharacterTable(im);;
    irr := Irr(tbl);;
    ccl := ConjugacyClasses(tbl);;
    pivals := List(ccl, C -> Number([1..n], x -> x^Representative(C) = x));;
    pi := ClassFunction(tbl, pivals);;

    mults := [];;
    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m > 0 then Add(mults, [chi[1], m]); fi;
    od;

    Print("q = ", q, "  |G| = ", Size(G), "  n = ", n,
          "  rank = ", Sum(mults, p -> p[2]^2), "\n");;
    Print("  Constituents (deg, mult): ", mults, "\n");;
    Print("  Wedderburn: ");;
    has2 := false;;
    multset := List(mults, p -> p[2]);;
    for p in mults do
        if p[2] = 1 then Print("C ");
        else Print("M_", p[2], " "); fi;
        if p[2] = 2 then has2 := true; fi;
    od;
    Print("\n");;
    Print("  Mult set: ", Set(multset), "\n");;
    if has2 then Print("  *** HAS 2-DIM BLOCK ***\n"); fi;
    Print("\n");;
od;

Print("════════ DONE ════════\n");;
