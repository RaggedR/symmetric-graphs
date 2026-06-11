# hecke_spherical.g — spherical function tables (1-dim Hecke irreps) for the
# Gelfand graphs/pairs of CATALOGUE.md.  For a Gelfand pair every irrep of
# H(G,H) is 1-dimensional: omega_chi(T_i) = tr(A_i P_chi)/chi(1).
#
# Usage: gap -q hecke_spherical.g > hecke_spherical.out 2>&1

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

spherical := function(name, im, n)
    local tbl, irr, ccl, csums, C, M, g, x, st, orbs, prs, o, L, lens,
          pivals, pi, idx, chi, m, P, vals, i, p;
    Print("############ ", name, " ############\n");;
    st := Stabilizer(im, 1);;
    orbs := ShallowCopy(Orbits(st, [1..n]));;
    SortBy(orbs, Length);;
    prs := [];;
    for o in orbs do Add(prs, Orbit(im, [1, o[1]], OnPairs)); od;
    lens := List(orbs, Length);;
    Print("rank ", Length(orbs), ",  suborbit lengths ", lens, "\n");;
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
    pivals := List(ccl, C -> Number([1..n], x -> x^Representative(C) = x));;
    pi := ClassFunction(tbl, pivals);;
    for idx in [1..Length(irr)] do
        chi := irr[idx];;
        m := ScalarProduct(pi, chi);;
        if m = 0 then continue; fi;
        if m > 1 then
            Print("  chi_", idx, " (deg ", chi[1], "): mult ", m,
                  " — NOT GELFAND HERE, see hecke_irreps*.g\n");;
            continue;;
        fi;
        P := (chi[1] / Size(im)) * Sum([1..Length(ccl)],
                 i -> ComplexConjugate(chi[i]) * csums[i]);;
        vals := List(prs, L -> Sum(L, p -> P[p[2]][p[1]]) / chi[1]);;
        Print("  omega on chi_", idx, " (deg ", chi[1], "): ", vals, "\n");;
    od;
    Print("\n");;
end;;

graphspherical := function(name, n, edges)
    local dir, gamma, A, o, im;
    dir := Concatenation(edges, List(edges, e -> [e[2], e[1]]));;
    gamma := EdgeOrbitsGraph(Group(()), dir, n);;
    A := AutGroupGraph(gamma);;
    if IsTransitive(A, [1..n]) then
        spherical(Concatenation(name, "  [Aut, vertex stab]"), A, n);;
    else
        for o in Orbits(A, [1..n]) do
            im := Image(FactorCosetAction(A, Stabilizer(A, o[1])));;
            spherical(Concatenation(name, "  [orbit size ",
                      String(Length(o)), "]"), im, Length(o));;
        od;
    fi;
end;;

# ── Foster graphs ≤ 30 (F26A excluded: non-Gelfand, done in hecke_irreps.g) ──
graphspherical("K4", 4, [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);;
graphspherical("K33", 6, [[1,4],[1,5],[1,6],[2,4],[2,5],[2,6],[3,4],[3,5],[3,6]]);;
graphspherical("Cube", 8, lcfedges([3,-3], 4));;
pairs := Combinations([1..5], 2);;
pedges := [];;
for i in [1..10] do for j in [i+1..10] do
    if Intersection(pairs[i], pairs[j]) = [] then Add(pedges, [i,j]); fi;
od; od;
graphspherical("Petersen", 10, pedges);;
graphspherical("Heawood", 14, lcfedges([5,-5], 7));;
graphspherical("MoebiusKantor", 16, lcfedges([5,-5], 8));;
graphspherical("Pappus", 18, lcfedges([5,7,-7,7,-7,-5], 3));;
graphspherical("Desargues", 20, lcfedges([5,-5,9,-9], 5));;
graphspherical("Dodecahedron", 20, lcfedges([10,7,4,-4,-7,10,-4,7,-7,4], 2));;
graphspherical("Nauru", 24, lcfedges([5,-9,7,-7,9,-5], 4));;
graphspherical("TutteCoxeter", 30, lcfedges([-13,-9,7,-7,9,13], 5));;

# ── Tutte 12-cage: Levi graph of GH(2,2), two orbits of 63 ──
G63 := PrimitiveGroup(63, 2);;
stab1 := Stabilizer(G63, 1);;
six := First(Orbits(stab1, [1..63]), o -> Length(o) = 6);;
hexedges := Orbit(G63, Set([1, six[1]]), OnSets);;
nbrs := List([1..63], i -> []);;
for e in hexedges do
    AddSet(nbrs[e[1]], e[2]);; AddSet(nbrs[e[2]], e[1]);;
od;
hexlines := [];;
for e in hexedges do
    for c in Intersection(nbrs[e[1]], nbrs[e[2]]) do
        AddSet(hexlines, Set([e[1], e[2], c]));;
    od;
od;
levi := [];;
for i in [1..Length(hexlines)] do
    for p in hexlines[i] do Add(levi, [p, 63 + i]); od;
od;
graphspherical("Tutte12Cage (Levi of GH(2,2))", 126, levi);;

# ── coset pairs not represented by a graph run above ──
G := PSL(2,5);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(3)) do
    im := Image(FactorCosetAction(G, Image(emb)));;
    spherical("PSL(2,5) / S3   (Petersen as coset graph)", im, 10);;
od;
G := PGL(2,7);;
cls := List(ConjugacyClassesSubgroups(G), Representative);;
cls := Filtered(cls, s -> Size(s) = 12);;
k := 0;;
for H in cls do
    k := k + 1;;
    im := Image(FactorCosetAction(G, H));;
    spherical(Concatenation("PGL(2,7) / ", StructureDescription(H),
              " #", String(k), "   (28, Coxeter candidate)"), im, 28);;
od;
G := PSL(2,17);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(4)) do
    im := Image(FactorCosetAction(G, Image(emb)));;
    spherical("PSL(2,17) / S4   (Biggs-Smith, n = 102)", im, 102);;
od;

# ── geometry point actions ──
G := PSL(3,3);;
im := Image(FactorCosetAction(G, Stabilizer(G, 1)));;
spherical("PSL(3,3) / P1   (13 points PG(2,3))", im, 13);;
G := PSL(3,4);;
im := Image(FactorCosetAction(G, Stabilizer(G, 1)));;
spherical("PSL(3,4) / P1   (21 points PG(2,4))", im, 21);;
G := PSp(4,3);;
im := Image(FactorCosetAction(G, Stabilizer(G, 1)));;
spherical("PSp(4,3) / P1   (40 points GQ(3,3))", im, 40);;
G := PrimitiveGroup(63, 2);;
spherical("G2(2) / P1   (63 points GH(2,2))", G, 63);;
G := PSU(3,3);;
im := Image(FactorCosetAction(G, Stabilizer(G, 1)));;
spherical("PSU(3,3) / H28   (Hermitian unital)", im, 28);;
cand := AllPrimitiveGroups(NrMovedPoints, 65, Size, 29120);;
if Length(cand) > 0 then
    spherical("Sz(8) / H65   (Suzuki ovoid)", cand[1], 65);;
fi;

Print("════════ DONE ════════\n");;
