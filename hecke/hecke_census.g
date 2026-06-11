# hecke_census.g — Hecke algebras H(G,H) for every (G,H) pair / graph in
# ~/git/graphviz/docs/CATALOGUE.md.
#
# For each pair: rank, suborbit lengths, Wedderburn decomposition
#   H(G,H) = ⊕_χ M_{m_χ}(C),  m_χ = <ind_H^G 1, χ>,
# commutativity (Gelfand pair iff multiplicity-free).
# For Foster-census graphs: build via LCF/explicit edges in GRAPE,
# G = Aut(graph), H = vertex stabilizer.
#
# Usage: gap -q hecke_census.g > hecke_census.out 2>&1

LoadPackage("grape");;

IRR_SIZE_LIMIT := 100000;;   # skip Irr() for groups bigger than this

report := function(name, G, H)
    local n, pi, irr, mults, wedd, rank, hom, im, st, sub, i;
    n := Index(G, H);;
    Print("=== ", name, " ===\n");;
    Print("  |G| = ", Size(G), ",  n = [G:H] = ", n, "\n");;
    hom := FactorCosetAction(G, H);;
    im := Image(hom);;
    st := Stabilizer(im, 1);;
    sub := ShallowCopy(OrbitLengths(st, [1..n]));;
    Sort(sub);;
    Print("  suborbits (", Length(sub), "): ", sub, "\n");;
    if Size(G) <= IRR_SIZE_LIMIT then
        pi := PermutationCharacter(G, H);;
        irr := Irr(G);;
        mults := List(irr, chi -> ScalarProduct(pi, chi));;
        wedd := [];;
        for i in [1..Length(irr)] do
            if mults[i] > 0 then Add(wedd, [irr[i][1], mults[i]]); fi;
        od;
        rank := Sum(mults, m -> m^2);;
        Print("  rank = ", rank, "   (check vs #suborbits: ",
              Length(sub), ")\n");;
        Print("  Wedderburn (chi-degree, mult): ", Collected(wedd), "\n");;
        Print("  commutative (Gelfand pair): ",
              ForAll(mults, m -> m <= 1), "\n");;
    else
        Print("  rank = ", Length(sub),
              "  [group too large for Irr; suborbits only]\n");;
    fi;
    Print("\n");;
end;;

graphreport := function(name, n, edges, girthexp)
    local dir, gamma, A, H, g;
    dir := Concatenation(edges, List(edges, e -> [e[2], e[1]]));;
    gamma := EdgeOrbitsGraph(Group(()), dir, n);;
    A := AutGroupGraph(gamma);;
    g := Girth(gamma);;
    Print("--- graph ", name, ": |V| = ", n, ", girth = ", g,
          " (expect ", girthexp, "), |Aut| = ", Size(A), " ---\n");;
    if IsTransitive(A, [1..n]) then
        H := Stabilizer(A, 1);;
        report(Concatenation(name, "  [G = Aut, H = vertex stab]"), A, H);;
    else
        Print("  NOT vertex-transitive; orbit lengths: ",
              List(Orbits(A, [1..n]), Length), "\n");;
        for g in Orbits(A, [1..n]) do
            report(Concatenation(name, "  [orbit of size ",
                   String(Length(g)), "]"),
                   A, Stabilizer(A, g[1]));;
        od;
    fi;
end;;

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

# ════════════════ Section A: group pairs (Sabidussi / coset) ════════════════
Print("════════ SECTION A: (G,H) coset pairs ════════\n\n");;

# A1 family
G := PSL(2,5);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
report("PSL(2,5) / C3   (20 = double cover of Petersen)", G, H);;

for emb in IsomorphicSubgroups(G, SymmetricGroup(3)) do
    report("PSL(2,5) / S3   (10 = Petersen)", G, Image(emb));;
od;

G := PSL(2,7);;
H := Group(Representative(First(ConjugacyClasses(G),
        c -> Order(Representative(c)) = 3)));;
report("PSL(2,7) / C3   (56)", G, H);;

G := PGL(2,7);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(4)) do
    report("PGL(2,7) / S4   (14 = Heawood)", G, Image(emb));;
od;
# Coxeter graph candidates: order-12 subgroups (D12)
cls := List(ConjugacyClassesSubgroups(G), Representative);;
cls := Filtered(cls, s -> Size(s) = 12);;
for H in cls do
    report(Concatenation("PGL(2,7) / ", StructureDescription(H),
           "   (28, Coxeter candidate)"), G, H);;
od;

G := PSL(2,13);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(3)) do
    report("PSL(2,13) / S3   (182 = Foster 182)", G, Image(emb));;
od;

G := PSL(2,17);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(4)) do
    report("PSL(2,17) / S4   (102 = Biggs-Smith)", G, Image(emb));;
od;

G := PSL(2,19);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(3)) do
    report("PSL(2,19) / S3   (570)", G, Image(emb));;
od;

G := PSL(2,31);;
for emb in IsomorphicSubgroups(G, SymmetricGroup(4)) do
    report("PSL(2,31) / S4   (620)", G, Image(emb));;
od;

# A2 family (projective planes; 2-transitive, rank 2)
G := PSL(3,3);;  report("PSL(3,3) / P1   (13 points of PG(2,3))",
                        G, Stabilizer(G, 1));;
G := PSL(3,4);;  report("PSL(3,4) / P1   (21 points of PG(2,4))",
                        G, Stabilizer(G, 1));;

# B2 family
G := PSp(4,3);;  report("PSp(4,3) / P1   (40 points of GQ(3,3))",
                        G, Stabilizer(G, 1));;

# G2 family
G := PrimitiveGroup(63, 2);;
report("G2(2) / P1   (63 points of GH(2,2))", G, Stabilizer(G, 1));;

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
report("G2(2) / PSL(3,2)   (72 = Meinhold)", G, K);;

# Future-direction groups (CATALOGUE §9.2)
G := PSU(3,3);;  report("PSU(3,3) / H28   (28 points of Hermitian unital)",
                        G, Stabilizer(G, 1));;
# Sz(8) is returned as a matrix group; use its primitive action on 65 points
cand := AllPrimitiveGroups(NrMovedPoints, 65, Size, 29120);;
if Length(cand) > 0 then
    G := cand[1];;
    report("Sz(8) / H65   (65 points of Suzuki ovoid)", G, Stabilizer(G, 1));;
fi;
cand := AllPrimitiveGroups(NrMovedPoints, 100, Size, 604800);;
if Length(cand) > 0 then
    G := cand[1];;
    report("J2 / U3(3)   (100 points, rank-3 SRG(100,36,14,12))",
           G, Stabilizer(G, 1));;
fi;

# ════════════════ Section B: Foster graphs via GRAPE ════════════════
Print("════════ SECTION B: graphs (G = Aut, H = vertex stabilizer) ════════\n\n");;

graphreport("K4", 4, [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]], 3);;
graphreport("K33", 6, [[1,4],[1,5],[1,6],[2,4],[2,5],[2,6],[3,4],[3,5],[3,6]], 4);;
graphreport("Cube", 8, lcfedges([3,-3], 4), 4);;

# Petersen as Kneser K(5,2)
pairs := Combinations([1..5], 2);;
pedges := [];;
for i in [1..10] do for j in [i+1..10] do
    if Intersection(pairs[i], pairs[j]) = [] then Add(pedges, [i,j]); fi;
od; od;
graphreport("Petersen", 10, pedges, 5);;

graphreport("Heawood", 14, lcfedges([5,-5], 7), 6);;
graphreport("MoebiusKantor", 16, lcfedges([5,-5], 8), 6);;
graphreport("Pappus", 18, lcfedges([5,7,-7,7,-7,-5], 3), 6);;
graphreport("Desargues", 20, lcfedges([5,-5,9,-9], 5), 6);;
graphreport("Dodecahedron", 20, lcfedges([10,7,4,-4,-7,10,-4,7,-7,4], 2), 5);;
graphreport("Nauru", 24, lcfedges([5,-9,7,-7,9,-5], 4), 6);;
graphreport("F26A", 26, lcfedges([-7,7], 13), 6);;
graphreport("TutteCoxeter", 30, lcfedges([-13,-9,7,-7,9,13], 5), 8);;

# Tutte 12-cage = Levi graph of GH(2,2): rebuild lines as in why13_fano.g
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
Print("hexagon rebuilt: ", Length(hexlines), " lines (expect 63)\n");;
levi := [];;
for i in [1..Length(hexlines)] do
    for p in hexlines[i] do Add(levi, [p, 63 + i]); od;
od;
graphreport("Tutte12Cage (Levi of GH(2,2))", 126, levi, 12);;

Print("════════ DONE ════════\n");;
