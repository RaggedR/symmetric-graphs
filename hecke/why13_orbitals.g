# why13_orbitals.g — Geometric meaning of the 8 orbitals of H(G₂(2), PSL(3,2))
#
# Follow-up to why13_fano.g. Computes:
#   1. N_G(K) and whether the normalizing involution swaps the two spreads.
#   2. The 72-point action on cosets of K; the 8 orbitals with valencies.
#   3. Symmetry of each orbital (which D_i have D_i⁻¹ = D_i).
#   4. Geometric invariants of each orbital: for (K, gK), compare the
#      21-point sets P, P' and spreads — |P ∩ P'|, spread overlaps.
#   5. The graph Γ for each orbital: bipartite? spectrum (char poly factors)?
#   6. Re-derive T_a·T_a for the 7-valent orbitals by path counting,
#      to match against the canonical multiplication table.
#
# Usage: gap -q why13_orbitals.g

G := PrimitiveGroup(63, 2);;

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

Print("|K| = ", Size(K), "\n");;

# ─── 1. Normalizer ───
N := Normalizer(G, K);;
Print("|N_G(K)| = ", Size(N), "  (index in G: ", Size(G)/Size(N), ")\n");;

# the two spreads (from why13_fano.g)
SA := [ [1,12,58],[33,37,38],[14,26,42],[22,34,45],[18,29,46],[8,11,32],[21,40,50] ];;
SB := [ [1,26,32],[8,18,33],[12,37,45],[14,21,38],[34,42,46],[29,40,58],[11,22,50] ];;
SA := Set(List(SA, Set));;
SB := Set(List(SB, Set));;

if Size(N) > Size(K) then
    t := First(Elements(N), x -> not x in K);;
    SAimg := Set(List(SA, l -> OnSets(l, t)));;
    Print("Normalizing element maps spread A to: ");;
    if SAimg = SB then Print("spread B  (SWAPS the spreads)\n");
    elif SAimg = SA then Print("spread A  (fixes each spread setwise)\n");
    else Print("something else?!\n");
    fi;
fi;

# ─── 2. Action on 72 cosets ───
cosets := RightCosets(G, K);;
n := Length(cosets);;
Print("\n|G/K| = ", n, "\n");;
act := ActionHomomorphism(G, cosets, OnRight);;
G72 := Image(act);;
stab := Stabilizer(G72, 1);;
sub72 := Orbits(stab, [1..n]);;
sub72 := List(sub72, Set);;
SortBy(sub72, o -> [Length(o), Minimum(o)]);;
Print("Suborbit lengths on 72 cosets: ", List(sub72, Length), "\n\n");;

# ─── 3. Orbital symmetry + adjacency ───
# Orbital of suborbit O = G-orbit of pairs (1, x), x ∈ O.
# Build adjacency lists; check symmetry.
reps := List(cosets, Representative);;
# For pair (i,j): double coset of reps[i]*reps[j]^-1 ... we instead use the
# group-theoretic pairing: orbital containing (1,x) is paired with the one
# containing (x,1); transport (x,1) -> (1, x') by g with 1^g... use G72:
PairedSuborbit := function(x)
    local g, y;
    # find g in G72 with x^g = 1; then the pair (x,1) maps to (1, 1^g)
    g := RepresentativeAction(G72, x, 1);;
    return 1^g;;
end;;

Print("Suborbit pairing (D vs D⁻¹):\n");;
pairing := [];;
for o in sub72 do
    x := o[1];;
    y := PairedSuborbit(x);;
    # which suborbit contains y?
    j := PositionProperty(sub72, s -> y in s);;
    i := PositionProperty(sub72, s -> x in s);;
    Add(pairing, [i, j]);;
    if i = j then
        Print("  suborbit ", i, " (val ", Length(o), "): SYMMETRIC\n");;
    else
        Print("  suborbit ", i, " (val ", Length(o), "): paired with suborbit ", j, "\n");;
    fi;
od;

# ─── 4. Geometric invariants per suborbit ───
# Each coset Kg ↔ conjugate Fano structure: P^g, SA^g, SB^g where g = rep.
# Coset i corresponds to K * reps[i]; the structure is the image under reps[i]⁻¹?
# Careful: coset Kg acts: the stabilizer of coset Kg in G is g⁻¹Kg.
# The Fano structure of coset Kg is (P, SA, SB) transported by g⁻¹... we want
# the structure stabilized by g⁻¹Kg, which is the image of K's structure under g⁻¹.
P0 := Union(SA);;
Print("\nGeometric invariants (relative to coset 1 = K itself):\n");;
Print("suborbit | val | |P∩P'| | |SA∩SA'| | |SA∩SB'| | |SB∩SA'| | |SB∩SB'|\n");;
for idx in [1..Length(sub72)] do
    o := sub72[idx];;
    x := o[1];;
    g := reps[x];;   # coset K*reps[x]
    gi := g^-1;;
    P1 := OnSets(P0, gi);;
    SA1 := Set(List(SA, l -> OnSets(l, gi)));;
    SB1 := Set(List(SB, l -> OnSets(l, gi)));;
    Print("   ", idx, "     |  ", Length(o), "  |   ",
          Length(Intersection(P0, P1)), "   |    ",
          Length(Intersection(SA, SA1)), "    |    ",
          Length(Intersection(SA, SB1)), "    |    ",
          Length(Intersection(SB, SA1)), "    |    ",
          Length(Intersection(SB, SB1)), "\n");;
od;

# Also: is coset x's point set P' EQUAL to P0 for the valency-1 suborbit?
# (the chirality partner should have the same 21 points, spreads swapped)

# ─── 5. Orbital graphs: path counts T_a·T_a and bipartiteness of Γ₃ ───
# adjacency lists for each suborbit
adjlist := List(sub72, o -> List([1..n], i -> []));;
for idx in [1..Length(sub72)] do
    o := sub72[idx];;
    for x in o do
        orb := Orbit(G72, [1, x], OnPairs);;
        for p in orb do
            AddSet(adjlist[idx][p[1]], p[2]);;
        od;
    od;
od;
# NB: this loops over all x in o but the orbit of one pair already covers the
# whole orbital; harmless duplication via AddSet.

# Path counting: (T_a²)_{1,y} for each suborbit a: count z with 1→z→y both in orbital a
Print("\nPath counts: for each suborbit a, T_a·T_a evaluated at representative pairs:\n");;
for a in [1..Length(sub72)] do
    if Length(sub72[a]) = 7 then
        counts := [];;
        for y in [1..n] do
            c := Number(adjlist[a][1], z -> y in adjlist[a][z]);;
            if c > 0 then
                # which suborbit is y in (y=1 -> identity)
                if y = 1 then j := 0;; else j := PositionProperty(sub72, s -> y in s);; fi;
                Add(counts, [j, c]);;
            fi;
        od;
        Print("  T_", a, "²  hits suborbits (suborbit, coeff): ", Collected(counts), "\n");;
    fi;
od;

# Bipartiteness of each 7-valent orbital graph (underlying undirected)
Print("\nBipartiteness check (BFS 2-colouring of underlying undirected graph):\n");;
for a in [1..Length(sub72)] do
    if Length(sub72[a]) = 7 then
        und := List([1..n], i -> Set(adjlist[a][i]));;
        for i in [1..n] do for j in und[i] do AddSet(und[j], i); od; od;
        col := List([1..n], i -> -1);;
        col[1] := 0;;
        frontier := [1];;
        ok := true;;
        while frontier <> [] and ok do
            newf := [];;
            for x in frontier do
                for y in und[x] do
                    if col[y] = -1 then
                        col[y] := 1 - col[x];;
                        Add(newf, y);;
                    elif col[y] = col[x] then
                        ok := false;;
                    fi;
                od;
            od;
            frontier := newf;;
        od;
        Print("  suborbit ", a, ": bipartite = ", ok, "\n");;
    fi;
od;

# ─── 6. Spectrum of Γ₃-candidates: char poly of adjacency matrices ───
Print("\nCharacteristic polynomials (factored) of 7-valent orbital graphs:\n");;
x_ := Indeterminate(Rationals, "x");;
for a in [1..Length(sub72)] do
    if Length(sub72[a]) = 7 then
        M := NullMat(n, n);;
        for i in [1..n] do for j in adjlist[a][i] do M[i][j] := 1; od; od;
        cp := CharacteristicPolynomial(M);;
        Print("  suborbit ", a, ": ", Factors(cp), "\n");;
    fi;
od;

# And the 21-valent ones
Print("\nCharacteristic polynomials of 21-valent orbital graphs:\n");;
for a in [1..Length(sub72)] do
    if Length(sub72[a]) = 21 then
        M := NullMat(n, n);;
        for i in [1..n] do for j in adjlist[a][i] do M[i][j] := 1; od; od;
        cp := CharacteristicPolynomial(M);;
        Print("  suborbit ", a, ": ", Factors(cp), "\n");;
    fi;
od;

Print("\n═══ DONE ═══\n");;
