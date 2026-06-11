# why13_fano.g — What does PSL(3,2) stabilize in the split Cayley hexagon?
#
# Goal: make "72 Fano planes inside GH(2,2)" concrete.
#   1. Reconstruct the 63 lines of GH(2,2) as the triangles of the
#      valency-6 orbital (collinearity) graph of G = PrimitiveGroup(63,2).
#   2. Find K = PSL(3,2) <= G (order 168, index 72).
#   3. Compute K-orbits on points and on lines.
#   4. For the candidate Fano orbits, print the induced incidence structure
#      and pairwise collinearity-graph distances.
#
# Usage: gap -q why13_fano.g

G := PrimitiveGroup(63, 2);;
Print("G = ", StructureDescription(G), ", |G| = ", Size(G), "\n");;

# ─── Collinearity graph: the valency-6 suborbit ───
stab1 := Stabilizer(G, 1);;
suborbs := Orbits(stab1, [1..63]);;
Print("Suborbit lengths at point 1: ", List(suborbs, Length), "\n");;

six := First(suborbs, o -> Length(o) = 6);;
edges := Orbit(G, Set([1, six[1]]), OnSets);;
Print("Number of edges: ", Length(edges), " (expect 189)\n");;

nbrs := List([1..63], i -> []);;
for e in edges do
    AddSet(nbrs[e[1]], e[2]);;
    AddSet(nbrs[e[2]], e[1]);;
od;

# ─── Lines = triangles of the collinearity graph ───
lines := [];;
for e in edges do
    common := Intersection(nbrs[e[1]], nbrs[e[2]]);;
    if Length(common) <> 1 then
        Print("WARNING: edge ", e, " has ", Length(common), " common nbrs\n");;
    fi;
    for c in common do
        AddSet(lines, Set([e[1], e[2], c]));;
    od;
od;
Print("Number of lines (triangles): ", Length(lines), " (expect 63)\n");;

# sanity: each point on 3 lines
plcount := List([1..63], i -> Number(lines, l -> i in l));;
Print("Lines per point: ", Set(plcount), " (expect [3])\n\n");;

# ─── Distance function in the collinearity graph ───
# BFS from each point
distmat := [];;
for i in [1..63] do
    d := List([1..63], x -> -1);;
    d[i] := 0;;
    frontier := [i];;
    while frontier <> [] do
        newf := [];;
        for x in frontier do
            for y in nbrs[x] do
                if d[y] = -1 then
                    d[y] := d[x] + 1;;
                    Add(newf, y);;
                fi;
            od;
        od;
        frontier := newf;;
    od;
    Add(distmat, d);;
od;
Print("Collinearity graph diameter: ", Maximum(Flat(distmat)), " (expect 3)\n\n");;

# ─── Find K = PSL(3,2) ───
Print("Searching for PSL(3,2) (order 168)...\n");;
subs := ConjugacyClassesSubgroups(G);;
k168 := Filtered(List(subs, Representative), H -> Size(H) = 168);;
Print("Classes of order-168 subgroups: ", Length(k168), "\n");;
K := k168[1];;
Print("K = ", StructureDescription(K), "\n");;
Print("K generators (for reuse):\n", GeneratorsOfGroup(K), "\n\n");;

# ─── K-orbits on points ───
porb := Orbits(K, [1..63]);;
Print("K-orbits on points: ", List(porb, Length), "\n");;

# ─── K-orbits on lines ───
lorb := Orbits(K, lines, OnSets);;
Print("K-orbits on lines:  ", List(lorb, Length), "\n\n");;

# ─── Examine each small orbit ───
for o in porb do
    if Length(o) <= 21 then
        Print("Point orbit of size ", Length(o), ": ", Set(o), "\n");;
        # pairwise distances
        dd := [];;
        for i in [1..Length(o)] do
            for j in [i+1..Length(o)] do
                Add(dd, distmat[o[i]][o[j]]);;
            od;
        od;
        Print("  pairwise distances: ", Collected(dd), "\n");;
    fi;
od;
Print("\n");;
for o in lorb do
    if Length(o) <= 21 then
        Print("Line orbit of size ", Length(o), ":\n");;
        for l in o do Print("  ", l, "\n"); od;
    fi;
od;

# ─── Candidate Fano incidence: 7-point orbit vs 7-line orbit? ───
p7 := Filtered(porb, o -> Length(o) = 7);;
l7 := Filtered(lorb, o -> Length(o) = 7);;
if Length(p7) > 0 and Length(l7) > 0 then
    P := Set(p7[1]);;
    L := List(l7[1], l -> Intersection(l, P));;
    Print("\nIncidence of 7-point orbit with 7-line orbit:\n");;
    Print("  point set: ", P, "\n");;
    for l in l7[1] do
        Print("  line ", l, "  meets P in ", Intersection(l, P), "\n");;
    od;
fi;

# ─── Alternative: Fano structure on a 21-orbit? Heawood-style? ───
# Also report, for each point orbit O and line orbit M, the incidence counts:
Print("\nIncidence counts between point orbits and line orbits:\n");;
for o in porb do
    for m in lorb do
        # how many points of o on a line of m (constant by transitivity? check)
        counts := Set(List(m, l -> Length(Intersection(l, Set(o)))));;
        Print("  |", Length(o), "-pt orbit ∩ line of ", Length(m), "-orbit| ∈ ", counts, "\n");;
    od;
od;

Print("\n═══ DONE ═══\n");;
