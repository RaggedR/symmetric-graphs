LoadPackage("grape");;

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
stab192 := Stabilizer(G, 1);;
orbs63 := Orbits(stab192, [1..63]);;
orb6 := First(orbs63, o -> Length(o) = 6);;
langer := EdgeOrbitsGraph(G, [[1, orb6[1]]], 63);;

# Build dual Langer
maxclasses := MaximalSubgroupClassReps(G);;
M4 := First(maxclasses, m -> Size(m) = 192 and not IsConjugate(G, m, stab192));;
act4 := FactorCosetAction(G, M4);;
im4 := Image(act4);;
stab4 := Stabilizer(im4, 1);;
orbs4 := Orbits(stab4, [1..63]);;
orb6_4 := First(orbs4, o -> Length(o) = 6);;
dualLanger := EdgeOrbitsGraph(im4, [[1, orb6_4[1]]], 63);;

# Build 126-vertex graph
subsOfStab := ConjugacyClassesSubgroups(stab192);;
K96 := First(List(subsOfStab, Representative),
    k -> Size(k) = 96 and StructureDescription(k) = "SL(2,3) : C4");;
act := FactorCosetAction(G, K96);;
im := Image(act);;
stab1 := Stabilizer(im, 1);;
orbs := Orbits(stab1, [1..126]);;
orb6list := Filtered(orbs, o -> Length(o) = 6 and not (1 in o));;
gamma126 := EdgeOrbitsGraph(im, [[1, orb6list[2][1]]], 126);;

Print("gamma126: V=126, bipartite, 6-regular, girth 6, diam 4\n\n");;

# Since it's bipartite with parts of size 63,
# the halved graphs (induced on each part) tell us what it is.
# Build the bipartition
color := List([1..126], i -> -1);;
color[1] := 0;;
queue := [1];;
while Length(queue) > 0 do
    v := Remove(queue, 1);;
    for w in Adjacency(gamma126, v) do
        if color[w] = -1 then
            color[w] := 1 - color[v];;
            Add(queue, w);;
        fi;
    od;
od;
part0 := Filtered([1..126], v -> color[v] = 0);;
part1 := Filtered([1..126], v -> color[v] = 1);;
Print("Part 0: ", Length(part0), " vertices\n");;
Print("Part 1: ", Length(part1), " vertices\n\n");;

# Build the halved graphs (distance-2 within each part)
half0 := Graph(Group(()), part0, OnPoints,
    function(x, y) return Distance(gamma126, x, y) = 2; end, true);;
half1 := Graph(Group(()), part1, OnPoints,
    function(x, y) return Distance(gamma126, x, y) = 2; end, true);;

Print("Half-0: deg=", Length(Adjacency(half0, part0[1])),
      ", connected=", IsConnectedGraph(half0), "\n");;
Print("Half-1: deg=", Length(Adjacency(half1, part1[1])),
      ", connected=", IsConnectedGraph(half1), "\n\n");;

# Are the halved graphs isomorphic to Langer and/or dual Langer?
Print("Half-0 iso to Langer: ", IsIsomorphicGraph(half0, langer), "\n");;
Print("Half-0 iso to dual Langer: ", IsIsomorphicGraph(half0, dualLanger), "\n");;
Print("Half-1 iso to Langer: ", IsIsomorphicGraph(half1, langer), "\n");;
Print("Half-1 iso to dual Langer: ", IsIsomorphicGraph(half1, dualLanger), "\n\n");;

# KEY INSIGHT: If this is a bipartite graph where one part = points,
# other part = lines, and adjacency = some relation between points and lines
# in GH(2,2), then it might be a DISTANCE graph of the Tutte 12-cage.

# Recover lines from Langer triangles
lines := [];;
for v in [1..63] do
    nbrs := Adjacency(langer, v);;
    for i in [1..Length(nbrs)] do
        for j in [i+1..Length(nbrs)] do
            if nbrs[j] in Adjacency(langer, nbrs[i]) then
                line := Set([v, nbrs[i], nbrs[j]]);;
                if not line in lines then Add(lines, line); fi;
            fi;
        od;
    od;
od;

# Build Tutte 12-cage
tutte := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        if x <= 63 and y > 63 then return x in lines[y-63];
        elif x > 63 and y <= 63 then return y in lines[x-63];
        else return false; fi;
    end, true);;

# What about: build a bipartite graph where point p ~ line l iff
# they are NOT incident (anti-flag graph)?
Print("=== ANTI-FLAG / NON-INCIDENCE GRAPH ===\n");;
antiflag := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        if x <= 63 and y > 63 then return not (x in lines[y-63]);
        elif x > 63 and y <= 63 then return not (y in lines[x-63]);
        else return false; fi;
    end, true);;
Print("Anti-flag: deg=", Length(Adjacency(antiflag, 1)),
      ", girth=", Girth(antiflag), "\n");;
# Anti-flag would be bipartite with deg = 63-3 = 60. Too high.

# What about: point p ~ line l iff they are at distance 3 in Tutte?
# (distance 1 = incident, distance 3 = "almost incident")
Print("\n=== TUTTE DISTANCE-3 POINT-LINE PAIRS ===\n");;
tutteD3pl := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        if x <= 63 and y > 63 then return Distance(tutte, x, y) = 3;
        elif x > 63 and y <= 63 then return Distance(tutte, x, y) = 3;
        else return false; fi;
    end, true);;
Print("Tutte d3 (point-line): deg=", Length(Adjacency(tutteD3pl, 1)),
      ", girth=", Girth(tutteD3pl),
      ", connected=", IsConnectedGraph(tutteD3pl), "\n");;
Print("Iso to gamma126: ", IsIsomorphicGraph(tutteD3pl, gamma126), "\n\n");;

# Distance 5 in Tutte (point-line)
Print("=== TUTTE DISTANCE-5 POINT-LINE PAIRS ===\n");;
tutteD5pl := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        if x <= 63 and y > 63 then return Distance(tutte, x, y) = 5;
        elif x > 63 and y <= 63 then return Distance(tutte, x, y) = 5;
        else return false; fi;
    end, true);;
Print("Tutte d5 (point-line): deg=", Length(Adjacency(tutteD5pl, 1)),
      ", connected=", IsConnectedGraph(tutteD5pl), "\n");;
Print("Iso to gamma126: ", IsIsomorphicGraph(tutteD5pl, gamma126), "\n\n");;

# What if it's the "local Langer" graph?
# p ~ l iff the (unique) point of l closest to p is a Langer-neighbor of p?
# This is distance-3 in Tutte restricted to point-line pairs.
# Already checked above.

# Let's try a completely different approach: just identify the bipartite
# adjacency matrix. For each point v in part0, what lines (part1) is it
# connected to in gamma126?
Print("=== DIRECT INSPECTION ===\n");;
# Map part0 and part1 vertices to 1..63
map0 := List([1..126], v -> Position(part0, v));;
map1 := List([1..126], v -> Position(part1, v));;

# For vertex 1 (which is in part0), find its neighbors (all in part1)
v := part0[1];;
nbrs_v := Adjacency(gamma126, v);;
Print("Vertex ", v, " (part0[1]) neighbors: ", nbrs_v, "\n");;
Print("These map to part1 positions: ",
      List(nbrs_v, w -> Position(part1, w)), "\n");;

# Now check: in the block quotient gamma126 → Langer,
# vertex v maps to block b(v). What are the blocks?
allBl := AllBlocks(im);;
block2 := First(allBl, bs -> Length(bs) = 2 and 1 in bs);;
allBlocks := Orbit(im, Set(block2), OnSets);;
Print("\nBlocks (first 5): ", allBlocks{[1..5]}, "\n");;

# Check: do the blocks mix parts or stay within parts?
blockMixes := 0;;
for bl in allBlocks do
    cols := Set(List(bl, v -> color[v]));;
    if Length(cols) > 1 then blockMixes := blockMixes + 1; fi;
od;
Print("Blocks that mix bipartite parts: ", blockMixes, " / ", Length(allBlocks), "\n");;

# So: either blocks stay within parts (graph is 2 copies of Langer),
# or blocks cross parts (each block = {point, line}).
# If blocks cross: the quotient graph has edges between blocks,
# and each block pairs a point with a line.

if blockMixes = Length(allBlocks) then
    Print("ALL blocks cross parts — each block = {point in part0, vertex in part1}\n");;
    Print("This means gamma126 pairs points and lines via the C2 of H_192.\n\n");;

    # What pairing? For each block, the two vertices are in different parts.
    # This gives a bijection part0 ↔ part1, hence points ↔ lines.
    Print("The C2 element of H_192/K_96 sends each GH(2,2) point to a line.\n");;
    Print("gamma126 connects point p to line l iff they are 'associated'\n");;
    Print("under this pairing AND adjacent in the block quotient (Langer).\n");;
elif blockMixes = 0 then
    Print("Blocks stay within parts — this is a genuine bipartite graph.\n");;
    Print("Each part quotients independently to 63/2... wait, 63 is odd.\n");;
    Print("Blocks must cross parts since 63 is odd and blocks have size 2.\n");;
fi;

Print("\n=== DONE ===\n");;
QuitGap(0);
