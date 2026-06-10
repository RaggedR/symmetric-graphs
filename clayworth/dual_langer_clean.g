LoadPackage("grape");;

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
stab192 := Stabilizer(G, 1);;
orbs63 := Orbits(stab192, [1..63]);;
orb6 := First(orbs63, o -> Length(o) = 6);;
langer := EdgeOrbitsGraph(G, [[1, orb6[1]]], 63);;

Print("=== RECOVERING GH(2,2) LINES FROM LANGER TRIANGLES ===\n\n");;

# Each line of GH(2,2) has 3 collinear points forming a triangle in Langer.
# Lambda = 1 means each edge lies in exactly 1 triangle.
lines := [];;
for v in [1..63] do
    nbrs := Adjacency(langer, v);;
    for i in [1..Length(nbrs)] do
        for j in [i+1..Length(nbrs)] do
            a := nbrs[i]; b := nbrs[j];;
            if b in Adjacency(langer, a) then
                line := Set([v, a, b]);;
                if not line in lines then
                    Add(lines, line);;
                fi;
            fi;
        od;
    od;
od;
Print("Lines recovered: ", Length(lines), "\n");;
Print("Expected: 63 (one per GH(2,2) line)\n");;

# Verify: each point on exactly 3 lines
for v in [1..63] do
    count := Number(lines, l -> v in l);;
    if count <> 3 then
        Print("ERROR: vertex ", v, " on ", count, " lines\n");;
    fi;
od;
Print("Each point on exactly 3 lines: verified\n\n");;

# === BUILD TUTTE 12-CAGE ===
# Vertices: 1..63 = points, 64..126 = lines
# Edge: point p ~ line l iff p in lines[l-63]
Print("=== BUILDING TUTTE 12-CAGE ===\n");;
tutte := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        if x <= 63 and y > 63 then
            return x in lines[y-63];
        elif x > 63 and y <= 63 then
            return y in lines[x-63];
        else
            return false;
        fi;
    end, true);;
Print("Tutte 12-cage: V=126, deg=", Length(Adjacency(tutte, 1)),
      ", girth=", Girth(tutte),
      ", diam=", Diameter(tutte), "\n");;

# Distance-2 restricted to points = Langer?
langerFromTutte := Graph(Group(()), [1..63], OnPoints,
    function(x, y)
        return Distance(tutte, x, y) = 2;
    end, true);;
Print("Tutte dist-2 on points iso to Langer: ",
      IsIsomorphicGraph(langerFromTutte, langer), "\n");;

# Distance-2 restricted to lines = dual Langer?
dualFromTutte := Graph(Group(()), [1..63], OnPoints,
    function(x, y)
        return Distance(tutte, x+63, y+63) = 2;
    end, true);;

# Build dual Langer from the other maximal subgroup
maxclasses := MaximalSubgroupClassReps(G);;
M4 := First(maxclasses, m -> Size(m) = 192 and not IsConjugate(G, m, stab192));;
act4 := FactorCosetAction(G, M4);;
im4 := Image(act4);;
stab4 := Stabilizer(im4, 1);;
orbs4 := Orbits(stab4, [1..63]);;
orb6_4 := First(orbs4, o -> Length(o) = 6);;
dualLanger := EdgeOrbitsGraph(im4, [[1, orb6_4[1]]], 63);;

Print("Tutte dist-2 on lines iso to dual Langer: ",
      IsIsomorphicGraph(dualFromTutte, dualLanger), "\n\n");;

# === CHARACTERIZE THE DUAL LANGER ===
Print("=== DUAL LANGER CHARACTERIZATION ===\n");;
Print("Same intersection array {6,4,4;1,1,3}\n");;

# Key distinguishing feature: subgraph induced on distance-3 vertices
# (the "opposite" vertices)
d3_langer := DistanceSet(langer, 3, [1]);;
d3_dual := DistanceSet(dualLanger, 3, [1]);;
Print("Langer: |d3(v)| = ", Length(d3_langer), "\n");;
Print("Dual:   |d3(v)| = ", Length(d3_dual), "\n");;

# Induced subgraph on distance-3 set
sub_langer := InducedSubgraph(langer, d3_langer);;
sub_dual := InducedSubgraph(dualLanger, d3_dual);;
Print("Langer d3-subgraph: V=", Length(Vertices(sub_langer)),
      ", edges per vertex=", Length(Adjacency(sub_langer, 1)),
      ", connected=", IsConnectedGraph(sub_langer), "\n");;
Print("Dual d3-subgraph:   V=", Length(Vertices(sub_dual)),
      ", edges per vertex=", Length(Adjacency(sub_dual, 1)),
      ", connected=", IsConnectedGraph(sub_dual), "\n\n");;

# === THE 126-VERTEX GRAPH: CAREFUL ANALYSIS ===
Print("=== 126-VERTEX GRAPH: DISTANCE-REGULARITY CHECK ===\n\n");;

subsOfStab := ConjugacyClassesSubgroups(stab192);;
K96 := First(List(subsOfStab, Representative),
    k -> Size(k) = 96 and StructureDescription(k) = "SL(2,3) : C4");;
act := FactorCosetAction(G, K96);;
im := Image(act);;
stab1 := Stabilizer(im, 1);;
orbs := Orbits(stab1, [1..126]);;
orb6list := Filtered(orbs, o -> Length(o) = 6 and not (1 in o));;
gamma126 := EdgeOrbitsGraph(im, [[1, orb6list[2][1]]], 126);;
diam := Diameter(gamma126);;

Print("V=126, deg=6, girth=", Girth(gamma126), ", diam=", diam, "\n");;

# Thorough distance-regularity: check ALL intersection numbers
Print("Checking full intersection number constancy...\n");;
# For each distance d, compute a_d, b_d, c_d from multiple vertices
isStrictlyDR := true;;
for d in [0..diam] do
    bvals := []; cvals := []; avals := [];;
    # Sample 20 vertex pairs at distance d
    count := 0;;
    for v in [1..126] do
        for w in DistanceSet(gamma126, d, [v]) do
            if count >= 20 then break; fi;
            nbrs_w := Adjacency(gamma126, w);;
            if d > 0 then
                cd := Number(nbrs_w, u -> Distance(gamma126, v, u) = d-1);;
                AddSet(cvals, cd);;
            fi;
            ad := Number(nbrs_w, u -> Distance(gamma126, v, u) = d);;
            AddSet(avals, ad);;
            if d < diam then
                bd := Number(nbrs_w, u -> Distance(gamma126, v, u) = d+1);;
                AddSet(bvals, bd);;
            fi;
            count := count + 1;;
        od;
        if count >= 20 then break; fi;
    od;
    Print("  d=", d, ": a=", avals);;
    if d < diam then Print(", b=", bvals); fi;
    if d > 0 then Print(", c=", cvals); fi;
    Print("\n");;
    if Length(avals) > 1 or Length(bvals) > 1 or Length(cvals) > 1 then
        isStrictlyDR := false;;
    fi;
od;
Print("Strictly distance-regular: ", isStrictlyDR, "\n\n");;

# Is it the bipartite double cover of Langer?
Print("Checking bipartiteness of 126-vertex graph...\n");;
# A bipartite double cover would be bipartite by definition
# Check: does a 2-coloring exist?
isBip := true;;
color := List([1..126], i -> -1);;
color[1] := 0;;
queue := [1];;
while Length(queue) > 0 do
    v := Remove(queue, 1);;
    for w in Adjacency(gamma126, v) do
        if color[w] = -1 then
            color[w] := 1 - color[v];;
            Add(queue, w);;
        elif color[w] = color[v] then
            isBip := false;;
        fi;
    od;
od;
Print("Bipartite: ", isBip, "\n");;
if isBip then
    part0 := Filtered([1..126], v -> color[v] = 0);;
    part1 := Filtered([1..126], v -> color[v] = 1);;
    Print("Parts: ", Length(part0), " + ", Length(part1), "\n");;
fi;

# If not bipartite, check odd girth
if not isBip then
    Print("(Not bipartite — girth ", Girth(gamma126), " confirms odd cycles)\n\n");;
fi;

# Is it the canonical double cover (tensor product with K2)?
# CDC of Langer would have adjacency: (i,0)~(j,1) iff i~j in Langer
# This is always bipartite. Since our graph has girth 6 (even) it COULD be bipartite.
# But if it's not bipartite, it's NOT the CDC.

# Let's also check: is this the "extended bipartite double" or some other construction?
# Key: |Aut| = 24192 = 2 * |G_2(2)|. The extra C2 must be meaningful.
Print("The 126-vertex graph has |Aut| = 24192 = 2 * |G_2(2)|\n");;
Print("The extra C2 in Aut = C2 x G_2(2) likely relates to the\n");;
Print("point-line duality of GH(2,2).\n\n");;

# Check: is the 126-vertex graph isomorphic to distance-3 of Tutte 12-cage?
# Distance 3 in bipartite Tutte connects points to lines
Print("=== DISTANCE-ODD GRAPHS OF TUTTE 12-CAGE ===\n");;
tutteD3 := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        return Distance(tutte, x, y) = 3;
    end, true);;
Print("Tutte dist-3: deg=", Length(Adjacency(tutteD3, 1)),
      ", connected=", IsConnectedGraph(tutteD3), "\n");;

tutteD1or3 := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        local d;
        d := Distance(tutte, x, y);
        return d = 1 or d = 3;
    end, true);;
Print("Tutte dist-{1,3}: deg=", Length(Adjacency(tutteD1or3, 1)),
      ", connected=", IsConnectedGraph(tutteD1or3), "\n");;

tutteD5 := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        return Distance(tutte, x, y) = 5;
    end, true);;
Print("Tutte dist-5: deg=", Length(Adjacency(tutteD5, 1)),
      ", connected=", IsConnectedGraph(tutteD5), "\n");;

# The halved graphs of Tutte 12-cage
# Halved = distance-2 restricted to one bipartite class
# We already checked: halved on points = Langer, halved on lines = dual Langer

# What about distance-{1,3,5} (all odd distances)?
tutteOdd := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        local d;
        d := Distance(tutte, x, y);
        return d in [1, 3, 5];
    end, true);;
Print("Tutte odd-dist: deg=", Length(Adjacency(tutteOdd, 1)),
      ", connected=", IsConnectedGraph(tutteOdd), "\n");;

# Check isomorphism with gamma126
Print("\ngamma126 iso to Tutte dist-3: ",
      IsIsomorphicGraph(gamma126, tutteD3), "\n");;

# What about the complement-in-bipartite?
# Or distance-{2,6} in Tutte?
tutteD2 := Graph(Group(()), [1..126], OnPoints,
    function(x, y)
        return Distance(tutte, x, y) = 2;
    end, true);;
Print("Tutte dist-2 (full): deg=", Length(Adjacency(tutteD2, 1)),
      ", connected=", IsConnectedGraph(tutteD2), "\n");;
Print("gamma126 iso to Tutte dist-2: ",
      IsIsomorphicGraph(gamma126, tutteD2), "\n");;

Print("\n=== DONE ===\n");;
QuitGap(0);
