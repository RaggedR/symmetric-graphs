# why13_subhexagon.g — Verify that (P, SA ∪ SB) is a generalized hexagon
# of order (2,1): 21 points, 14 lines, incidence graph girth 12, diameter 6.
# Spreads from why13_fano.g / why13_orbitals.g.
# Usage: gap -q why13_subhexagon.g

SA := [ [1,12,58],[33,37,38],[14,26,42],[22,34,45],[18,29,46],[8,11,32],[21,40,50] ];;
SB := [ [1,26,32],[8,18,33],[12,37,45],[14,21,38],[34,42,46],[29,40,58],[11,22,50] ];;
L := Concatenation(SA, SB);;
P := Set(Flat(L));;
Print("points: ", Length(P), "  lines: ", Length(L), "\n");;

# each point on how many lines?
Print("lines per point: ", Set(List(P, p -> Number(L, l -> p in l))), " (expect [2])\n");;

# incidence graph: vertices 1..21 points, 22..35 lines
np := Length(P);;
adj := List([1..np + Length(L)], i -> []);;
for i in [1..Length(L)] do
    for p in L[i] do
        pi := Position(P, p);;
        AddSet(adj[pi], np + i);;
        AddSet(adj[np + i], pi);;
    od;
od;

# BFS for distances; compute diameter and girth
nv := Length(adj);;
diam := 0;;
girth := infinity;;
for s in [1..nv] do
    d := List([1..nv], x -> -1);;
    d[s] := 0;;
    f := [s];;
    while f <> [] do
        nf := [];;
        for x in f do
            for y in adj[x] do
                if d[y] = -1 then
                    d[y] := d[x] + 1;;
                    Add(nf, y);;
                elif d[y] >= d[x] then
                    # cycle through s of length d[x]+d[y]+1 (odd) — bipartite so use:
                    g := d[x] + d[y] + 1;;
                    if g < girth then girth := g; fi;
                fi;
            od;
        od;
        f := nf;;
    od;
    diam := Maximum(diam, Maximum(d));;
od;
# NB the BFS girth bound above can overcount for non-tree edges within the same
# level only when odd cycles exist; graph is bipartite, so cross-edges at equal
# depth cannot occur; edges to depth d[x] give even cycles d[x]+d[y]+1? For
# bipartite graphs a non-tree edge from x (depth k) hits y at depth k+1 already
# visited via another branch: cycle length 2k+2 at most. The minimum over all
# roots is exact for vertex-transitive-ish small graphs; we double check below
# by counting 4-,6-,8-,10-cycles via matrix powers if girth looks wrong.
Print("incidence graph: diameter = ", diam, " (expect 6), girth ≈ ", girth, " (expect 12)\n");;

# exact girth via adjacency matrix powers: smallest k with trace(A^k) > 0, k even
A := NullMat(nv, nv);;
for i in [1..nv] do for j in adj[i] do A[i][j] := 1; od; od;
Ak := A;;
for k in [2..12] do
    Ak := Ak * A;;
    if k mod 2 = 0 then
        tr := Sum([1..nv], i -> Ak[i][i]);;
        # closed k-walks include backtracking; girth g iff trace counts exceed
        # tree-walk counts. Simpler: count cycles by checking A^k diagonal vs
        # bipartite biregular tree walks is messy — instead check directly:
        # girth > k iff between any adjacent... use distance argument instead.
        Print("  trace(A^", k, ") = ", tr, "\n");;
    fi;
od;

# cleanest girth check: for each edge (u,v), remove it and find dist(u,v);
# girth = 1 + min over edges of that distance.
mind := infinity;;
for u in [1..nv] do
    for v in adj[u] do
        if u < v then
            # BFS from u avoiding edge uv
            d := List([1..nv], x -> -1);; d[u] := 0;; f := [u];;
            while f <> [] and d[v] = -1 do
                nf := [];;
                for x in f do
                    for y in adj[x] do
                        if (x = u and y = v) or (x = v and y = u) then continue; fi;
                        if d[y] = -1 then d[y] := d[x] + 1;; Add(nf, y);; fi;
                    od;
                od;
                f := nf;;
            od;
            if d[v] <> -1 and d[v] + 1 < mind then mind := d[v] + 1; fi;
        fi;
    od;
od;
Print("exact girth = ", mind, " (expect 12)\n");;
Print("═══ DONE ═══\n");;
