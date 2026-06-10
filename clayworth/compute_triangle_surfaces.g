# Triangle group surface construction
# Given G with presentation <R, S | R^r = S^3 = (RS)^2 = 1>,
# compute the regular representation on darts (= group elements)
# and coset projection maps for vertices, edges, faces.
#
# Output: JSON-like data for generate_triangle_lean.py

LoadPackage("grape");;

# ============================================================
# Choose the group: PSL(2,7) for Klein, or G2(2) for Clayworth
# ============================================================

mode := "clayworth";;  # "klein" for Phase 1, "clayworth" for Phase 2

if mode = "klein" then
    G := PSL(2, 7);;
    r := 7;;
    name := "Klein";;
elif mode = "clayworth" then
    # G2(2) as normalizer of PSU(3,3) in Sym(63), then coset action on 4032 pts
    psu33 := PrimitiveGroup(63, 1);;
    g2raw := Normalizer(SymmetricGroup(63), psu33);;
    # Find C3 with 0 fixed points on 63 points
    c3cls := Filtered(ConjugacyClassesSubgroups(g2raw),
        c -> Size(Representative(c)) = 3);;
    cw3 := First(c3cls, cl ->
        Number([1..63], p -> p^GeneratorsOfGroup(Representative(cl))[1] = p) = 0);;
    cwsub := Representative(cw3);;
    # Coset action gives G2(2) on 4032 points (but we need the abstract group)
    G := g2raw;;
    r := 12;;
    name := "Clayworth";;
fi;;

nDarts := Size(G);;
Print("# ", name, ": |G| = ", nDarts, ", r = ", r, "\n");;

# Find R (order r) and S (order 3) with (RS)^2 = 1 and <R,S> = G
# For large groups, iterate conjugacy class representatives instead of all elements
found := false;;
cc := ConjugacyClasses(G);;
eltsR := [];;
eltsS := [];;
for c in cc do
    rep := Representative(c);;
    if Order(rep) = r then Append(eltsR, Elements(c)); fi;
    if Order(rep) = 3 then Append(eltsS, Elements(c)); fi;
od;;
Print("# Candidates: ", Length(eltsR), " of order ", r, ", ", Length(eltsS), " of order 3\n");;

for R in eltsR do
    for S in eltsS do
        if Order(R*S) = 2 then
            if Size(Group(R, S)) = nDarts then
                found := true;
                break;
            fi;
        fi;
    od;
    if found then break; fi;
od;;

if not found then
    Error("No (2,3,", r, ") generators found!");
fi;;

T := R * S;;
Print("# R order=", Order(R), " S order=", Order(S), " T=RS order=", Order(T), "\n");;
Print("# <R,S> = G? ", Size(Group(R, S)) = nDarts, "\n");;

# Enumerate group elements: darts = [1..nDarts] <-> Elements(G)
elts := Elements(G);;
Assert(0, Length(elts) = nDarts);;

# Build element -> dart index lookup (1-based)
elt_idx := NewDictionary(Identity(G), true);;
for i in [1..nDarts] do
    AddDictionary(elt_idx, elts[i], i);
od;;

# R, S, T as permutations of [1..nDarts] (right multiplication)
# dart g -> dart g*R
R_perm := List([1..nDarts], i -> LookupDictionary(elt_idx, elts[i] * R));;
S_perm := List([1..nDarts], i -> LookupDictionary(elt_idx, elts[i] * S));;
T_perm := List([1..nDarts], i -> LookupDictionary(elt_idx, elts[i] * T));;

# Verify T = RS on darts
for i in [1..nDarts] do
    Assert(0, T_perm[i] = S_perm[R_perm[i]]);
od;;
Print("# Verified: T = RS on all darts\n");;

# Inverse permutations
R_inv := List([1..nDarts], i -> 0);;
S_inv := List([1..nDarts], i -> 0);;
for i in [1..nDarts] do
    R_inv[R_perm[i]] := i;
    S_inv[S_perm[i]] := i;
od;;

# Coset maps: vertexOf (S-cosets), edgeOf (T-cosets), faceOf (R-cosets)
# Must use LEFT cosets g<H> so right multiplication by H preserves them:
#   g·S ∈ g<S> (vertex coset preserved), g·T ∈ g<T> (edge), g·R ∈ g<R> (face)
H_S := Group(S);;
H_T := Group(T);;
H_R := Group(R);;

# Compute left cosets manually: g<H> for representative g
# Use right transversal to get coset representatives
computeLeftCosets := function(G, H)
    local reps, cosets, g, found, c;
    reps := [];
    cosets := [];
    for g in Elements(G) do
        found := false;
        for c in cosets do
            if g in c then found := true; break; fi;
        od;
        if not found then
            Add(reps, g);
            Add(cosets, Set(List(Elements(H), h -> g * h)));
        fi;
    od;
    return rec(reps := reps, cosets := cosets);
end;;

lc_S := computeLeftCosets(G, H_S);;
lc_T := computeLeftCosets(G, H_T);;
lc_R := computeLeftCosets(G, H_R);;

nV := Length(lc_S.cosets);;
nE := Length(lc_T.cosets);;
nF := Length(lc_R.cosets);;

Print("# V=", nV, " E=", nE, " F=", nF, "\n");;
Print("# chi = ", nV - nE + nF, " = 2 - 2*", (2 - (nV - nE + nF))/2, "\n");;

# Build coset index lookups
vertex_of := List([1..nDarts], i -> 0);;
edge_of := List([1..nDarts], i -> 0);;
face_of := List([1..nDarts], i -> 0);;

for ci in [1..nV] do
    for g in lc_S.cosets[ci] do
        idx := LookupDictionary(elt_idx, g);;
        vertex_of[idx] := ci;
    od;
od;;

for ci in [1..nE] do
    for g in lc_T.cosets[ci] do
        idx := LookupDictionary(elt_idx, g);;
        edge_of[idx] := ci;
    od;
od;;

for ci in [1..nF] do
    for g in lc_R.cosets[ci] do
        idx := LookupDictionary(elt_idx, g);;
        face_of[idx] := ci;
    od;
od;;

# Verify all assigned
Assert(0, ForAll(vertex_of, x -> x > 0));;
Assert(0, ForAll(edge_of, x -> x > 0));;
Assert(0, ForAll(face_of, x -> x > 0));;

# Verify coset axioms: S preserves vertex cosets, T preserves edge cosets, R preserves face cosets
for i in [1..nDarts] do
    Assert(0, vertex_of[S_perm[i]] = vertex_of[i]);
    Assert(0, edge_of[T_perm[i]] = edge_of[i]);
    Assert(0, face_of[R_perm[i]] = face_of[i]);
od;;
Print("# Verified: coset axioms hold\n");;

# Coset representatives: smallest-index dart in each coset
vertex_rep := List([1..nV], v -> Minimum(Filtered([1..nDarts], i -> vertex_of[i] = v)));;
edge_rep := List([1..nE], e -> Minimum(Filtered([1..nDarts], i -> edge_of[i] = e)));;
face_rep := List([1..nF], f -> Minimum(Filtered([1..nDarts], i -> face_of[i] = f)));;

# Verify round-trip
for v in [1..nV] do Assert(0, vertex_of[vertex_rep[v]] = v); od;;
for e in [1..nE] do Assert(0, edge_of[edge_rep[e]] = e); od;;
for f in [1..nF] do Assert(0, face_of[face_rep[f]] = f); od;;

# Verify edge_ne: for each edge, src and tgt vertices differ
for e in [1..nE] do
    d := edge_rep[e];;
    Assert(0, vertex_of[d] <> vertex_of[T_perm[d]]);
od;;
Print("# Verified: edge_ne (all edges connect distinct vertices)\n");;

# Verify face_inj: for each face, the r darts hit r distinct edges
for f in [1..nF] do
    d := face_rep[f];;
    edge_set := [];;
    dd := d;;
    for k in [1..r] do
        e := edge_of[dd];;
        Assert(0, not e in edge_set);
        Add(edge_set, e);
        dd := R_perm[dd];
    od;
    Assert(0, dd = d);  # R^r returns to start
od;;
Print("# Verified: face_inj (each face has ", r, " distinct edges)\n");;

# Left action of R and S on vertices (cosets of <S>)
# L_h(v) = index of the coset containing h * rep(v)
R_vertex := List([1..nV], v -> vertex_of[
    LookupDictionary(elt_idx, R * lc_S.reps[v])]);;
S_vertex := List([1..nV], v -> vertex_of[
    LookupDictionary(elt_idx, S * lc_S.reps[v])]);;

# Verify they're permutations
Assert(0, Set(R_vertex) = [1..nV]);;
Assert(0, Set(S_vertex) = [1..nV]);;

# Inverse permutations
Ri := R^(-1);; Si := S^(-1);;
R_vertex_inv := List([1..nV], v -> vertex_of[
    LookupDictionary(elt_idx, Ri * lc_S.reps[v])]);;
S_vertex_inv := List([1..nV], v -> vertex_of[
    LookupDictionary(elt_idx, Si * lc_S.reps[v])]);;

# Verify inverses
for v in [1..nV] do
    Assert(0, R_vertex[R_vertex_inv[v]] = v);
    Assert(0, S_vertex[S_vertex_inv[v]] = v);
od;;
Print("# Verified: vertex-action generators are permutations with correct inverses\n");;

# BFS witness words: for each vertex v, a word in {R, R^-1, S, S^-1}
# mapping vertex 1 (= 0 in Lean) to vertex v
vertex_gens := [R_vertex, R_vertex_inv, S_vertex, S_vertex_inv];;
bfs_words := List([1..nV], v -> fail);;
bfs_words[1] := [];;
bfs_queue := [1];;
while Length(bfs_queue) > 0 do
    v := Remove(bfs_queue, 1);;
    for gi in [1..4] do
        w := vertex_gens[gi][v];;
        if bfs_words[w] = fail then
            bfs_words[w] := Concatenation(bfs_words[v], [gi-1]);
            Add(bfs_queue, w);
        fi;
    od;
od;;
Assert(0, ForAll(bfs_words, w -> w <> fail));;
max_word_len := Maximum(List(bfs_words, Length));;
Print("# BFS: all ", nV, " vertices reached, max word length = ", max_word_len, "\n");;

# Left multiplication on DARTS by R and S
# L_R(d) = index of R * elts[d], L_S(d) = index of S * elts[d]
LR_perm := List([1..nDarts], i -> LookupDictionary(elt_idx, R * elts[i]));;
LS_perm := List([1..nDarts], i -> LookupDictionary(elt_idx, S * elts[i]));;

Assert(0, Set(LR_perm) = [1..nDarts]);;
Assert(0, Set(LS_perm) = [1..nDarts]);;

# Inverse permutations
LR_inv := List([1..nDarts], i -> 0);;
LS_inv := List([1..nDarts], i -> 0);;
for i in [1..nDarts] do
    LR_inv[LR_perm[i]] := i;
    LS_inv[LS_perm[i]] := i;
od;;

# Verify: left mult commutes with right mult (and hence with T)
for i in [1..nDarts] do
    Assert(0, LR_perm[R_perm[i]] = R_perm[LR_perm[i]]);
    Assert(0, LR_perm[S_perm[i]] = S_perm[LR_perm[i]]);
    Assert(0, LS_perm[R_perm[i]] = R_perm[LS_perm[i]]);
    Assert(0, LS_perm[S_perm[i]] = S_perm[LS_perm[i]]);
od;;

# Verify vertex induction: vertexOf(L_R(d)) = R_vertex[vertexOf(d)]
for i in [1..nDarts] do
    Assert(0, vertex_of[LR_perm[i]] = R_vertex[vertex_of[i]]);
    Assert(0, vertex_of[LS_perm[i]] = S_vertex[vertex_of[i]]);
od;;
Print("# Verified: dart left-mult commutes with right-mult and induces vertex action\n");;

# Output everything (0-indexed for Lean)
Print("\n");;
Print("DATA_START\n");;
Print("name=", name, "\n");;
Print("nDarts=", nDarts, "\n");;
Print("nV=", nV, "\n");;
Print("nE=", nE, "\n");;
Print("nF=", nF, "\n");;
Print("r=", r, "\n");;

Print("R_fwd=", List(R_perm, x -> x-1), "\n");;
Print("R_inv=", List(R_inv, x -> x-1), "\n");;
Print("S_fwd=", List(S_perm, x -> x-1), "\n");;
Print("S_inv=", List(S_inv, x -> x-1), "\n");;
Print("vertex_of=", List(vertex_of, x -> x-1), "\n");;
Print("edge_of=", List(edge_of, x -> x-1), "\n");;
Print("face_of=", List(face_of, x -> x-1), "\n");;
Print("vertex_rep=", List(vertex_rep, x -> x-1), "\n");;
Print("edge_rep=", List(edge_rep, x -> x-1), "\n");;
Print("face_rep=", List(face_rep, x -> x-1), "\n");;
Print("R_vertex_fwd=", List(R_vertex, x -> x-1), "\n");;
Print("R_vertex_inv=", List(R_vertex_inv, x -> x-1), "\n");;
Print("S_vertex_fwd=", List(S_vertex, x -> x-1), "\n");;
Print("S_vertex_inv=", List(S_vertex_inv, x -> x-1), "\n");;
Print("bfs_words=", bfs_words, "\n");;
Print("LR_fwd=", List(LR_perm, x -> x-1), "\n");;
Print("LR_inv=", List(LR_inv, x -> x-1), "\n");;
Print("LS_fwd=", List(LS_perm, x -> x-1), "\n");;
Print("LS_inv=", List(LS_inv, x -> x-1), "\n");;

Print("DATA_END\n");;

QUIT;
