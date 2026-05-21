# Klein quartic via (2,3,7) triangle group
# R (order 7) = face rotation, S (order 3) = vertex rotation, T=RS (order 2) = edge flip
# Vertices = cosets of <S>, neighbors via T, rotation via S
LoadPackage("grape");;

G := PSL(2, 7);;
Print("# |PSL(2,7)| = ", Size(G), "\n");;

# Find R (order 7) and S (order 3) with (RS)^2 = 1
found := false;;
for R in Elements(G) do
    if Order(R) = 7 then
        for S in Elements(G) do
            if Order(S) = 3 and Order(R*S) = 2 then
                found := true;
                break;
            fi;
        od;
        if found then break; fi;
    fi;
od;;

if not found then
    Print("ERROR: No (2,3,7) generators found!\n");
fi;;

T := R * S;;
Print("# R order=", Order(R), " S order=", Order(S), " T=RS order=", Order(T), "\n");;
Print("# <R,S> = G? ", Size(Group(R, S)) = Size(G), "\n");;

H := Group(S);;
cosets := RightCosets(G, H);;
n := Length(cosets);;
Print("# n = ", n, " vertices\n");;

# Build coset lookup: element -> coset index (1-based)
coset_idx := NewDictionary(Identity(G), true);;
for i in [1..n] do
    for h in Elements(H) do
        AddDictionary(coset_idx, h * Representative(cosets[i]), i);
    od;
od;;

# Neighbors of coset c_i: c_i * T gives the edge-flip neighbor
# Under the S-orbit: neighbors are c_i*T, c_i*S*T, c_i*S^2*T
# These are the 3 neighbors in rotation order

# Build edges and rotation system
edges := [];;
edgeset := [];;
rotation := List([1..n], i -> []);;

for i in [1..n] do
    g := Representative(cosets[i]);;
    rot := [];;
    for k in [0..2] do
        nb_elem := T * S^k * g;;
        j := LookupDictionary(coset_idx, nb_elem);;
        if j = fail then
            # Try all coset reps
            for ci in [1..n] do
                if nb_elem in cosets[ci] then j := ci; break; fi;
            od;
        fi;
        Add(rot, j);;
        # Add edge (undirected)
        e := Set([i, j]);;
        if not e in edgeset then
            Add(edgeset, e);
            Add(edges, [i-1, j-1]);  # 0-indexed
        fi;
    od;
    rotation[i] := List(rot, x -> x-1);;  # 0-indexed
od;;

Print("# ", n, " vertices, ", Length(edges), " edges\n");;

# Output rotation
Print("ROTATION_START\n");;
for i in [1..n] do
    Print(rotation[i], "\n");;
od;;
Print("ROTATION_END\n");;

# Output edges
Print("EDGES_START\n");;
for e in edges do
    Print(e[1], " ", e[2], "\n");;
od;;
Print("EDGES_END\n");;

QUIT;
