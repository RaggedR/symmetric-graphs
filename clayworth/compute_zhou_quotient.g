# Compute the Z₂ quotient of the Zhou graph
LoadPackage("GRAPE");;

G := PSL(2,13);;
subs := ConjugacyClassesSubgroups(G);;
s3_classes := Filtered(subs, c -> Size(Representative(c)) = 6 and
    (not IsAbelian(Representative(c))));;

H := Representative(s3_classes[1]);;
cosets := RightCosets(G, H);;
act := ActionHomomorphism(G, cosets, OnRight);;
P := Image(act);;

# Get the block system
blocks := AllBlocks(P);;
Print("Block: ", blocks[1], "\n");;
Print("Block size: ", Length(blocks[1]), "\n");;

# Compute full block system from one block
bs := Orbit(P, Set(blocks[1]), OnSets);;
Print("Number of blocks: ", Length(bs), "\n");;

# Build the quotient: map each vertex to its block index
blockMap := List([1..182], function(v)
    local i;
    for i in [1..Length(bs)] do
        if v in bs[i] then return i; fi;
    od;
    return fail;
end);;

# Build the Zhou graph first (from suborbit of size 3)
stab := Stabilizer(P, 1);;
orbs_of_1 := Orbits(stab, [2..182]);;
size3_orbs := Filtered(orbs_of_1, o -> Length(o) = 3);;

# Try each size-3 suborbit to find the cubic one
for oi in [1..Length(size3_orbs)] do
    gamma := EdgeOrbitsGraph(P, List(size3_orbs[oi], x -> [1, x]));;
    if IsSimpleGraph(gamma) then
        if ForAll([1..182], function(v)
            return Length(Adjacency(gamma, v)) = 3; end) then

            Print("\n=== Zhou graph from suborbit ", oi, " ===\n");;

            # Build quotient graph
            qn := Length(bs);;
            qadj := List([1..qn], i -> []);;
            for i in [1..qn] do
                for j in [i+1..qn] do
                    # Check if any vertex in block i is adjacent to any in block j
                    found := false;;
                    for u in bs[i] do
                        if found then break; fi;
                        for v in bs[j] do
                            if v in Adjacency(gamma, u) then
                                found := true;;
                                break;
                            fi;
                        od;
                    od;
                    if found then
                        Add(qadj[i], j);;
                        Add(qadj[j], i);;
                    fi;
                od;
            od;

            # Quotient graph properties
            qdeg := List([1..qn], i -> Length(qadj[i]));;
            Print("Quotient vertices: ", qn, "\n");;
            Print("Quotient degree set: ", Set(qdeg), "\n");;
            qedges := Sum(qdeg) / 2;;
            Print("Quotient edges: ", qedges, "\n");;

            # Check if quotient is vertex-transitive
            # Find the overgroup K containing H with [K:H] = 2
            overgroups := Filtered(
                ConjugacyClassesSubgroups(G),
                c -> Size(Representative(c)) = 12 and IsSubgroup(Representative(c), H));;
            Print("Overgroups of order 12 containing H: ", Length(overgroups), "\n");;

            for og in overgroups do
                K := Representative(og);;
                Print("  Overgroup structure: ", StructureDescription(K), "\n");;
                Print("  Index in PSL(2,13): ", Index(G, K), "\n");;
            od;

            # Also check: all overgroups of H
            allover := IntermediateSubgroups(G, H);;
            Print("Intermediate subgroups H < K < G:\n");;
            for K in allover.subgroups do
                Print("  |K| = ", Size(K), ", [G:K] = ", Index(G,K),
                      ", structure: ", StructureDescription(K), "\n");;
            od;
        fi;
    fi;
od;
