# Compute Klein graph quotients
LoadPackage("GRAPE");;

G := PSL(2,7);;
c3_classes := Filtered(ConjugacyClassesSubgroups(G),
    c -> Size(Representative(c)) = 3);;
H := Representative(c3_classes[1]);;
cosets := RightCosets(G, H);;
act := ActionHomomorphism(G, cosets, OnRight);;
P := Image(act);;

# Build the Klein graph (cubic orbital)
stab := Stabilizer(P, 1);;
orbs := Orbits(stab, [2..56]);;
size3 := Filtered(orbs, o -> Length(o) = 3);;

for oi in [1..Length(size3)] do
    gamma := EdgeOrbitsGraph(P, List(size3[oi], x -> [1, x]));;
    if IsSimpleGraph(gamma) and ForAll([1..56], function(v)
        return Length(Adjacency(gamma, v)) = 3; end) then
        Print("Klein graph found from suborbit ", oi, "\n");;

        # Get all block systems
        blocks_list := AllBlocks(P);;
        Print("Block systems:\n");;
        for b in blocks_list do
            bsize := Length(b);;
            nblocks := 56 / bsize;;
            Print("  Block size ", bsize, " -> ", nblocks, " blocks\n");;

            # Compute the full partition
            bs := Orbit(P, Set(b), OnSets);;

            # Build quotient
            blockMap := List([1..56], function(v)
                local i;
                for i in [1..Length(bs)] do
                    if v in bs[i] then return i; fi;
                od;
            end);;

            qn := Length(bs);;
            qadj := List([1..qn], i -> []);;
            for i in [1..qn] do
                for j in [i+1..qn] do
                    found := false;;
                    for u in bs[i] do
                        if found then break; fi;
                        for v in bs[j] do
                            if v in Adjacency(gamma, u) then
                                found := true; break;
                            fi;
                        od;
                    od;
                    if found then
                        Add(qadj[i], j); Add(qadj[j], i);
                    fi;
                od;
            od;

            qdeg := Set(List([1..qn], i -> Length(qadj[i])));;
            qedges := Sum(List([1..qn], i -> Length(qadj[i]))) / 2;;
            Print("    Quotient: ", qn, "v, deg ", qdeg, ", ", qedges, " edges\n");;

            # Print edge list
            Print("    Edges: ");;
            for i in [1..qn] do
                for j in qadj[i] do
                    if i < j then Print("(", i-1, ",", j-1, ") "); fi;
                od;
            od;
            Print("\n");;
        od;
        break;  # Only need one cubic graph
    fi;
od;
