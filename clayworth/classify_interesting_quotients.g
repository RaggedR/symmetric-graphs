# =============================================================================
# classify_interesting_quotients.g
#
# Identify the most interesting quotients of the Clayworth graph.
# Focus on: sparse (low degree/vertex ratio), high girth, full Aut = G₂(2).
# =============================================================================

LoadPackage("grape");;

# --- Setup ---
psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
cwC3 := First(
    Filtered(ConjugacyClassesSubgroups(G), c -> Size(Representative(c)) = 3),
    function(cl)
        return Number([1..63], function(p)
            return p^GeneratorsOfGroup(Representative(cl))[1] = p;
        end) = 0;
    end
);;
cwC3 := Representative(cwC3);;

cwAct := FactorCosetAction(G, cwC3);;
cwPerm := Image(cwAct);;
n := 4032;;

# Find connected self-paired orbital
stab_cw1 := Stabilizer(cwPerm, 1);;
orbsCW := Orbits(stab_cw1, [2..n]);;
gamma := fail;;
for orb in Filtered(orbsCW, function(o) return Length(o) = 3; end) do
    v := orb[1];;
    stabv := Stabilizer(cwPerm, v);;
    orbOf1 := First(Orbits(stabv, [1..n]), function(o) return 1 in o; end);;
    if orbOf1 = fail or Length(orbOf1) <> 3 then continue; fi;
    testG := EdgeOrbitsGraph(cwPerm, List(orb, function(v) return [1,v]; end), n);;
    if IsConnectedGraph(testG) then gamma := testG; break; fi;
od;
Print("Clayworth: ", n, "v, deg ", Length(Adjacency(gamma, 1)), "\n\n");

# --- Classify all quotients ---
ints := IntermediateSubgroups(G, cwC3);;

# Collect results
results := [];;
for ki in [1..Length(ints.subgroups)] do
    K := ints.subgroups[ki];;
    ksize := Size(K);;
    nblocks := 12096 / ksize;;

    Kimage := Image(cwAct, K);;
    block := Set(Orbit(Kimage, 1));;
    blocks := Orbit(cwPerm, block, OnSets);;
    blockMap := ListWithIdenticalEntries(n, 0);;
    for bi in [1..nblocks] do
        for v in blocks[bi] do blockMap[v] := bi; od;
    od;

    # Build quotient adjacency
    qAdj := List([1..nblocks], function(i) return []; end);;
    for v in [1..n] do
        bv := blockMap[v];
        for nbr in Adjacency(gamma, v) do
            bnbr := blockMap[nbr];
            if bv <> bnbr then AddSet(qAdj[bv], bnbr); fi;
        od;
    od;
    qDeg := Length(qAdj[1]);;

    # Build GRAPE quotient
    qGamma := Graph(Group(()), [1..nblocks],
        function(x,g) return x; end,
        function(x,y) return y in qAdj[x]; end, true);;
    comps := ConnectedComponents(qGamma);;
    nComp := Length(comps);;
    compSize := Length(comps[1]);;

    # Girth (skip very large)
    girth := -1;;
    if nblocks <= 1500 then girth := Girth(qGamma); fi;

    # Aut (skip very large)
    autSize := -1;;
    if nblocks <= 700 then
        autSize := Size(AutGroupGraph(qGamma));
    fi;

    # Edge multiplicity
    multi := 0;;
    if qDeg > 0 then
        target := qAdj[1][1];;
        for v in blocks[1] do
            for nbr in Adjacency(gamma, v) do
                if blockMap[nbr] = target then multi := multi + 1; fi;
            od;
        od;
    fi;

    # Density = deg / (V-1)
    density := 0;;
    if nblocks > 1 then density := Float(qDeg) / Float(nblocks - 1); fi;

    Add(results, rec(
        ki := ki, ksize := ksize, kstruct := StructureDescription(K),
        nblocks := nblocks, blockSize := Length(block),
        qDeg := qDeg, nComp := nComp, compSize := compSize,
        girth := girth, autSize := autSize, multi := multi,
        density := density
    ));
od;

# --- Sort by interest: low density, high girth, |Aut| = 12096 ---
Print("=======================================================================\n");
Print("QUOTIENT CLASSIFICATION (sorted by sparsity)\n");
Print("=======================================================================\n\n");

# Filter interesting: connected, not complete, deg > 0
interesting := Filtered(results, function(r)
    return r.nComp = 1 and r.qDeg > 0 and r.qDeg < r.nblocks - 1;
end);;

# Sort by density (sparsest first)
SortBy(interesting, function(r) return r.density; end);;

Print(String("V", -7), String("Deg", -6), String("Girth", -7),
      String("Mult", -6), String("Dens%", -8),
      String("|Aut|", -15), String("BlkSz", -7),
      String("K structure", -35), "\n");
Print(ListWithIdenticalEntries(91, '-'), "\n");

for r in interesting do
    densStr := String(Int(r.density * 1000) / 10);;

    autStr := "";;
    if r.autSize = 12096 then autStr := "G₂(2)";
    elif r.autSize = 24192 then autStr := "2×G₂(2)";
    elif r.autSize > 0 then autStr := String(r.autSize);
    else autStr := "?";
    fi;

    girthStr := "";;
    if r.girth > 0 then girthStr := String(r.girth);
    else girthStr := "?";
    fi;

    Print(String(r.nblocks, -7), String(r.qDeg, -6),
          String(girthStr, -7), String(r.multi, -6),
          String(densStr, -8), String(autStr, -15),
          String(r.blockSize, -7), r.kstruct, "\n");
od;

# --- Highlight the best ---
Print("\n=======================================================================\n");
Print("HIGHLIGHTS\n");
Print("=======================================================================\n\n");

# Sparsest with |Aut| = G₂(2)
sparse_g2 := Filtered(interesting, function(r)
    return r.autSize = 12096 and r.density < 0.15;
end);;
SortBy(sparse_g2, function(r) return r.density; end);;

Print("Sparsest quotients with Aut = G₂(2):\n\n");
for r in sparse_g2 do
    Print("  ", r.nblocks, "v, ", r.qDeg, "-regular, girth ",
          r.girth, ", density ", Int(r.density*1000)/10, "%\n");
    Print("    K = ", r.kstruct, " (order ", r.ksize,
          ", blocks of ", r.blockSize, ")\n");
    Print("    Edge multiplicity: ", r.multi, "\n\n");
od;

# Highest girth
highGirth := Filtered(interesting, function(r) return r.girth >= 6; end);;
SortBy(highGirth, function(r) return -r.girth; end);;

Print("Highest girth quotients:\n\n");
for r in highGirth{[1..Minimum(8, Length(highGirth))]} do
    Print("  ", r.nblocks, "v, ", r.qDeg, "-regular, girth ",
          r.girth, ", |Aut| = ");
    if r.autSize = 12096 then Print("G₂(2)");
    elif r.autSize = 24192 then Print("2×G₂(2)");
    elif r.autSize > 0 then Print(r.autSize);
    else Print("?");
    fi;
    Print("\n    K = ", r.kstruct, "\n\n");
od;

# Degenerate (complete or edgeless)
degenerate := Filtered(results, function(r)
    return r.qDeg = 0 or r.qDeg = r.nblocks - 1 or r.nComp > 1;
end);;
Print("Degenerate quotients: ", Length(degenerate), "\n");
for r in degenerate do
    if r.qDeg = 0 then
        Print("  ", r.nblocks, " isolated vertices (K = ", r.kstruct, ")\n");
    elif r.qDeg = r.nblocks - 1 then
        Print("  K_", r.nblocks, " (K = ", r.kstruct, ")\n");
    else
        Print("  ", r.nComp, " × ", r.compSize, "v (K = ", r.kstruct, ")\n");
    fi;
od;

QuitGap(0);
