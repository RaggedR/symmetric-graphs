# =============================================================================
# compute_meinhold_graphs.g
#
# Compute the Meinhold family of quotient graphs of the Clayworth graph.
# Each Meinhold graph is a connected vertex-transitive quotient of
# Sab(G₂(2), C₃) via a specific intermediate subgroup.
#
#   Meinhold-72:   72v,  7-reg, girth 4, via K = PSL(3,2)
#   Meinhold-504:  504v, 6-reg, girth 6, via K = S₄
#   Meinhold-576:  576v, 7-reg, girth 6, via K = C₇:C₃ (= F₂₁)
#   Meinhold-1344: 1344v, 9-reg, girth 6, via K = C₃×C₃
#   Meinhold-2016: 2016v, 6-reg, girth 8, via K = C₆
#
# Named after Nick Meinhold, co-founder of Imagineering Melbourne.
# =============================================================================

LoadPackage("grape");;

Print("=======================================================================\n");
Print("COMPUTING THE MEINHOLD FAMILY\n");
Print("=======================================================================\n\n");

# --- Build G₂(2) and Clayworth ---

psu33 := PrimitiveGroup(63, 1);;
G := Normalizer(SymmetricGroup(63), psu33);;
Assert(0, Size(G) = 12096);

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
chosenOrb := fail;;
for orb in Filtered(orbsCW, function(o) return Length(o) = 3; end) do
    v := orb[1];;
    stabv := Stabilizer(cwPerm, v);;
    orbOf1 := First(Orbits(stabv, [1..n]), function(o) return 1 in o; end);;
    if orbOf1 = fail or Length(orbOf1) <> 3 then continue; fi;
    testG := EdgeOrbitsGraph(cwPerm, List(orb, function(vv) return [1,vv]; end), n);;
    if IsConnectedGraph(testG) then gamma := testG; chosenOrb := orb; break; fi;
od;
Print("Clayworth orbital: ", chosenOrb, "\n");
Print("Clayworth: ", n, "v, ", Length(Adjacency(gamma, 1)), "-regular, ",
      Sum(List([1..n], function(v) return Length(Adjacency(gamma, v)); end))/2,
      " edges\n\n");

# --- Find intermediate subgroups ---
ints := IntermediateSubgroups(G, cwC3);;

# Helper: compute quotient graph and export
ComputeMeinhold := function(name, K, targetDeg, targetGirth)
    local ksize, nblocks, Kimage, block, blocks, blockMap, bi, v, qAdj,
          nbr, bv, bnbr, qDeg, qGamma, girth, qAut, edges, filename,
          comps, diam;

    ksize := Size(K);;
    nblocks := 12096 / ksize;;
    Print("--- ", name, " ---\n");
    Print("  K = ", StructureDescription(K), " (order ", ksize, ")\n");

    Kimage := Image(cwAct, K);;
    block := Set(Orbit(Kimage, 1));;
    blocks := Orbit(cwPerm, block, OnSets);;
    Assert(0, Length(blocks) = nblocks);

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
    Print("  Quotient: ", nblocks, "v, ", qDeg, "-regular\n");

    # Build GRAPE graph
    qGamma := Graph(Group(()), [1..nblocks],
        function(x,g) return x; end,
        function(x,y) return y in qAdj[x]; end, true);;

    comps := ConnectedComponents(qGamma);;
    Assert(0, Length(comps) = 1);

    girth := Girth(qGamma);;
    diam := Diameter(qGamma);;
    qAut := AutGroupGraph(qGamma);;

    Print("  Girth: ", girth, ", Diameter: ", diam, "\n");
    Print("  |Aut| = ", Size(qAut), "\n");
    Print("  Aut = ", StructureDescription(qAut), "\n");

    # Check if we got the target
    if qDeg = targetDeg and girth = targetGirth then
        Print("  ✓ Matches target (", targetDeg, "-reg, girth ", targetGirth, ")\n");
    else
        Print("  ✗ Does NOT match target (wanted ", targetDeg, "-reg girth ",
              targetGirth, ", got ", qDeg, "-reg girth ", girth, ")\n");
        return fail;
    fi;

    # Export edges (0-indexed for Python compatibility)
    edges := [];;
    for v in [1..nblocks] do
        for nbr in qAdj[v] do
            if v < nbr then Add(edges, [v-1, nbr-1]); fi;
        od;
    od;
    Print("  Edges: ", Length(edges), "\n");

    # Save to file
    filename := Concatenation("named_graphs/", name, ".g");;
    PrintTo(filename, "# ", name, "\n");
    AppendTo(filename, "# Quotient of Clayworth graph Sab(G₂(2), C₃)\n");
    AppendTo(filename, "# via K = ", StructureDescription(K),
             " (order ", String(ksize), ")\n");
    AppendTo(filename, "# ", String(nblocks), " vertices, ",
             String(qDeg), "-regular, girth ", String(girth),
             ", diameter ", String(diam), "\n");
    AppendTo(filename, "# |Aut| = ", String(Size(qAut)), "\n\n");
    AppendTo(filename, name, "_n := ", nblocks, ";\n");
    AppendTo(filename, name, "_deg := ", qDeg, ";\n");
    AppendTo(filename, name, "_girth := ", girth, ";\n");
    AppendTo(filename, name, "_edges := ", edges, ";\n");
    AppendTo(filename, name, "_adj := ", qAdj, ";\n");

    # Also save block map (quotient projection from Clayworth)
    AppendTo(filename, "\n# Block map: Clayworth vertex v → Meinhold vertex blockMap[v]\n");
    AppendTo(filename, name, "_blockMap := ", blockMap, ";\n");

    # Save automorphism group generators (on the quotient vertices)
    AppendTo(filename, "\n# Automorphism group generators\n");
    AppendTo(filename, name, "_autGens := ",
             GeneratorsOfGroup(qAut), ";\n");

    Print("  Saved to ", filename, "\n\n");
    return qGamma;
end;;

# --- Find the right subgroup for each Meinhold graph ---

# We need to select the specific subgroup that gives the target degree/girth.
# Different conjugates of the same abstract group can give different quotients.

# Meinhold-72: K = PSL(3,2), order 168, target 7-regular girth 4
Print("=== Meinhold-72 ===\n");
k168list := Filtered(ints.subgroups, function(K) return Size(K) = 168; end);;
Print("PSL(3,2) overgroups: ", Length(k168list), "\n");
for K in k168list do
    result := ComputeMeinhold("meinhold72", K, 7, 4);;
    if result <> fail then break; fi;
od;

# Meinhold-504: K = S₄, order 24, target 6-regular girth 6
Print("=== Meinhold-504 ===\n");
k24list := Filtered(ints.subgroups, function(K) return Size(K) = 24; end);;
Print("Order-24 overgroups: ", Length(k24list), "\n");
for K in k24list do
    result := ComputeMeinhold("meinhold504", K, 6, 6);;
    if result <> fail then break; fi;
od;

# Meinhold-576: K = C₇:C₃ = F₂₁, order 21, target 7-regular girth 6
Print("=== Meinhold-576 ===\n");
k21list := Filtered(ints.subgroups, function(K) return Size(K) = 21; end);;
Print("Order-21 overgroups: ", Length(k21list), "\n");
for K in k21list do
    result := ComputeMeinhold("meinhold576", K, 7, 6);;
    if result <> fail then break; fi;
od;

# Meinhold-1344: K = C₃×C₃, order 9, target 9-regular girth 6
Print("=== Meinhold-1344 ===\n");
k9list := Filtered(ints.subgroups, function(K) return Size(K) = 9; end);;
Print("Order-9 overgroups: ", Length(k9list), "\n");
for K in k9list do
    result := ComputeMeinhold("meinhold1344", K, 9, 6);;
    if result <> fail then break; fi;
od;

# Meinhold-2016: K = C₆, order 6, target 6-regular girth 8
Print("=== Meinhold-2016 ===\n");
k6list := Filtered(ints.subgroups, function(K) return Size(K) = 6; end);;
Print("Order-6 overgroups: ", Length(k6list), "\n");
for K in k6list do
    result := ComputeMeinhold("meinhold2016", K, 6, 8);;
    if result <> fail then break; fi;
od;

Print("=======================================================================\n");
Print("MEINHOLD FAMILY COMPLETE\n");
Print("=======================================================================\n");

QuitGap(0);
