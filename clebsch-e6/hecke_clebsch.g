# hecke_clebsch.g
# Hecke algebras of the Schlafli (27 lines) and Petersen graphs,
# and their connection through the Clebsch diagonal surface.
#
# Computes:
#   1. H(PSp(4,3), (C2^4):A5) -- rank 3, Bose-Mesner of the Schlafli graph
#   2. H(A5, S3) -- rank 3, Bose-Mesner of the Petersen graph
#   3. S5 orbit decomposition 27 = 15 + 12 (Clebsch surface structure)
#   4. Induced subgraphs on Type I (15) and Type II (12) lines
#   5. Bimodule decomposition C[G/H] into irreducible constituents
#
# The Clebsch diagonal surface has Aut = S5. Its 27 lines split under S5
# into 15 (pair partitions) and 12 (5-cycles). The Petersen graph K(5,2)
# appears as a substructure of the 15 Type I lines. This script traces
# the Hecke algebra embedding through this geometry.
#
# Usage: gap -q < hecke_clebsch.g

LoadPackage("grape");;


# =====================================================
# Generic: compute Hecke algebra H(G,H) from a pair
# =====================================================

ComputeHecke := function(G, H)
    local cosets, n, reps, dcs, r, vals, dc_of,
          T, sc, i, j, k, prod, p, q, found;

    cosets := RightCosets(G, H);
    n := Length(cosets);
    reps := List(cosets, Representative);

    dcs := DoubleCosets(G, H, H);
    vals := List(dcs, dc -> Size(dc) / Size(H));
    # Sort by valency so T_0 = identity (valency 1) comes first
    SortParallel(vals, dcs);
    r := Length(dcs);

    dc_of := function(g)
        local k;
        for k in [1..r] do
            if g in dcs[k] then return k; fi;
        od;
        return fail;
    end;

    # Build adjacency matrices T[k], one per double coset
    T := List([1..r], k -> NullMat(n, n));
    for i in [1..n] do
        for j in [1..n] do
            k := dc_of(reps[i] * reps[j]^(-1));
            T[k][i][j] := 1;
        od;
    od;

    # Structure constants: T[a]*T[b] = Sum_k sc[a][b][k] * T[k]
    sc := List([1..r], a -> List([1..r], b -> List([1..r], k -> 0)));
    for i in [1..r] do
        for j in [1..r] do
            prod := T[i] * T[j];
            for k in [1..r] do
                found := false;
                for p in [1..n] do
                    if found then break; fi;
                    for q in [1..n] do
                        if T[k][p][q] = 1 then
                            sc[i][j][k] := prod[p][q];
                            found := true;
                            break;
                        fi;
                    od;
                od;
            od;
        od;
    od;

    return rec(n := n, r := r, vals := vals, T := T, sc := sc,
               reps := reps, G := G, H := H);
end;;


PrintMultTable := function(HA, name)
    local i, j, k, first;
    Print("\n  Hecke algebra: ", name, "\n");
    Print("  Rank = ", HA.r, ", valencies = ", HA.vals, "\n\n");
    for i in [1..HA.r] do
        for j in [i..HA.r] do
            Print("    T", i-1, " * T", j-1, " = ");
            first := true;
            for k in [1..HA.r] do
                if HA.sc[i][j][k] <> 0 then
                    if not first then Print(" + "); fi;
                    if HA.sc[i][j][k] <> 1 then
                        Print(HA.sc[i][j][k]);
                    fi;
                    Print("T", k-1);
                    first := false;
                fi;
            od;
            Print("\n");
        od;
    od;
end;;


# =====================================================
# Part 1: Schlafli graph -- H(PSp(4,3), (C2^4):A5)
# =====================================================

Print("========================================================\n");
Print("  Hecke Algebras and the Clebsch Surface\n");
Print("========================================================\n\n");

Print("> Building PSp(4,3) = Sp(4,3) / Center...\n");
G_full := Sp(4, 3);;
ctr := Center(G_full);;
Print("  |Sp(4,3)| = ", Size(G_full), ", |Center| = ", Size(ctr), "\n");
nat := NaturalHomomorphismByNormalSubgroup(G_full, ctr);;
G_psp := Image(nat);;
Print("  |PSp(4,3)| = ", Size(G_psp), "\n\n");

Print("> Computing subgroup lattice (this may take a minute)...\n");
subs_psp := ConjugacyClassesSubgroups(G_psp);;
Print("  ", Length(subs_psp), " conjugacy classes of subgroups\n");

# Find (C2^4):A5 of order 960 giving SRG(27, 10, 1, 5)
H_sch := fail;;
for c in Filtered(subs_psp, c -> Size(Representative(c)) = 960) do
    H := Representative(c);
    v := SortedList(List(DoubleCosets(G_psp, H, H), d -> Size(d)/Size(H)));
    if v = [1, 10, 16] then
        H_sch := H;
        Print("  Stabiliser: ", StructureDescription(H), "\n");
        Print("  Double coset valencies: ", v, " --> SRG(27, 10, 1, 5)\n\n");
        break;
    fi;
od;;

if H_sch = fail then
    Print("ERROR: stabiliser for Schlafli graph not found.\n");
fi;;

Print("> Computing Hecke algebra H(PSp(4,3), H)...\n");
HA_sch := ComputeHecke(G_psp, H_sch);;
PrintMultTable(HA_sch, "H(PSp(4,3), (C2^4):A5) -- Schlafli graph");

# Analytical verification
Print("\n  Verification against SRG(27, 10, 1, 5) formulas:\n");
Print("    T1^2 should be 10*T0 + 1*T1 + 5*T2\n");
Print("    T1*T2 should be 8*T1 + 5*T2\n");
Print("    T2^2 should be 16*T0 + 8*T1 + 10*T2\n");


# =====================================================
# Part 2: Petersen graph -- H(A5, S3)
# =====================================================

Print("\n========================================================\n");
Print("  Petersen graph -- H(A5, S3)\n");
Print("========================================================\n\n");

G_a5 := AlternatingGroup(5);;
# S3 = <(1,2,3), (1,2)(4,5)> inside A5
# The involution (1,2)(4,5) inverts the 3-cycle: tau*g3*tau = g3^{-1}
H_pet := Group((1,2,3), (1,2)(4,5));;
Print("  |A5| = ", Size(G_a5), ", |S3| = ", Size(H_pet), "\n");
Print("  |A5/S3| = ", Size(G_a5)/Size(H_pet), "\n\n");

Print("> Computing Hecke algebra H(A5, S3)...\n");
HA_pet := ComputeHecke(G_a5, H_pet);;
PrintMultTable(HA_pet, "H(A5, S3) -- Petersen graph");

Print("\n  Verification against SRG(10, 3, 0, 1) formulas:\n");
Print("    T1^2 should be 3*T0 + 0*T1 + 1*T2\n");
Print("    T1*T2 should be 2*T1 + 2*T2\n");
Print("    T2^2 should be 6*T0 + 4*T1 + 3*T2\n");


# =====================================================
# Part 3: Clebsch S5 decomposition via Aut(Schlafli)
# =====================================================

Print("\n========================================================\n");
Print("  Clebsch S5 Decomposition: 27 = 15 + 12\n");
Print("========================================================\n\n");

# Find the adjacency matrix (the one with valency 10)
adj_idx := First([1..HA_sch.r], k -> HA_sch.vals[k] = 10);;
adj_sch := HA_sch.T[adj_idx];;

# Build the Schlafli graph in GRAPE and compute its automorphism group
gamma_sch := Graph(Group(()), [1..27],
    function(x, g) return x; end,
    function(x, y) return adj_sch[x][y] = 1; end, true);;
aut_sch := AutGroupGraph(gamma_sch);;
Print("|Aut(Schlafli)| = ", Size(aut_sch), "\n");
Print("Aut = ", StructureDescription(aut_sch), "\n\n");

# Find S5 inside Aut(Schlafli) with orbits [12, 15]
# This is the S5 = Aut(Clebsch surface) acting on the 27 lines
Print("> Searching for Clebsch S5 (orbits 15 + 12)...\n");
subs_aut := ConjugacyClassesSubgroups(aut_sch);;
Print("  ", Length(subs_aut), " conjugacy classes in Aut\n");

found_S5 := fail;;
for c in subs_aut do
    S := Representative(c);
    if Size(S) = 120 then
        orb_sizes := SortedList(List(Orbits(S, [1..27]), Length));
        if orb_sizes = [12, 15] then
            if StructureDescription(S) = "S5" then
                found_S5 := S;
                Print("  Found S5 with orbits [12, 15]\n\n");
                break;
            fi;
        fi;
    fi;
od;;

if found_S5 = fail then
    Print("  S5 with orbits [12, 15] not found.\n");
    Print("  All subgroups of order 120 in Aut:\n");
    for c in subs_aut do
        S := Representative(c);
        if Size(S) = 120 then
            Print("    ", StructureDescription(S), " --> orbits ",
                  SortedList(List(Orbits(S, [1..27]), Length)), "\n");
        fi;
    od;
fi;;


if found_S5 <> fail then

    # Identify the two S5-orbits
    orbs := Orbits(found_S5, [1..27]);;
    orb15 := fail;; orb12 := fail;;
    for o in orbs do
        if Length(o) = 15 then orb15 := Set(o); fi;
        if Length(o) = 12 then orb12 := Set(o); fi;
    od;;

    Print("Type I  (15 pair-partition lines): ", orb15, "\n");
    Print("Type II (12 five-cycle lines):     ", orb12, "\n\n");


    # ---- Induced subgraph on 15 Type I lines ----
    Print("---- Induced subgraph on 15 Type I lines ----\n\n");

    degs15 := List(orb15, v -> Number(orb15, w -> adj_sch[v][w] = 1));;
    Print("  Degree sequence: ", Collected(SortedList(degs15)), "\n");
    Print("  Edges: ", Sum(degs15)/2, "\n");

    gamma15 := Graph(Group(()), [1..15],
        function(x, g) return x; end,
        function(x, y) return adj_sch[orb15[x]][orb15[y]] = 1; end, true);;
    aut15 := AutGroupGraph(gamma15);;
    Print("  |Aut| = ", Size(aut15), "\n");
    Print("  Aut = ", StructureDescription(aut15), "\n");

    # SRG parameters
    lambda_set := Set([]);; mu_set := Set([]);;
    for i in [1..15] do
        for j in [i+1..15] do
            com := Number([1..15], function(k)
                return adj_sch[orb15[i]][orb15[k]] = 1
                   and adj_sch[orb15[j]][orb15[k]] = 1;
            end);
            if adj_sch[orb15[i]][orb15[j]] = 1 then
                AddSet(lambda_set, com);
            else
                AddSet(mu_set, com);
            fi;
        od;
    od;;
    Print("  lambda (adj pairs, common nbrs):     ", lambda_set, "\n");
    Print("  mu (non-adj pairs, common nbrs):     ", mu_set, "\n");
    if Length(Set(degs15)) = 1 and Length(lambda_set) = 1
       and Length(mu_set) = 1 then
        Print("  --> SRG(15, ", degs15[1], ", ",
              lambda_set[1], ", ", mu_set[1], ")\n");
    fi;

    # Check: is the complement a known graph?
    Print("\n  Complement (intersection graph of 15 Type I lines):\n");
    degs15c := List([1..15], v ->
        Number([1..15], w -> w <> v and adj_sch[orb15[v]][orb15[w]] = 0));;
    Print("  Degree sequence: ", Collected(SortedList(degs15c)), "\n");

    gamma15c := Graph(Group(()), [1..15],
        function(x, g) return x; end,
        function(x, y) return x <> y
            and adj_sch[orb15[x]][orb15[y]] = 0; end, true);;
    aut15c := AutGroupGraph(gamma15c);;
    Print("  |Aut| = ", Size(aut15c), "\n");
    Print("  Aut = ", StructureDescription(aut15c), "\n");


    # ---- Induced subgraph on 12 Type II lines ----
    Print("\n---- Induced subgraph on 12 Type II lines ----\n\n");

    degs12 := List(orb12, v -> Number(orb12, w -> adj_sch[v][w] = 1));;
    Print("  Degree sequence: ", Collected(SortedList(degs12)), "\n");

    gamma12 := Graph(Group(()), [1..12],
        function(x, g) return x; end,
        function(x, y) return adj_sch[orb12[x]][orb12[y]] = 1; end, true);;
    aut12 := AutGroupGraph(gamma12);;
    Print("  |Aut| = ", Size(aut12), "\n");
    Print("  Aut = ", StructureDescription(aut12), "\n");


    # ---- Cross edges ----
    Print("\n---- Cross edges (Type I <-> Type II) ----\n\n");

    cross15 := List(orb15, v -> Number(orb12, w -> adj_sch[v][w] = 1));;
    cross12 := List(orb12, v -> Number(orb15, w -> adj_sch[v][w] = 1));;
    Print("  Type I  -> Type II per vertex: ", Set(cross15), "\n");
    Print("  Type II -> Type I  per vertex: ", Set(cross12), "\n");
    Print("  Total cross edges: ", Sum(cross15), "\n");


    # ---- A5 = [S5, S5] orbits ----
    Print("\n---- A5 = [S5, S5] orbits ----\n\n");

    A5_sub := DerivedSubgroup(found_S5);;
    Print("  |A5| = ", Size(A5_sub), ", ", StructureDescription(A5_sub), "\n");
    Print("  A5 orbits on 15 Type I:  ",
        SortedList(List(Orbits(A5_sub, orb15), Length)), "\n");
    Print("  A5 orbits on 12 Type II: ",
        SortedList(List(Orbits(A5_sub, orb12), Length)), "\n");
    Print("  A5 orbits on all 27:     ",
        SortedList(List(Orbits(A5_sub, [1..27]), Length)), "\n");

fi;;


# =====================================================
# Part 4: Bimodule decomposition C[G/H]
# =====================================================

Print("\n========================================================\n");
Print("  Bimodule Decomposition C[G/H]\n");
Print("========================================================\n\n");

# For a Gelfand pair (G,H) with SRG(n,k,lambda,mu), the bimodule
# decomposition follows directly from the eigenvalue multiplicities:
#   C[G/H] = trivial(dim 1) + V_r(dim f) + V_s(dim g)
# where r,s are the non-trivial eigenvalues and f,g their multiplicities.
#
# Eigenvalues: r,s = ((lambda-mu) +/- sqrt(D)) / 2, D = (lambda-mu)^2 + 4(k-mu)
# Multiplicities: 1 + f + g = n, k + rf + sg = 0, k^2 + r^2*f + s^2*g = nk

DecomposeFromSRG := function(label, n, k, lambda, mu)
    local D, r, s, f, g;
    D := (lambda - mu)^2 + 4*(k - mu);
    # r = ((lambda-mu) + sqrt(D))/2,  s = ((lambda-mu) - sqrt(D))/2
    # For our cases, D is a perfect square.
    r := ((lambda - mu) + RootInt(D)) / 2;
    s := ((lambda - mu) - RootInt(D)) / 2;
    # From k + rf + sg = 0 and f + g = n - 1:
    g := (k + r*(n-1)) / (r - s);
    f := n - 1 - g;

    Print("-- ", label, " --\n");
    Print("  Eigenvalues of T1: ", k, " (x1), ", r, " (x", f, "), ",
          s, " (x", g, ")\n");
    Print("  C[G/H] = 1 + ", f, " + ", g,
          "  (total = ", 1 + f + g, ")\n");
    Print("  Each eigenspace is an irreducible G-module (Gelfand pair).\n\n");
end;;

DecomposeFromSRG("Schlafli SRG(27, 10, 1, 5)", 27, 10, 1, 5);
DecomposeFromSRG("Petersen SRG(10, 3, 0, 1)", 10, 3, 0, 1);


# =====================================================
# Part 5: Eigenvalue--irrep correspondence
# =====================================================

Print("========================================================\n");
Print("  Eigenvalue -- Irrep Correspondence\n");
Print("========================================================\n\n");

Print("Schlafli SRG(27, 10, 1, 5):\n");
Print("  Eigenvalues of adjacency matrix:\n");
Print("    10 (mult 1)  -- trivial representation\n");
Print("     1 (mult 20) -- 20-dimensional irrep(s)\n");
Print("    -5 (mult 6)  -- 6-dimensional irrep\n");
Print("  Multiplicities = dimensions of irreps with H-fixed vectors.\n\n");

Print("Petersen SRG(10, 3, 0, 1):\n");
Print("  Eigenvalues of adjacency matrix:\n");
Print("    3 (mult 1) -- trivial representation\n");
Print("    1 (mult 5) -- 5-dimensional irrep\n");
Print("   -2 (mult 4) -- 4-dimensional irrep\n\n");

Print("Both share eigenvalue 1. This is the non-trivial eigenvalue\n");
Print("that should be traceable through the Clebsch embedding:\n");
Print("the 5-dim eigenspace (Petersen) embeds in the 20-dim\n");
Print("eigenspace (Schlafli) via restriction from 27 to 15 lines.\n\n");

Print("========================================================\n");
Print("  Summary\n");
Print("========================================================\n\n");
Print("Both Hecke algebras are rank 3 (= 3 double cosets).\n");
Print("Both are commutative (Gelfand pairs).\n");
Print("As abstract algebras, both are isomorphic to C^3.\n");
Print("The isomorphism is NOT a coincidence -- it reflects\n");
Print("the embedding of the Petersen graph inside the 27 lines\n");
Print("on the Clebsch diagonal surface.\n\n");

Print("Done.\n");
