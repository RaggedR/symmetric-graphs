# symmetric-layout.g
# General tool: find subgroups of Aut(Γ) whose orbit quotients have low degree
# Output JSON for the Python layout engine
#
# Usage: gap -q symmetric-layout.g -- graphname [max_degree]
# e.g.:  gap -q symmetric-layout.g -- dodecahedron 3

LoadPackage("grape");
LoadPackage("json", false);  # may not be available, we'll format manually

# ============================================================
# Graph constructors
# ============================================================

# Generalized Petersen graph GP(n,k)
# Outer cycle: 1-2-...-n-1, inner star: (n+i) adjacent to n+((i+k-1) mod n)+1
# Spokes: i adjacent to n+i
GPGraph := function(n, k)
    local adj, i, ip, im;
    adj := List([1..2*n], i -> []);
    for i in [1..n] do
        # Outer cycle
        ip := (i mod n) + 1;
        Add(adj[i], ip); Add(adj[ip], i);
        # Spoke
        Add(adj[i], n+i); Add(adj[n+i], i);
        # Inner star: step k
        ip := n + ((i-1+k) mod n) + 1;
        if not ip in adj[n+i] then
            Add(adj[n+i], ip); Add(adj[ip], n+i);
        fi;
    od;
    # Deduplicate
    for i in [1..2*n] do adj[i] := Set(adj[i]); od;
    return adj;
end;

# LCF graph: Hamiltonian cycle 1-2-...-n-1 plus chords from shift pattern
# shifts is a list, repeated to fill n vertices
LCFGraph := function(n, shifts)
    local adj, i, s, j, period;
    adj := List([1..n], i -> []);
    period := Length(shifts);
    for i in [1..n] do
        # Hamiltonian cycle
        j := (i mod n) + 1;
        if not j in adj[i] then Add(adj[i], j); Add(adj[j], i); fi;
        # LCF chord
        s := shifts[((i-1) mod period) + 1];
        j := ((i-1+s) mod n) + 1;
        if not j in adj[i] then Add(adj[i], j); Add(adj[j], i); fi;
    od;
    for i in [1..n] do adj[i] := Set(adj[i]); od;
    return adj;
end;

# Voltage graph: derived graph from a voltage assignment on a base multigraph.
# base_n: number of base vertices
# base_edges: list of [i, j, v] where i,j are base vertices (1-indexed),
#   v is the voltage (integer mod group_order) on the directed edge i→j.
#   For i=j this is a self-loop. Reverse edge j→i gets voltage -v.
#   Each edge should appear only once (i ≤ j, or i=j for loops).
# group_order: size of the cyclic voltage group Z_m.
# Returns: adjacency list on base_n * group_order vertices.
# Vertex (i, g) maps to index (i-1)*group_order + g + 1, where g ∈ {0,..,m-1}.
VoltageGraph := function(base_n, base_edges, group_order)
    local n, adj, e, i, j, v, g, u, w, m;
    m := group_order;
    n := base_n * m;
    adj := List([1..n], x -> []);
    for e in base_edges do
        i := e[1]; j := e[2]; v := e[3];
        for g in [0..m-1] do
            u := (i-1)*m + g + 1;
            w := (j-1)*m + ((g + v) mod m) + 1;
            if u <> w then
                Add(adj[u], w);
                Add(adj[w], u);
            fi;
        od;
    od;
    for i in [1..n] do adj[i] := Set(adj[i]); od;
    return adj;
end;

# Complete graph K_n
CompleteAdj := function(n)
    return List([1..n], i -> Filtered([1..n], j -> j <> i));
end;

# Complete bipartite K_{a,b}
CompleteBipartiteAdj := function(a, b)
    local adj, i;
    adj := [];
    for i in [1..a] do Add(adj, List([1..b], j -> a+j)); od;
    for i in [1..b] do Add(adj, List([1..a], j -> j)); od;
    return adj;
end;

# ============================================================
# Graph library
# ============================================================

GraphLibrary := rec();

# --- Foster census F004: K4 ---
GraphLibrary.k4 := function()
    return rec(name := "K4", n := 4, adj := CompleteAdj(4));
end;

# --- Foster census F006: K_{3,3} ---
GraphLibrary.k33 := function()
    return rec(name := "K33", n := 6, adj := CompleteBipartiteAdj(3, 3));
end;

# --- Foster census F008: Cube = GP(4,1) ---
GraphLibrary.cube := function()
    return rec(name := "Cube", n := 8, adj := GPGraph(4, 1));
end;

# --- Foster census F010: Petersen = GP(5,2) ---
GraphLibrary.petersen := function()
    return rec(name := "Petersen", n := 10, adj := GPGraph(5, 2));
end;

# --- Foster census F014: Heawood, LCF [5,-5]^7 ---
GraphLibrary.heawood := function()
    return rec(name := "Heawood", n := 14, adj := LCFGraph(14, [5, -5]));
end;

# --- Foster census F016: Mobius-Kantor = GP(8,3) ---
GraphLibrary.moebiuskantor := function()
    return rec(name := "Moebius-Kantor", n := 16, adj := GPGraph(8, 3));
end;

# --- Foster census F018: Pappus, LCF [5,7,-7,7,-7,-5]^3 ---
GraphLibrary.pappus := function()
    return rec(name := "Pappus", n := 18, adj := LCFGraph(18, [5, 7, -7, 7, -7, -5]));
end;

# --- Foster census F020A: Dodecahedron ---
# LCF [10,7,4,-4,-7,10,-4,7,-7,4]^2, but we use a labeling where the C5
# orbit quotient is a path — this gives a Schlegel diagram with 0 crossings
GraphLibrary.dodecahedron := function()
    return rec(name := "Dodecahedron", n := 20, adj := [
        [2,5,6], [1,3,8], [2,4,10], [3,5,12], [1,4,14],
        [1,7,15], [6,8,16], [2,7,9], [8,10,17], [3,9,11],
        [10,12,18], [4,11,13], [12,14,19], [5,13,15], [6,14,20],
        [7,17,20], [9,16,18], [11,17,19], [13,18,20], [15,16,19]
    ]);
end;

# --- Foster census F020B: Desargues = GP(10,3) ---
GraphLibrary.desargues := function()
    return rec(name := "Desargues", n := 20, adj := GPGraph(10, 3));
end;

# --- Foster census F024: Nauru = GP(12,5) ---
GraphLibrary.nauru := function()
    return rec(name := "Nauru", n := 24, adj := GPGraph(12, 5));
end;

# --- Foster census F026A, LCF [-7,7]^13 ---
GraphLibrary.f26a := function()
    return rec(name := "F26A", n := 26, adj := LCFGraph(26, [-7, 7]));
end;

# --- Foster census F028: Coxeter graph ---
# Hub (1-7) connected to three 7-cycle rings with steps 1, 2, 3
GraphLibrary.coxeter := function()
    local adj, i, j, steps;
    adj := List([1..28], i -> []);
    steps := [1, 2, 3];
    for j in [0..6] do
        # Spokes from hub vertex j+1 to each ring
        for i in [1..3] do
            Add(adj[j+1], 7*i + j + 1);
            Add(adj[7*i + j + 1], j+1);
        od;
        # Ring edges: ring i has step steps[i]
        for i in [1..3] do
            Add(adj[7*i + j + 1], 7*i + ((j + steps[i]) mod 7) + 1);
            Add(adj[7*i + ((j + steps[i]) mod 7) + 1], 7*i + j + 1);
        od;
    od;
    for i in [1..28] do adj[i] := Set(adj[i]); od;
    return rec(name := "Coxeter", n := 28, adj := adj);
end;

# --- Foster census F030: Tutte-Coxeter (Tutte 8-cage), LCF [-13,-9,7,-7,9,13]^5 ---
GraphLibrary.tuttecoxeter := function()
    return rec(name := "Tutte-Coxeter", n := 30,
               adj := LCFGraph(30, [-13, -9, 7, -7, 9, 13]));
end;

# ============================================================
# Voltage graph examples
# ============================================================

# Heawood as Z_7-covering of K_2 (3 parallel edges, voltages {0,4,6})
# Base: 2 vertices, edges 1→2 with voltages 0, 4, 6
# Recovers the incidence graph of PG(2,2) = Fano plane Levi graph
GraphLibrary.heawood_voltage := function()
    return rec(name := "Heawood (voltage)", n := 14,
               adj := VoltageGraph(2, [[1,2,0], [1,2,4], [1,2,6]], 7));
end;

# Petersen as Z_5-covering of a "theta" multigraph
# Base: 2 vertices, self-loop on 1 (voltage 1 = outer cycle step),
#   self-loop on 2 (voltage 2 = inner star step), edge 1-2 (voltage 0 = spoke)
GraphLibrary.petersen_voltage := function()
    return rec(name := "Petersen (voltage)", n := 10,
               adj := VoltageGraph(2, [[1,1,1], [2,2,2], [1,2,0]], 5));
end;

# GP(n,k) as Z_n-covering: outer cycle step 1, inner star step k, spoke step 0
GPVoltage := function(n, k)
    return VoltageGraph(2, [[1,1,1], [2,2,k], [1,2,0]], n);
end;

# Cube = GP(4,1) via voltage
GraphLibrary.cube_voltage := function()
    return rec(name := "Cube (voltage)", n := 8, adj := GPVoltage(4, 1));
end;

# Moebius-Kantor = GP(8,3) via voltage
GraphLibrary.mk_voltage := function()
    return rec(name := "Moebius-Kantor (voltage)", n := 16, adj := GPVoltage(8, 3));
end;

# Nauru = GP(12,5) via voltage
GraphLibrary.nauru_voltage := function()
    return rec(name := "Nauru (voltage)", n := 24, adj := GPVoltage(12, 5));
end;

# Desargues = GP(10,3) via voltage
GraphLibrary.desargues_voltage := function()
    return rec(name := "Desargues (voltage)", n := 20, adj := GPVoltage(10, 3));
end;

# PG(2,3) with identity matching removed → cubic graph on 26 vertices
# Singer cycle: points = Z_13, lines = translates of {0,1,3,9}
# Remove "identity" edges (point i on line i) → voltages {1,3,9} mod 13
GraphLibrary.pg23 := function()
    return rec(name := "PG(2,3) projection", n := 26,
               adj := VoltageGraph(2, [[1,2,1], [1,2,3], [1,2,9]], 13));
end;

# Cubic residue graph CR(p): for prime p ≡ 1 mod 3
# Voltage graph on K₂ with voltages = cubic residues {1, g, g²} mod p
CRGraph := function(p)
    local g, i, cr;
    for g in [2..p-1] do
        if PowerMod(g, 3, p) = 1 and g <> 1 then
            cr := [1, g, PowerMod(g, 2, p)];
            return VoltageGraph(2, List(cr, v -> [1,2,v]), p);
        fi;
    od;
    return fail;
end;

GraphLibrary.cr19 := function()
    return rec(name := "CR(19)", n := 38, adj := CRGraph(19));
end;

GraphLibrary.cr31 := function()
    return rec(name := "CR(31)", n := 62, adj := CRGraph(31));
end;

# ============================================================
# PSL(2,q) Sabidussi graphs
# ============================================================

# Biggs-Smith graph: Sab(PSL(2,17), S4, a)
# 102 vertices, cubic, girth 9, Aut = PSL(2,17), distance-transitive
GraphLibrary.biggssmith := function()
    return rec(name := "Biggs-Smith", n := 102, adj :=
    [ [ 32, 40, 71 ], [ 11, 44, 72 ], [ 15, 27, 73 ], [ 30, 39, 74 ],
      [ 41, 49, 75 ], [ 50, 58, 76 ], [ 13, 56, 69 ], [ 17, 29, 70 ],
      [ 46, 61, 93 ], [ 29, 35, 76 ], [ 2, 38, 91 ], [ 22, 58, 74 ],
      [ 7, 28, 86 ], [ 30, 34, 84 ], [ 3, 68, 87 ], [ 23, 47, 92 ],
      [ 8, 52, 77 ], [ 53, 65, 83 ], [ 59, 68, 95 ], [ 28, 54, 93 ],
      [ 41, 46, 87 ], [ 12, 45, 96 ], [ 16, 33, 97 ], [ 55, 66, 98 ],
      [ 31, 51, 80 ], [ 40, 47, 69 ], [ 3, 48, 92 ], [ 13, 20, 75 ],
      [ 8, 10, 89 ], [ 4, 14, 78 ], [ 25, 45, 101 ], [ 1, 61, 81 ],
      [ 23, 67, 86 ], [ 14, 64, 94 ], [ 10, 62, 98 ], [ 50, 55, 78 ],
      [ 43, 60, 79 ], [ 11, 51, 70 ], [ 4, 57, 80 ], [ 1, 26, 90 ],
      [ 5, 21, 82 ], [ 56, 59, 82 ], [ 37, 54, 97 ], [ 2, 66, 84 ],
      [ 22, 31, 89 ], [ 9, 21, 99 ], [ 16, 26, 95 ], [ 27, 60, 71 ],
      [ 5, 63, 79 ], [ 6, 36, 85 ], [ 25, 38, 94 ], [ 17, 64, 85 ],
      [ 18, 62, 101 ], [ 20, 43, 90 ], [ 24, 36, 100 ], [ 7, 42, 88 ],
      [ 39, 65, 72 ], [ 6, 12, 83 ], [ 19, 42, 102 ], [ 37, 48, 99 ],
      [ 9, 32, 88 ], [ 35, 53, 91 ], [ 49, 67, 73 ], [ 34, 52, 96 ],
      [ 18, 57, 100 ], [ 24, 44, 77 ], [ 33, 63, 102 ], [ 15, 19, 81 ],
      [ 7, 26, 84 ], [ 8, 38, 87 ], [ 1, 48, 78 ], [ 2, 57, 82 ],
      [ 3, 63, 85 ], [ 4, 12, 88 ], [ 5, 28, 77 ], [ 6, 10, 81 ],
      [ 17, 66, 75 ], [ 30, 36, 71 ], [ 37, 49, 100 ], [ 25, 39, 99 ],
      [ 32, 68, 76 ], [ 41, 42, 72 ], [ 18, 58, 102 ], [ 14, 44, 69 ],
      [ 50, 52, 73 ], [ 13, 33, 96 ], [ 15, 21, 70 ], [ 56, 61, 74 ],
      [ 29, 45, 93 ], [ 40, 54, 98 ], [ 11, 62, 95 ], [ 16, 27, 94 ],
      [ 9, 20, 89 ], [ 34, 51, 92 ], [ 19, 47, 91 ], [ 22, 64, 86 ],
      [ 23, 43, 101 ], [ 24, 35, 90 ], [ 46, 60, 80 ], [ 55, 65, 79 ],
      [ 31, 53, 97 ], [ 59, 67, 83 ] ]);
end;

# PSL(2,13)/S3: cubic arc-transitive, 182 vertices, girth 7, Aut = PSL(2,13)
GraphLibrary.psl2_13 := function()
    return rec(name := "PSL2(13)", n := 182, adj :=
    [ [ 30, 115, 160 ], [ 37, 120, 130 ], [ 51, 85, 170 ], [ 65, 93, 168 ],
      [ 77, 100, 138 ], [ 18, 107, 153 ], [ 42, 120, 169 ], [ 56, 150, 171 ],
      [ 48, 135, 176 ], [ 34, 133, 178 ], [ 75, 112, 121 ], [ 26, 83, 115 ],
      [ 74, 112, 116 ], [ 71, 97, 163 ], [ 25, 124, 144 ], [ 24, 154, 170 ],
      [ 55, 105, 150 ], [ 6, 84, 159 ], [ 50, 85, 157 ], [ 39, 84, 120 ],
      [ 69, 164, 181 ], [ 60, 140, 147 ], [ 47, 127, 135 ], [ 16, 86, 119 ],
      [ 15, 119, 121 ], [ 12, 105, 179 ], [ 38, 89, 132 ], [ 32, 142, 168 ],
      [ 67, 110, 164 ], [ 1, 79, 166 ], [ 64, 93, 176 ], [ 28, 94, 124 ],
      [ 53, 79, 85 ], [ 10, 129, 180 ], [ 68, 152, 163 ], [ 59, 140, 161 ],
      [ 2, 80, 159 ], [ 27, 86, 124 ], [ 20, 110, 167 ], [ 52, 96, 171 ],
      [ 45, 138, 146 ], [ 7, 117, 180 ], [ 76, 100, 147 ], [ 72, 163, 178 ],
      [ 41, 89, 101 ], [ 61, 80, 93 ], [ 23, 143, 154 ], [ 9, 128, 179 ],
      [ 57, 150, 153 ], [ 19, 125, 154 ], [ 3, 81, 166 ], [ 40, 89, 94 ],
      [ 33, 117, 133 ], [ 66, 104, 181 ], [ 17, 107, 152 ], [ 8, 127, 179 ],
      [ 49, 96, 108 ], [ 73, 81, 100 ], [ 36, 142, 144 ], [ 22, 167, 169 ],
      [ 46, 125, 135 ], [ 78, 112, 129 ], [ 70, 160, 164 ], [ 31, 90, 142 ],
      [ 4, 82, 159 ], [ 54, 96, 101 ], [ 29, 115, 128 ], [ 35, 133, 157 ],
      [ 21, 161, 167 ], [ 63, 104, 116 ], [ 14, 82, 107 ], [ 44, 132, 146 ],
      [ 58, 90, 140 ], [ 13, 119, 143 ], [ 11, 130, 180 ], [ 43, 97, 146 ],
      [ 5, 83, 166 ], [ 62, 104, 108 ], [ 30, 33, 162 ], [ 37, 46, 165 ],
      [ 51, 58, 149 ], [ 65, 71, 173 ], [ 12, 77, 156 ], [ 18, 20, 131 ],
      [ 3, 19, 33 ], [ 24, 38, 91 ], [ 114, 128, 158 ], [ 102, 137, 181 ],
      [ 27, 45, 52 ], [ 64, 73, 141 ], [ 86, 131, 155 ], [ 109, 129, 141 ],
      [ 4, 31, 46 ], [ 32, 52, 98 ], [ 123, 169, 177 ], [ 40, 57, 66 ],
      [ 14, 76, 145 ], [ 94, 162, 175 ], [ 113, 143, 145 ], [ 5, 43, 58 ],
      [ 45, 66, 103 ], [ 88, 134, 157 ], [ 101, 137, 165 ], [ 54, 70, 78 ],
      [ 17, 26, 174 ], [ 122, 144, 174 ], [ 6, 55, 71 ], [ 57, 78, 111 ],
      [ 92, 158, 176 ], [ 29, 39, 155 ], [ 108, 141, 149 ], [ 11, 13, 62 ],
      [ 99, 147, 177 ], [ 87, 132, 155 ], [ 1, 12, 67 ], [ 13, 70, 118 ],
      [ 42, 53, 175 ], [ 116, 145, 173 ], [ 24, 25, 74 ], [ 2, 7, 20 ],
      [ 11, 25, 126 ], [ 106, 134, 152 ], [ 95, 171, 175 ], [ 15, 32, 38 ],
      [ 50, 61, 137 ], [ 121, 156, 174 ], [ 23, 56, 177 ], [ 48, 67, 87 ],
      [ 34, 62, 92 ], [ 2, 75, 172 ], [ 84, 91, 182 ], [ 27, 72, 114 ],
      [ 10, 53, 68 ], [ 102, 122, 161 ], [ 9, 23, 61 ], [ 149, 153, 182 ],
      [ 88, 103, 125 ], [ 5, 41, 151 ], [ 148, 162, 168 ], [ 22, 36, 73 ],
      [ 90, 92, 111 ], [ 28, 59, 64 ], [ 47, 74, 99 ], [ 15, 59, 106 ],
      [ 97, 99, 118 ], [ 41, 72, 76 ], [ 22, 43, 113 ], [ 139, 160, 173 ],
      [ 81, 111, 136 ], [ 8, 17, 49 ], [ 138, 165, 172 ], [ 35, 55, 122 ],
      [ 6, 49, 136 ], [ 16, 47, 50 ], [ 91, 110, 114 ], [ 83, 126, 172 ],
      [ 19, 68, 102 ], [ 87, 109, 178 ], [ 18, 37, 65 ], [ 1, 63, 148 ],
      [ 36, 69, 134 ], [ 79, 98, 139 ], [ 14, 35, 44 ], [ 21, 29, 63 ],
      [ 80, 103, 151 ], [ 30, 51, 77 ], [ 39, 60, 69 ], [ 4, 28, 139 ],
      [ 7, 60, 95 ], [ 3, 16, 182 ], [ 8, 40, 123 ], [ 130, 151, 156 ],
      [ 82, 118, 148 ], [ 105, 106, 126 ], [ 98, 117, 123 ], [ 9, 31, 109 ],
      [ 95, 113, 127 ], [ 10, 44, 158 ], [ 26, 48, 56 ], [ 34, 42, 75 ],
      [ 21, 54, 88 ], [ 131, 136, 170 ] ]);
end;

# PSL(2,5)/C3: double cover of Petersen, 20 vertices, girth 5, Aut = C2 x A5
GraphLibrary.psl2_5_c3 := function()
    return rec(name := "PSL2(5)/C3", n := 20, adj :=
    [ [ 2, 5, 11 ], [ 1, 10, 14 ], [ 4, 7, 15 ], [ 3, 9, 17 ],
      [ 1, 6, 18 ], [ 5, 8, 13 ], [ 3, 8, 12 ], [ 6, 7, 19 ], [ 4, 10, 16 ],
      [ 2, 9, 20 ], [ 1, 13, 20 ], [ 7, 13, 17 ], [ 6, 11, 12 ], [ 2, 16, 18 ],
      [ 3, 16, 19 ], [ 9, 14, 15 ], [ 4, 12, 20 ], [ 5, 14, 19 ], [ 8, 15, 18 ],
      [ 10, 11, 17 ] ]);
end;

# PSL(2,7)/C3: 56 vertices, girth 7, Aut = PGL(2,7) = PSL(3,2) : C2
GraphLibrary.psl2_7_c3 := function()
    return rec(name := "PSL2(7)/C3", n := 56, adj :=
    [ [ 41, 48, 55 ], [ 17, 35, 46 ], [ 19, 27, 56 ], [ 14, 31, 53 ],
      [ 28, 33, 42 ], [ 10, 21, 39 ], [ 12, 25, 49 ], [ 17, 48, 56 ],
      [ 34, 42, 45 ], [ 6, 29, 50 ], [ 28, 38, 51 ], [ 7, 20, 40 ],
      [ 16, 24, 33 ], [ 4, 26, 44 ], [ 31, 42, 55 ], [ 13, 49, 52 ],
      [ 2, 8, 36 ], [ 25, 37, 45 ], [ 3, 34, 47 ], [ 12, 23, 30 ], [ 6, 26, 51 ],
      [ 34, 48, 50 ], [ 20, 41, 43 ], [ 13, 36, 55 ], [ 7, 18, 54 ],
      [ 14, 21, 35 ], [ 3, 32, 40 ], [ 5, 11, 47 ], [ 10, 41, 49 ],
      [ 20, 38, 56 ], [ 4, 15, 43 ], [ 27, 44, 52 ], [ 5, 13, 54 ],
      [ 9, 19, 22 ], [ 2, 26, 37 ], [ 17, 24, 44 ], [ 18, 35, 50 ],
      [ 11, 30, 46 ], [ 6, 47, 52 ], [ 12, 27, 53 ], [ 1, 23, 29 ], [ 5, 9, 15 ],
      [ 23, 31, 51 ], [ 14, 32, 36 ], [ 9, 18, 53 ], [ 2, 38, 54 ],
      [ 19, 28, 39 ], [ 1, 8, 22 ], [ 7, 16, 29 ], [ 10, 22, 37 ],
      [ 11, 21, 43 ], [ 16, 32, 39 ], [ 4, 40, 45 ], [ 25, 33, 46 ],
      [ 1, 15, 24 ], [ 3, 8, 30 ] ]);
end;

# --- Split Cayley hexagon GH(2,2) incidence graph = Tutte 12-cage ---
# 126 vertices, cubic, girth 12, Aut = G_2(2) (order 12096)
# Constructed from the Zorn vector matrix model of split octonions over F_2
# Points: vertices 1-63 (on parabolic quadric Q(6,2))
# Lines: vertices 64-126 (totally isotropic lines with p·q ∈ ⟨p,q⟩)
GraphLibrary.tuttecage12 := function()
    return rec(name := "Tutte-12-cage", n := 126, adj := [
    [ 64, 65, 66 ], [ 67, 68, 69 ], [ 70, 71, 72 ], [ 73, 74, 75 ],
    [ 76, 77, 78 ], [ 79, 80, 81 ], [ 82, 83, 84 ], [ 67, 73, 79 ],
    [ 67, 85, 86 ], [ 73, 87, 88 ], [ 79, 89, 90 ], [ 91, 92, 93 ],
    [ 94, 95, 96 ], [ 97, 98, 99 ], [ 100, 101, 102 ], [ 64, 74, 76 ],
    [ 64, 103, 104 ], [ 74, 91, 97 ], [ 76, 105, 106 ], [ 87, 107, 108 ],
    [ 100, 109, 110 ], [ 88, 111, 112 ], [ 94, 113, 114 ], [ 70, 75, 82 ],
    [ 70, 115, 116 ], [ 75, 94, 100 ], [ 82, 117, 118 ], [ 97, 119, 120 ],
    [ 88, 121, 122 ], [ 91, 123, 124 ], [ 87, 125, 126 ], [ 65, 68, 71 ],
    [ 65, 107, 109 ], [ 68, 92, 95 ], [ 71, 119, 121 ], [ 85, 103, 115 ],
    [ 101, 104, 125 ], [ 86, 113, 123 ], [ 98, 111, 116 ], [ 69, 77, 83 ],
    [ 69, 98, 101 ], [ 77, 108, 114 ], [ 83, 122, 124 ], [ 95, 105, 126 ],
    [ 92, 112, 117 ], [ 86, 106, 118 ], [ 85, 110, 120 ], [ 66, 80, 84 ],
    [ 66, 111, 113 ], [ 80, 93, 102 ], [ 84, 120, 126 ], [ 89, 109, 124 ],
    [ 99, 107, 118 ], [ 90, 104, 117 ], [ 96, 103, 122 ], [ 72, 78, 81 ],
    [ 72, 123, 125 ], [ 78, 110, 112 ], [ 81, 96, 99 ], [ 102, 106, 121 ],
    [ 90, 114, 119 ], [ 89, 105, 116 ], [ 93, 108, 115 ],
    [ 1, 16, 17 ], [ 1, 32, 33 ], [ 1, 48, 49 ], [ 2, 8, 9 ],
    [ 2, 32, 34 ], [ 2, 40, 41 ], [ 3, 24, 25 ], [ 3, 32, 35 ],
    [ 3, 56, 57 ], [ 4, 8, 10 ], [ 4, 16, 18 ], [ 4, 24, 26 ],
    [ 5, 16, 19 ], [ 5, 40, 42 ], [ 5, 56, 58 ], [ 6, 8, 11 ],
    [ 6, 48, 50 ], [ 6, 56, 59 ], [ 7, 24, 27 ], [ 7, 40, 43 ],
    [ 7, 48, 51 ], [ 9, 36, 47 ], [ 9, 38, 46 ], [ 10, 20, 31 ],
    [ 10, 22, 29 ], [ 11, 52, 62 ], [ 11, 54, 61 ], [ 12, 18, 30 ],
    [ 12, 34, 45 ], [ 12, 50, 63 ], [ 13, 23, 26 ], [ 13, 34, 44 ],
    [ 13, 55, 59 ], [ 14, 18, 28 ], [ 14, 39, 41 ], [ 14, 53, 59 ],
    [ 15, 21, 26 ], [ 15, 37, 41 ], [ 15, 50, 60 ], [ 17, 36, 55 ],
    [ 17, 37, 54 ], [ 19, 44, 62 ], [ 19, 46, 60 ], [ 20, 33, 53 ],
    [ 20, 42, 63 ], [ 21, 33, 52 ], [ 21, 47, 58 ], [ 22, 39, 49 ],
    [ 22, 45, 58 ], [ 23, 38, 49 ], [ 23, 42, 61 ], [ 25, 36, 63 ],
    [ 25, 39, 62 ], [ 27, 45, 54 ], [ 27, 46, 53 ], [ 28, 35, 61 ],
    [ 28, 47, 51 ], [ 29, 35, 60 ], [ 29, 43, 55 ], [ 30, 38, 57 ],
    [ 30, 43, 52 ], [ 31, 37, 57 ], [ 31, 44, 51 ]
    ]);
end;

# GH(2,2) collinearity graph: point graph of the generalized hexagon of order (2,2)
# 63 vertices, 6-regular, arc-transitive, distance-regular {6,4,4;1,1,3}
# Aut = G_2(2), order 12096
GraphLibrary.gh22 := function()
    return rec(name := "GH(2,2) collinearity", n := 63, adj := [
    [ 19, 20, 31, 41, 44, 60 ], [ 15, 16, 35, 38, 49, 51 ],
    [ 23, 24, 27, 36, 51, 60 ], [ 5, 21, 40, 44, 45, 57 ],
    [ 4, 25, 38, 43, 47, 57 ], [ 7, 17, 48, 51, 52, 63 ],
    [ 6, 28, 29, 39, 45, 63 ], [ 10, 13, 55, 60, 61, 62 ],
    [ 11, 12, 33, 44, 56, 58 ], [ 8, 26, 32, 34, 43, 62 ],
    [ 9, 25, 32, 46, 52, 58 ], [ 9, 14, 30, 38, 50, 56 ],
    [ 8, 22, 29, 30, 37, 55 ], [ 12, 21, 28, 50, 53, 61 ],
    [ 2, 28, 33, 37, 49, 54 ], [ 2, 18, 34, 35, 45, 59 ],
    [ 6, 18, 26, 33, 42, 48 ], [ 16, 17, 25, 42, 59, 61 ],
    [ 1, 30, 32, 41, 42, 54 ], [ 1, 22, 31, 39, 43, 59 ],
    [ 4, 14, 22, 40, 52, 53 ], [ 13, 20, 21, 37, 52, 59 ],
    [ 3, 27, 34, 37, 47, 50 ], [ 3, 26, 29, 36, 46, 53 ],
    [ 5, 11, 18, 46, 47, 61 ], [ 10, 17, 24, 33, 43, 53 ],
    [ 3, 23, 39, 40, 42, 58 ], [ 7, 14, 15, 39, 54, 61 ],
    [ 7, 13, 24, 30, 45, 46 ], [ 12, 13, 19, 29, 38, 42 ],
    [ 1, 20, 46, 48, 49, 50 ], [ 10, 11, 19, 34, 52, 54 ],
    [ 9, 15, 17, 26, 37, 44 ], [ 10, 16, 23, 32, 45, 50 ],
    [ 2, 16, 41, 53, 55, 58 ], [ 3, 24, 54, 56, 57, 59 ],
    [ 13, 15, 22, 23, 33, 47 ], [ 2, 5, 12, 30, 43, 51 ],
    [ 7, 20, 27, 28, 43, 58 ], [ 4, 21, 27, 42, 49, 62 ],
    [ 1, 19, 35, 47, 53, 63 ], [ 17, 18, 19, 27, 30, 40 ],
    [ 5, 10, 20, 26, 38, 39 ], [ 1, 4, 9, 33, 45, 60 ],
    [ 4, 7, 16, 29, 34, 44 ], [ 11, 24, 25, 29, 31, 49 ],
    [ 5, 23, 25, 37, 41, 63 ], [ 6, 17, 31, 50, 55, 57 ],
    [ 2, 15, 31, 40, 46, 62 ], [ 12, 14, 23, 31, 34, 48 ],
    [ 2, 3, 6, 38, 52, 60 ], [ 6, 11, 21, 22, 32, 51 ],
    [ 14, 21, 24, 26, 35, 41 ], [ 15, 19, 28, 32, 36, 57 ],
    [ 8, 13, 35, 48, 57, 58 ], [ 9, 12, 36, 59, 62, 63 ],
    [ 4, 5, 36, 48, 54, 55 ], [ 9, 11, 27, 35, 39, 55 ],
    [ 16, 18, 20, 22, 36, 56 ], [ 1, 3, 8, 44, 51, 61 ],
    [ 8, 14, 18, 25, 28, 60 ], [ 8, 10, 40, 49, 56, 63 ],
    [ 6, 7, 41, 47, 56, 62 ]
    ]);
end;

# --- Meinhold-72: Quotient of Clayworth via K = PSL(3,2) ---
# 72 vertices, 7-regular, girth 4, diameter 5
# Aut = G_2(2), order 12096
GraphLibrary.meinhold72 := function()
    return rec(name := "Meinhold-72", n := 72, adj := [
    [ 4, 24, 32, 38, 43, 52, 59 ],
    [ 10, 15, 17, 20, 33, 44, 69 ],
    [ 10, 11, 20, 29, 34, 37, 47 ],
    [ 1, 14, 18, 28, 48, 65, 71 ],
    [ 11, 22, 25, 35, 37, 56, 66 ],
    [ 7, 15, 22, 33, 38, 59, 66 ],
    [ 6, 8, 23, 39, 41, 54, 57 ],
    [ 7, 9, 22, 24, 25, 32, 53 ],
    [ 8, 12, 23, 26, 58, 61, 64 ],
    [ 2, 3, 31, 39, 40, 51, 57 ],
    [ 3, 5, 14, 36, 40, 65, 72 ],
    [ 9, 25, 30, 44, 56, 60, 69 ],
    [ 26, 45, 49, 50, 61, 67, 70 ],
    [ 4, 11, 16, 20, 43, 56, 69 ],
    [ 2, 6, 18, 45, 57, 70, 71 ],
    [ 14, 18, 26, 36, 45, 46, 64 ],
    [ 2, 18, 19, 23, 39, 48, 64 ],
    [ 4, 15, 16, 17, 21, 38, 63 ],
    [ 17, 21, 30, 34, 44, 47, 68 ],
    [ 2, 3, 14, 45, 48, 49, 62 ],
    [ 18, 19, 27, 28, 46, 50, 70 ],
    [ 5, 6, 8, 36, 45, 61, 62 ],
    [ 7, 9, 17, 30, 33, 55, 63 ],
    [ 1, 8, 31, 39, 48, 58, 62 ],
    [ 5, 8, 12, 31, 57, 65, 71 ],
    [ 9, 13, 16, 32, 43, 60, 63 ],
    [ 21, 35, 38, 42, 52, 66, 68 ],
    [ 4, 21, 30, 55, 56, 59, 66 ],
    [ 3, 41, 49, 51, 54, 67, 72 ],
    [ 12, 19, 23, 28, 54, 65, 72 ],
    [ 10, 24, 25, 34, 35, 44, 52 ],
    [ 1, 8, 26, 36, 54, 65, 67 ],
    [ 2, 6, 23, 40, 49, 61, 72 ],
    [ 3, 19, 31, 48, 50, 65, 67 ],
    [ 5, 27, 31, 40, 50, 58, 61 ],
    [ 11, 16, 22, 32, 38, 47, 68 ],
    [ 3, 5, 41, 45, 46, 50, 57 ],
    [ 1, 6, 18, 27, 36, 39, 40 ],
    [ 7, 10, 17, 24, 38, 42, 47 ],
    [ 10, 11, 33, 35, 38, 43, 63 ],
    [ 7, 29, 37, 42, 53, 55, 66 ],
    [ 27, 39, 41, 46, 51, 58, 64 ],
    [ 1, 14, 26, 40, 49, 51, 58 ],
    [ 2, 12, 19, 31, 61, 62, 70 ],
    [ 13, 15, 16, 20, 22, 37, 53 ],
    [ 16, 21, 37, 42, 47, 56, 60 ],
    [ 3, 19, 36, 39, 46, 54, 62 ],
    [ 4, 17, 20, 24, 34, 53, 55 ],
    [ 13, 20, 29, 33, 43, 55, 59 ],
    [ 13, 21, 34, 35, 37, 55, 63 ],
    [ 10, 29, 42, 43, 52, 60, 69 ],
    [ 1, 27, 31, 51, 67, 70, 71 ],
    [ 8, 41, 45, 48, 64, 67, 71 ],
    [ 7, 29, 30, 32, 47, 59, 60 ],
    [ 23, 28, 41, 48, 49, 50, 58 ],
    [ 5, 12, 14, 28, 46, 58, 62 ],
    [ 7, 10, 15, 25, 37, 60, 63 ],
    [ 9, 24, 35, 42, 43, 55, 56 ],
    [ 1, 6, 28, 49, 54, 62, 70 ],
    [ 12, 26, 46, 51, 54, 57, 70 ],
    [ 9, 13, 22, 33, 35, 44, 68 ],
    [ 20, 22, 24, 44, 47, 56, 59 ],
    [ 18, 23, 26, 40, 50, 57, 65 ],
    [ 9, 16, 17, 42, 53, 68, 69 ],
    [ 4, 11, 25, 30, 32, 34, 63 ],
    [ 5, 6, 27, 28, 41, 71, 72 ],
    [ 13, 29, 32, 34, 52, 53, 68 ],
    [ 19, 27, 36, 61, 64, 67, 72 ],
    [ 2, 12, 14, 51, 64, 71, 72 ],
    [ 13, 15, 21, 44, 52, 59, 60 ],
    [ 4, 15, 25, 52, 53, 66, 69 ],
    [ 11, 29, 30, 33, 66, 68, 69 ]
    ]);
end;

# ============================================================
# Core functions
# ============================================================

# Build GRAPE graph from adjacency list
BuildGraph := function(data)
    return Graph(Group(()), [1..data.n], function(x,g) return x; end,
        function(x,y) return y in data.adj[x]; end, true);
end;

# Compute quotient from a list of orbits
ComputeQuotient := function(orbits, adj, n)
    local n_orbs, quot_adj, i, j, found, u, v, valencies, max_val,
          edge_counts, internal_counts;

    n_orbs := Length(orbits);
    quot_adj := List([1..n_orbs], i -> []);
    edge_counts := List([1..n_orbs], i -> List([1..n_orbs], j -> 0));
    internal_counts := List([1..n_orbs], i -> 0);

    for i in [1..n_orbs] do
        # Internal edges
        for u in orbits[i] do
            for v in orbits[i] do
                if u < v and v in adj[u] then
                    internal_counts[i] := internal_counts[i] + 1;
                fi;
            od;
        od;
        # Cross edges
        for j in [i+1..n_orbs] do
            found := 0;
            for u in orbits[i] do
                for v in orbits[j] do
                    if v in adj[u] then
                        found := found + 1;
                    fi;
                od;
            od;
            edge_counts[i][j] := found;
            edge_counts[j][i] := found;
            if found > 0 then
                Add(quot_adj[i], j);
                Add(quot_adj[j], i);
            fi;
        od;
    od;

    valencies := List(quot_adj, Length);
    max_val := Maximum(valencies);

    return rec(
        n_orbits := n_orbs,
        adj := quot_adj,
        valencies := valencies,
        max_valency := max_val,
        edge_counts := edge_counts,
        internal_counts := internal_counts
    );
end;

# Analyse angular offsets between orbits
# For each pair of adjacent orbits, determine how cross-edges connect positions
AnalyseOffsets := function(orbits, adj, generator, n)
    local n_orbs, offsets, i, j, orbit_order, pos_of, k, v, seq, start,
          u, w, pos_u, pos_w, diff, diffs, mode_diff;

    n_orbs := Length(orbits);

    # For each orbit, determine the ordering under the generator
    orbit_order := [];
    pos_of := List([1..n], x -> rec(orbit := 0, pos := 0));

    for i in [1..n_orbs] do
        start := Minimum(orbits[i]);
        seq := [start];
        v := start ^ generator;
        while v <> start do
            Add(seq, v);
            v := v ^ generator;
        od;
        Add(orbit_order, seq);
        for k in [1..Length(seq)] do
            pos_of[seq[k]] := rec(orbit := i, pos := k - 1);
        od;
    od;

    # For each pair of adjacent orbits, find the position offsets
    offsets := List([1..n_orbs], i -> List([1..n_orbs], j -> 0));
    for i in [1..n_orbs] do
        for j in [i+1..n_orbs] do
            diffs := [];
            for u in orbits[i] do
                for w in orbits[j] do
                    if w in adj[u] then
                        pos_u := pos_of[u].pos;
                        pos_w := pos_of[w].pos;
                        diff := (pos_w - pos_u) mod Length(orbits[i]);
                        Add(diffs, diff);
                    fi;
                od;
            od;
            if Length(diffs) > 0 then
                # Find most common offset (mode)
                Sort(diffs);
                mode_diff := diffs[1];  # simple: take smallest
                offsets[i][j] := mode_diff;
                offsets[j][i] := (-mode_diff) mod Length(orbits[i]);
            fi;
        od;
    od;

    return rec(
        orbit_order := orbit_order,
        pos_of := pos_of,
        offsets := offsets
    );
end;

# ============================================================
# Search for good quotients
# ============================================================

SearchQuotients := function(data, max_degree)
    local graph, aut, adj, n, results, seen, ProcessSubgroup,
          norms, conj_classes, cc, rep, H, g, ord, orders_seen, cyc;

    n := data.n;
    adj := data.adj;
    graph := BuildGraph(data);
    aut := AutGroupGraph(graph);

    data.aut_order := Size(aut);
    data.aut_structure := StructureDescription(aut);

    Print("Graph: ", data.name, " (", n, " vertices)\n");
    Print("|Aut| = ", data.aut_order, "  ", data.aut_structure, "\n\n");

    results := [];
    seen := [];  # track orbit sets we've already processed

    ProcessSubgroup := function(H, label)
        local orbs, orb_key, quot, offset_data, orbit_sizes, gen,
              has_pentagon, i, seq, start, v, is_cycle;

        if Size(H) = 1 or Size(H) = Size(aut) then return; fi;

        orbs := Orbits(H, [1..n]);
        orbs := List(orbs, Set);
        Sort(orbs, function(a,b) return Minimum(a) < Minimum(b); end);

        orb_key := orbs;
        if orb_key in seen then return; fi;
        Add(seen, orb_key);

        orbit_sizes := List(orbs, Length);

        quot := ComputeQuotient(orbs, adj, n);

        if quot.max_valency > max_degree then return; fi;
        # Allow 2 orbits if they're equal-sized (generalized Petersen pattern)
        if quot.n_orbits < 2 then return; fi;
        if quot.n_orbits = 2 and Length(Set(orbit_sizes)) > 1 then return; fi;

        # Check if quotient is a cycle
        is_cycle := ForAll(quot.valencies, v -> v = 2);

        # Find the highest-order element of H that preserves each orbit
        # (needed for cycle detection and offset analysis)
        gen := fail;
        for g in H do
            if Order(g) > 1 and ForAll(orbs, o -> Set(List(o, x -> x^g)) = o) then
                if gen = fail or Order(g) > Order(gen) then
                    gen := g;
                fi;
            fi;
        od;

        # Check which orbits form cycles under any generator of H
        has_pentagon := List([1..Length(orbs)], i -> false);
        if gen <> fail then
            for i in [1..Length(orbs)] do
                start := Minimum(orbs[i]);
                seq := [start];
                v := start ^ gen;
                while v <> start do
                    Add(seq, v);
                    v := v ^ gen;
                od;
                # Check if consecutive elements are adjacent
                if Length(seq) = Length(orbs[i]) then
                    has_pentagon[i] := ForAll([1..Length(seq)], function(k)
                        local next;
                        next := seq[(k mod Length(seq)) + 1];
                        return next in adj[seq[k]];
                    end);
                fi;
            od;
        fi;

        Add(results, rec(
            label := label,
            subgroup_order := Size(H),
            subgroup_structure := StructureDescription(H),
            n_orbits := quot.n_orbits,
            orbit_sizes := orbit_sizes,
            orbits := orbs,
            quotient_adj := quot.adj,
            quotient_valencies := quot.valencies,
            max_valency := quot.max_valency,
            edge_counts := quot.edge_counts,
            internal_counts := quot.internal_counts,
            is_cycle := is_cycle,
            has_internal_cycle := has_pentagon,
            generator := gen
        ));

        Print("  ", label, ": |H|=", Size(H), " (", StructureDescription(H), ")");
        Print(", ", quot.n_orbits, " orbits ", orbit_sizes);
        Print(", max_deg=", quot.max_valency);
        if is_cycle then Print(" [CYCLE]"); fi;
        Print("\n");
    end;

    # 1. Normal subgroups
    Print("--- Normal subgroups ---\n");
    norms := NormalSubgroups(aut);
    for H in norms do
        ProcessSubgroup(H, Concatenation("Normal |", String(Size(H)), "|"));
    od;

    # 2. Conjugacy classes of subgroups (for small groups)
    if Size(aut) <= 1000 then
        Print("\n--- Conjugacy classes of subgroups ---\n");
        conj_classes := ConjugacyClassesSubgroups(aut);
        for cc in conj_classes do
            rep := Representative(cc);
            ProcessSubgroup(rep, Concatenation("Conj |", String(Size(rep)), "| ",
                StructureDescription(rep)));
        od;
    else
        # For large groups, just try cyclic subgroups
        Print("\n--- Cyclic subgroups (group too large for full search) ---\n");
        orders_seen := [];
        for g in aut do
            ord := Order(g);
            if ord >= 2 and not ord in orders_seen then
                Add(orders_seen, ord);
                cyc := Group(g);
                ProcessSubgroup(cyc, Concatenation("Cyclic Z", String(ord)));
            fi;
        od;
    fi;

    return results;
end;

# ============================================================
# Output JSON
# ============================================================

OutputJSON := function(data, results)
    local r, i, j, offset_data;

    Print("\n\n=== JSON OUTPUT ===\n");
    Print("{\n");
    Print("  \"graph\": \"", data.name, "\",\n");
    Print("  \"n\": ", data.n, ",\n");
    Print("  \"aut_order\": ", data.aut_order, ",\n");
    Print("  \"aut_structure\": \"", data.aut_structure, "\",\n");
    Print("  \"adj\": ", data.adj, ",\n");
    Print("  \"results\": [\n");

    for i in [1..Length(results)] do
        r := results[i];
        Print("    {\n");
        Print("      \"label\": \"", r.label, "\",\n");
        Print("      \"subgroup_order\": ", r.subgroup_order, ",\n");
        Print("      \"subgroup_structure\": \"", r.subgroup_structure, "\",\n");
        Print("      \"n_orbits\": ", r.n_orbits, ",\n");
        Print("      \"orbit_sizes\": ", r.orbit_sizes, ",\n");
        Print("      \"orbits\": ", r.orbits, ",\n");
        Print("      \"quotient_adj\": ", r.quotient_adj, ",\n");
        Print("      \"quotient_valencies\": ", r.quotient_valencies, ",\n");
        Print("      \"max_valency\": ", r.max_valency, ",\n");
        Print("      \"edge_counts\": ", r.edge_counts, ",\n");
        Print("      \"internal_counts\": ", r.internal_counts, ",\n");
        Print("      \"is_cycle\": ", r.is_cycle, ",\n");
        Print("      \"has_internal_cycle\": ", r.has_internal_cycle, "\n");

        # Compute offsets if we have a generator
        if r.generator <> fail then
            offset_data := AnalyseOffsets(r.orbits, data.adj, r.generator, data.n);
            Print("      ,\"orbit_order\": ", offset_data.orbit_order, "\n");
            Print("      ,\"offsets\": ", offset_data.offsets, "\n");
        fi;

        if i < Length(results) then
            Print("    },\n");
        else
            Print("    }\n");
        fi;
    od;

    Print("  ]\n");
    Print("}\n");
    Print("=== END JSON ===\n");
end;

# ============================================================
# Main
# ============================================================

# Parse arguments — set these before running
# Override by defining graph_name and max_deg before Read()ing this file
# Set skip_main := true before Read()ing to import definitions only
if not (IsBound(skip_main) and skip_main = true) then

if not IsBound(graph_name) then
    graph_name := "dodecahedron";
fi;
if not IsBound(max_deg) then
    max_deg := 4;
fi;

if not IsBound(GraphLibrary.(graph_name)) then
    Print("Unknown graph: ", graph_name, "\n");
    Print("Available: ");
    for graph_name in RecNames(GraphLibrary) do
        Print(graph_name, " ");
    od;
    Print("\n");
else

data := GraphLibrary.(graph_name)();
Print("Searching for quotients of ", data.name, " with max degree <= ", max_deg, "\n\n");
results := SearchQuotients(data, max_deg);

Print("\n=== SUMMARY: ", Length(results), " quotients found with max degree <= ", max_deg, " ===\n\n");
for i in [1..Length(results)] do
    r := results[i];
    Print(i, ". ", r.label, "\n");
    Print("   ", r.n_orbits, " orbits of sizes ", r.orbit_sizes, "\n");
    Print("   Quotient degree sequence: ", r.quotient_valencies, "\n");
    Print("   Internal cycles: ", r.has_internal_cycle, "\n");
    if r.is_cycle then Print("   ** Quotient is a cycle! **\n"); fi;
    Print("\n");
od;

OutputJSON(data, results);

fi; # end if graph found
fi; # end skip_main guard
