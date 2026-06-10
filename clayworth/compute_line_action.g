LoadPackage("grape");;

# Compute the G₂(2) generator action on LINES of GH(2,2)
# Given the generators acting on points (from G2Action.lean),
# compute the induced action on lines (via the lineToTriangle data).

# Generator 1 (order 7): point permutation from G2Action.lean
gen1 := PermList([
  2, 4, 6, 3, 1, 7, 5, 24, 26, 25, 27, 29, 31, 28, 30, 32, 34, 35, 33, 36,
  38, 39, 37, 56, 59, 57, 58, 61, 62, 60, 63, 8, 9, 10, 11, 13, 12, 15, 14,
  16, 18, 17, 19, 20, 22, 21, 23, 40, 41, 43, 42, 46, 47, 45, 44, 48, 50, 49,
  51, 52, 54, 53, 55
]);;

# Generator 2 (order 6): point permutation from G2Action.lean
gen2 := PermList([
  11, 21, 30, 50, 56, 37, 45, 15, 26, 60, 41, 1, 20, 51, 36, 6, 8, 48, 59,
  19, 25, 35, 42, 12, 18, 63, 34, 7, 29, 49, 46, 52, 62, 33, 43, 4, 9, 23,
  28, 58, 47, 5, 22, 53, 32, 13, 24, 54, 61, 17, 27, 39, 44, 2, 10, 57, 38,
  3, 31, 55, 40, 14, 16
]);;

# The 63 lines as sorted triples of point indices (1-indexed for GAP)
# From compute_triangle_data.g output, converted to 1-indexed
lines := [
  [1, 16, 17], [1, 32, 33], [1, 48, 49],
  [2, 8, 9], [2, 32, 34], [2, 40, 41],
  [3, 24, 25], [3, 32, 35], [3, 56, 57],
  [4, 8, 10], [4, 16, 18], [4, 24, 26],
  [5, 16, 19], [5, 40, 42], [5, 56, 58],
  [6, 8, 11], [6, 48, 50], [6, 56, 59],
  [7, 24, 27], [7, 40, 43], [7, 48, 51],
  [9, 36, 47], [9, 38, 46], [10, 20, 31],
  [10, 22, 29], [11, 52, 62], [11, 54, 61],
  [12, 18, 30], [12, 34, 45], [12, 50, 63],
  [13, 23, 26], [13, 34, 44], [13, 55, 59],
  [14, 18, 28], [14, 39, 41], [14, 53, 59],
  [15, 21, 26], [15, 37, 41], [15, 50, 60],
  [17, 36, 55], [17, 37, 54], [19, 44, 62],
  [19, 46, 60], [20, 33, 53], [20, 42, 63],
  [21, 33, 52], [21, 47, 58], [22, 39, 49],
  [22, 45, 58], [23, 38, 49], [23, 42, 61],
  [25, 36, 63], [25, 39, 62], [27, 45, 54],
  [27, 46, 53], [28, 35, 61], [28, 47, 51],
  [29, 35, 60], [29, 43, 55], [30, 38, 57],
  [30, 43, 52], [31, 37, 57], [31, 44, 51]
];;

# For each generator, compute the induced line permutation
# A generator σ maps line {p1,p2,p3} to {σ(p1),σ(p2),σ(p3)}, sorted
# Then find which line index this corresponds to

lineIndex := function(triple)
  local sorted, i;
  sorted := ShallowCopy(triple);
  Sort(sorted);
  for i in [1..63] do
    if lines[i] = sorted then return i; fi;
  od;
  return fail;
end;;

for genName in ["gen1", "gen2"] do
  if genName = "gen1" then gen := gen1; else gen := gen2; fi;
  perm := [];
  for i in [1..63] do
    img := List(lines[i], p -> p^gen);
    idx := lineIndex(img);
    if idx = fail then
      Print("ERROR: line ", i, " maps to non-line ", img, "\n");
    fi;
    Add(perm, idx);
  od;
  Print(genName, "_on_lines (0-indexed): [");
  for i in [1..63] do
    Print(perm[i]-1);
    if i < 63 then Print(", "); fi;
  od;
  Print("]\n\n");
od;
