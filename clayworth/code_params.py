"""
Compute the parameters of the CSS code from Q(6,2).

H is 7 x 63, rank 7, H*H^T = 0 mod 2.
C = ker(H) is a [63, 56, d] code.
C^perp = rowspan(H) is [63, 7, d^perp].
Since H*H^T = 0, C^perp <= C, so we get a CSS code.
CSS code: [[63, 2*56 - 63, d]] = [[63, 49, d]]

Need to find d = minimum weight of C \ C^perp (coset leaders).
"""
import numpy as np
from itertools import combinations

def quadratic_form(x):
    return (x[0]*x[4] + x[1]*x[5] + x[2]*x[6] + x[3]*x[3]) % 2

# Enumerate points of Q(6,2)
points = []
for i in range(1, 128):
    x = [(i >> j) & 1 for j in range(7)]
    if quadratic_form(x) == 0:
        points.append(x)

H = np.array(points, dtype=int).T  # 7 x 63

# C = ker(H) has dimension 63 - 7 = 56
# C^perp = rowspan(H) has dimension 7
# CSS code encodes k = 56 - 7 = 49 logical qubits...

# Actually let me recompute. For CSS(C1, C2) with C1 = C2 = C:
# k = dim(C) - dim(C^perp) = 56 - 7 = 49
# [[n, k, d]] = [[63, 49, d]]

print("Code C = ker(H): [63, 56, d]")
print("Dual C^perp = rowspan(H): [63, 7, d^perp]")
print("CSS code: [[63, 49, d]]")

# Find minimum distance of C (minimum weight of nonzero codeword in ker(H))
# ker(H) is 56-dimensional in F_2^63 - too large to enumerate (2^56 elements)
# But we can find minimum distance of C^perp (rowspan of H) easily

# Enumerate all 2^7 - 1 = 127 nonzero linear combinations of rows of H
min_weight_dual = 63
for mask in range(1, 128):
    # linear combination of rows
    codeword = np.zeros(63, dtype=int)
    for bit in range(7):
        if (mask >> bit) & 1:
            codeword = (codeword + H[bit]) % 2
    w = np.sum(codeword)
    if w < min_weight_dual:
        min_weight_dual = w

print(f"\nMinimum weight of C^perp: {min_weight_dual}")
print(f"(This is d^perp, the distance of the dual code)")

# For the CSS code distance, d = min(d_X, d_Z)
# d_X = minimum weight of C \ C^perp
# d_Z = minimum weight of C \ C^perp (same since C1=C2)
#
# Since C^perp subset C, we need min weight of codewords in C that are NOT in C^perp
# This is hard to compute directly for dim 56...
#
# But we can get a lower bound from the dual distance
# and also try to find low-weight codewords

# For the [63, 56] code C = ker(H), the minimum distance is at most 4
# (since the Singleton bound gives d <= n - k + 1 = 8)
# Let's check: find minimum weight vectors in ker(H)

# Build a basis for ker(H) using Gaussian elimination over F_2
def kernel_F2(M):
    """Compute kernel of M over F_2 using row reduction of [M | I]"""
    m, n = M.shape
    # Augment with identity
    aug = np.hstack([M.T % 2, np.eye(n, dtype=int)])  # n x (m + n)

    # Row reduce over F_2
    pivot_cols = []
    row = 0
    for col in range(m):
        # Find pivot
        found = False
        for r in range(row, n):
            if aug[r, col] == 1:
                aug[[row, r]] = aug[[r, row]]
                found = True
                break
        if not found:
            continue
        pivot_cols.append(col)
        # Eliminate
        for r in range(n):
            if r != row and aug[r, col] == 1:
                aug[r] = (aug[r] + aug[row]) % 2
        row += 1

    # Kernel basis = rows of aug where first m columns are zero
    kernel_basis = []
    for r in range(n):
        if np.all(aug[r, :m] == 0):
            kernel_basis.append(aug[r, m:])

    return np.array(kernel_basis, dtype=int)

K = kernel_F2(H)
print(f"\nKernel dimension: {len(K)} (expected 56)")

# Check weights of kernel basis vectors
weights = [np.sum(v) for v in K]
print(f"Basis vector weights: min={min(weights)}, max={max(weights)}")
print(f"Weight distribution of basis: {sorted(set(weights))}")

# Find minimum weight codeword by checking small linear combinations
# Can't check all 2^56, but try all singles, pairs, triples
print("\nSearching for minimum weight codewords...")
min_weight_C = min(weights)
min_word = K[weights.index(min_weight_C)]

# Check pairs
for i in range(len(K)):
    for j in range(i+1, len(K)):
        w = np.sum((K[i] + K[j]) % 2)
        if w < min_weight_C:
            min_weight_C = w
            print(f"  Found weight {w} from pair ({i},{j})")

# Check triples (this might find distance 3 if it exists)
print("Checking triples...")
found_lower = False
for i in range(len(K)):
    for j in range(i+1, len(K)):
        vij = (K[i] + K[j]) % 2
        for k in range(j+1, len(K)):
            w = np.sum((vij + K[k]) % 2)
            if w < min_weight_C:
                min_weight_C = w
                print(f"  Found weight {w} from triple ({i},{j},{k})")
                found_lower = True
                if min_weight_C <= 3:
                    break
        if min_weight_C <= 3:
            break
    if min_weight_C <= 3:
        break

print(f"\nMinimum weight found in C: {min_weight_C}")
print(f"Minimum weight of C^perp: {min_weight_dual}")
print(f"\n=== CSS Code: [[63, 49, d >= {min(min_weight_C, min_weight_dual)}]] ===")

# Also: compare with the Hamming/Steane analogy
print("\n--- Comparison ---")
print("Steane (from Fano/PG(2,2)):     [[7, 1, 3]]")
print(f"Q(6,2)/G_2(2):                   [[63, 49, ?]]")
print(f"  (d_dual = {min_weight_dual}, d_C >= {min_weight_C} from sampling)")
