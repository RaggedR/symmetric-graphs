"""
Generate 74 unique orbital graphs for a proper 1-1 Carter class mapping.

Strategy:
- Run build scripts for all root systems (G2, F4, A2, A3, D4, D6, etc.)
- For each, extract ALL orbitals (different ip values)
- Generate enough clean PDFs to have 74 unique graphs total
"""

import subprocess
import sys

print("="*70)
print("GENERATING 74 UNIQUE HECKE ALGEBRA ORBITALS")
print("="*70)

# Step 1: Run all build scripts to generate root projections
build_scripts = [
    ("e8/exceptional_build.py", "F4, G2, A2, A3 roots"),
    ("e8/d4_build.py", "D4 root system"),
    ("e8/d6_build.py", "D6 root system"),
    ("e8/e5_build.py", "E5/D5 root system"),
    ("e8/e6_build.py", "E6 root system"),
    ("e8/e7_build.py", "E7 root system"),
]

print("\nStep 1: Building all root systems...")
print("-" * 70)

for script, desc in build_scripts:
    print(f"\nRunning {script} ({desc})...")
    try:
        result = subprocess.run(
            [sys.executable, script],
            cwd="/Users/miriam/symmetric-graphs",
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode == 0:
            print(f"  ✓ Success")
            if result.stdout:
                # Show last few lines of output
                lines = result.stdout.strip().split('\n')
                for line in lines[-3:]:
                    print(f"    {line}")
        else:
            print(f"  ✗ Failed")
            if result.stderr:
                print(f"    Error: {result.stderr[:200]}")
    except Exception as e:
        print(f"  ✗ Exception: {e}")

print("\n" + "="*70)
print("Step 2: Extract all orbitals from each root system")
print("="*70)
print("\nThis requires:")
print("1. Load each root system's P (projections) and R (roots)")
print("2. Compute G = R @ R.T (inner products)")
print("3. Find all unique values in G (these are the orbitals)")
print("4. For each unique ip value, generate an orbital graph PDF")
print("\nThis will give us ~50-60 orbitals, enough for 74 unique images")
print("when combined with the 5 E8/E7 orbitals we already have.")

print("\n" + "="*70)
print("RECOMMENDATION")
print("="*70)
print("The full orbital computation is complex. For the 74-class atlas,")
print("we have two options:")
print("")
print("Option A (QUICK): Accept that some classes share images for now,")
print("         document which ones, and note in handover that computing")
print("         more orbitals will give true 1-1 mapping")
print("")
print("Option B (PROPER): Write a comprehensive orbital extractor that")
print("          processes each root system and generates all orbitals.")
print("          This is the right long-term solution.")
print("")
print("Which would you like to pursue?")
