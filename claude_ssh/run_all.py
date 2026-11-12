"""
run_all.py -- reproduce the entire study from scratch.

Runs the four analysis parts in order; each one regenerates its figures
(into figures/), data files (into data/) and LTSpice netlists (netlists/).
Afterwards, compile the report with:

    cd report && tectonic report.tex

Everything is deterministic (fixed random seeds), so re-running reproduces
the same numbers exactly.
"""
import subprocess, sys, time

PARTS = [
    "part1_pristine.py",
    "part2_single_impurity.py",
    "part3_domain_wall.py",
    "part4_monte_carlo.py",
]

t0 = time.time()
for p in PARTS:
    print(f"\n>>> running {p}")
    r = subprocess.run([sys.executable, p])
    if r.returncode != 0:
        sys.exit(f"FAILED: {p}")
print(f"\nAll parts complete in {time.time()-t0:.1f} s.")
print("Now compile the report:  cd report && tectonic report.tex")
