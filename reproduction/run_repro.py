"""
run_repro.py — simulate every impurity configuration with real ngspice and
record the normalised edge-state voltage profiles needed to recreate the two
target JPGs.

For each configuration we:
  * run an ngspice AC sweep 35-55 kHz (4001 pts),
  * locate the edge-state resonance (the |V(n1)| peak),
  * read |V| on the odd nodes n1,n3,...,n15 at that resonance and normalise to n1.

Output: data/profiles.csv  (figure, panel, series, node, voltage)
        data/resonances.csv (figure, panel, series, f_peak_kHz)
"""
import csv
import numpy as np
import ssh_spice as S

NPTS = 4001
GAMMA_FIXED = -0.1                      # cell-3 strength sweep fixes gamma = -0.1
DELTAS = [-0.3, -0.1, 0.1, 0.3]         # cell-3 strength sweep
POSITIONS = [2, 4, 6]                   # "impurity" position sweep (unit cell)
# Position ("impurity") panels: fixed strength, cell swept. Matching the target:
#   AB impurity dips DOWN at the impurity node (delta>0) then returns to the line;
#   BA impurity bumps UP  at the impurity node (delta<0), nearly flat;
#   Fig.2 single A bumps up (delta<0), single B stays flat.  gamma fixed at -0.1.
AB_POS = (0.3, -0.1)                    # (delta, gamma)  -> downward dip
BA_POS = (-0.1, -0.1)                   # (delta, gamma)  -> small upward bump (near-flat, like target)
POS_SINGLE_DELTA = -0.1
NODES = list(range(1, 16, 2))           # 1,3,5,...,15


def profile_for(caps):
    freqs, V, _ = S.run_ngspice(caps, npts=NPTS)
    idx, fpk = S.find_resonance(freqs, V)
    return S.voltage_profile_odd(V, idx), fpk


def main():
    prof_rows = []
    res_rows = []

    def record(figure, panel, series, caps):
        prof, fpk = profile_for(caps)
        for node, v in zip(NODES, prof):
            prof_rows.append((figure, panel, series, node, f"{v:.6f}"))
        res_rows.append((figure, panel, series, f"{fpk/1e3:.4f}"))
        print(f"{figure:8s} {panel:14s} {series:14s} f_peak={fpk/1e3:6.3f} kHz  "
              f"n15={prof[-1]:.4f}")

    # ---- Figure 1 : AB / BA pair impurity --------------------------------
    for pos in POSITIONS:
        record("fig1", "AB impurity", f"AB_{pos}", S.impurity_caps("AB", pos, *AB_POS))
    for pos in POSITIONS:
        record("fig1", "BA impurity", f"BA_{pos}", S.impurity_caps("BA", pos, *BA_POS))
    # AB cell3: SYMMETRIC strength (delta=gamma) keeps cell-3's intra/inter ratio fixed,
    # so after the node-5 bump the decay stays parallel to pristine (a clean domain-wall offset).
    for d in DELTAS:
        record("fig1", "AB cell3", f"AB_{d:g}/{d:g}", S.impurity_caps("AB", 3, d, d))
    # BA cell3: gamma=-0.1 fixed, delta swept -> node-5 spread (delta<0 up, delta>0 down) that reconverges.
    for d in DELTAS:
        record("fig1", "BA cell3", f"BA_{GAMMA_FIXED:g}/{d:g}", S.impurity_caps("BA", 3, d, GAMMA_FIXED))

    # ---- Figure 2 : single A / B impurity --------------------------------
    for pos in POSITIONS:
        record("fig2", "A impurity", f"N={pos}", S.single_caps("A", pos, POS_SINGLE_DELTA))
    for pos in POSITIONS:
        record("fig2", "B impurity", f"N={pos}", S.single_caps("B", pos, POS_SINGLE_DELTA))
    for d in DELTAS:
        record("fig2", "A cell3", f"A_{d:g}", S.single_caps("A", 3, d))
    for d in DELTAS:
        record("fig2", "B cell3", f"B_{d:g}", S.single_caps("B", 3, d))

    # pristine reference line
    ref = S.pristine_reference()
    for node, v in zip(NODES, ref):
        prof_rows.append(("ref", "pristine", "theory", node, f"{v:.6f}"))

    with open("data/profiles.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["figure", "panel", "series", "node", "voltage"])
        w.writerows(prof_rows)
    with open("data/resonances.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["figure", "panel", "series", "f_peak_kHz"])
        w.writerows(res_rows)
    print(f"\nwrote data/profiles.csv ({len(prof_rows)} rows) and data/resonances.csv")


if __name__ == "__main__":
    main()
