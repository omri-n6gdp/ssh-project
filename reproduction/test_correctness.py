"""
test_correctness.py — verify the ngspice netlist is correct.

Two independent checks:
  1. ngspice node voltages vs. a from-scratch numpy nodal-admittance solve of the
     IDENTICAL topology, across the whole 35-55 kHz sweep, for several circuits.
  2. Physical sanity: pristine resonance ~= analytic f_R, odd-node decay ~= 0.667/cell.
"""
import numpy as np
import ssh_spice as S

def max_rel_err(a, b):
    denom = np.maximum(np.abs(b), 1e-12)
    return float(np.max(np.abs(a - b) / denom))

def check(label, caps):
    freqs, Vsp, _ = S.run_ngspice(caps, npts=801)
    Vnp = S.solve_numpy(caps, freqs)
    # compare magnitudes (phase reference is identical, but mag is what we plot)
    err = max_rel_err(np.abs(Vsp), np.abs(Vnp))
    # compare at resonance specifically
    idx, fpk = S.find_resonance(freqs, Vsp)
    err_res = max_rel_err(np.abs(Vsp[idx]), np.abs(Vnp[idx]))
    print(f"{label:28s} f_peak={fpk/1e3:6.3f} kHz   "
          f"max|.| rel-err(sweep)={err:.2e}   (at resonance)={err_res:.2e}")
    return err

if __name__ == "__main__":
    print(f"Analytic f_R = {S.F_R/1e3:.4f} kHz\n")
    print("=== ngspice  vs  independent numpy nodal solve (same circuit) ===")
    errs = []
    errs.append(check("pristine", S.pristine_caps()))
    errs.append(check("AB cell3 d=-0.3 g=-0.1", S.impurity_caps("AB", 3, -0.3, -0.1)))
    errs.append(check("AB cell3 d=+0.3 g=-0.1", S.impurity_caps("AB", 3, +0.3, -0.1)))
    errs.append(check("BA cell3 d=-0.3 g=-0.1", S.impurity_caps("BA", 3, -0.3, -0.1)))
    errs.append(check("AB pos N=2 d=g=0.3",     S.impurity_caps("AB", 2, 0.3, 0.3)))
    errs.append(check("A single cell3 d=-0.3",  S.single_caps("A", 3, -0.3)))
    errs.append(check("B single cell3 d=+0.3",  S.single_caps("B", 3, +0.3)))

    print("\n=== pristine physical sanity ===")
    caps = S.pristine_caps()
    freqs, V, _ = S.run_ngspice(caps, npts=4001)
    idx, fpk = S.find_resonance(freqs, V)
    prof = S.voltage_profile_odd(V, idx)
    ref = S.pristine_reference()
    # decay ratio per cell from a log-linear fit on odd nodes
    cells = np.arange(8)
    slope = np.polyfit(cells, np.log(prof), 1)[0]
    print(f"resonance peak       = {fpk/1e3:.4f} kHz  (analytic {S.F_R/1e3:.4f})")
    print(f"measured decay ratio = {np.exp(slope):.4f} /cell  (theory {S.C_WEAK/S.C_STRONG:.4f})")
    print(f"node-15 amplitude    = {prof[-1]:.4f}  (theory {ref[-1]:.4f})")

    worst = max(errs)
    print(f"\nWORST ngspice-vs-numpy rel error across all circuits = {worst:.2e}")
    print("PASS" if worst < 5e-3 else "FAIL — investigate")
