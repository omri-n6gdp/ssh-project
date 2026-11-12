"""
Automated single-impurity sweep over impurity cell position and strength delta.
Replicates Anzurez et al. (2025) Figs 4-6 (A-impurity case).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyLTSpice import SimRunner, SpiceEditor, RawRead

# ---- Configuration ----
LTSPICE_PATH = "/Applications/LTspice.app/Contents/MacOS/LTspice"
NETLIST = os.path.expanduser("~/code/ssh_project/netlists/ssh_pristine.cir")
OUTDIR = os.path.expanduser("~/code/ssh_project/raw")
N_CELLS = 8
DELTAS = [-0.3, -0.1, 0.0, 0.1, 0.3]
CELLS = list(range(2, 8))  # cells 2 through 7
F_R = 189685.0  # pristine resonance, Hz
C1_BASE = 0.1e-9
C2_BASE = 0.22e-9

# Map cell index -> capacitor designators in the netlist
# (you must edit the pristine .cir to use consistent names: C1_1...C1_8, C2_1...C2_7)
def c1_name(cell):
    return f"C1_{cell}"
def c2_name(cell):
    return f"C2_{cell}"

results = []

for cell in CELLS:
    for delta in DELTAS:
        editor = SpiceEditor(NETLIST)
        # Replace the two capacitors of the impurity cell
        editor.set_component_value(c1_name(cell), C1_BASE * (1 + delta))
        editor.set_component_value(c2_name(cell), C2_BASE * (1 + delta))
        # Save modified netlist to a uniquely-named file
        out_cir = os.path.join(OUTDIR, f"impA_cell{cell}_d{delta:+.2f}.cir")
        editor.save_netlist(out_cir)

        # Run LTSpice in batch
        runner = SimRunner(simulator=LTSPICE_PATH, output_folder=OUTDIR)
        raw, log = runner.run_now(out_cir)

        # Read raw, find peak frequency
        r = RawRead(raw)
        freqs = r.get_trace('frequency').get_wave()
        v_n1 = r.get_trace('V(n1)').get_wave()
        i_v1 = r.get_trace('I(V1)').get_wave()
        Z = np.abs(v_n1 / (-i_v1))
        idx = np.argmax(Z)
        f_peak = freqs[idx].real  # LTSpice AC gives complex frequencies
        results.append({'cell': cell, 'delta': delta,
                        'f_peak': f_peak, 'f_ratio': f_peak / F_R})
        print(f"Cell {cell}, delta={delta:+.2f}: f_peak = {f_peak/1e3:.2f} kHz, f/fR = {f_peak/F_R:.4f}")

df = pd.DataFrame(results)
df.to_csv(os.path.expanduser("~/code/ssh_project/figures/single_impurity_A_sweep.csv"), index=False)

# Plot like Anzurez Fig 4 upper-left
fig, ax = plt.subplots(figsize=(7, 5))
for delta in DELTAS:
    sub = df[df.delta == delta]
    ax.plot(sub.cell, sub.f_ratio, marker='o', label=f"δ = {delta:+.1f}")
ax.axhline(1.0, color='red', linestyle='--', alpha=0.6, label='pristine')
ax.set_xlabel('Impurity cell n_i')
ax.set_ylabel(r'$f / f_R$')
ax.set_title('A-impurity resonance displacement (replication of Anzurez 2025 Fig 4 UL)')
ax.set_ylim(0.85, 1.15)
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.expanduser("~/code/ssh_project/figures/single_impurity_A_sweep.png"), dpi=120)
plt.show()