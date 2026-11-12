"""
Part 2 -- Single-impurity parametric sweeps.

A single unit cell n_i is made "impure": both of its capacitors are scaled
by (1 + delta).  We sweep the impurity over cells 2..7 and the strength over
delta = +-0.1, +-0.3, track the displacement of the TBR, and compare against
first-order perturbation theory.  Produces figures 4-5 of the report.
"""
import numpy as np
import ssh_circuit as ssh
from common import newfig, save, update_results, CLR, DATADIR, NETDIR
import os

print("=" * 70)
print("PART 2  --  SINGLE-IMPURITY SWEEPS")
print("=" * 70)

c1_0, c2_0, cb, lvals = ssh.pristine_components()
F_R = ssh.F_R
DELTAS = [-0.3, -0.1, 0.1, 0.3]
CELLS = list(range(2, 8))                       # cells 2..7

# --------------------------------------------------------------------------
# First-order perturbation theory
# Scaling cell n_i's capacitors by (1+delta) changes the TBR eigenvalue mu
# of the capacitance Laplacian by  delta * sum_{cell caps} C (psi_a-psi_b)^2.
# Since f ~ 1/sqrt(mu),  f/f_R ~= 1 - (delta/2) * (d mu / mu).
# --------------------------------------------------------------------------
W0 = ssh.cap_laplacian(c1_0, c2_0, cb)
mu, vecs = np.linalg.eigh(W0)
mid = int(np.argmin(np.abs(mu - (ssh.C1 + ssh.C2))))
psi = vecs[:, mid]
mu_tbr = mu[mid]

def cell_response(cell):
    """Fractional eigenvalue response r = d mu / mu for unit scaling of cell."""
    dmu = 0.0
    a, b = 2*cell - 2, 2*cell - 1               # C1 of cell: nodes (2c-1,2c)
    dmu += ssh.C1 * (psi[a] - psi[b])**2
    if cell <= 7:
        a, b = 2*cell - 1, 2*cell               # C2 of cell: nodes (2c,2c+1)
        dmu += ssh.C2 * (psi[a] - psi[b])**2
    return dmu / mu_tbr

r_cell = {c: cell_response(c) for c in CELLS}

# --------------------------------------------------------------------------
# Simulation sweep: cell position x delta
# --------------------------------------------------------------------------
ratio = {d: [] for d in DELTAS}
for d in DELTAS:
    for cell in CELLS:
        c1 = c1_0.copy(); c2 = c2_0.copy()
        c1[cell-1] *= (1 + d)
        if cell <= 7:
            c2[cell-1] *= (1 + d)
        f_peak, _ = ssh.find_tbr(c1, c2, cb, lvals)
        ratio[d].append(f_peak / F_R)
    print(f"  delta = {d:+.1f}:  f/f_R over cells 2..7 = "
          + ", ".join(f"{x:.4f}" for x in ratio[d]))

# theory check: sign of the shift
allpos = all(np.array(ratio[d])[0] < 1 for d in DELTAS if d > 0)
allneg = all(np.array(ratio[d])[0] > 1 for d in DELTAS if d < 0)
print(f"\nTheory check -- positive delta lowers f_R: {allpos}")
print(f"Theory check -- negative delta raises f_R: {allneg}")
print(f"Theory check -- cell 2 shifts most, cell 7 least: "
      f"{r_cell[2] > r_cell[7]}  (r2={r_cell[2]:.3f}, r7={r_cell[7]:.4f})")

# --------------------------------------------------------------------------
# Figure 4 -- TBR displacement vs impurity position
# --------------------------------------------------------------------------
fig, ax = newfig()
colors = {-0.3: CLR["theory"], -0.1: CLR["s5"], 0.1: CLR["even"], 0.3: CLR["tbr"]}
for d in DELTAS:
    ax.plot(CELLS, ratio[d], "o-", color=colors[d], ms=7,
            label=fr"simulation, $\delta$ = {d:+.1f}")
    theory = [1 - 0.5 * d * r_cell[c] for c in CELLS]
    ax.plot(CELLS, theory, "--", color=colors[d], lw=1.2, alpha=0.8)
ax.axhline(1.0, color="k", lw=1.0, alpha=0.6)
ax.plot([], [], "k--", lw=1.2, label="first-order perturbation theory")
ax.set_xlabel("impurity unit cell  $n_i$")
ax.set_ylabel(r"resonance ratio  $f_{\rm peak}\,/\,f_R$")
ax.set_title("Single impurity  --  TBR displacement vs position and strength")
ax.set_xticks(CELLS)
ax.legend(fontsize=8.5, ncol=2, loc="best")
save(fig, "fig04_single_impurity_shift.png")

with open(os.path.join(DATADIR, "part2_impurity_shift.csv"), "w") as fh:
    fh.write("delta,cell,f_over_fR,theory_f_over_fR\n")
    for d in DELTAS:
        for i, c in enumerate(CELLS):
            fh.write(f"{d},{c},{ratio[d][i]:.6f},{1-0.5*d*r_cell[c]:.6f}\n")

# --------------------------------------------------------------------------
# Quantitative analysis of Figure 4
# --------------------------------------------------------------------------
shift_neg = (ratio[-0.3][0] - 1) * 100          # cell 2, delta = -0.3
shift_pos = (1 - ratio[0.3][0]) * 100           # cell 2, delta = +0.3
asym = shift_neg / shift_pos
leverage = r_cell[2] / r_cell[7]
print(f"\nFig 4 analysis:")
print(f"  largest shifts (cell 2): delta=-0.3 -> +{shift_neg:.2f}%, "
      f"delta=+0.3 -> -{shift_pos:.2f}%   (asymmetry x{asym:.3f})")
print(f"  position leverage r(cell2)/r(cell7) = {leverage:.0f}")

# --------------------------------------------------------------------------
# Figure 5 -- voltage profiles for an impurity in cell 3
# Odd (A) and even (B) sublattices drawn as SEPARATE curves -- no zig-zag.
# --------------------------------------------------------------------------
nodes = np.arange(1, 17)
odd, even = nodes[::2], nodes[1::2]
fp0, _ = ssh.find_tbr(c1_0, c2_0, cb, lvals)
prof_pris = ssh.voltage_profile(fp0, c1_0, c2_0, cb, lvals)

configs = [("pristine", 0.0, CLR["bulk"]),
           (r"impurity $\delta=+0.3$", 0.3, CLR["tbr"]),
           (r"impurity $\delta=-0.3$", -0.3, CLR["theory"])]
fig, ax = newfig()
prof_store, devs, prof_by_d = {}, {}, {}
for name, d, col in configs:
    c1, c2 = c1_0.copy(), c2_0.copy()
    if d != 0:
        c1[2] *= (1 + d); c2[2] *= (1 + d)              # cell 3
    fp = fp0 if d == 0 else ssh.find_tbr(c1, c2, cb, lvals)[0]
    prof = ssh.voltage_profile(fp, c1, c2, cb, lvals)
    ax.semilogy(odd, prof[::2], "o-", color=col, ms=6, label=name)
    ax.semilogy(even, prof[1::2], "s--", color=col, ms=5, alpha=0.8,
                label="_nolegend_")
    prof_store[name] = prof.tolist()
    prof_by_d[d] = prof
    if d != 0:
        devs[d] = np.abs(prof - prof_pris) / prof_pris * 100

ax.axvspan(4.5, 7.5, color=CLR["warn"], alpha=0.15, label="impurity cell 3")
ax.text(2.4, 0.62, "A sublattice (odd nodes) -- solid", fontsize=8.5)
ax.text(8.6, 4.5e-3, "B sublattice (even nodes) -- dashed", fontsize=8.5)
ax.set_xlabel("node index  $i$")
ax.set_ylabel(r"normalised voltage  $|V(N_i)|\,/\,V_0$")
ax.set_title("Single impurity in cell 3  --  voltage profile (sublattices separated)")
ax.set_xticks(nodes)
ax.set_ylim(2e-3, 2.4)
ax.legend(fontsize=8.5, loc="upper right")
save(fig, "fig05_single_impurity_profile.png")

# A-sublattice (dominant edge mode): upstream nodes 1,3 vs downstream 5..15
odd_i = np.arange(0, 16, 2)
up_A = {d: float(devs[d][odd_i[:2]].max()) for d in devs}
down_A = {d: float(devs[d][odd_i[2:]].max()) for d in devs}
even_pris = float(prof_pris[1::2].max())
even_imp = {d: float(prof_by_d[d][1::2].max()) for d in devs}
print("Fig 5 analysis (A-sublattice = dominant edge mode):")
for d in devs:
    print(f"  delta={d:+.1f}:  A upstream (N1,N3) <= {up_A[d]:.2f}%,  "
          f"A downstream (N5-N15) up to {down_A[d]:.1f}%")
print(f"  B-sublattice: pristine max |V|={even_pris:.3f}, "
      f"lifted near defect to {max(even_imp.values()):.3f}  (still << A)")

update_results("part2", dict(
    deltas=DELTAS, cells=CELLS,
    ratio={str(d): ratio[d] for d in DELTAS},
    r_cell={str(c): r_cell[c] for c in CELLS},
    sign_check_pos=bool(allpos), sign_check_neg=bool(allneg),
    shift_neg=shift_neg, shift_pos=shift_pos, asym=asym, leverage=leverage,
    dev_up_A={str(d): up_A[d] for d in devs},
    dev_down_A={str(d): down_A[d] for d in devs},
    even_pris_max=even_pris, even_imp_max=max(even_imp.values()),
    profiles=prof_store,
))

# netlist of one representative impurity configuration
c1 = c1_0.copy(); c2 = c2_0.copy(); c1[2] *= 1.3; c2[2] *= 1.3
ssh.write_netlist(os.path.join(NETDIR, "ssh_impurity_cell3.cir"),
                  c1, c2, cb, lvals, title="SSH single impurity cell 3, delta=+0.3")
print("\nPart 2 complete.\n")
