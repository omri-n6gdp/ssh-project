"""
Part 3 -- Domain-wall regime (AB and BA impurity pairs).

Two adjacent impurity sites are introduced at cell 3.  Following Anzurez
et al. (2025, Table 1):

  AB pair:  C2_left  -> C2(1+delta)
            C1_mid   -> C1(1 + (delta+gamma)/2)
            C2_right -> C2(1+gamma)
  BA pair:  C1_left  -> C1(1+delta)
            C2_mid   -> C2(1 + (delta+gamma)/2)
            C1_right -> C1(1+gamma)

We sweep gamma at fixed delta, track the TBR, and inspect the voltage
profile for the emergent domain-wall state.  Produces figures 6-8.
"""
import numpy as np
import ssh_circuit as ssh
from common import newfig, save, update_results, CLR, DATADIR, NETDIR
import os

print("=" * 70)
print("PART 3  --  DOMAIN-WALL REGIME (AB / BA PAIRS)")
print("=" * 70)

c1_0, c2_0, cb, lvals = ssh.pristine_components()
F_R = ssh.F_R

# pristine TBR eigenvector for first-order theory
W0 = ssh.cap_laplacian(c1_0, c2_0, cb)
mu, vecs = np.linalg.eigh(W0)
mid = int(np.argmin(np.abs(mu - (ssh.C1 + ssh.C2))))
psi, mu_tbr = vecs[:, mid], mu[mid]


def make_AB(delta, gamma):
    c1, c2 = c1_0.copy(), c2_0.copy()
    c2[1] *= (1 + delta)                          # C2 left  (nodes 4-5)
    c1[2] *= (1 + (delta + gamma) / 2)            # C1 mid   (nodes 5-6)
    c2[2] *= (1 + gamma)                          # C2 right (nodes 6-7)
    return c1, c2


def make_BA(delta, gamma):
    c1, c2 = c1_0.copy(), c2_0.copy()
    c1[2] *= (1 + delta)                          # C1 left  (nodes 5-6)
    c2[2] *= (1 + (delta + gamma) / 2)            # C2 mid   (nodes 6-7)
    c1[3] *= (1 + gamma)                          # C1 right (nodes 7-8)
    return c1, c2


def theory_AB(delta, gamma):
    dmu = (delta * ssh.C2 * (psi[3] - psi[4])**2 +
           (delta + gamma) / 2 * ssh.C1 * (psi[4] - psi[5])**2 +
           gamma * ssh.C2 * (psi[5] - psi[6])**2)
    return 1.0 - dmu / (2 * mu_tbr)


def theory_BA(delta, gamma):
    dmu = (delta * ssh.C1 * (psi[4] - psi[5])**2 +
           (delta + gamma) / 2 * ssh.C2 * (psi[5] - psi[6])**2 +
           gamma * ssh.C1 * (psi[6] - psi[7])**2)
    return 1.0 - dmu / (2 * mu_tbr)


# --------------------------------------------------------------------------
# gamma sweeps at fixed delta  (Figures 6 and 7)
# --------------------------------------------------------------------------
GAMMAS = np.linspace(-0.3, 0.3, 25)
DELTAS = [-0.2, 0.0, 0.2]
dcolors = {-0.2: CLR["theory"], 0.0: CLR["bulk"], 0.2: CLR["tbr"]}

def exact_ratio(c1, c2):
    """Resonance ratio from the exact lossless capacitance-Laplacian eigenvalue."""
    fr = ssh.resonances_from_W(ssh.cap_laplacian(c1, c2, cb))
    return fr[np.argmin(np.abs(fr - F_R))] / F_R

p3 = {}
for tag, maker, theory in [("AB", make_AB, theory_AB),
                           ("BA", make_BA, theory_BA)]:
    fig, ax = newfig()
    store, dev_th, dev_exact = {}, [], []
    for d in DELTAS:
        sim, th, exa = [], [], []
        for g in GAMMAS:
            c1, c2 = maker(d, g)
            fp, _ = ssh.find_tbr(c1, c2, cb, lvals)
            sim.append(fp / F_R)
            th.append(theory(d, g))
            exa.append(exact_ratio(c1, c2))
        store[str(d)] = sim
        dev_th += [abs(s - t) for s, t in zip(sim, th)]
        dev_exact += [abs(s - e) for s, e in zip(sim, exa)]
        ax.plot(GAMMAS, sim, "o-", color=dcolors[d], ms=5,
                label=fr"simulation, $\delta$ = {d:+.1f}")
        ax.plot(GAMMAS, th, "--", color=dcolors[d], lw=1.2, alpha=0.85)
        print(f"  {tag} delta={d:+.1f}:  f/f_R range "
              f"{min(sim):.4f} .. {max(sim):.4f}")
    ax.axhline(1.0, color="k", lw=1.0, alpha=0.6)
    ax.plot([], [], "k--", lw=1.2, label="first-order perturbation theory")
    ax.set_xlabel(r"second impurity strength  $\gamma$")
    ax.set_ylabel(r"resonance ratio  $f_{\rm peak}\,/\,f_R$")
    ax.set_title(f"{tag} impurity pair at cell 3  --  TBR displacement vs "
                 r"$\gamma$")
    ax.legend(fontsize=9)
    fno = "fig06" if tag == "AB" else "fig07"
    save(fig, f"{fno}_{tag}_pair_sweep.png")

    # ---- quantitative analysis ----
    s0 = np.array(store["0.0"])
    curv = float(np.polyfit(GAMMAS, s0, 2)[0])          # curvature of delta=0 curve
    ip = int(np.argmin(abs(GAMMAS - 0.3)))
    inj = int(np.argmin(abs(GAMMAS + 0.3)))
    asym_g = abs(s0[inj] - 1) / abs(s0[ip] - 1)         # |-gamma| vs |+gamma|
    span = (max(max(store[str(d)]) for d in DELTAS) -
            min(min(store[str(d)]) for d in DELTAS))
    p3[tag] = dict(curv=curv, asym_g=asym_g, span=span,
                   dev_th=max(dev_th), dev_exact=max(dev_exact))
    print(f"  {tag} analysis: curvature={curv:+.4f}, gamma-asymmetry "
          f"x{asym_g:.2f}, span={span*100:.2f}%; max deviation vs first-order "
          f"theory={max(dev_th)*100:.2f}%, vs exact eigenvalue="
          f"{max(dev_exact)*100:.3f}%")

    with open(os.path.join(DATADIR, f"part3_{tag}_sweep.csv"), "w") as fh:
        fh.write("gamma," + ",".join(f"delta_{d}" for d in DELTAS) + "\n")
        for i, g in enumerate(GAMMAS):
            fh.write(f"{g:.4f}," + ",".join(f"{store[str(d)][i]:.6f}"
                                            for d in DELTAS) + "\n")
    update_results(f"part3_{tag}", dict(gammas=GAMMAS.tolist(),
                                        deltas=DELTAS, ratio=store, **p3[tag]))

# --------------------------------------------------------------------------
# Domain-wall voltage profiles  (Figure 8)
# --------------------------------------------------------------------------
nodes = np.arange(1, 17)
prof_pris = ssh.voltage_profile(ssh.find_tbr(c1_0, c2_0, cb, lvals)[0],
                                c1_0, c2_0, cb, lvals)

dw_delta, dw_gamma = 0.3, 0.3
c1_ab, c2_ab = make_AB(dw_delta, dw_gamma)
c1_ba, c2_ba = make_BA(dw_delta, dw_gamma)
fp_ab, _ = ssh.find_tbr(c1_ab, c2_ab, cb, lvals)
fp_ba, _ = ssh.find_tbr(c1_ba, c2_ba, cb, lvals)
prof_ab = ssh.voltage_profile(fp_ab, c1_ab, c2_ab, cb, lvals)
prof_ba = ssh.voltage_profile(fp_ba, c1_ba, c2_ba, cb, lvals)

# even-node (right-edge channel) IPR -- a "uniformity" check AB vs BA
even = np.arange(1, 16, 2)
ipr_ab = ssh.ipr(prof_ab, nodes=even)
ipr_ba = ssh.ipr(prof_ba, nodes=even)
print(f"\nDomain wall (delta=gamma={dw_delta}):")
print(f"  even-node IPR  AB = {ipr_ab:.4f}   BA = {ipr_ba:.4f}")
print(f"  -> BA more uniform than AB: {ipr_ba < ipr_ab}")

odd_n, even_n = nodes[::2], nodes[1::2]
ev_pris = float(prof_pris[1::2].max())
ev_ab = float(prof_ab[1::2].max())
ev_ba = float(prof_ba[1::2].max())
print(f"  B-sublattice (even-node) max |V|: pristine {ev_pris:.3f}, "
      f"AB {ev_ab:.3f}  (x{ev_ab/ev_pris:.0f}), BA {ev_ba:.3f}  (x{ev_ba/ev_pris:.0f})")

fig, ax = newfig()
for prof, name, col in [(prof_pris, "pristine", CLR["bulk"]),
                        (prof_ab, fr"AB pair ($\delta=\gamma=0.3$)", CLR["theory"]),
                        (prof_ba, fr"BA pair ($\delta=\gamma=0.3$)", CLR["tbr"])]:
    ax.semilogy(odd_n, prof[::2], "o-", color=col, ms=6, label=name)
    ax.semilogy(even_n, prof[1::2], "s--", color=col, ms=5, alpha=0.8,
                label="_nolegend_")
ax.axvspan(4.5, 8.5, color=CLR["warn"], alpha=0.15, label="impurity-pair region")
ax.text(2.2, 0.60, "A sublattice (odd nodes) -- solid", fontsize=8.5)
ax.text(2.2, 1.4e-3, "B sublattice (even nodes) -- dashed", fontsize=8.5)
ax.set_xlabel("node index  $i$")
ax.set_ylabel(r"normalised voltage  $|V(N_i)|\,/\,V_0$")
ax.set_title("Domain wall  --  AB and BA impurity pairs (sublattices separated)")
ax.set_xticks(nodes)
ax.set_ylim(8e-4, 2.4)
ax.legend(fontsize=8.5, loc="upper right")
save(fig, "fig08_domain_wall_profile.png")

update_results("part3_dw", dict(delta=dw_delta, gamma=dw_gamma,
                                ipr_even_AB=ipr_ab, ipr_even_BA=ipr_ba,
                                ev_pris=ev_pris, ev_ab=ev_ab, ev_ba=ev_ba,
                                profile_pristine=prof_pris.tolist(),
                                profile_AB=prof_ab.tolist(),
                                profile_BA=prof_ba.tolist()))

ssh.write_netlist(os.path.join(NETDIR, "ssh_AB_pair_cell3.cir"),
                  c1_ab, c2_ab, cb, lvals, title="SSH AB pair cell 3, delta=gamma=0.3")
ssh.write_netlist(os.path.join(NETDIR, "ssh_BA_pair_cell3.cir"),
                  c1_ba, c2_ba, cb, lvals, title="SSH BA pair cell 3, delta=gamma=0.3")
print("\nPart 3 complete.\n")
