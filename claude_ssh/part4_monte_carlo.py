"""
Part 4 -- Monte Carlo disorder analysis.

Random, uncorrelated Gaussian disorder is applied to every component, many
times over, to turn the single-shot impurity studies of Parts 2-3 into
statistical distributions.  Three studies:

  (A) capacitor disorder   -- the realistic component-tolerance question
  (B) inductor  disorder   -- a chiral-symmetry-breaking control
  (C) trivial-phase circuit -- a no-protection control (C1 <-> C2 swapped)

Produces figures 9-16 of the report.
"""
import numpy as np
import ssh_circuit as ssh
from common import newfig, save, update_results, CLR, DATADIR, NETDIR
import os

print("=" * 70)
print("PART 4  --  MONTE CARLO DISORDER ANALYSIS")
print("=" * 70)

F_R = ssh.F_R
SIGMAS = [0.01, 0.05, 0.10]
M = 8000
MC_WIN = dict(f_lo=34e3, f_hi=62e3)
DELOC = 1.0 / 15.0                       # IPR of a state spread over N2..N16

c1_0, c2_0, cb0, lv0 = ssh.pristine_components()

# pristine reference (deterministic)
fp0, zp0 = ssh.find_tbr(c1_0, c2_0, cb0, lv0, **MC_WIN)
prof0 = ssh.voltage_profile(fp0, c1_0, c2_0, cb0, lv0)
ipr0 = ssh.ipr(prof0)
xi0 = ssh.localization_length(prof0)
xi_theory = 1.0 / np.log(ssh.C2 / ssh.C1)
print(f"\nPristine: f_peak = {fp0/1e3:.4f} kHz, |Z| = {zp0:,.0f}, "
      f"IPR = {ipr0:.4f}, xi = {xi0:.4f}")


def draw(rng, sigma, n):
    """Gaussian disorder multipliers, clipped to keep components positive."""
    return np.clip(rng.normal(0.0, sigma, n), -0.85, None)


def boot_std(x, nboot=500, seed=0):
    """Bootstrap uncertainty on the standard deviation of a sample."""
    rng = np.random.default_rng(seed)
    s = [np.std(rng.choice(x, size=len(x), replace=True)) for _ in range(nboot)]
    return float(np.std(s))


def run_mc(sigma, M, seed0, target="cap", trivial=False):
    """One Monte Carlo ensemble. target = 'cap' or 'ind'."""
    base_c1 = np.full(ssh.N_C1, ssh.C2 if trivial else ssh.C1)
    base_c2 = np.full(ssh.N_C2, ssh.C1 if trivial else ssh.C2)
    base_cb = ssh.C1 if trivial else ssh.C2
    fpk = np.empty(M); zpk = np.empty(M)
    iprv = np.empty(M); xiv = np.empty(M)
    prof = np.empty((M, 16))
    for r in range(M):
        rng = np.random.default_rng(seed0 + r)
        c1, c2, cb, lv = base_c1.copy(), base_c2.copy(), base_cb, lv0.copy()
        if target == "cap":
            c1 *= 1 + draw(rng, sigma, ssh.N_C1)
            c2 *= 1 + draw(rng, sigma, ssh.N_C2)
            cb *= 1 + draw(rng, sigma, 1)[0]
        else:                                    # inductor disorder
            lv = lv * (1 + draw(rng, sigma, ssh.N_NODES))
        fp, zp = ssh.find_tbr(c1, c2, cb, lv, **MC_WIN)
        p = ssh.voltage_profile(fp, c1, c2, cb, lv)
        fpk[r], zpk[r] = fp, zp
        iprv[r], xiv[r] = ssh.ipr(p), ssh.localization_length(p)
        prof[r] = p
    return dict(f=fpk, z=zpk, ipr=iprv, xi=xiv, prof=prof)


# ==========================================================================
# (A) Capacitor-disorder ensembles
# ==========================================================================
print("\n(A) capacitor disorder ...")
cap = {}
for s in SIGMAS:
    cap[s] = run_mc(s, M, seed0=10_000 + int(s * 1e5), target="cap")
    d = cap[s]
    print(f"  sigma={s*100:4.0f}%:  <f>={d['f'].mean()/1e3:.3f} kHz "
          f"(std {d['f'].std()/1e3:.3f}),  <IPR>={d['ipr'].mean():.4f},  "
          f"<xi>={d['xi'].mean():.3f}")

# theory check: standard deviation of f_peak grows linearly with sigma
std1 = cap[0.01]['f'].std()
lin_5 = cap[0.05]['f'].std() / (5 * std1)
lin_10 = cap[0.10]['f'].std() / (10 * std1)
print(f"  theory check -- std(f) linear in sigma:  "
      f"ratio@5% = {lin_5:.2f},  ratio@10% = {lin_10:.2f}  (expect ~1)")

# small systematic biases of the ensemble mean away from the pristine value
f_bias = {s: (cap[s]['f'].mean() - fp0) / fp0 * 100 for s in SIGMAS}
z_bias = {s: (cap[s]['z'].mean() - zp0) / zp0 * 100 for s in SIGMAS}
xi_bias = {s: (cap[s]['xi'].mean() - xi0) / xi0 * 100 for s in SIGMAS}
print(f"  ensemble-mean bias @ 10%:  <f> {f_bias[0.10]:+.2f}%,  "
      f"<|Z|> {z_bias[0.10]:+.2f}%,  <xi> {xi_bias[0.10]:+.2f}%")

# ==========================================================================
# (B) Inductor-disorder ensembles  (chiral-symmetry-breaking control)
# ==========================================================================
print("\n(B) inductor disorder ...")
ind = {}
for s in SIGMAS:
    ind[s] = run_mc(s, M, seed0=20_000 + int(s * 1e5), target="ind")
    d = ind[s]
    print(f"  sigma={s*100:4.0f}%:  <f>={d['f'].mean()/1e3:.3f} kHz "
          f"(std {d['f'].std()/1e3:.3f}),  <IPR>={d['ipr'].mean():.4f}")

# ==========================================================================
# Figure 9 -- f_peak distributions (capacitor disorder)
# ==========================================================================
allf = np.concatenate([cap[s]['f'] for s in SIGMAS]) / 1e3
bins = np.linspace(allf.min(), allf.max(), 48)
scol = {0.01: CLR["s1"], 0.05: CLR["s5"], 0.10: CLR["s10"]}
fig, ax = newfig()
for s in SIGMAS:
    fwhm = 2.355 * cap[s]['f'].std() / 1e3
    ax.hist(cap[s]['f'] / 1e3, bins=bins, alpha=0.6, color=scol[s],
            label=fr"$\sigma$ = {s*100:.0f}%   (FWHM {fwhm:.2f} kHz)")
ax.axvline(fp0 / 1e3, color="k", ls="--", lw=1.5,
           label=f"pristine $f_R$ = {fp0/1e3:.2f} kHz")
ax.set_xlabel("TBR peak frequency  $f_{\\rm peak}$  (kHz)")
ax.set_ylabel(f"count  (M = {M} per $\\sigma$)")
ax.set_title("Monte Carlo  --  TBR frequency distribution vs capacitor disorder")
ax.legend(fontsize=9)
save(fig, "fig09_mc_fpeak_hist.png")

# ==========================================================================
# Figure 10 -- f_peak mean +/- std vs sigma
# ==========================================================================
xs = np.array([0] + [s * 100 for s in SIGMAS])
fmean = np.array([fp0] + [cap[s]['f'].mean() for s in SIGMAS]) / 1e3
fstd = np.array([0] + [cap[s]['f'].std() for s in SIGMAS]) / 1e3
fig, ax = newfig()
ax.errorbar(xs, fmean, yerr=fstd, fmt="o-", color=CLR["tbr"], ms=8,
            capsize=6, lw=1.8, label=r"$\langle f_{\rm peak}\rangle \pm 1$ s.d.")
ax.axhline(fp0 / 1e3, color=CLR["theory"], ls="--", lw=1.5,
           label=f"pristine $f_R$ = {fp0/1e3:.3f} kHz")
ax.set_xlabel("capacitor disorder strength  $\\sigma$  (%)")
ax.set_ylabel("TBR peak frequency  (kHz)")
ax.set_title("Monte Carlo  --  TBR frequency: ensemble mean and spread vs $\\sigma$")
ax.legend(fontsize=9)
save(fig, "fig10_mc_fpeak_vs_sigma.png")

# ==========================================================================
# Figure 11 -- peak impedance vs sigma
# ==========================================================================
zmean = np.array([zp0] + [cap[s]['z'].mean() for s in SIGMAS]) / 1e3
zstd = np.array([0] + [cap[s]['z'].std() for s in SIGMAS]) / 1e3
fig, ax = newfig()
ax.errorbar(xs, zmean, yerr=zstd, fmt="s-", color=CLR["even"], ms=8,
            capsize=6, lw=1.8, label=r"$\langle |Z|_{\rm peak}\rangle \pm 1$ s.d.")
ax.axhline(zp0 / 1e3, color=CLR["theory"], ls="--", lw=1.5,
           label=f"pristine $|Z|_{{\\rm peak}}$ = {zp0/1e3:.1f} k$\\Omega$")
ax.set_xlabel("capacitor disorder strength  $\\sigma$  (%)")
ax.set_ylabel(r"peak impedance  $|Z|_{\rm peak}$  (k$\Omega$)")
ax.set_title("Monte Carlo  --  resonance strength (peak impedance) vs $\\sigma$")
ax.legend(fontsize=9)
save(fig, "fig11_mc_zpeak_vs_sigma.png")

# ==========================================================================
# Figure 12 -- IPR distributions (localization quality)
# ==========================================================================
alli = np.concatenate([cap[s]['ipr'] for s in SIGMAS])
ibins = np.linspace(alli.min(), max(alli.max(), ipr0) * 1.02, 55)
fig, ax = newfig()
for s in SIGMAS:
    ax.hist(cap[s]['ipr'], bins=ibins, alpha=0.6, color=scol[s],
            label=fr"$\sigma$ = {s*100:.0f}%")
ax.axvline(ipr0, color="k", ls="--", lw=1.5, label=f"pristine IPR = {ipr0:.3f}")
ax.axvline(DELOC, color=CLR["warn"], ls=":", lw=1.8,
           label=f"fully delocalised = 1/15 = {DELOC:.3f}")
ax.set_xlabel("inverse participation ratio  (nodes N2..N16)")
ax.set_ylabel(f"count  (M = {M} per $\\sigma$)")
ax.set_title("Monte Carlo  --  edge-state localization (IPR) vs capacitor disorder")
ax.legend(fontsize=9)
save(fig, "fig12_mc_ipr_hist.png")

# ==========================================================================
# Figure 13 -- localization length vs sigma
# ==========================================================================
xmean = np.array([xi0] + [cap[s]['xi'].mean() for s in SIGMAS])
xstd = np.array([0] + [cap[s]['xi'].std() for s in SIGMAS])
fig, ax = newfig()
ax.errorbar(xs, xmean, yerr=xstd, fmt="D-", color=CLR["tbr"], ms=7,
            capsize=6, lw=1.8, label=r"$\langle \xi\rangle \pm 1$ s.d.")
ax.axhline(xi_theory, color=CLR["theory"], ls="--", lw=1.5,
           label=fr"theory  $\xi = 1/\ln(C_2/C_1)$ = {xi_theory:.3f} cells")
ax.set_xlabel("capacitor disorder strength  $\\sigma$  (%)")
ax.set_ylabel(r"localization length  $\xi$  (unit cells)")
ax.set_title("Monte Carlo  --  edge-state decay length vs $\\sigma$")
ax.legend(fontsize=9)
save(fig, "fig13_mc_xi_vs_sigma.png")

# ==========================================================================
# Figure 14 -- voltage-profile envelope  (A sublattice, no zig-zag)
# ==========================================================================
nodes = np.arange(1, 17)
odd_n = nodes[::2]
oi = np.arange(0, 16, 2)                       # odd-node (A-sublattice) indices
fig, ax = newfig()
for s in SIGMAS:
    P = cap[s]['prof'][:, oi]
    med = np.median(P, axis=0)
    p5, p95 = np.percentile(P, 5, axis=0), np.percentile(P, 95, axis=0)
    ax.fill_between(odd_n, p5, p95, color=scol[s], alpha=0.22)
    ax.semilogy(odd_n, med, "-", color=scol[s], lw=1.7,
                label=fr"$\sigma$ = {s*100:.0f}%  (median + 5-95% band)")
ax.semilogy(odd_n, prof0[oi], "o-", color="k", ms=6, lw=2.0,
            label="pristine", zorder=6)
# relative band half-width at a representative mid-chain node (N9)
bw10 = (np.percentile(cap[0.10]['prof'][:, 8], 95) -
        np.percentile(cap[0.10]['prof'][:, 8], 5)) / \
       (2 * np.median(cap[0.10]['prof'][:, 8])) * 100
print(f"\nFig 14: sigma=10% 5-95% band half-width at node N9 = +/-{bw10:.0f}%")
ax.set_xlabel("node index  $i$   (A sublattice = odd nodes)")
ax.set_ylabel(r"normalised voltage  $|V(N_i)|\,/\,V_0$")
ax.set_title("Monte Carlo  --  edge-state (A-sublattice) voltage envelope")
ax.set_xticks(odd_n)
ax.legend(fontsize=8.5)
save(fig, "fig14_mc_voltage_envelope.png")

# ==========================================================================
# Figure 15 -- capacitor vs inductor disorder (chiral-symmetry control)
# ==========================================================================
# resonance-frequency spread (% of f_R) -- the clearest chiral signature
spr_cap = [0.0] + [cap[s]['f'].std() / F_R * 100 for s in SIGMAS]
spr_cap_e = [0.0] + [boot_std(cap[s]['f'], seed=1) / F_R * 100 for s in SIGMAS]
spr_ind = [0.0] + [ind[s]['f'].std() / F_R * 100 for s in SIGMAS]
spr_ind_e = [0.0] + [boot_std(ind[s]['f'], seed=2) / F_R * 100 for s in SIGMAS]
ipr_cap = [ipr0] + [cap[s]['ipr'].mean() for s in SIGMAS]
ipr_ind = [ipr0] + [ind[s]['ipr'].mean() for s in SIGMAS]
ratio_10 = spr_ind[-1] / spr_cap[-1]
print(f"\n(B) chiral control -- f-spread at 10%: capacitor {spr_cap[-1]:.2f}% "
      f"vs inductor {spr_ind[-1]:.2f}%  (inductor x{ratio_10:.2f} wider)")
print(f"    mean IPR at 10%: capacitor {ipr_cap[-1]:.4f} vs inductor {ipr_ind[-1]:.4f}")
fig, ax = newfig()
ax.errorbar(xs, spr_cap, yerr=spr_cap_e, fmt="o-", color=CLR["tbr"], ms=8,
            capsize=6, lw=1.8,
            label="capacitor disorder (off-diagonal, chiral-symmetric)")
ax.errorbar(xs, spr_ind, yerr=spr_ind_e, fmt="s-", color=CLR["warn"], ms=8,
            capsize=6, lw=1.8,
            label="inductor disorder (on-site, chiral-breaking)")
ax.set_xlabel("disorder strength  $\\sigma$  (%)")
ax.set_ylabel(r"resonance-frequency spread  $\sigma_f / f_R$  (%)")
ax.set_title("Monte Carlo  --  chiral-symmetry test: capacitor vs inductor disorder")
ax.legend(fontsize=9)
save(fig, "fig15_mc_chiral_control.png")

# ==========================================================================
# (C) Trivial-phase control -- impedance ensembles at sigma = 10%
# ==========================================================================
print("\n(C) trivial-phase control ...")
fgrid = np.linspace(30e3, 70e3, 700)
M_C = 400


def imp_ensemble(trivial, sigma, M, seed0):
    base_c1 = np.full(ssh.N_C1, ssh.C2 if trivial else ssh.C1)
    base_c2 = np.full(ssh.N_C2, ssh.C1 if trivial else ssh.C2)
    base_cb = ssh.C1 if trivial else ssh.C2
    Z = np.empty((M, fgrid.size))
    for r in range(M):
        rng = np.random.default_rng(seed0 + r)
        c1 = base_c1 * (1 + draw(rng, sigma, ssh.N_C1))
        c2 = base_c2 * (1 + draw(rng, sigma, ssh.N_C2))
        cb = base_cb * (1 + draw(rng, sigma, 1)[0])
        Z[r] = ssh.impedance_curve(fgrid, c1, c2, cb, lv0)
    return Z

z_topo0 = ssh.impedance_curve(fgrid, c1_0, c2_0, cb0, lv0)
z_triv0 = ssh.impedance_curve(fgrid, np.full(8, ssh.C2), np.full(7, ssh.C1),
                              ssh.C1, lv0)
Z_topo = imp_ensemble(False, 0.10, M_C, seed0=30_000)
Z_triv = imp_ensemble(True, 0.10, M_C, seed0=40_000)
print(f"  pristine peak |Z|:  topological {z_topo0.max():,.0f}  "
      f"vs  trivial {z_triv0.max():,.0f}")

fig, ax = newfig()
for Z, z0, col, name in [(Z_topo, z_topo0, CLR["tbr"], "topological  ($C_1<C_2$)"),
                         (Z_triv, z_triv0, CLR["s5"], "trivial  ($C_1>C_2$)")]:
    med = np.median(Z, axis=0)
    p5, p95 = np.percentile(Z, 5, axis=0), np.percentile(Z, 95, axis=0)
    ax.fill_between(fgrid / 1e3, p5, p95, color=col, alpha=0.14)
    ax.semilogy(fgrid / 1e3, z0, ":", color=col, lw=1.3, alpha=0.9)
    ax.semilogy(fgrid / 1e3, med, "-", color=col, lw=2.1,
                label=f"{name}  --  10% median")
ax.plot([], [], "k:", lw=1.3, label="pristine spectra (no disorder)")
ax.axvline(F_R / 1e3, color="k", ls="--", lw=1.2,
           label=f"gap centre = {F_R/1e3:.2f} kHz")
ax.set_xlabel("frequency  (kHz)")
ax.set_ylabel(r"input impedance  $|Z_{\rm in}|$  ($\Omega$)")
ax.set_title("Monte Carlo  --  topological vs trivial circuit at 10% disorder")
ax.legend(fontsize=9)
save(fig, "fig16_mc_topo_vs_trivial.png")

# --------------------------------------------------------------------------
# Data files
# --------------------------------------------------------------------------
with open(os.path.join(DATADIR, "part4_mc_results.csv"), "w") as fh:
    fh.write("disorder,sigma,realization,f_peak_Hz,Z_peak_ohm,IPR,xi_cells\n")
    for label, ens in [("capacitor", cap), ("inductor", ind)]:
        for s in SIGMAS:
            d = ens[s]
            for r in range(M):
                fh.write(f"{label},{s},{r},{d['f'][r]:.3f},{d['z'][r]:.3f},"
                         f"{d['ipr'][r]:.6f},{d['xi'][r]:.6f}\n")

with open(os.path.join(DATADIR, "part4_mc_summary.csv"), "w") as fh:
    fh.write("disorder,sigma,f_mean_kHz,f_std_kHz,Z_mean,Z_std,"
             "IPR_mean,IPR_std,xi_mean,xi_std\n")
    fh.write(f"pristine,0,{fp0/1e3:.4f},0,{zp0:.1f},0,{ipr0:.4f},0,{xi0:.4f},0\n")
    for label, ens in [("capacitor", cap), ("inductor", ind)]:
        for s in SIGMAS:
            d = ens[s]
            fh.write(f"{label},{s},{d['f'].mean()/1e3:.4f},{d['f'].std()/1e3:.4f},"
                     f"{d['z'].mean():.1f},{d['z'].std():.1f},"
                     f"{d['ipr'].mean():.4f},{d['ipr'].std():.4f},"
                     f"{d['xi'].mean():.4f},{d['xi'].std():.4f}\n")

update_results("part4", dict(
    M=M, sigmas=SIGMAS,
    pristine=dict(f_peak=fp0, z_peak=zp0, ipr=ipr0, xi=xi0),
    xi_theory=xi_theory, deloc=DELOC,
    std_linearity=dict(ratio5=lin_5, ratio10=lin_10),
    bias={str(s): dict(f=f_bias[s], z=z_bias[s], xi=xi_bias[s]) for s in SIGMAS},
    band_halfwidth_N9=float(bw10),
    cap_summary={str(s): dict(f_mean=cap[s]['f'].mean(), f_std=cap[s]['f'].std(),
                              z_mean=cap[s]['z'].mean(), z_std=cap[s]['z'].std(),
                              ipr_mean=cap[s]['ipr'].mean(), ipr_std=cap[s]['ipr'].std(),
                              xi_mean=cap[s]['xi'].mean(), xi_std=cap[s]['xi'].std(),
                              f_min=cap[s]['f'].min(), f_max=cap[s]['f'].max())
                 for s in SIGMAS},
    ind_summary={str(s): dict(f_mean=ind[s]['f'].mean(), f_std=ind[s]['f'].std(),
                              ipr_mean=ind[s]['ipr'].mean(), ipr_std=ind[s]['ipr'].std())
                 for s in SIGMAS},
    ipr_cap=ipr_cap, ipr_ind=ipr_ind,
    chiral=dict(spr_cap=spr_cap, spr_ind=spr_ind, ratio10=ratio_10),
    trivial=dict(topo_peak=float(z_topo0.max()), triv_peak=float(z_triv0.max())),
))
print("\nPart 4 complete.\n")
