"""
Part 1 -- Pristine baseline reproduction.

Builds the pristine SSH topolectrical circuit with the updated component
values, locates the topological boundary resonance (TBR), measures the
edge-state voltage profile, and cross-checks every number against analytic
theory.  Produces figures 1-3 of the report.
"""
import numpy as np
import ssh_circuit as ssh
from common import newfig, save, update_results, CLR, FIGDIR, DATADIR, NETDIR
import os

print("=" * 70)
print("PART 1  --  PRISTINE BASELINE")
print("=" * 70)

c1, c2, cb, lvals = ssh.pristine_components()

# --------------------------------------------------------------------------
# Theory anchors
# --------------------------------------------------------------------------
F_R = ssh.F_R
xi_theory = 1.0 / np.log(ssh.C2 / ssh.C1)          # 1/ln(w/v) cells
ratio_theory = ssh.C1 / ssh.C2
print(f"\nAnalytic TBR frequency   f_R = {F_R/1e3:.4f} kHz")
print(f"Theory localization len  xi  = {xi_theory:.4f} cells   (= 1/ln(C2/C1))")
print(f"Theory decay ratio           = {ratio_theory:.4f}      (= C1/C2)")

# --------------------------------------------------------------------------
# Impedance sweep  (Figure 1)
# --------------------------------------------------------------------------
freqs = np.logspace(np.log10(20e3), np.log10(120e3), 1600)
Z = ssh.impedance_curve(freqs, c1, c2, cb, lvals)
f_peak, z_peak = ssh.find_tbr(c1, c2, cb, lvals)
print(f"\nSimulated TBR peak           = {f_peak/1e3:.4f} kHz   |Z|_peak = {z_peak:,.0f} ohm")
print(f"  agreement with analytic    = {abs(f_peak-F_R)/F_R*100:.3f} %")

# quality factor from the -3 dB (half-power in |Z|) width
half = z_peak / np.sqrt(2.0)
above = freqs[Z >= half]
if len(above) > 1:
    fwhm = above.max() - above.min()
    Q = f_peak / fwhm
else:
    fwhm, Q = np.nan, np.nan
print(f"  resonance FWHM             = {fwhm:.1f} Hz   ->  Q = {Q:.0f}")

# resonances predicted by the capacitance-Laplacian eigenvalues
W = ssh.cap_laplacian(c1, c2, cb)
f_modes = ssh.resonances_from_W(W)               # 16 resonances, sorted
f_mid = f_modes[np.argmin(np.abs(f_modes - F_R))]
print(f"  W-eigenvalue midgap reson  = {f_mid/1e3:.4f} kHz")

# analytic bulk-band edges (frequency)
v, w = ssh.C1, ssh.C2
band_lo = [1.0/(2*np.pi*np.sqrt(ssh.L*2*w*( (w+v)/w ))), 0]  # placeholder
f_of_mu = lambda mu: 1.0/(2*np.pi*np.sqrt(ssh.L*mu))
band_hi_mu = [2*w, 2*(v+w)]          # high-mu band -> low-freq bulk band
band_lo_mu = 2*v                     # low-mu band edge -> high-freq bulk band
bandA = sorted([f_of_mu(band_hi_mu[0]), f_of_mu(band_hi_mu[1])])   # low-freq band
bandB_edge = f_of_mu(band_lo_mu)                                   # high-freq band starts here

overlap = (ssh.C1 / ssh.C2) ** ssh.N_CELLS          # edge-mode overlap
gap_width = bandB_edge - bandA[1]
tbr_in_gap = (f_peak - bandA[1]) / gap_width
print(f"  band gap = {bandA[1]/1e3:.2f} - {bandB_edge/1e3:.2f} kHz "
      f"(width {gap_width/1e3:.2f} kHz); TBR at {tbr_in_gap*100:.0f}% into gap")
print(f"  finite-size edge-mode overlap (C1/C2)^N = {overlap*100:.2f}%")

fig, ax = newfig()
ax.semilogy(freqs/1e3, Z, color=CLR["tbr"], lw=1.6)
ax.axvspan(bandA[0]/1e3, bandA[1]/1e3, color=CLR["bulk"], alpha=0.20)
ax.axvspan(bandB_edge/1e3, 120, color=CLR["bulk"], alpha=0.20,
           label="SSH bulk bands (theory)")
ax.axvline(F_R/1e3, color=CLR["theory"], ls="--", lw=1.5,
           label=f"analytic $f_R$ = {F_R/1e3:.2f} kHz")
ax.plot(f_peak/1e3, z_peak, "o", color="k", ms=8, zorder=5,
        label=f"simulated TBR peak = {f_peak/1e3:.2f} kHz")
ax.annotate("TBR", xy=(f_peak/1e3, z_peak), xytext=(f_peak/1e3+9, z_peak*0.9),
            fontsize=10, ha="left")
ax.set_xlabel("frequency  (kHz)")
ax.set_ylabel(r"input impedance  $|Z_{\rm in}|$  ($\Omega$)")
ax.set_title("Pristine circuit  --  impedance spectrum and the TBR")
ax.set_xlim(20, 120)
ax.legend(loc="lower left", fontsize=9)
save(fig, "fig01_pristine_impedance.png")

np.savetxt(os.path.join(DATADIR, "part1_impedance.csv"),
           np.column_stack([freqs, Z]), delimiter=",",
           header="frequency_Hz,abs_Z_ohm", comments="")

# --------------------------------------------------------------------------
# Voltage profile at the TBR  (Figure 2)
# --------------------------------------------------------------------------
prof = ssh.voltage_profile(f_peak, c1, c2, cb, lvals)
nodes = np.arange(1, 17)
odd, even = nodes[::2], nodes[1::2]
v_odd, v_even = prof[::2], prof[1::2]

xi_meas = ssh.localization_length(prof)
ratio_meas = ssh.decay_ratio(prof)
ipr_meas = ssh.ipr(prof)
print(f"\nVoltage profile at the TBR:")
print(f"  measured decay ratio       = {ratio_meas:.4f}   (theory {ratio_theory:.4f})")
print(f"  measured xi                = {xi_meas:.4f} cells (theory {xi_theory:.4f})")
print(f"  measured IPR (N2..N16)     = {ipr_meas:.4f}")

# exponential fit on the odd (A-sublattice) nodes, cells 2..8
cells_odd = np.array([2, 3, 4, 5, 6, 7, 8])
logv = np.log(prof[2:15:2])
slope, intercept = np.polyfit(cells_odd, logv, 1)
fitline = intercept + slope * cells_odd
r2 = 1 - np.sum((logv - fitline)**2) / np.sum((logv - logv.mean())**2)
print(f"  exponential-fit quality   R^2 = {r2:.5f}")
fit_cells = np.linspace(2, 8, 50)
fit_nodes = 2 * fit_cells - 1
fit_v = np.exp(intercept + slope * fit_cells)

fig, ax = newfig()
ax.semilogy(odd, v_odd, "o-", color=CLR["odd"], ms=7,
            label="odd nodes  (A sublattice, left-edge mode)")
ax.semilogy(even, v_even, "s-", color=CLR["even"], ms=7,
            label="even nodes (B sublattice, right-edge mode)")
ax.semilogy(fit_nodes, fit_v, "--", color=CLR["theory"], lw=1.6,
            label=fr"exp. fit: $\xi$ = {xi_meas:.2f} cells (theory {xi_theory:.2f})")
ax.set_xlabel("node index  $i$")
ax.set_ylabel(r"normalised voltage  $|V(N_i)|\,/\,V_0$")
ax.set_title("Pristine circuit  --  edge-state voltage profile at the TBR")
ax.set_xticks(nodes)
ax.legend(loc="upper center", fontsize=9)
save(fig, "fig02_pristine_voltage.png")

np.savetxt(os.path.join(DATADIR, "part1_voltage_profile.csv"),
           np.column_stack([nodes, prof]), delimiter=",",
           header="node,abs_V_over_V0", comments="")

# --------------------------------------------------------------------------
# Resonance spectrum cross-check  (Figure 3)
# --------------------------------------------------------------------------
midgap_idx = int(np.argmin(np.abs(f_modes - F_R)))
mode_idx = np.arange(1, 17)
is_mid = np.zeros(16, bool); is_mid[midgap_idx] = True

fig, ax = newfig()
ax.axhspan(bandA[0]/1e3, bandA[1]/1e3, color=CLR["bulk"], alpha=0.20)
ax.axhspan(bandB_edge/1e3, 130, color=CLR["bulk"], alpha=0.20,
           label="SSH bulk bands (theory)")
ax.plot(mode_idx[~is_mid], f_modes[~is_mid]/1e3, "o", color=CLR["bulk"],
        ms=8, label="bulk resonances (sim.)")
ax.plot(mode_idx[is_mid], f_modes[is_mid]/1e3, "*", color=CLR["tbr"],
        ms=20, label="topological boundary resonance (sim.)")
ax.axhline(F_R/1e3, color=CLR["theory"], ls="--", lw=1.5,
           label=f"analytic $f_R$ = {F_R/1e3:.2f} kHz")
ax.set_xlabel("mode index  (resonances sorted by frequency)")
ax.set_ylabel("resonance frequency  (kHz)")
ax.set_title("Pristine resonance spectrum  --  TBR isolated in the band gap")
ax.set_xticks(mode_idx)
ax.set_ylim(28, 95)
ax.legend(loc="upper left", fontsize=9)
save(fig, "fig03_resonance_spectrum.png")

np.savetxt(os.path.join(DATADIR, "part1_resonance_spectrum.csv"),
           np.column_stack([mode_idx, f_modes]), delimiter=",",
           header="mode_index,resonance_freq_Hz", comments="")

# --------------------------------------------------------------------------
# LTSpice netlist (for completeness / lab reproduction)
# --------------------------------------------------------------------------
ssh.write_netlist(os.path.join(NETDIR, "ssh_pristine.cir"),
                  c1, c2, cb, lvals, title="SSH pristine -- updated values 2.2nF/3.3nF/2200uH")

# --------------------------------------------------------------------------
# Bookkeeping
# --------------------------------------------------------------------------
update_results("part1", dict(
    f_R_analytic=F_R, f_peak_sim=f_peak, z_peak=z_peak,
    f_midgap_eig=f_mid, agreement_pct=abs(f_peak-F_R)/F_R*100,
    fwhm=fwhm, Q=Q,
    xi_theory=xi_theory, xi_meas=xi_meas, fit_r2=float(r2),
    ratio_theory=ratio_theory, ratio_meas=ratio_meas,
    ipr_pristine=ipr_meas,
    bandA=[bandA[0], bandA[1]], bandB_edge=bandB_edge,
    gap_width=gap_width, tbr_in_gap=tbr_in_gap, overlap=overlap,
    profile=prof.tolist(),
))
print("\nPart 1 complete.\n")
