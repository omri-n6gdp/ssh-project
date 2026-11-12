"""Generate Monte Carlo plots from the saved CSVs (no re-simulation)."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FIG_DIR = os.path.expanduser("~/code/ssh_project/figures")
RESULTS = os.path.join(FIG_DIR, "monte_carlo_results.csv")
VOLTAGES = os.path.join(FIG_DIR, "monte_carlo_voltages.csv")
SUMMARY = os.path.join(FIG_DIR, "mc_summary.csv")

df = pd.read_csv(RESULTS)
df_v = pd.read_csv(VOLTAGES)

# Pristine reference from a deterministic re-derivation: use mean of the lowest-disorder sub-bin
# as a stand-in (the actual pristine values were printed by monte_carlo.py).
PRISTINE = dict(
    f_peak=189689.5,      # Hz, refined sub-grid value
    z_peak=81436.1,
    IPR_full=0.6578,
    IPR_bulk=0.6555,
    xi=2.532,
)
N_NODES = 16
sigmas = sorted(df.sigma.unique())


def _hist(ax, data, pristine_val, color, label_fmt):
    ax.hist(data, bins=40, alpha=0.75, edgecolor="black", color=color)
    ax.axvline(pristine_val, color="red", ls="--",
               label=f"pristine = {label_fmt.format(pristine_val)}")


# ---- f_peak ----
fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
for ax, sig in zip(axes, sigmas):
    sub = df[df.sigma == sig]
    _hist(ax, sub.f_peak / 1e3, PRISTINE["f_peak"] / 1e3,
          color="C0", label_fmt="{:.2f} kHz")
    fwhm = 2.355 * sub.f_peak.std()
    ax.set_title(f"sigma = {sig*100:.0f}%  (FWHM = {fwhm/1e3:.2f} kHz)")
    ax.set_xlabel("f_peak [kHz]")
    ax.legend(fontsize=8)
axes[0].set_ylabel("count")
fig.suptitle("Monte Carlo distribution of resonance frequency under capacitor disorder")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "mc_fpeak_histograms.png"), dpi=130)
plt.close(fig)
print("wrote mc_fpeak_histograms.png")

# ---- IPR (bulk) ----
fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
for ax, sig in zip(axes, sigmas):
    sub = df[df.sigma == sig]
    ax.hist(sub.IPR_bulk, bins=40, alpha=0.75, edgecolor="black", color="C2")
    ax.axvline(PRISTINE["IPR_bulk"], color="red", ls="--",
               label=f"pristine = {PRISTINE['IPR_bulk']:.3f}")
    ax.axvline(1 / (N_NODES - 1), color="gray", ls=":",
               label=f"uniform = {1/(N_NODES-1):.3f}")
    ax.set_title(f"sigma = {sig*100:.0f}%")
    ax.set_xlabel("IPR (bulk, nodes 2..16)")
    ax.legend(fontsize=8)
axes[0].set_ylabel("count")
fig.suptitle("Inverse Participation Ratio -- chiral-symmetry protection of edge state")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "mc_ipr_histograms.png"), dpi=130)
plt.close(fig)
print("wrote mc_ipr_histograms.png")

# ---- Z_peak ----
fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
for ax, sig in zip(axes, sigmas):
    sub = df[df.sigma == sig]
    ax.hist(sub.z_peak / 1e3, bins=40, alpha=0.75, edgecolor="black", color="C1")
    ax.axvline(PRISTINE["z_peak"] / 1e3, color="red", ls="--",
               label=f"pristine = {PRISTINE['z_peak']/1e3:.1f} k")
    ax.set_title(f"sigma = {sig*100:.0f}%")
    ax.set_xlabel("|Z|_peak [k$\\Omega$ × (V_src=1V)]")
    ax.legend(fontsize=8)
axes[0].set_ylabel("count")
fig.suptitle("Peak impedance at the resonance under disorder")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "mc_zpeak_histograms.png"), dpi=130)
plt.close(fig)
print("wrote mc_zpeak_histograms.png")

# ---- xi ----
fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
for ax, sig in zip(axes, sigmas):
    sub = df[df.sigma == sig].dropna(subset=["xi"])
    xi_arr = np.clip(sub.xi.values, 0, 6)
    ax.hist(xi_arr, bins=40, alpha=0.75, edgecolor="black", color="C3")
    ax.axvline(PRISTINE["xi"], color="red", ls="--",
               label=f"pristine = {PRISTINE['xi']:.2f}")
    ax.set_title(f"sigma = {sig*100:.0f}%")
    ax.set_xlabel("localization length $\\xi$ (sites)")
    ax.legend(fontsize=8)
axes[0].set_ylabel("count")
fig.suptitle("Edge-state localization length from log-linear fit on odd nodes N3..N15")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "mc_xi_histograms.png"), dpi=130)
plt.close(fig)
print("wrote mc_xi_histograms.png")

# ---- Voltage envelope ----
# For each sigma, show the median and 5%/95% percentiles of |V_i|/V_0
fig, ax = plt.subplots(figsize=(8, 5.5))
node_idx = np.arange(1, N_NODES + 1)
V_cols = [f"V{i}" for i in node_idx]
pristine_voltages = np.array([1.0, 0.0213, 0.4539, 0.0096, 0.2060, 0.0044,
                              0.0935, 0.0020, 0.0424, 0.0009, 0.0193, 0.0005,
                              0.0087, 0.0008, 0.0040, 0.0018])
ax.semilogy(node_idx, pristine_voltages, "ko-", lw=1.5, label="pristine")
colors = ["C0", "C2", "C3"]
for sig, color in zip(sigmas, colors):
    sub = df_v[df_v.sigma == sig][V_cols].values
    med = np.median(sub, axis=0)
    p5 = np.percentile(sub, 5, axis=0)
    p95 = np.percentile(sub, 95, axis=0)
    ax.semilogy(node_idx, med, "o-", color=color, alpha=0.9,
                label=f"$\\sigma$ = {sig*100:.0f}%  (median)")
    ax.fill_between(node_idx, p5, p95, color=color, alpha=0.15)
ax.set_xlabel("node index $i$")
ax.set_ylabel("|V$_i$| / V$_0$")
ax.set_xticks(node_idx)
ax.set_title("Voltage profile at resonance: pristine vs disorder ensemble (5-95 pctile band)")
ax.legend(fontsize=9)
ax.grid(True, which="both", ls=":", alpha=0.4)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "mc_voltage_envelope.png"), dpi=130)
plt.close(fig)
print("wrote mc_voltage_envelope.png")

# ---- f_peak vs sigma summary (one figure) ----
fig, ax = plt.subplots(figsize=(7, 4.5))
means = df.groupby("sigma").f_peak.mean()
stds = df.groupby("sigma").f_peak.std()
sigs_pct = np.array(sigmas) * 100
ax.errorbar(sigs_pct, means / 1e3, yerr=stds / 1e3, fmt="o-", lw=1.5,
            capsize=4, color="C0", label="MC ensemble (mean ± std)")
ax.axhline(PRISTINE["f_peak"] / 1e3, color="red", ls="--",
           label=f"pristine = {PRISTINE['f_peak']/1e3:.3f} kHz")
ax.set_xlabel("disorder strength $\\sigma$ (%)")
ax.set_ylabel("f_peak [kHz]")
ax.set_title("Resonance frequency under capacitor disorder")
ax.legend()
ax.grid(True, ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "mc_fpeak_vs_sigma.png"), dpi=130)
plt.close(fig)
print("wrote mc_fpeak_vs_sigma.png")

# ---- Concise summary text ----
print()
print("SUMMARY")
print("=======")
summary = pd.read_csv(SUMMARY, index_col=0)
print(summary.to_string(float_format="{:.4f}".format))
