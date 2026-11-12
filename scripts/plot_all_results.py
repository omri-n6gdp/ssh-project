"""
plot_all_results.py
===================
Generates all robustness-analysis plots for the SSH topolectrical project.

Reads data from:  ~/code/ssh_project/scripts/*.xlsx
Saves plots to:   ~/code/ssh_project/figures/

Data files expected:
  impurity_sweep.xlsx    — single-impurity f_peak vs cell, sweep over delta
  AB_delta_sweep.xlsx    — AB pair, delta sweep at fixed gamma
  AB_gamma_sweep.xlsx    — AB pair, gamma sweep at fixed delta
  BA_delta_sweep.xlsx    — BA pair, delta sweep at fixed gamma
  BA_gamma_sweep.xlsx    — BA pair, gamma sweep at fixed delta
  AB_voltage_domain.xlsx — AB pair voltage profile vs node
  BA_voltage_domain.xlsx — BA pair voltage profile vs node

Error bars: derived from the AC sweep frequency resolution
(.ac dec 1000 over the swept range gives ~0.23% spacing between adjacent
points, so the peak-frequency uncertainty is ~half a bin width).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─── User-specific paths ──────────────────────────────────────────────────────
HOME       = os.path.expanduser("~")
DATA_DIR   = os.path.join(HOME, "ssh_project", "scripts")
OUT_DIR    = os.path.join(HOME, "ssh_project", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Physics constants ────────────────────────────────────────────────────────
F_R        = 189685.0      # pristine TBR frequency in Hz
PTS_PER_DECADE = 1000      # matches your .ac dec 1000 directive
# Frequency uncertainty from finite sweep resolution:
#   delta_log10_f = 1/N_per_decade  →  σ_f ≈ ln(10)/N_per_decade × f / 2
SIGMA_F_REL = np.log(10) / PTS_PER_DECADE / 2.0   # ≈ 0.115%
SIGMA_F     = SIGMA_F_REL * F_R                   # in Hz

# Consistent colour map for delta/gamma values
DELTA_COLORS = {
     0.3: '#d73027',     # red
     0.1: '#fc8d59',     # orange
     0.0: 'black',
    -0.1: '#4575b4',     # blue
    -0.3: '#313695',     # dark blue
}
MARKER_FOR = {0.3: 's', 0.1: '^', 0.0: 'D', -0.1: 'v', -0.3: 'o'}


# ─── Generic plot builder for f_peak vs cell sweeps ───────────────────────────
def plot_frequency_sweep(filename, title, outname,
                         param_symbol=r'\delta',
                         add_pristine_line=True):
    """Plot f/f_R vs impurity cell for several values of an impurity parameter.

    filename : str  — name of the xlsx in DATA_DIR
    title    : str  — figure title
    outname  : str  — output file name (PNG)
    param_symbol : LaTeX string for the swept parameter (delta or gamma)
    """
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"  [skip] {filename} not found")
        return

    df = pd.read_excel(path)
    cells = df['cell'].values
    param_cols = [c for c in df.columns if c != 'cell']

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for p in param_cols:
        try:
            p_float = float(p)
        except (TypeError, ValueError):
            continue
        freqs = df[p].values.astype(float)
        ratios = freqs / F_R
        # propagate σ_f to σ(f/fR) (f_R itself is taken as exact)
        sigma_ratio = SIGMA_F / F_R * np.ones_like(ratios)

        color  = DELTA_COLORS.get(round(p_float, 3), None)
        marker = MARKER_FOR.get(round(p_float, 3), 'o')

        ax.errorbar(cells, ratios, yerr=sigma_ratio,
                    marker=marker, ms=8, lw=1.8,
                    capsize=4, capthick=1.2, elinewidth=1.2,
                    color=color,
                    label=fr'${param_symbol} = {p_float:+.1f}$')

    if add_pristine_line:
        ax.axhline(1.0, color='red', linestyle='--', lw=1.2, alpha=0.7,
                   label='Pristine ($f / f_R = 1$)')

    ax.set_xlabel('Impurity starting cell', fontsize=12)
    ax.set_ylabel(r'$f_{\rm peak} / f_R$', fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xticks(cells)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10, loc='best')

    out = os.path.join(OUT_DIR, outname)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  saved → {out}")


# ─── Voltage-profile plotter (for domain-wall files) ──────────────────────────
def plot_voltage_profile(filename, title, outname,
                         param_symbol=r'\gamma',
                         pristine_db=None,
                         impurity_region=None):
    """Plot |V|/V0 vs node for several values of an impurity parameter.

    filename       : xlsx file containing 'node' and one column per param value
    pristine_db    : optional list of 16 pristine dB values for reference
    impurity_region: tuple (left, right) for the shaded impurity region
    """
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"  [skip] {filename} not found")
        return

    df = pd.read_excel(path)
    nodes = df['node'].values
    param_cols = [c for c in df.columns if c != 'node']

    # Check whether data is filled in
    data_block = df[param_cols].values
    if np.all(pd.isna(data_block)):
        print(f"  [skip] {filename} contains no data yet")
        return

    odd_nodes  = nodes[nodes % 2 != 0]
    even_nodes = nodes[nodes % 2 == 0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, sublattice, node_list, marker in zip(
            axes,
            ['Odd nodes (A-sublattice)', 'Even nodes (B-sublattice)'],
            [odd_nodes, even_nodes],
            ['o', 's']):

        idx = node_list - 1

        # Pristine reference (if supplied)
        if pristine_db is not None:
            V_pristine = 10**(np.array(pristine_db) / 20.0)  # V/V0
            ax.semilogy(node_list, V_pristine[idx],
                        marker+'--', color='gray', lw=1.5, ms=6, alpha=0.6,
                        label='Pristine', zorder=5)

        for p in param_cols:
            try:
                p_float = float(p)
            except (TypeError, ValueError):
                continue
            db_vals = df[p].values.astype(float)
            if np.all(np.isnan(db_vals)):
                continue
            V = 10**(db_vals / 20.0)
            color = DELTA_COLORS.get(round(p_float, 3), None)
            # Error in dB measurement is small; we propagate ~0.05 dB → 0.6%
            sigma_db = 0.05
            sigma_V  = V * np.log(10) / 20.0 * sigma_db
            ax.errorbar(node_list, V[idx], yerr=sigma_V[idx],
                        marker=marker, ms=7, lw=2,
                        capsize=3, capthick=1, elinewidth=1,
                        color=color,
                        label=fr'${param_symbol} = {p_float:+.1f}$')

        if impurity_region is not None:
            ax.axvspan(*impurity_region, alpha=0.12, color='purple',
                       label='Impurity region')

        ax.set_title(sublattice, fontsize=12)
        ax.set_xlabel('Node index', fontsize=12)
        ax.set_ylabel(r'$|V|/V_0$', fontsize=12)
        ax.set_xticks(nodes)
        ax.set_xlim(0.5, 16.5)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3, which='both')

    plt.suptitle(title, fontsize=13, fontweight='bold')
    plt.tight_layout()
    out = os.path.join(OUT_DIR, outname)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  saved → {out}")


# ─── Pristine voltage profile (your earlier baseline measurement) ─────────────
PRISTINE_DB = [0, -32.6779, -6.85828, -39.5375, -13.7166, -46.4004,
               -20.5748, -53.2718, -27.433, -60.0825, -34.2911, -65.0451,
               -41.1485, -61.5394, -48.0043, -54.8527]


# ─── Generate all plots ───────────────────────────────────────────────────────
def main():
    print("Generating plots...")

    # 1) Single-impurity frequency sweep (Anzurez Fig. 4 analog)
    plot_frequency_sweep(
        filename='impurity_sweep.xlsx',
        title='Single impurity: resonance displacement vs cell position\n'
              r'(Replication of Anzurez et al. 2025, Fig. 4)',
        outname='fig_impurity_sweep.png',
        param_symbol=r'\delta',
    )

    # 2) AB-pair delta sweep (Anzurez Fig. 7 analog)
    plot_frequency_sweep(
        filename='AB_delta_sweep.xlsx',
        title='AB pair: resonance displacement vs cell position ($\\delta$ swept)\n'
              r'(Replication of Anzurez et al. 2025, Fig. 7)',
        outname='fig_AB_delta_sweep.png',
        param_symbol=r'\delta',
    )

    # 3) AB-pair gamma sweep
    plot_frequency_sweep(
        filename='AB_gamma_sweep.xlsx',
        title='AB pair: resonance displacement vs cell position ($\\gamma$ swept)\n'
              r'(Replication of Anzurez et al. 2025, Fig. 7)',
        outname='fig_AB_gamma_sweep.png',
        param_symbol=r'\gamma',
    )

    # 4) BA-pair delta sweep (Anzurez Fig. 8 analog)
    plot_frequency_sweep(
        filename='BA_delta_sweep.xlsx',
        title='BA pair: resonance displacement vs cell position ($\\delta$ swept)\n'
              r'(Replication of Anzurez et al. 2025, Fig. 8)',
        outname='fig_BA_delta_sweep.png',
        param_symbol=r'\delta',
    )

    # 5) BA-pair gamma sweep
    plot_frequency_sweep(
        filename='BA_gamma_sweep.xlsx',
        title='BA pair: resonance displacement vs cell position ($\\gamma$ swept)\n'
              r'(Replication of Anzurez et al. 2025, Fig. 8)',
        outname='fig_BA_gamma_sweep.png',
        param_symbol=r'\gamma',
    )

    # 6) AB-pair voltage profile (will skip if all-NaN)
    plot_voltage_profile(
        filename='AB_voltage_domain.xlsx',
        title='AB-pair domain wall — voltage profile\n'
              r'(Replication of Anzurez et al. 2025, Fig. 14)',
        outname='fig_AB_voltage_domain.png',
        param_symbol=r'\gamma',
        pristine_db=PRISTINE_DB,
        impurity_region=(2.5, 5.5),
    )

    # 7) BA-pair voltage profile (will skip if all-NaN)
    plot_voltage_profile(
        filename='BA_voltage_domain.xlsx',
        title='BA-pair domain wall — voltage profile\n'
              r'(Replication of Anzurez et al. 2025, Fig. 14)',
        outname='fig_BA_voltage_domain.png',
        param_symbol=r'\gamma',
        pristine_db=PRISTINE_DB,
        impurity_region=(2.5, 5.5),
    )

    print("Done.")


if __name__ == "__main__":
    main()
