"""
ssh_circuit.py
==============
From-scratch AC nodal-admittance simulator for the SSH topolectrical circuit.

This module recreates, from first principles, the circuit described in
``ssh_project/robustness_manual.md`` (after Lee et al. 2018 and
Anzurez et al. 2025).  No data from the existing project is used: every
number this study reports is produced by the solver below.

Circuit topology (Su-Schrieffer-Heeger topolectrical chain)
-----------------------------------------------------------
  * N = 8 unit cells  ->  16 nodes, numbered 1..16
  * a grounded inductor L at every node                       (16 inductors)
  * a chain of capacitors connecting adjacent nodes, alternating
        C1 (intracell, analogue of hopping "v") : (1,2),(3,4),...,(15,16)
        C2 (intercell, analogue of hopping "w") : (2,3),(4,5),...,(14,15)
  * a grounded boundary capacitor Cb = C2 attached to node 1
  * an AC voltage source V0 driving node 1

Simulation engine
-----------------
Build the 16x16 complex nodal-admittance matrix Y(omega), force the driven
node to V0, and solve  Y_ff v_f = -Y_fd V0  for the other 15 node voltages.
This is *exactly* the linear system a SPICE ``.ac`` analysis solves at each
frequency point -- so the results are SPICE-accurate by construction, while
remaining fully reproducible and theory-checkable.  The assembly and solve
are vectorised over frequency (batched LAPACK) so a full sweep, or a Monte
Carlo of thousands of circuits, runs in seconds.

UPDATED component values for this study
---------------------------------------
    C1 = 2.2 nF   intracell capacitor
    C2 = 3.3 nF   intercell capacitor
    L  = 2200 uH  grounded inductor
    series resistance = 1 ohm on every capacitor and every inductor
"""

import numpy as np

# --------------------------------------------------------------------------
# Nominal component values  (the "updated circuit" of this study)
# --------------------------------------------------------------------------
C1 = 2.2e-9      # F  -- intracell capacitor  (analogue of v)
C2 = 3.3e-9      # F  -- intercell capacitor  (analogue of w)
L  = 2200e-6     # H  -- grounded inductor at every node
R_C = 1.0        # ohm -- series resistance of every capacitor
R_L = 1.0        # ohm -- series resistance of every inductor
V0  = 1.0        # V  -- AC source amplitude at node 1

N_CELLS = 8
N_NODES = 16
N_C1 = 8         # intracell capacitors  c1[0..7]
N_C2 = 7         # intercell capacitors  c2[0..6]

# Fixed topology: node-index endpoints of the 15 inter-node capacitors.
# Order is [C1 x8, C2 x7] so a value vector concatenates as cvals = [c1, c2].
_IA = np.array([0, 2, 4, 6, 8, 10, 12, 14,  1, 3, 5, 7, 9, 11, 13])
_IB = np.array([1, 3, 5, 7, 9, 11, 13, 15,  2, 4, 6, 8, 10, 12, 14])


# --------------------------------------------------------------------------
# Analytic reference (textbook theory)
# --------------------------------------------------------------------------
def f_resonance(L=L, c1=C1, c2=C2):
    """Topological boundary resonance, f_R = 1/(2 pi sqrt(L (C1+C2)))."""
    return 1.0 / (2.0 * np.pi * np.sqrt(L * (c1 + c2)))

F_R = f_resonance()          # ~= 45.75 kHz for the updated values


def pristine_components():
    """Return (c1, c2, cb, lvals) for the perfectly pristine circuit."""
    c1 = np.full(N_C1, C1)
    c2 = np.full(N_C2, C2)
    cb = C2
    lvals = np.full(N_NODES, L)
    return c1, c2, cb, lvals


# --------------------------------------------------------------------------
# Core: vectorised nodal-admittance assembly and AC solve
# --------------------------------------------------------------------------
def _Y_stack(freqs, cvals, cb, lvals, rc, rl):
    """Assemble the nodal-admittance matrix at every frequency.

    Returns Y with shape (Nf, 16, 16), complex.  Node indices 0..15 are the
    physical nodes 1..16; ground is the reference.
    """
    omega = 2.0 * np.pi * np.atleast_1d(np.asarray(freqs, dtype=float))
    Nf = omega.size

    yL  = 1.0 / (rl + 1j * omega[:, None] * lvals[None, :])          # (Nf,16)
    yC  = 1.0 / (rc + 1.0 / (1j * omega[:, None] * cvals[None, :]))  # (Nf,15)
    yCb = 1.0 / (rc + 1.0 / (1j * omega * cb))                       # (Nf,)

    Y = np.zeros((Nf, N_NODES, N_NODES), dtype=complex)
    di = np.arange(N_NODES)
    Y[:, di, di] += yL                       # grounded inductor at every node
    Y[:, 0, 0]  += yCb                       # grounded boundary capacitor

    fr  = np.broadcast_to(np.arange(Nf)[:, None], (Nf, 15))
    iam = np.broadcast_to(_IA[None, :], (Nf, 15))
    ibm = np.broadcast_to(_IB[None, :], (Nf, 15))
    np.add.at(Y, (fr, iam, iam),  yC)        # capacitor degree (diagonal)
    np.add.at(Y, (fr, ibm, ibm),  yC)
    np.add.at(Y, (fr, iam, ibm), -yC)        # capacitor coupling (off-diag)
    np.add.at(Y, (fr, ibm, iam), -yC)
    return Y


def solve_stack(freqs, c1, c2, cb, lvals, rc=R_C, rl=R_L, v0=V0):
    """Solve the driven circuit at an array of frequencies.

    Returns
    -------
    V    : complex ndarray (Nf, 16)  -- node voltages, V[:,0] = v0
    Zin  : complex ndarray (Nf,)     -- input impedance V0 / I_source
    """
    cvals = np.concatenate([np.asarray(c1, float), np.asarray(c2, float)])
    Y = _Y_stack(freqs, cvals, cb, lvals, rc, rl)

    # node 1 (index 0) is held at v0 by the ideal source; solve the other 15
    Yff = Y[:, 1:, 1:]
    Yfd = Y[:, 1:, 0]
    vf = np.linalg.solve(Yff, (-Yfd * v0)[:, :, None])[:, :, 0]   # (Nf,15)

    Nf = vf.shape[0]
    V = np.empty((Nf, N_NODES), dtype=complex)
    V[:, 0] = v0
    V[:, 1:] = vf

    Isrc = np.einsum('fj,fj->f', Y[:, 0, :], V)          # source current
    Zin = v0 / Isrc
    return V, Zin


def solve_ac(freq, c1, c2, cb, lvals, rc=R_C, rl=R_L, v0=V0):
    """Single-frequency solve. Returns (V[16], Zin, Isrc)."""
    V, Zin = solve_stack([freq], c1, c2, cb, lvals, rc, rl, v0)
    return V[0], Zin[0], v0 / Zin[0]


def impedance_curve(freqs, c1, c2, cb, lvals, rc=R_C, rl=R_L, v0=V0):
    """Return |Zin| over an array of frequencies (Hz)."""
    _, Zin = solve_stack(freqs, c1, c2, cb, lvals, rc, rl, v0)
    return np.abs(Zin)


def build_Y(omega, c1, c2, cb, lvals, rc=R_C, rl=R_L):
    """Single-frequency nodal-admittance matrix (readable reference form)."""
    cvals = np.concatenate([np.asarray(c1, float), np.asarray(c2, float)])
    return _Y_stack([omega / (2 * np.pi)], cvals, cb, lvals, rc, rl)[0]


# --------------------------------------------------------------------------
# Lossless theory: capacitance Laplacian and its resonances
# --------------------------------------------------------------------------
def cap_laplacian(c1, c2, cb):
    """Lossless capacitance (degree-)Laplacian W.

    The lossless nodal matrix is  Y = j omega W - (j/omega) diag(1/L).
    A resonance (zero eigenvalue) for uniform L occurs at
        f_n = 1 / (2 pi sqrt(L * mu_n))
    where mu_n are the eigenvalues of W.
    """
    W = np.zeros((N_NODES, N_NODES))
    W[0, 0] += cb
    for k in range(N_C1):
        a, b = 2 * k, 2 * k + 1
        W[a, a] += c1[k]; W[b, b] += c1[k]
        W[a, b] -= c1[k]; W[b, a] -= c1[k]
    for k in range(N_C2):
        a, b = 2 * k + 1, 2 * k + 2
        W[a, a] += c2[k]; W[b, b] += c2[k]
        W[a, b] -= c2[k]; W[b, a] -= c2[k]
    return W


def resonances_from_W(W, Lval=L):
    """Resonance frequencies (Hz) from the eigenvalues of W (uniform Lval)."""
    mu = np.linalg.eigvalsh(W)
    mu = mu[mu > 0]
    return np.sort(1.0 / (2.0 * np.pi * np.sqrt(Lval * mu)))


# --------------------------------------------------------------------------
# Peak finding for the topological boundary resonance (TBR)
# --------------------------------------------------------------------------
def find_tbr(c1, c2, cb, lvals, rc=R_C, rl=R_L, v0=V0,
             f_lo=38e3, f_hi=56e3, n_coarse=320, n_fine=160):
    """Locate the dominant impedance peak (the TBR) in [f_lo, f_hi].

    Two-stage search: a coarse sweep brackets the peak, a fine local sweep
    resolves it below the resonance line-width, and a three-point parabolic
    fit refines it to sub-grid accuracy.

    Returns (f_peak [Hz], z_peak [ohm]).
    """
    fc = np.linspace(f_lo, f_hi, n_coarse)
    zc = impedance_curve(fc, c1, c2, cb, lvals, rc, rl, v0)
    i = int(np.argmax(zc))

    lo = fc[max(0, i - 1)]
    hi = fc[min(n_coarse - 1, i + 1)]
    ff = np.linspace(lo, hi, n_fine)
    zf = impedance_curve(ff, c1, c2, cb, lvals, rc, rl, v0)
    j = int(np.argmax(zf))

    if 0 < j < n_fine - 1:
        x = ff[j - 1:j + 2]
        y = zf[j - 1:j + 2]
        d = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
        A = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) +
             x[0] * (y[2] - y[1])) / d
        B = (x[2]**2 * (y[0] - y[1]) + x[1]**2 * (y[2] - y[0]) +
             x[0]**2 * (y[1] - y[2])) / d
        if A < 0:
            f_peak = float(-B / (2 * A))
            _, Zin = solve_stack([f_peak], c1, c2, cb, lvals, rc, rl, v0)
            return f_peak, float(abs(Zin[0]))
    return float(ff[j]), float(zf[j])


def voltage_profile(freq, c1, c2, cb, lvals, rc=R_C, rl=R_L, v0=V0):
    """Return |V(N_i)| / v0 for i = 1..16 at a single frequency."""
    V, _ = solve_stack([freq], c1, c2, cb, lvals, rc, rl, v0)
    return np.abs(V[0]) / v0


# --------------------------------------------------------------------------
# Localization metrics
# --------------------------------------------------------------------------
def ipr(profile, nodes=None):
    """Inverse Participation Ratio of a voltage profile.

    profile : |V| at the 16 nodes (any positive scaling).
    nodes   : 0-based indices to include; default = nodes 2..16 (index 1..15),
              excluding the driven node 1 whose voltage is pinned to V0.

    A perfectly delocalised state over m nodes gives IPR = 1/m; a state
    localised on a single node gives IPR = 1.
    """
    if nodes is None:
        nodes = np.arange(1, N_NODES)
    p = profile[nodes]**2
    p = p / p.sum()
    return float(np.sum(p**2))


def localization_length(profile):
    """Edge-state decay length xi, in unit-cell units.

    The left-edge topological mode lives on the odd (A-sublattice) nodes.
    Fit  log|V| = a - (1/xi) * cell_index  on nodes N3, N5, ..., N15
    (cell indices 2..8); xi is the e-folding length per unit cell.
    """
    odd_idx = np.array([2, 4, 6, 8, 10, 12, 14])
    cells = np.array([2, 3, 4, 5, 6, 7, 8])
    slope, _ = np.polyfit(cells, np.log(profile[odd_idx]), 1)
    return float(1.0 / abs(slope)) if slope != 0 else np.inf


def decay_ratio(profile):
    """Geometric-mean voltage ratio between successive odd (A-sublattice)
    nodes -- equals C1/C2 for the pristine topological edge mode."""
    odd_idx = np.array([2, 4, 6, 8, 10, 12, 14])
    v = profile[odd_idx]
    return float(np.exp(np.mean(np.log(v[1:] / v[:-1]))))


# --------------------------------------------------------------------------
# LTSpice netlist export (for completeness / lab reproduction)
# --------------------------------------------------------------------------
def write_netlist(path, c1, c2, cb, lvals, rc=R_C, rl=R_L, v0=V0,
                  ac_directive=".ac dec 400 20k 120k",
                  title="SSH topolectrical circuit"):
    """Write an LTSpice-readable .cir netlist of the current circuit."""
    lines = [f"* {title}"]
    lines.append(f"V1 N001 0 AC {v0}")
    lines.append(f"Cb N001 0 {cb:.6e} Rser={rc}")
    for i in range(N_NODES):
        lines.append(f"L{i+1} 0 N{i+1:03d} {lvals[i]:.6e} Rser={rl}")
    cn = 1
    for k in range(N_C1):
        a, b = 2 * k + 1, 2 * k + 2
        lines.append(f"C{cn} N{a:03d} N{b:03d} {c1[k]:.6e} Rser={rc}  ; C1 cell {k+1}")
        cn += 1
    for k in range(N_C2):
        a, b = 2 * k + 2, 2 * k + 3
        lines.append(f"C{cn} N{a:03d} N{b:03d} {c2[k]:.6e} Rser={rc}  ; C2 cell {k+1}")
        cn += 1
    lines.append(ac_directive)
    lines.append(".end")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------
if __name__ == "__main__":
    c1, c2, cb, lvals = pristine_components()
    print(f"Analytic f_R              = {F_R/1e3:.4f} kHz")
    f_peak, z_peak = find_tbr(c1, c2, cb, lvals)
    print(f"Impedance-sweep TBR peak  = {f_peak/1e3:.4f} kHz   |Z| = {z_peak:,.1f} ohm")
    W = cap_laplacian(c1, c2, cb)
    fr = resonances_from_W(W)
    f_mid = fr[np.argmin(np.abs(fr - F_R))]
    print(f"W-eigenvalue midgap reson = {f_mid/1e3:.4f} kHz")
    prof = voltage_profile(f_peak, c1, c2, cb, lvals)
    print(f"IPR (N2..N16)             = {ipr(prof):.4f}")
    print(f"localization length xi    = {localization_length(prof):.4f} cells")
    print(f"decay ratio (meas / C1/C2)= {decay_ratio(prof):.4f} / {C1/C2:.4f}")
