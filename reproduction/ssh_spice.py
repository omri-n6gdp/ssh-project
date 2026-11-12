"""
ssh_spice.py
============
Build, simulate (real ngspice) and cross-check the SSH topolectrical circuit
taken from the most-recent schematic  netlists/Impurities/AB_3.3_2.2.asc .

Topology (verified against the .asc):

    V1(AC 1) --[C1=3.3n]-- n1 --[C2=2.2n]-- n2 --[C3=3.3n]-- n3 -- ... -- n16
                            |               |               |
                           L1              L2              L3   ... (to ground)

  * 16 nodes, each grounded through L = 2200 uH (Rser = 8 ohm)
  * 16 series capacitors C1..C16 chaining source->n1->n2->...->n16
        odd-index caps  = 3.3 nF (Rser = 30)   <- "strong"  (intercell / boundary)
        even-index caps = 2.2 nF (Rser = 40)   <- "weak"    (intracell)
  * AC drive 35 kHz .. 55 kHz

Unit cell k = (n[2k-1], n[2k]).  Intracell bond C[2k]=2.2 (weak), intercell
bond C[2k+1]=3.3 (strong).  Because the bond leaving the left edge is weak, a
topological edge mode lives on the odd (A-sublattice) nodes n1,n3,...,n15 and
decays by C_intra/C_inter = 2.2/3.3 = 0.667 per cell.
"""

import os
import subprocess
import tempfile
import numpy as np

# ---------------------------------------------------------------------------
# Nominal component values (exactly as in AB_3.3_2.2.asc / BA_3.3_2.2.asc)
# ---------------------------------------------------------------------------
C_STRONG = 3.3e-9        # F   odd caps  C1,C3,...  (Rser 30)
C_WEAK   = 2.2e-9        # F   even caps C2,C4,...  (Rser 40)
RSER_STRONG = 30.0       # ohm
RSER_WEAK   = 40.0       # ohm
L_GND   = 2200e-6        # H   grounding inductor at every node (Rser 8)
RSER_L  = 8.0            # ohm
V0      = 1.0            # V   AC source amplitude

N_NODES = 16
N_CAPS  = 16             # C1..C16 : C1 is source->n1, Ck is n(k-1)->nk

NGSPICE = "/opt/homebrew/bin/ngspice"

# Analytic boundary resonance: f = 1/(2 pi sqrt(L (C_strong + C_weak)))
F_R = 1.0 / (2.0 * np.pi * np.sqrt(L_GND * (C_STRONG + C_WEAK)))


# ---------------------------------------------------------------------------
# Capacitor-bank construction
# ---------------------------------------------------------------------------
def pristine_caps():
    """16 series-cap values, alternating strong/weak starting with strong (C1)."""
    return np.array([C_STRONG if (k % 2 == 0) else C_WEAK for k in range(N_CAPS)])


def impurity_caps(kind, cell, delta, gamma):
    """Return the 16-cap array with an AB or BA *pair* impurity at `cell`
    (Anzurez et al. 2025, Table 1, Models AB / BA).

    A pair impurity is two adjacent on-site impurities. In BOTH cases the
    A-site of `cell` (odd node 2*cell-1, where the left-edge mode lives) carries
    the swept strength `delta`; the partner B-site carries `gamma`. They differ
    only in which side the B-partner sits:

      AB : partner = B-site of the SAME cell (node 2*cell), to the RIGHT.
           A-site caps C(2c-1),C(2c) *= (1+delta); B caps C(2c),C(2c+1) *= (1+gamma).
           Net (strong,weak,strong) triple:  C(2c-1)(1+delta), C(2c)(1+avg), C(2c+1)(1+gamma).
      BA : partner = B-site of the PREVIOUS cell (node 2*cell-2), to the LEFT.
           A-site caps C(2c-1),C(2c) *= (1+delta); B caps C(2c-2),C(2c-1) *= (1+gamma).
           Net (weak,strong,weak) triple:    C(2c-2)(1+gamma), C(2c-1)(1+avg), C(2c)(1+delta).

    avg = (delta+gamma)/2 on the shared bond.  Cap Ck (1-based) -> 0-based k-1.
    """
    c = pristine_caps()
    avg = (delta + gamma) / 2.0
    if kind == "AB":
        iL, iM, iR = 2 * cell - 2, 2 * cell - 1, 2 * cell      # C(2c-1),C(2c),C(2c+1)
        c[iL] *= (1.0 + delta)
        c[iM] *= (1.0 + avg)
        c[iR] *= (1.0 + gamma)
    else:  # 'BA' : triple shifted one bond left, delta on the A-site (right) cap
        iL, iM, iR = 2 * cell - 3, 2 * cell - 2, 2 * cell - 1  # C(2c-2),C(2c-1),C(2c)
        c[iL] *= (1.0 + gamma)
        c[iM] *= (1.0 + avg)
        c[iR] *= (1.0 + delta)
    return c


def single_caps(sublattice, cell, delta):
    """Single-sublattice on-site impurity (Anzurez et al. 2025, Table 1, Models 1A/1B).

    A single impurity perturbs the on-site capacitance of ONE node by scaling the
    TWO capacitors that touch it by (1+delta) -- not a single bond.

      'A' impurity in cell n -> perturb the A-site = odd node (2n-1).
          Caps touching node (2n-1): C(2n-1) [intercell/strong] and C(2n) [intracell/weak].
          The left-edge zero-mode lives on the odd (A) sublattice, so this perturbs it
          directly -> a visible, localised bump at node (2n-1).
      'B' impurity in cell n -> perturb the B-site = even node (2n).
          Caps touching node (2n): C(2n) [intracell/weak] and C(2n+1) [intercell/strong].
          The mode has ~zero weight on even nodes, so the odd-node profile stays ~flat.

    Cap index Ck (1-based) connects node(k-1)->node k; 0-based array index = k-1.
    """
    c = pristine_caps()
    node = (2 * cell - 1) if sublattice == "A" else (2 * cell)   # 1-based node index
    left = node - 1                 # 0-based idx of cap C(node)   touching node from the left
    right = node                    # 0-based idx of cap C(node+1) touching node from the right
    c[left] *= (1.0 + delta)
    if right < N_CAPS:              # last node (16) has no cap to its right
        c[right] *= (1.0 + delta)
    return c


def _rser_for(cval):
    """Match each cap to the Rser used in the schematic by its nominal family."""
    return RSER_STRONG if abs(cval - C_STRONG) < abs(cval - C_WEAK) else RSER_WEAK


# ---------------------------------------------------------------------------
# ngspice netlist generation
# ---------------------------------------------------------------------------
def build_netlist(caps, fstart=35e3, fstop=55e3, npts=2001, raw_path="out.raw",
                  title="SSH topolectrical chain"):
    """Emit an ngspice .cir string for the given 16-cap bank."""
    lines = [f"* {title}", f"V1 nsrc 0 DC 0 AC {V0}"]
    # 16 grounding inductors, each with explicit series resistor (ngspice has no Rser=)
    for i in range(1, N_NODES + 1):
        lines.append(f"L{i} n{i} nl{i} {L_GND:.6e}")
        lines.append(f"RL{i} nl{i} 0 {RSER_L}")
    # 16 series caps: C1 = nsrc->n1, Ck = n(k-1)->nk, each with explicit series resistor
    for k in range(1, N_CAPS + 1):
        a = "nsrc" if k == 1 else f"n{k-1}"
        b = f"n{k}"
        cval = caps[k - 1]
        lines.append(f"C{k} {a} nc{k} {cval:.8e}")
        lines.append(f"RC{k} nc{k} {b} {_rser_for(cval)}")
    lines.append(".control")
    lines.append("set filetype=ascii")
    lines.append(f"ac lin {npts} {fstart:g} {fstop:g}")
    lines.append(f"write {raw_path} " + " ".join(f"v(n{i})" for i in range(1, N_NODES + 1)))
    lines.append("quit")
    lines.append(".endc")
    lines.append(".end")
    return "\n".join(lines) + "\n"


def _parse_ascii_raw(path):
    """Parse an ngspice ASCII rawfile -> (freqs[Nf], V[Nf, n_vars] complex, names)."""
    with open(path) as fh:
        text = fh.read()
    head, _, body = text.partition("Values:")
    names = []
    nvars = npoints = None
    in_vars = False
    for line in head.splitlines():
        s = line.strip()
        if s.startswith("No. Variables:"):
            nvars = int(s.split(":")[1])
        elif s.startswith("No. Points:"):
            npoints = int(s.split(":")[1])
        elif s.startswith("Variables:"):
            in_vars = True
            rest = s[len("Variables:"):].strip()
            if rest:
                names.append(rest.split()[1])
        elif in_vars:
            parts = s.split()
            if len(parts) >= 2 and parts[0].isdigit():
                names.append(parts[1])
            else:
                in_vars = False
    # Values block: blocks of nvars lines, each "re,im" (first line prefixed by index)
    toks = body.split()
    vals = []
    for t in toks:
        if "," in t:
            re, im = t.split(",")
            vals.append(complex(float(re), float(im)))
    data = np.array(vals, dtype=complex).reshape(npoints, nvars)
    freqs = np.abs(data[:, 0])                # first var is frequency (real)
    return freqs, data[:, 1:], names[1:]


def run_ngspice(caps, **kw):
    """Run a real ngspice AC sweep. Returns (freqs, V[Nf,16], node_names)."""
    with tempfile.TemporaryDirectory() as d:
        raw = os.path.join(d, "out.raw")
        cir = os.path.join(d, "ckt.cir")
        with open(cir, "w") as fh:
            fh.write(build_netlist(caps, raw_path=raw, **kw))
        r = subprocess.run([NGSPICE, "-b", cir], capture_output=True, text=True, timeout=120)
        if not os.path.exists(raw):
            raise RuntimeError(f"ngspice produced no rawfile.\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}")
        freqs, V, names = _parse_ascii_raw(raw)
    order = np.argsort(freqs)
    return freqs[order], V[order], names


# ---------------------------------------------------------------------------
# Independent numpy nodal-admittance solver (cross-check of the SAME circuit)
# ---------------------------------------------------------------------------
def solve_numpy(caps, freqs):
    """Solve the identical topology with a from-scratch nodal-admittance matrix.

    Returns V[Nf, 16] complex node voltages (n1..n16), source held at V0 through C1.
    """
    freqs = np.atleast_1d(np.asarray(freqs, float))
    omega = 2 * np.pi * freqs
    Nf = omega.size
    rser = np.array([_rser_for(c) for c in caps])
    # admittances (Nf, .)
    yL = 1.0 / (RSER_L + 1j * omega[:, None] * L_GND)              # (Nf,16)
    yC = 1.0 / (rser[None, :] + 1.0 / (1j * omega[:, None] * caps[None, :]))  # (Nf,16)

    Y = np.zeros((Nf, N_NODES, N_NODES), dtype=complex)
    di = np.arange(N_NODES)
    Y[:, di, di] += yL
    # cap C1 connects source(known V0) -> n1 : adds to diagonal of n1, RHS gets yC1*V0
    Y[:, 0, 0] += yC[:, 0]
    # caps C2..C16 connect n(k-1)-nk  (k=2..16 -> nodes idx k-2, k-1)
    for k in range(2, N_CAPS + 1):
        a, b = k - 2, k - 1
        Y[:, a, a] += yC[:, k - 1]
        Y[:, b, b] += yC[:, k - 1]
        Y[:, a, b] -= yC[:, k - 1]
        Y[:, b, a] -= yC[:, k - 1]
    rhs = np.zeros((Nf, N_NODES), dtype=complex)
    rhs[:, 0] = yC[:, 0] * V0
    V = np.linalg.solve(Y, rhs[:, :, None])[:, :, 0]
    return V


def input_impedance(caps, freqs):
    """|Zin| seen by the source: Zin = V0 / I, I = yC1*(V0 - V_n1)."""
    freqs = np.atleast_1d(np.asarray(freqs, float))
    omega = 2 * np.pi * freqs
    yC1 = 1.0 / (_rser_for(caps[0]) + 1.0 / (1j * omega * caps[0]))
    V = solve_numpy(caps, freqs)
    I = yC1 * (V0 - V[:, 0])
    return np.abs(V0 / I)


# ---------------------------------------------------------------------------
# Edge-state resonance + voltage profile
# ---------------------------------------------------------------------------
def find_resonance(freqs, V):
    """Edge-state resonance = frequency that maximises |V(n1)| (the driven edge)."""
    i = int(np.argmax(np.abs(V[:, 0])))
    return i, freqs[i]


def voltage_profile_odd(V, idx):
    """Normalised |V| on the odd (A-sublattice) nodes n1,n3,...,n15 at freq index idx."""
    odd = np.arange(0, N_NODES, 2)        # 0-based -> nodes n1,n3,...,n15
    amp = np.abs(V[idx, odd])
    return amp / amp[0]


def pristine_reference():
    """Theoretical pristine decay line on odd nodes: 0.667^(cell-1)."""
    ratio = C_WEAK / C_STRONG
    cells = np.arange(8)                  # 8 odd nodes -> 0..7 cells
    return ratio ** cells


if __name__ == "__main__":
    print(f"Analytic boundary resonance f_R = {F_R/1e3:.4f} kHz")
    caps = pristine_caps()
    freqs, V, names = run_ngspice(caps)
    idx, fpk = find_resonance(freqs, V)
    print(f"ngspice |V(n1)| peak           = {fpk/1e3:.4f} kHz")
    prof = voltage_profile_odd(V, idx)
    ref = pristine_reference()
    print("node  ngspice   theory")
    for j, (p, r) in enumerate(zip(prof, ref)):
        print(f" n{2*j+1:<3d} {p:8.4f} {r:8.4f}")
