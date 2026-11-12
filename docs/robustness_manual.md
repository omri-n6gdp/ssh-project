# Robustness Analysis of the SSH Topolectrical Circuit
## A Complete Step-by-Step Manual

**Scope.** This manual takes you from an empty LTSpice and Python environment to a finished robustness study comprising four parts: (1) pristine baseline reproduction, (2) single-impurity parametric sweeps, (3) domain-wall (AB / BA) configurations, and (4) Monte Carlo disorder analysis. Every step assumes no prior knowledge of the simulation tooling beyond having LTSpice and Python 3 installed on a Mac.

**Reference papers.** The two papers anchoring this work are:
- Lee, C. H., Imhof, S., Berger, C., Bayer, F., Brehm, J., Molenkamp, L. W., Kiessling, T., & Thomale, R. (2018). Topolectrical Circuits. *Communications Physics*, **1**, 39. DOI: 10.1038/s42005-018-0035-2.
- Anzurez, J., Luna Reyes, I., Pérez-Pedraza, J. C., Raya, A., Sáenz-Trujillo, A. U., Tapia Sánchez, R., & Tinoco-Santillán, G. (2025). Edge state robustness in the SSH model with impurities: Topoelectric circuit simulations toward topological heterostructures. *Physica B: Condensed Matter*, **715**, 417609. DOI: 10.1016/j.physb.2025.417609.

**Pedagogical references** for the SSH model physics, used throughout: Asbóth, Oroszlány, & Pályi, *A Short Course on Topological Insulators*, Lecture Notes in Physics **919**, Springer (2016), arXiv:1509.02295; Su, Schrieffer, & Heeger, "Solitons in polyacetylene," *Phys. Rev. Lett.* **42**, 1698 (1979).

**Numerical tooling references**: Analog Devices LTspice documentation, in particular the AC analysis and `.step`/`.meas` guidance (Analog Devices, 2024–2025); PyLTSpice (Brum, 2024) for batch automation; NumPy and Matplotlib for post-processing.

---

# Part 0 — Software Setup

You should already have LTSpice installed. If not, download it from analog.com. Everything else is Python.

## 0.1 Install Python dependencies

Open the macOS Terminal and run:

```bash
python3 -m venv ~/ssh_robustness_env
source ~/ssh_robustness_env/bin/activate
pip install --upgrade pip
pip install numpy scipy matplotlib pandas PyLTSpice
```

This installs:
- **numpy / scipy** — numerical arrays and the eigensolver used for analytic cross-checks.
- **matplotlib** — plotting.
- **pandas** — small convenience for tabulating sweep results.
- **PyLTSpice** — programmatically edits LTSpice netlists and runs LTSpice in batch mode, returning the raw waveform data to Python (Brum, 2024). This is what makes parametric sweeps and the Monte Carlo loop tractable.

Verify the install:

```bash
python3 -c "import PyLTSpice; import numpy; import matplotlib; print('All OK')"
```

## 0.2 Locate the LTSpice executable

PyLTSpice needs to know where LTSpice lives. On macOS the binary is usually at:

```
/Applications/LTspice.app/Contents/MacOS/LTspice
```

Confirm by running it in Terminal:

```bash
/Applications/LTspice.app/Contents/MacOS/LTspice -h
```

If that fails, find it with `mdfind -name LTspice.app`. Note the path; you'll need it later.

## 0.3 Create the project folder

```bash
mkdir -p ~/code/ssh_project/{netlists,raw,figures,scripts}
cd ~/code/ssh_project
```

This gives you four working folders:
- `netlists/` for `.cir` and `.asc` files,
- `raw/` for LTSpice's raw simulation output,
- `figures/` for output plots,
- `scripts/` for the Python files you'll write.

---

# Part 1 — Pristine Baseline Reproduction

**Goal.** Build a working SSH topolectrical circuit in LTSpice with N = 8 unit cells (16 nodes), reproduce the topological boundary resonance (TBR), and confirm the voltage profile predicted by Lee et al. (2018) and Anzurez et al. (2025). This baseline becomes the reference all later results compare against.

## 1.1 Why N = 8 unit cells (and not 10 or 20)

Anzurez et al. (2025, §3) explicitly use N = 8 unit cells (16 sites), noting "A larger number of sites has a very mild effect on the results." Lee et al. (2018, Fig. 2) used N = 10. Eight unit cells is the smallest size where (a) edge cells 1 and 8 can be left pristine, (b) impurities can be placed in cells 2–7, and (c) the exponential edge-state decay is clearly resolvable on a log plot. We use N = 8 for direct comparison with the Anzurez paper.

## 1.2 The circuit topology

Refer back to Fig. 2a of Lee et al. (2018) and Fig. 2 of Anzurez et al. (2025). The topology is:

- A horizontal chain of **alternating capacitors** C_1 and C_2 (16 capacitors total for N = 8). The pattern is C_1, C_2, C_1, C_2, ... starting from the left. C_1 represents the intracell hopping v; C_2 represents the intercell hopping w.
- Between every pair of adjacent capacitors there is a **node** (numbered 1 through 16). Lee et al. (2018) confirm there are 16 such nodes for an N = 8 chain (their Fig. 2c shows nodes 1 through 15, which is the 16-node case minus the input terminal).
- Each node is connected to ground through an **inductor L = 10 μH**. So 16 inductors in total.
- A **sinusoidal AC voltage source V_0** is connected to node 1 with the other terminal grounded. This is the impedance probe.
- **Boundary condition (critical).** Lee et al. (2018, Fig. 2a) show that the **leftmost capacitor is grounded** — i.e. one terminal of the chain's first capacitor C_1 connects to node 1, the other to ground. This grounded boundary capacitor is what your earlier verification of the Lee circuit identified as essential for the correct TBR. Equivalently, in the Anzurez convention, node 1 sits between two capacitors of the same value (C_2 on the ground side, then C_1 going right), so the circuit Laplacian has equal diagonal entries at every node. Without this grounded boundary capacitor, the diagonal of the Laplacian at node 1 differs from the rest and the resonance shifts substantially.

## 1.3 Component values (Anzurez et al. parameter set)

From Anzurez et al. (2025, §3):

| Component | Value | Role |
|-----------|-------|------|
| L | 10 μH | All grounded inductors |
| C_1 | 0.22 μF | Intracell capacitor (analog of v) |
| C_2 | 0.33 μF | Intercell capacitor (analog of w) |
| V_0 | 10 V | AC source amplitude |

Note that C_1 < C_2, which places the circuit in the **topologically non-trivial phase** (Lee et al. 2018: the SSH midgap mode exists for t = C_1/C_2 < 1). The dimensionless ratio is t = 0.22 / 0.33 ≈ 0.667.

## 1.4 Predicted resonant frequency

From Lee et al. (2018, Methods) and Anzurez et al. (2025, Eq. 12), the topological boundary resonance sits at:

f_R = 1 / (2π √(L (C_1 + C_2)))

Plugging in:
- L (C_1 + C_2) = 10×10⁻⁶ × 0.55×10⁻⁶ = 5.5×10⁻¹²
- √ = 2.345×10⁻⁶
- 2π √ = 1.473×10⁻⁵
- **f_R ≈ 67.9 kHz**

This is the number you should see in your AC sweep. Write it on a sticky note.

## 1.5 Open LTSpice and start a new schematic

1. Launch LTSpice from Applications.
2. **File → New Schematic**.
3. Save it immediately to `~/code/ssh_project/netlists/ssh_pristine.asc`.

The Mac version of LTSpice has the same functionality as the Windows version but a different toolbar layout. The hotkeys differ slightly; the most useful ones are listed in the LTSpice on macOS help menu.

## 1.6 Build the schematic — slow walkthrough

You'll place 16 inductors, 17 capacitors (16 in the chain + 1 grounded boundary), and one AC source. Working systematically left to right prevents you from miscounting.

### 1.6.1 Place the AC source

1. **Edit → Component** (or press **F2**).
2. Search for `voltage`, double-click. Place it on the left side of the schematic.
3. Right-click the source. Set:
   - **DC value**: 0
   - **Small signal AC analysis**: tick "AC Amplitude" and enter `10` (volts).
4. Wire the **bottom terminal to ground** (**Edit → Place GND**, or **F4** then click).

### 1.6.2 Place the grounded boundary capacitor (the critical one)

This is the capacitor that ensures the Laplacian boundary condition is correct (your earlier finding confirmed by Lee et al. 2018 Fig. 2a). It has value **C_2 = 0.33 μF** and sits between **node 1** and **ground**.

Why C_2 here, not C_1? Lee et al. (2018, Fig. 2a) place the grounded capacitor before the first inductor; in the convention where the first capacitor in the chain going right is C_1, the grounded boundary capacitor must be C_2 so that the unit-cell repeating pattern (C_1, C_2, C_1, C_2, ...) starts cleanly at node 1. Equivalently: every node in the chain must have one C_1 and one C_2 connecting to it.

1. **F2** → search `cap` → place a capacitor.
2. Wire the top terminal where node 1 will be (right of the AC source).
3. Wire the bottom terminal to ground.
4. Right-click → set value to `0.33u` (LTSpice convention: `u` means micro).
5. Right-click → rename to `Cb` (b for "boundary").

### 1.6.3 Place the first inductor at node 1

1. **F2** → search `ind` → place. Inductors in LTSpice point top-to-bottom by default; that's what you want for the grounded inductors.
2. Place it just to the right of `Cb`. Wire its top terminal to the same node as `Cb`'s top and the AC source's top (this is node 1).
3. Wire the bottom terminal to ground.
4. Set value to `10u`. Rename to `L1`.

### 1.6.4 Place capacitor C_1 between node 1 and node 2

1. **F2** → cap → place. Rotate it 90° (**Ctrl+R**) so it's horizontal.
2. Wire its left terminal to node 1. Its right terminal becomes node 2.
3. Set value to `0.22u`. Rename to `C1_1`.

### 1.6.5 Place inductor L_2 at node 2

Same as L_1. Top terminal at node 2, bottom to ground. Value `10u`. Rename `L2`.

### 1.6.6 Place capacitor C_2 between node 2 and node 3

Horizontal cap, value `0.33u`, rename `C2_1`.

### 1.6.7 Repeat the unit cell

The repeating unit-cell pattern from here is:
- Inductor to ground at each new node,
- C_1 (0.22 μF) between odd→even nodes (i.e. node 3→4, 5→6, ...),
- C_2 (0.33 μF) between even→odd nodes (i.e. node 4→5, 6→7, ...).

Continue placing components until you have:
- 16 inductors: L1, L2, ..., L16, each `10u`, each going from its node to ground.
- Capacitors between adjacent nodes following the alternating C_1, C_2 pattern. With nodes 1 through 16, that's 15 inter-node capacitors named C1_1, C2_1, C1_2, C2_2, ..., C1_8 (or C2_7 depending on parity — see next note).

**Counting check.** For N = 8 unit cells with 16 nodes, indexed 1–16:
- Capacitors C_1 (value 0.22 μF) go between nodes (1,2), (3,4), (5,6), (7,8), (9,10), (11,12), (13,14), (15,16) → 8 of them.
- Capacitors C_2 (value 0.33 μF) go between nodes (2,3), (4,5), (6,7), (8,9), (10,11), (12,13), (14,15) → 7 of them.
- Plus the **grounded boundary** C_b (0.33 μF) attached to node 1.
- Total capacitors: 16. Total inductors: 16. 

Don't terminate the right end with an extra grounded capacitor — Lee et al. (2018) and Anzurez et al. (2025) use **open boundary conditions** on the right. The right edge is just the last node 16 with its grounded inductor and nothing further to the right.

### 1.6.8 Label the nodes

LTSpice auto-numbers nodes, but for plotting and analysis you want named nodes.

1. **Edit → Label Net** (or **F4** on the wire).
2. Click on each node in turn and assign a label: `N1`, `N2`, ..., `N16`.

This labeling is essential because you'll later read the node voltages by name from the simulation output.

## 1.7 Add a series resistance to each inductor (optional but realistic)

A perfectly ideal LC circuit has infinite Q and the resonance peak becomes a delta function. Real inductors have a small DC resistance. Lee et al. (2018, Methods and Fig. 5) explicitly note that "one has to take into account ... capacitive and resistive losses." Adding a small ESR makes the peak finite and the simulation more numerically stable.

For each inductor, right-click → set **Series Resistance** to `0.1` (i.e. 0.1 Ω). Don't add resistors as separate components; use the inductor's built-in `Rser` field. With Rser = 0.1 Ω and L = 10 μH, the inductor's Q at f_R = 67.9 kHz is Q = 2π f L / R ≈ 42, which is realistic.

If you want a sharper ideal-limit peak for the baseline check, leave Rser = 0 for now and add it later.

## 1.8 Add the AC analysis directive

1. **Edit → Spice Directive** (or **S** key).
2. Type:

   ```
   .ac dec 200 1k 1Meg
   ```

   This runs an AC sweep from 1 kHz to 1 MHz with 200 points per decade on a logarithmic frequency axis. The Analog Devices LTspice AC analysis documentation confirms that decade sweeps with 100–1000 points per decade are standard for resolving sharp resonances. With 200 points/decade across three decades you get 600 frequency points, which resolves a Q ≈ 40 peak comfortably.

3. Click anywhere on the schematic to place the directive.

## 1.9 Run the simulation

1. **Simulate → Run** (or click the runner icon).
2. The waveform viewer opens.
3. Click on node N1 in the schematic. The plot pane shows V(N1) versus frequency.
4. The voltage at node 1 with the source drive is V_0; the *impedance* between node 1 and ground is Z = V_0 / I_source. Since LTSpice sources a fixed-amplitude AC voltage, the natural impedance probe is to measure the current through V1 (right-click the V_source → look for "current"). Plot **−I(V1)** (LTSpice's current convention is into the + terminal; flipping the sign gives the current you draw from the source). Then **Z = V_0 / I_source = 10 / I(V1)**.

   To plot impedance directly: in the plot pane, right-click → **Add Trace** → enter `V(N1)/(-I(V1))`. You should see a clear peak around 68 kHz.

5. **Sanity check.** Read the peak frequency. If it's at 67.9 kHz ± 0.1 kHz, your circuit is correct. If it's anywhere else (e.g. 89 kHz, 107 kHz), you have a boundary condition or component error — go back to Step 1.6.2 and check the grounded boundary capacitor C_b.

## 1.10 Voltage profile measurement (reproducing Anzurez et al. Fig. 9)

For the voltage profile you don't need a sweep — you just need the complex node voltages at exactly the resonant frequency.

1. Add a second SPICE directive (**S** key):

   ```
   .ac list 67864
   ```

   This tells LTSpice to perform AC analysis at the single frequency 67,864 Hz (your predicted f_R). The `list` keyword takes a list of frequencies; with one value it runs at exactly that frequency. This is much cleaner than reading the sweep at one point because there's no interpolation.

2. Comment out or delete the original `.ac dec ...` directive (you can only run one AC analysis at a time). Right-click → "Comment out" prefixes a `;` to disable it.

3. Run again. LTSpice now computes complex voltages at every node at 67,864 Hz.

4. Use the **Tools → Export data** menu, or simply hover over each node and read off |V(N_i)| from the cursor. For systematic analysis, export to a text file.

5. **Plot |V(N_i)| / V_0 versus i on log y-axis.** You should see:
   - Even nodes (2, 4, 6, ..., 16): voltage starts near 0 at the left and rises approximately exponentially toward the right edge. This is the right-edge localization of the topological mode (Anzurez et al. 2025, Fig. 9).
   - Odd nodes (1, 3, 5, ..., 15): voltage starts at 1.0 at node 1 and falls approximately exponentially to the right. This is the left-edge localization, one order of magnitude smaller scale (Anzurez et al. 2025, Fig. 9).

If the profile is qualitatively right (alternating, exponential, two-edge structure), your baseline is fully verified. Save this schematic — every later step is a modification of it.

## 1.11 Generate the corresponding netlist

Some later steps (Part 4 in particular) need the textual netlist, not the schematic.

1. With the schematic open, **View → SPICE Netlist** (Mac LTSpice: this is sometimes under the Simulate menu).
2. Save the netlist as `~/code/ssh_project/netlists/ssh_pristine.cir`.

The netlist is a plain text file you can read in a text editor; it lists every component, its terminals, and its value. This is what PyLTSpice manipulates programmatically.

---

# Part 2 — Single-Impurity Parametric Sweeps

**Goal.** Reproduce the resonance-displacement and voltage-localization results of Anzurez et al. (2025) Figs. 4–6 and 10–11 for cases 1 and 2 (single-sublattice impurities). Specifically: place an impurity in cells 2–7 with strength δ ∈ {−0.3, −0.1, +0.1, +0.3}, and for each case extract (a) the resonance frequency shift f/f_R and (b) the voltage profile.

## 2.1 How to encode an impurity in the circuit

From Anzurez et al. (2025, Table 1):

- **A-site impurity in unit cell n_i**: the two capacitors of cell n_i become C'_1 = C_1(1+δ), C'_2 = C_2(1+δ). Both capacitors in cell n_i scale together with parameter δ.
- **B-site impurity in unit cell n_i**: the two capacitors of cell n_i become C''_1 = C_1(1+γ), C''_2 = C_2(1+γ). Same form but with parameter γ (which we set equal to δ for cases 1 and 2).

In the schematic, **unit cell n_i** corresponds to:
- the C_1 capacitor between nodes (2n_i − 1) and (2n_i), and
- the C_2 capacitor between nodes (2n_i) and (2n_i + 1).

So for example, an impurity in cell 3 modifies C1_3 (between nodes 5 and 6) and C2_3 (between nodes 6 and 7).

**Mapping table (memorise this or print it):**

| Unit cell n_i | C_1 between nodes | C_2 between nodes |
|---------------|-------------------|-------------------|
| 1 | 1–2 | 2–3 |
| 2 | 3–4 | 4–5 |
| 3 | 5–6 | 6–7 |
| 4 | 7–8 | 8–9 |
| 5 | 9–10 | 10–11 |
| 6 | 11–12 | 12–13 |
| 7 | 13–14 | 14–15 |
| 8 | 15–16 | — |

## 2.2 Parameterize the schematic

Instead of editing capacitor values manually for every sweep point, you'll parameterize them so a single `.step` directive runs the whole sweep.

1. In `ssh_pristine.asc`, **right-click** on the C_1 capacitor in the cell you want to impurify (say cell 3, so right-click on capacitor C1_3 between nodes 5 and 6).
2. Change its value from `0.22u` to `{C1_imp}`. The curly braces tell LTSpice this is a parameter, as documented in the Analog Devices LTspice parameter-sweep guide.
3. Right-click on C2_3 (between nodes 6 and 7) and change its value to `{C2_imp}`.
4. Add a SPICE directive (**S** key):

   ```
   .param delta = 0
   .param C1_imp = 0.22u * (1 + delta)
   .param C2_imp = 0.33u * (1 + delta)
   ```

5. Add another SPICE directive for the sweep:

   ```
   .step param delta list -0.3 -0.1 0 0.1 0.3
   ```

   This runs the AC analysis five times, once per value of δ.

6. Save the schematic as `ssh_impurity_cell3_A.asc`.

## 2.3 Run the sweep and extract f/f_R

1. Re-enable the `.ac dec 200 1k 1Meg` directive (comment out the `.ac list` line for now).
2. **Simulate → Run.**
3. LTSpice will run five simulations. The waveform viewer shows one trace per δ value.
4. Plot the impedance `V(N1)/(-I(V1))` again. You should see five curves, the middle one (δ = 0) at 67.9 kHz, and the others displaced.
5. To extract the peak frequency programmatically per step, add `.meas` directives. The Analog Devices `.meas` documentation describes the FIND/WHEN syntax for extremum detection. Add this directive:

   ```
   .meas AC fpeak FIND frequency WHEN mag(V(N1)/(-I(V1)))=MAX
   ```

   After the run, **View → SPICE Error Log** (Cmd+L on Mac) shows the measured `fpeak` for each step. Copy these into a spreadsheet or read them with PyLTSpice in Part 2.5.

6. Compute f/f_R for each δ by dividing by 67864 Hz. Plot f/f_R vs δ — this reproduces one curve of Anzurez et al. (2025) Fig. 4.

## 2.4 Sweep over impurity position

To reproduce all four panels of Fig. 4 in Anzurez et al. (2025), you also need to sweep the impurity *position* (cells 2–7). LTSpice can't easily move a parameter "around" the schematic, so the cleanest approach is to make six schematics (one per impurity position) and run each one.

Save copies:

```
ssh_impurity_cell2_A.asc
ssh_impurity_cell3_A.asc
...
ssh_impurity_cell7_A.asc
```

In each, parameterize a different pair of capacitors (the C_1 and C_2 of that cell) and keep the rest pristine. Run each, extract f_peak for the five δ values via `.meas`, and accumulate into a table:

| Cell n_i | δ = −0.3 | δ = −0.1 | δ = +0.1 | δ = +0.3 |
|----------|----------|----------|----------|----------|
| 2 | ? | ? | ? | ? |
| 3 | ? | ? | ? | ? |
| ... | | | | |
| 7 | ? | ? | ? | ? |

The expected qualitative behaviour from Anzurez et al. (2025, §3.1):
- Positive δ → resonance shifts to **lower** frequency (f/f_R < 1).
- Negative δ → resonance shifts to **higher** frequency (f/f_R > 1).
- Cell 2 shows the largest shift; cell 7 shows almost no shift.
- A-impurities (cases 1) are more symmetric in sign than B-impurities (case 2).

## 2.5 Automate the cell-position sweep with PyLTSpice

Doing this by hand is tedious. Here's the automated version. Save as `~/code/ssh_project/scripts/single_impurity_sweep.py`:

```python
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
F_R = 67864.0  # pristine resonance, Hz
C1_BASE = 0.22e-6
C2_BASE = 0.33e-6

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
```

**A few practical notes about the script:**

- You'll need to make sure your pristine `.cir` has consistently named capacitors C1_1..C1_8 and C2_1..C2_7. The auto-generated netlist from LTSpice uses whatever names you gave components in the schematic. Open the `.cir` in a text editor and rename if needed before running the script.
- `SpiceEditor` is part of PyLTSpice; the API is documented in the PyLTSpice readthedocs (Brum, 2024).
- LTSpice on Mac runs the simulator non-interactively when launched with a `.cir` argument; PyLTSpice handles this.
- The peak detection uses `argmax`, which is crude — if your resonance is narrow and barely sampled, refine by fitting a Lorentzian around the peak. For 200 points/decade this is rarely needed.

## 2.6 Voltage localization sweeps (replicating Fig. 10–11)

For voltage profiles, instead of sweeping frequency you fix it at f_R and read the 16 node voltages. Modify the script above so that the SPICE simulation uses `.ac list 67864` instead of the decade sweep, then save and plot |V(N_i)|/V_0 versus i on a log y-axis.

A small subtlety: when you change a capacitor far from the source, the resonance frequency itself shifts (that's Part 2.3). For the voltage profile, you have a choice:
1. **Drive at fixed f_R**: this measures how the original mode responds when perturbed. The voltage at f_R will drop because the resonance has moved away. This is what is shown in Lee et al. (2018) Fig. 2c, where the voltage profile is sampled at the unperturbed ω̃.
2. **Drive at the shifted peak f_peak(δ)**: this measures the mode that *now* exists in the perturbed circuit, which is closer to what Anzurez et al. (2025) Figs. 10–11 show.

Use option 2 (drive at the shifted peak) for direct comparison with Anzurez. You'll need to chain two AC analyses per impurity configuration: first a decade sweep to find f_peak, then a `.ac list f_peak` to extract the node voltages. Or do both in one script: run the sweep, identify f_peak, then run a transient at that frequency. PyLTSpice handles the back-to-back runs without restarting LTSpice between them.

The plotting code:

```python
# After identifying f_peak for a given (cell, delta) configuration:
editor = SpiceEditor(NETLIST_TEMPLATE)
editor.set_component_value(c1_name(cell), C1_BASE * (1 + delta))
editor.set_component_value(c2_name(cell), C2_BASE * (1 + delta))
# Replace the .ac directive
editor.set_simulation_command(f".ac list {f_peak}")
# Run, then read all 16 node voltages
runner.run_now(...)
r = RawRead(raw)
V = np.array([np.abs(r.get_trace(f'V(n{i})').get_wave()[0]) for i in range(1, 17)])
V_normalized = V / 10.0  # divide by V_0 = 10 V

# Plot
fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogy(range(1, 17), V_normalized, 'o-', label=f'cell {cell}, δ={delta:+.2f}')
ax.semilogy(range(1, 17), V_pristine, 'r--', label='pristine', alpha=0.5)
ax.set_xlabel('Node index')
ax.set_ylabel(r'|V| / $V_0$')
ax.legend()
plt.savefig(...)
```

You should see the same qualitative features Anzurez reports: the **even-site exponential rise toward the right edge** remains nearly intact for A-impurities, while **B-impurities disrupt odd-site voltages** (Anzurez et al. 2025, §3.2 and Fig. 12).

---

# Part 3 — Domain-Wall Regime (AB and BA Configurations)

**Goal.** Reproduce Anzurez et al. (2025) Figs. 7–8 (frequency displacement) and Figs. 13–14 (voltage localization with two distinct contributions and emergent domain walls).

## 3.1 What an AB-ordered pair of impurities does

From Anzurez et al. (2025, Eq. 8 and Table 1), an **AB pair** in unit cell n_i means: the **B-site of cell n_i is impurified with strength δ**, and the **A-site of cell n_i + 1 is impurified with strength γ**. The hoppings around the pair change as:

- w' = w (1 + δ)   ← the C_2 capacitor immediately *left* of the pair, between nodes (2n_i − 1) and (2n_i)... wait, let's be careful here.

In Anzurez et al.'s convention (their Fig. 1d for AB and Fig. 1e for BA): the impurities replace two adjacent sites, and the hopping *between* the two impurity sites becomes:

- v''' = v [1 + (δ + γ)/2]   for the AB case (the bond between the impurities is type-v),

while the hoppings to the neighbours on either side become w' and w''.

Translating to the circuit (capacitors substituting for hoppings):

- The C_2 capacitor *to the left* of the AB pair (between cell n_i − 1's B-site and cell n_i's A-site, but with the A-site of cell n_i replaced) becomes C_2(1 + δ).
- The C_1 capacitor *between* the two impurities (which sits inside the supercell, between the B of cell n_i and the A of cell n_i+1 — actually no, between the A and B of cell n_i' where the impurities sit, depending on indexing) becomes C_1[1 + (δ + γ)/2].
- The C_2 capacitor *to the right* of the AB pair becomes C_2(1 + γ).

The cleanest way to avoid confusion is to read Table 1 of Anzurez et al. (2025) line by line and identify exactly which three capacitors change. The cases listed there are:

- **Model 1A (single A impurity):** C'_1 = C_1(1 + δ), C'_2 = C_2(1 + δ). Two capacitors change.
- **Model 1B (single B impurity):** C''_1 = C_1(1 + γ), C''_2 = C_2(1 + γ). Two capacitors change.
- **Model AB (two impurities, A then B):** C'_2 = C_2(1 + δ), C''_2 = C_2(1 + γ), C'''_1 = C_1[1 + (δ + γ)/2]. **Three capacitors change.**
- **Model BA (two impurities, B then A):** C'_1 = C_1(1 + δ), C''_1 = C_1(1 + γ), C'''_2 = C_2[1 + (δ + γ)/2]. **Three capacitors change.**

## 3.2 Practical procedure

Anzurez et al. (2025) place the AB and BA pair starting in **cell 3** (their Figs. 7 and 8 title says "site 3"). For consistency, follow them.

1. Open `ssh_pristine.asc`, save as `ssh_AB_pair_cell3.asc`.
2. Identify the three capacitors that need to change for an AB pair starting at cell 3:
   - The C_2 to the *left* of cell 3's A-site: this is C2_2 (between nodes 4 and 5).
   - The C_1 *inside* the AB pair: this is C1_3 (between nodes 5 and 6).
   - The C_2 to the *right* of cell 3's B-site: this is C2_3 (between nodes 6 and 7).
3. Parameterize these three:
   - C2_2 = `{C2_left}`
   - C1_3 = `{C1_mid}`
   - C2_3 = `{C2_right}`
4. Add SPICE directives:

   ```
   .param delta = 0.1
   .param gamma = -0.1
   .param C2_left  = 0.33u * (1 + delta)
   .param C1_mid   = 0.22u * (1 + (delta + gamma)/2)
   .param C2_right = 0.33u * (1 + gamma)
   ```

5. To sweep γ at fixed δ (Anzurez et al. Fig. 7 protocol):

   ```
   .step param gamma list -0.3 -0.1 0.1 0.3
   ```

   Run the simulation. You get four impedance curves. Extract f_peak / f_R for each.

6. Repeat with δ as the swept parameter at fixed γ, etc., to cover the full grid in Fig. 7 and Fig. 8.

7. For the **BA pair**, the three modified capacitors are different (see Table 1). Save as `ssh_BA_pair_cell3.asc` and adapt.

## 3.3 Voltage profile for the domain wall

Repeat the Part 2.6 procedure (find f_peak, then `.ac list f_peak`, then read all 16 node voltages, plot on log axis). The signature you're looking for (Anzurez et al. 2025, §3.2 closing paragraphs and Figs. 13–14):

- For **two consecutive impurities**, the voltage near the right edge shows **two distinct contributions**: one from the original right-edge localization and one localized near the impurity pair. This is the **emergent domain wall** — charge accumulation around the impurities.
- **BA ordering gives a more uniform distribution** than AB ordering (Anzurez et al. 2025, Fig. 13 comparison).

This is the most physically interesting bit of the Anzurez paper because it ties directly to the **domain-wall soliton states** of the original Su-Schrieffer-Heeger (1979) polyacetylene model. In SSH theory, a domain wall between two regions of opposite dimerization hosts a mid-gap localized state — the "kink" or soliton. Here the AB/BA impurity pair acts as a local dimerization flip, and you're seeing the circuit-level analog of that soliton (Asbóth, Oroszlány, Pályi 2016, Ch. 1; Su, Schrieffer, Heeger 1979).

## 3.4 Position-resolved domain-wall sweep

Anzurez et al. (2025) Figs. 13–14 also examine how the voltage profile changes when the AB/BA pair is placed in different cells (not just cell 3). The procedure is identical to Part 2.4: make six schematics, one per starting cell, parameterize the three capacitors each time, sweep δ and γ, and extract the voltage profiles.

If you want the full grid of (start cell, δ, γ) automated, extend the PyLTSpice script from Part 2.5. The only structural change is that now three capacitors per configuration get modified instead of two.

---

# Part 4 — Monte Carlo Disorder Analysis

**Goal.** Move beyond *structured* impurities (Parts 2 and 3) to **random, uncorrelated disorder** on every capacitor in the chain, and ask the quantitative robustness question: how does the TBR survive when many realistic component variations are present simultaneously?

This is the part that goes beyond Anzurez et al. and connects to Lee et al. (2018, Fig. 5), where they explicitly compare 1% and 10% disorder tolerances and argue that topological circuits with extensive boundary mode degeneracy are markedly more disorder-robust than non-topological alternatives.

## 4.1 Why this is the right question to ask

In a real laboratory circuit, every capacitor has a manufacturing tolerance. Class-1 ceramic C0G/NP0 capacitors come at ±1% tolerance; the more common Class-2 X7R types are ±5% to ±10% (PCBSync 2026; weishielectronics 2024). When you build the SSH circuit on a breadboard or PCB, every capacitor's actual value is drawn from a distribution around its nominal value. The structured A/B/AB/BA impurities of Anzurez et al. are a *worst-case structured* perturbation; random disorder is the *realistic* perturbation.

Theoretically, this is closely related to **Anderson localization** in 1D systems (Anderson 1958; Mott & Twose 1961). In 1D, *all* eigenstates of a disordered system are localized for any finite disorder strength (Anderson 1958; Mott & Twose 1961). For the topological SSH model, the question is whether the *protected* edge state remains localized at the edge under disorder, or whether disorder mixes it with bulk states and washes it out. The bulk–boundary correspondence theorem (Asbóth, Oroszlány, Pályi 2016, Ch. 1) guarantees survival as long as **(a) the bulk gap doesn't close and (b) chiral symmetry isn't broken**. Capacitance disorder preserves chiral symmetry (it only modifies off-diagonal terms of the Laplacian); inductance disorder breaks it (it modifies diagonal terms).

So the Monte Carlo gives two qualitatively different scenarios:
1. **Capacitor-only disorder (chiral-preserving):** edge state should be very robust.
2. **Inductor-only disorder (chiral-breaking):** edge state should be more fragile.

You can do both. The contrast is the most informative result.

## 4.2 The Monte Carlo recipe

For each realisation r = 1, ..., M (typically M = 500 to 5000):

1. Draw 15 independent random values for the chain capacitors:
   - C_1 capacitors: C_1^(r,k) = 0.22 μF × (1 + ε_k), ε_k ~ N(0, σ²), k = 1, ..., 8.
   - C_2 capacitors: C_2^(r,k) = 0.33 μF × (1 + ε'_k), ε'_k ~ N(0, σ²), k = 1, ..., 7.
   - And the grounded boundary C_b: same kind of perturbation.
2. (Optionally also perturb inductors.)
3. Run the AC simulation, extract:
   - **f_peak**: the resonance peak frequency,
   - **Z_peak**: the peak impedance height (proxy for the Q-factor and the edge-state quality),
   - **Voltage profile** at f_peak.
4. Compute **per-realisation localization metrics** (see 4.4 below).
5. Repeat M times.
6. Plot histograms / scatter of f_peak, Z_peak, and the localization metrics.

The Gaussian choice of perturbation is standard (Akhanjee 2008, arXiv:0810.4532 §I; Croy & Schreiber 2011, arXiv:1103.2866 §2). A uniform box distribution ε ~ U(−W, W) is equivalent up to a factor; the standard convention is Gaussian with standard deviation σ matching the nominal tolerance. So **σ = 0.05 corresponds to a ±5% (1-sigma) tolerance band**, which approximately matches common ceramic capacitor tolerance specs.

## 4.3 Choice of M (number of realisations)

The standard deviation of the mean of M realisations scales as 1/√M. Practical guidelines from the disorder literature:
- For **histograms of f_peak**, M = 500 is enough to resolve a smooth distribution.
- For **ensemble averages** of localization length, M = 2000–5000 gives 1–2% precision (Croy & Schreiber 2011 used L = 2¹⁹ chain sizes; for our small N = 8 chain, we need more realisations rather than longer chains).
- For your project, **start with M = 500** to verify everything runs, then scale up to M = 2000 for the final result.

At ~1 second per simulation (typical for this circuit size), M = 2000 takes ~30 minutes. Workable.

## 4.4 Localization metrics

You have the voltage profile |V(N_i)|/V_0 for i = 1, ..., 16. Two standard scalar summaries:

**Inverse Participation Ratio (IPR).** Defined as (Akhanjee 2008; Murphy et al. 2010, in Emergent Mind 2026; Asbóth, Oroszlány, Pályi 2016, problem set):

$$\text{IPR} = \frac{\sum_i |V_i|^4}{\left(\sum_i |V_i|^2\right)^2}$$

A perfectly delocalized state has IPR = 1/N (with N = 16, that's 0.0625). A perfectly localized state on one site has IPR = 1. The topological edge state typically has IPR ~ 0.1–0.3 in N = 16.

**Localization length ξ.** Fit |V_i| ∝ exp(−i / ξ) at even sites for i = N down to some cutoff. From Anzurez et al. (2025, §3.2): for the pristine SSH circuit, even-site voltages grow exponentially from left to right with a characteristic scale. ξ is the decay length of the topological edge state. Robustness of the topological phase means **ξ remains small (well-localized) under disorder**.

For each Monte Carlo realisation, compute IPR and ξ. The ensemble distributions tell you how disorder affects edge-state localization.

## 4.5 Full Monte Carlo script

Save as `~/code/ssh_project/scripts/monte_carlo.py`:

```python
"""
Monte Carlo disorder analysis of the SSH topolectrical circuit.
Generates M realisations of capacitor disorder, extracts f_peak,
Z_peak, voltage profile, IPR, and localization length.

References:
- Lee et al. 2018 (disordered topolectrical circuits, Fig. 5)
- Anderson 1958 (1D Anderson localization)
- Croy & Schreiber 2011 (numerical disorder methods)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyLTSpice import SimRunner, SpiceEditor, RawRead

# ---- Configuration ----
LTSPICE_PATH = "/Applications/LTspice.app/Contents/MacOS/LTspice"
TEMPLATE = os.path.expanduser("~/code/ssh_project/netlists/ssh_pristine.cir")
OUTDIR = os.path.expanduser("~/code/ssh_project/raw/mc")
os.makedirs(OUTDIR, exist_ok=True)

M = 500                    # number of realisations - start small
SIGMA_LIST = [0.01, 0.05, 0.10]   # 1%, 5%, 10% disorder
N_NODES = 16
F_R_NOMINAL = 67864.0

C1_BASE = 0.22e-6
C2_BASE = 0.33e-6
N_C1 = 8   # number of intracell capacitors
N_C2 = 7   # number of intercell capacitors

# Per-realisation result storage
def run_realisation(sigma, seed):
    rng = np.random.default_rng(seed)
    eps_C1 = rng.normal(0, sigma, N_C1)
    eps_C2 = rng.normal(0, sigma, N_C2)
    eps_Cb = rng.normal(0, sigma)

    editor = SpiceEditor(TEMPLATE)
    for k in range(1, N_C1 + 1):
        editor.set_component_value(f"C1_{k}", C1_BASE * (1 + eps_C1[k-1]))
    for k in range(1, N_C2 + 1):
        editor.set_component_value(f"C2_{k}", C2_BASE * (1 + eps_C2[k-1]))
    editor.set_component_value("Cb", C2_BASE * (1 + eps_Cb))

    # Use a sweep that brackets f_R; tighter than the 1k-1Meg pristine sweep
    editor.set_simulation_command(".ac dec 400 30k 200k")

    out_cir = os.path.join(OUTDIR, f"mc_sig{int(sigma*100):02d}_r{seed:04d}.cir")
    editor.save_netlist(out_cir)

    runner = SimRunner(simulator=LTSPICE_PATH, output_folder=OUTDIR)
    raw, log = runner.run_now(out_cir)
    r = RawRead(raw)

    freqs = np.real(r.get_trace('frequency').get_wave())
    v_n1 = r.get_trace('V(n1)').get_wave()
    i_v1 = r.get_trace('I(V1)').get_wave()
    Z = np.abs(v_n1 / (-i_v1))
    idx = np.argmax(Z)
    f_peak = freqs[idx]
    z_peak = Z[idx]

    # Voltage profile at peak: re-run a single-frequency AC
    editor.set_simulation_command(f".ac list {f_peak}")
    editor.save_netlist(out_cir)
    raw2, _ = runner.run_now(out_cir)
    r2 = RawRead(raw2)
    V = np.array([np.abs(r2.get_trace(f'V(n{i})').get_wave()[0])
                  for i in range(1, N_NODES + 1)])

    # IPR
    V2 = V**2
    IPR = np.sum(V**4) / (np.sum(V2))**2

    # Localization length: fit log|V_even| vs even-site index
    even_idx = np.arange(2, N_NODES + 1, 2)  # nodes 2, 4, ..., 16
    V_even = V[even_idx - 1]
    # Fit only the right half (where exponential growth is monotone)
    # site index measured from right edge: i' = 16 - i
    log_V = np.log(V_even)
    site = even_idx
    # Linear fit log|V| = a + b * site; localization length = 1/|b|
    coeffs = np.polyfit(site, log_V, 1)
    xi = 1.0 / np.abs(coeffs[0]) if coeffs[0] != 0 else np.inf

    return f_peak, z_peak, V, IPR, xi


# ---- Pristine baseline run ----
print("Running pristine baseline...")
f0, z0, V0_profile, IPR0, xi0 = run_realisation(0.0, seed=0)
print(f"  f_peak = {f0/1e3:.3f} kHz, Z_peak = {z0:.1f}, IPR = {IPR0:.4f}, xi = {xi0:.3f}")

# ---- Monte Carlo loop ----
all_results = []
for sigma in SIGMA_LIST:
    print(f"\nMonte Carlo at sigma = {sigma*100:.0f}%, M = {M} realisations:")
    for r in range(M):
        try:
            fp, zp, Vp, IPR, xi = run_realisation(sigma, seed=10000 + r)
            all_results.append({'sigma': sigma, 'realisation': r,
                                'f_peak': fp, 'Z_peak': zp,
                                'IPR': IPR, 'xi': xi})
        except Exception as e:
            print(f"  realisation {r} failed: {e}")
        if r % 50 == 0:
            print(f"  done {r}/{M}")

df = pd.DataFrame(all_results)
df.to_csv(os.path.expanduser("~/code/ssh_project/figures/monte_carlo_results.csv"),
          index=False)

# ---- Plotting ----
# 1) f_peak distribution per sigma
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
for ax, sigma in zip(axes, SIGMA_LIST):
    sub = df[df.sigma == sigma]
    ax.hist(sub.f_peak / 1e3, bins=40, alpha=0.7, edgecolor='black')
    ax.axvline(F_R_NOMINAL/1e3, color='red', linestyle='--',
               label=f'pristine f_R = {F_R_NOMINAL/1e3:.2f} kHz')
    ax.set_xlabel('f_peak (kHz)')
    ax.set_title(f'σ = {sigma*100:.0f}% disorder')
    ax.legend()
axes[0].set_ylabel('Count')
plt.tight_layout()
plt.savefig(os.path.expanduser("~/code/ssh_project/figures/mc_fpeak_histograms.png"),
            dpi=120)

# 2) IPR distribution
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
for ax, sigma in zip(axes, SIGMA_LIST):
    sub = df[df.sigma == sigma]
    ax.hist(sub.IPR, bins=40, alpha=0.7, edgecolor='black', color='C2')
    ax.axvline(IPR0, color='red', linestyle='--',
               label=f'pristine IPR = {IPR0:.3f}')
    ax.set_xlabel('IPR')
    ax.set_title(f'σ = {sigma*100:.0f}% disorder')
    ax.legend()
axes[0].set_ylabel('Count')
plt.tight_layout()
plt.savefig(os.path.expanduser("~/code/ssh_project/figures/mc_ipr_histograms.png"),
            dpi=120)

# 3) Summary statistics
summary = df.groupby('sigma').agg({
    'f_peak': ['mean', 'std', 'min', 'max'],
    'Z_peak': ['mean', 'std'],
    'IPR': ['mean', 'std'],
    'xi': ['mean', 'std']
})
print("\nSummary statistics:")
print(summary)
summary.to_csv(os.path.expanduser("~/code/ssh_project/figures/mc_summary.csv"))

plt.show()
```

## 4.6 Interpreting the Monte Carlo output

Three figures and a CSV table fall out of the script. Things to look for:

1. **f_peak histograms.** At σ = 1%, the distribution should be very narrow (FWHM well under 1% of f_R). At σ = 10%, the distribution broadens substantially (FWHM ~5–10% of f_R), but the peak should still be clearly identifiable. This directly verifies the Lee et al. (2018, Fig. 5b) claim that the SSH circuit's TBR is robust against 10% component disorder.

2. **IPR histograms.** The pristine IPR is around 0.1–0.2 (well localized). Under disorder, the IPR distribution can either stay tight (topological protection) or broaden and shift toward 1/N = 0.0625 (delocalization). For pure capacitor disorder (chiral-preserving), the IPR should stay tight. This is your direct numerical check of bulk-boundary correspondence.

3. **Z_peak.** The peak impedance height drops with disorder (because energy bleeds into multiple modes and the resonance broadens). The *rate* at which Z_peak falls with σ tells you how much loss is acceptable in practice.

4. **Comparison to Anzurez:** Anzurez et al. (2025) studied structured impurities only. Your Monte Carlo extends their analysis to randomly placed disorder of the same magnitude. If you find that random disorder of strength σ = 0.1 (i.e. δ_typical ≈ 0.1) gives a smaller resonance shift than a single structured impurity with δ = 0.1, that's the rigorous statement of "topological protection averages out worst-case structured perturbations." This is the genuinely novel contribution your project adds beyond the literature.

## 4.7 Two extensions to consider

**(a) Inductor disorder (chiral-breaking).** Repeat the Monte Carlo but perturb the inductors L_i instead of (or in addition to) the capacitors. Inductor disorder modifies the *diagonal* of the circuit Laplacian, which corresponds to **on-site disorder** in the TB model — and on-site disorder breaks chiral symmetry. You should observe that the f_peak distribution is much broader and the IPR distribution shifts toward delocalization much more rapidly with σ. The contrast with the capacitor-only case is the cleanest demonstration of chiral-symmetry protection (Asbóth, Oroszlány, Pályi 2016, Ch. 1 §1.5; Anderson 1958).

**(b) Trivial-phase comparison.** Swap C_1 and C_2 in your pristine circuit so that t = C_1/C_2 > 1, putting the circuit in the **topologically trivial phase** (Lee et al. 2018 Fig. 2b, dashed curves). Run the same Monte Carlo. The trivial phase has no protected edge state, so disorder should destroy any sharp impedance feature very quickly. This is exactly the comparison Lee et al. (2018, Fig. 5) made between their topological and non-topological circuits. Demonstrating it on your own setup makes the topological-protection argument complete.

---

# Part 5 — Putting It All Together

By the end of these four parts you'll have:

- A verified pristine LTSpice schematic reproducing the Anzurez et al. (2025) parameter set TBR at 67.9 kHz.
- A grid of single-impurity sweep results (6 cells × 5 deltas × 2 sublattices = 60 simulations) reproducing the resonance-displacement trends of Anzurez et al. Figs. 4–6.
- Voltage profiles for each configuration replicating Figs. 10–11.
- AB- and BA-pair sweeps reproducing Figs. 7–8 and the domain-wall signatures of Figs. 13–14.
- A Monte Carlo robustness analysis at three disorder strengths (1%, 5%, 10%) with histograms of f_peak and IPR — a result that goes beyond the published literature and directly tests the Lee et al. (2018, Fig. 5) protection argument with your own circuit.

## Suggested ordering of the work

| Week | Task | Output |
|------|------|--------|
| 1 | Part 1 (pristine) | `ssh_pristine.asc`, baseline TBR plot, voltage profile |
| 2 | Part 2 (single impurity) | Sweep CSV + replicated Fig. 4 |
| 3 | Part 3 (AB/BA pair) | Domain-wall voltage profiles + replicated Figs. 13–14 |
| 4 | Part 4 (Monte Carlo) | MC histograms + summary table |
| 5 | Cross-checks + write-up | Comparison plots, extensions if time |

Estimated wall-clock simulation time at M = 2000: ~30 min. Estimated total active work time across all parts: ~30–40 hours.

---

# Part 6 — Common Pitfalls and Sanity Checks

These are mistakes that cost real debugging time if not anticipated:

1. **The boundary capacitor.** If your pristine TBR is not at 67.9 kHz to within 1%, recheck the grounded boundary capacitor C_b at node 1. The choice between attaching it to the source side or omitting it shifts the TBR by tens of percent. Your earlier successful Lee-paper reproduction confirms this.

2. **Capacitor naming consistency.** PyLTSpice modifies components by their netlist designators (`C1_1`, `C2_3`, etc.). If you let LTSpice auto-name them as `C1`, `C2`, `C3`, ..., you'll lose the cell-index information. Always rename systematically before exporting the netlist.

3. **Frequency-domain vs single-point AC.** The two `.ac` directives — `.ac dec ...` and `.ac list ...` — cannot both be active simultaneously. Comment out the unused one each time you switch.

4. **"Current convention" gotcha.** LTSpice defines the current through a voltage source as flowing into the + terminal. So the source current you're measuring is `-I(V1)`, not `I(V1)`. Get this sign wrong and your impedance is negative.

5. **Drive frequency vs natural frequency for voltage profile.** When you perturb the circuit, the resonance frequency shifts. If you drive at the *old* (unperturbed) f_R, the voltage profile will be the response of a now-off-resonance circuit, which is small everywhere. Always re-find f_peak first, then drive there.

6. **Monte Carlo seeds.** Use distinct seeds for each realisation (`seed=10000+r` in the script). Re-using seeds gives correlated samples and underestimates the variance.

7. **Sample size in IPR.** IPR has known high-variance issues (Bauer et al. 2015, arXiv:1506.02048) at small N. With N = 16 nodes, expect IPR to fluctuate substantially across realisations. The ensemble distribution is what matters, not any single number.

8. **Negative ε draws.** If σ is large enough (say σ = 0.3) you can draw ε < −1, meaning a *negative capacitance*. LTSpice will complain or produce nonsense. Either clip the distribution (truncate to ε > −0.9), use a log-normal distribution, or stick to σ ≤ 0.15. Real capacitor tolerances never exceed ±50%, so σ ≤ 0.15 (1-sigma) covers all realistic cases.

9. **Series resistance.** Setting Rser = 0 on inductors gives an idealised infinite-Q peak. AC simulators handle this fine analytically, but the peak becomes infinitely narrow and harder to resolve numerically. Adding Rser = 0.1 Ω makes the peak finite and gives Q ≈ 40. Both choices are fine; report which you use.

10. **Don't forget the inductor at every node.** Easy to miss one when building 16 of them. After placing them all, count: there should be exactly N_nodes = 16 inductors. A missing inductor at node k creates a Laplacian asymmetry that shifts the resonance by ~1 kHz, big enough to mislead you.

---

# References Cited

1. **Anderson, P. W.** (1958). Absence of diffusion in certain random lattices. *Physical Review*, **109**, 1492.

2. **Anzurez, J., Luna Reyes, I., Pérez-Pedraza, J. C., Raya, A., Sáenz-Trujillo, A. U., Tapia Sánchez, R., & Tinoco-Santillán, G.** (2025). Edge state robustness in the SSH model with impurities: Topoelectric circuit simulations toward topological heterostructures. *Physica B: Condensed Matter*, **715**, 417609.

3. **Asbóth, J. K., Oroszlány, L., & Pályi, A.** (2016). *A Short Course on Topological Insulators: Band-structure Topology and Edge States in One and Two Dimensions*. Lecture Notes in Physics, vol. 919, Springer. arXiv:1509.02295.

4. **Akhanjee, S.** (2008). Perturbative criteria for Anderson localization in long-ranged 1D tight-binding models. arXiv:0810.4532.

5. **Bauer, B., et al.** (2015). Moments of the inverse participation ratio for the Laplacian on finite regular graphs. arXiv:1506.02048.

6. **Brum, N.** (2024). PyLTSpice: Python utilities to interact with LTSpice and NGSpice. Python package documentation, https://pyltspice.readthedocs.io.

7. **Croy, A., & Schreiber, M.** (2011). Anderson Localization in 1D Systems with Correlated Disorder. arXiv:1103.2866.

8. **Analog Devices** (2024). LTspice: Using the .STEP Command to Perform Repeated Analysis. Technical Article. https://www.analog.com/en/resources/technical-articles/ltspice-using-the-step-command-to-perform-repeated-analysis.html.

9. **Analog Devices** (2025). LTspice: AC Analysis Using The .STEP Command. Technical Article. https://www.analog.com/en/resources/technical-articles/ltspice-ac-analysis-using-the-step-command.html.

10. **Analog Devices** (2024). LTspice: Using .MEAS and .STEP Commands to Calculate Efficiency. Technical Article. https://www.analog.com/en/resources/technical-articles/ltspice-using-meas-and-step-commands-to-calculate-efficiency.html.

11. **Lee, C. H., Imhof, S., Berger, C., Bayer, F., Brehm, J., Molenkamp, L. W., Kiessling, T., & Thomale, R.** (2018). Topolectrical Circuits. *Communications Physics*, **1**, 39.

12. **Mott, N. F., & Twose, W. D.** (1961). The theory of impurity conduction. *Advances in Physics*, **10**, 107–163.

13. **PCBSync** (2026). Capacitor Tolerance: Accuracy & Application Selection. Technical reference, https://pcbsync.com/capacitor-tolerance/.

14. **Su, W. P., Schrieffer, J. R., & Heeger, A. J.** (1979). Solitons in polyacetylene. *Physical Review Letters*, **42**, 1698.

