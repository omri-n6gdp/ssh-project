# SSH Topolectrical Circuit — Robustness Study (independent reconstruction)

This folder rebuilds the SSH topolectrical-circuit robustness study described
in `~/code/ssh_project/robustness_manual.md` **entirely from scratch**. No data from
`ssh_project/` was reused — every circuit was re-simulated and every number
re-derived here.

## Updated circuit

| Component | Value | Series resistance |
|-----------|-------|-------------------|
| Intracell capacitor `C1` | 2.2 nF | 1 Ω |
| Intercell capacitor `C2` | 3.3 nF | 1 Ω |
| Grounded inductor `L`    | 2200 µH | 1 Ω |

8 unit cells → 16 nodes. Topological boundary resonance: **45.74 kHz**.

## The deliverable

**`SSH_Robustness_Report.pdf`** — a 23-page LaTeX report. One graph per page,
each with a short *Theory in brief* box and three analysis blocks (what the
plot shows / how the data was obtained / quantitative analysis with every
deviation explained). It ends with an appendix: sources, a full table of
numerical results vs. theory, and a glossary. Source: `report/report.tex`.

## How it works

The simulation engine (`ssh_circuit.py`) is a from-scratch AC nodal-admittance
solver: it assembles the 16×16 complex admittance matrix of the circuit at each
frequency and solves the linear system for the node voltages — exactly what a
SPICE `.ac` analysis does. Results are cross-checked against analytic theory at
every step.

## Files

| Path | Contents |
|------|----------|
| `ssh_circuit.py` | the circuit solver + theory functions |
| `common.py` | plot styling and bookkeeping helpers |
| `part1_pristine.py` | pristine baseline — Figures 1–3 |
| `part2_single_impurity.py` | single-impurity sweeps — Figures 4–5 |
| `part3_domain_wall.py` | AB/BA impurity pairs — Figures 6–8 |
| `part4_monte_carlo.py` | Monte Carlo disorder — Figures 9–16 |
| `run_all.py` | runs all four parts in order |
| `figures/` | the 16 figures (PNG) |
| `data/` | CSV data + `results.json` |
| `netlists/` | LTSpice `.cir` netlists of the key circuits |
| `report/` | `report.tex` and the compiled `report.pdf` |

## Reproduce everything

```bash
cd ~/claude_ssh
python3 run_all.py            # regenerates all figures + data (~80 s, 8000-sample Monte Carlo)
cd report && tectonic report.tex   # rebuilds the PDF
```

Random seeds are fixed, so a re-run reproduces every number exactly.
