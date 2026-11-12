"""
Monte Carlo disorder analysis of the SSH topolectrical circuit (Part 4).

For each realization r and disorder strength sigma:
  - Draw 16 Gaussian multipliers (1 + eps_k) for the 16 chain capacitors
    (1 grounded boundary C_b, 8 intracell C1_k, 7 intercell C2_k).
  - Apply to the pristine template netlist, run an AC sweep,
    extract the peak frequency f_peak and peak impedance Z_peak.
  - Read all 16 node voltages |V(N_i)| at the peak index and compute:
      IPR_full  -- sum|V|^4 / (sum|V|^2)^2 over all 16 nodes
      IPR_bulk  -- same restricted to N2..N16 (excludes the source node)
      xi        -- localization length from a log-linear fit on odd
                   nodes N3..N15 (the left-edge mode under N1 drive).

Outputs (all under ~/code/ssh_project/figures):
  monte_carlo_results.csv         per-realization records
  monte_carlo_voltages.csv        per-realization full voltage profile
  mc_summary.csv                  per-sigma mean/std summary
  mc_fpeak_histograms.png
  mc_ipr_histograms.png
  mc_zpeak_histograms.png
  mc_xi_histograms.png
  mc_voltage_envelope.png         pristine + bands of |V|/V0 vs node

References:
  Lee et al. (2018) Communications Physics 1, 39 -- Fig. 5 disorder argument.
  Anzurez et al. (2025) Physica B 715, 417609 -- structured impurity baseline.
  Anderson (1958); Croy & Schreiber (2011) -- 1D disorder methodology.
"""
import argparse
import os
import sys
import time
import warnings
from contextlib import redirect_stdout, redirect_stderr

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyLTSpice import SimRunner, SpiceEditor, RawRead


# ---------- Configuration ----------
LTSPICE = "/Applications/LTspice.app/Contents/MacOS/LTspice"
TEMPLATE = os.path.expanduser("~/code/ssh_project/netlists/mc_template.cir")
RAW_ROOT = os.path.expanduser("~/code/ssh_project/raw/mc")
FIG_DIR = os.path.expanduser("~/code/ssh_project/figures")
os.makedirs(RAW_ROOT, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

N_NODES = 16
N_C1 = 8                       # intracell capacitors C1_1..C1_8
N_C2 = 7                       # intercell capacitors C2_1..C2_7
C1_BASE = 0.1e-9               # 0.1 nF
C2_BASE = 0.22e-9              # 0.22 nF
CB_BASE = 0.22e-9              # 0.22 nF (= C2)


# ---------- Core simulation ----------
def run_realization(runner, sigma, seed):
    """Build, simulate, and post-process one disorder realization."""
    rng = np.random.default_rng(seed)
    if sigma > 0:
        eps_C1 = rng.normal(0.0, sigma, N_C1)
        eps_C2 = rng.normal(0.0, sigma, N_C2)
        eps_Cb = float(rng.normal(0.0, sigma))
    else:
        eps_C1 = np.zeros(N_C1)
        eps_C2 = np.zeros(N_C2)
        eps_Cb = 0.0

    # Clip to keep capacitor values positive (|eps| < 0.9 worst case)
    eps_C1 = np.clip(eps_C1, -0.9, 5.0)
    eps_C2 = np.clip(eps_C2, -0.9, 5.0)
    eps_Cb = float(np.clip(eps_Cb, -0.9, 5.0))

    editor = SpiceEditor(TEMPLATE)
    editor.set_component_value("Cb", CB_BASE * (1.0 + eps_Cb))
    for k, eps in enumerate(eps_C1, 1):
        editor.set_component_value(f"C1_{k}", C1_BASE * (1.0 + eps))
    for k, eps in enumerate(eps_C2, 1):
        editor.set_component_value(f"C2_{k}", C2_BASE * (1.0 + eps))

    tag = f"sig{int(round(sigma * 1000)):04d}_r{seed:08d}"
    out_cir = os.path.join(RAW_ROOT, f"{tag}.cir")
    editor.save_netlist(out_cir)

    # PyLTSpice writes lots of progress to stdout; swallow it
    with open(os.devnull, "w") as devnull, redirect_stdout(devnull), redirect_stderr(devnull):
        raw_path, log_path = runner.run_now(out_cir)

    if raw_path is None or not os.path.exists(raw_path):
        # Simulation failed; clean up and return NaNs
        _safe_remove(out_cir)
        return None

    try:
        r = RawRead(raw_path)
        freqs = np.real(np.asarray(r.get_trace('frequency').get_wave(), dtype=complex))
        v_n1 = np.asarray(r.get_trace('V(n1)').get_wave(), dtype=complex)
        i_v1 = np.asarray(r.get_trace('I(V1)').get_wave(), dtype=complex)
        Z = np.abs(v_n1 / (-i_v1))
        idx = int(np.argmax(Z))
        f_peak, z_peak = _parabolic_peak(freqs, Z, idx)

        V = np.empty(N_NODES)
        for i in range(1, N_NODES + 1):
            wave = np.asarray(r.get_trace(f'V(n{i})').get_wave(), dtype=complex)
            V[i - 1] = float(np.abs(wave[idx]))
    finally:
        _safe_remove(raw_path)
        _safe_remove(log_path)
        _safe_remove(out_cir)
        # PyLTSpice also writes a .net/.op.raw sometimes; clean broadly
        for ext in (".net", ".op.raw"):
            _safe_remove(os.path.join(RAW_ROOT, f"{tag}{ext}"))

    V_norm = V / V[0] if V[0] > 0 else V.copy()

    V2 = V_norm ** 2
    s2 = float(np.sum(V2))
    s4 = float(np.sum(V_norm ** 4))
    IPR_full = s4 / s2 ** 2 if s2 > 0 else float("nan")

    V_bulk = V_norm[1:]
    s2b = float(np.sum(V_bulk ** 2))
    s4b = float(np.sum(V_bulk ** 4))
    IPR_bulk = s4b / s2b ** 2 if s2b > 0 else float("nan")

    odd_idx = np.array([3, 5, 7, 9, 11, 13, 15])
    V_odd = V_norm[odd_idx - 1]
    pos = V_odd > 0
    if pos.sum() >= 3:
        slope, _ = np.polyfit(odd_idx[pos], np.log(V_odd[pos]), 1)
        xi = float(1.0 / abs(slope)) if slope != 0 else float("inf")
    else:
        xi = float("nan")

    return dict(
        f_peak=f_peak,
        z_peak=z_peak,
        IPR_full=IPR_full,
        IPR_bulk=IPR_bulk,
        xi=xi,
        voltages=V_norm,
    )


def _parabolic_peak(freqs, Z, idx):
    """Sub-grid peak via quadratic fit in (log f, log|Z|) around idx.
    Returns (f_peak, |Z|_peak). Falls back to the grid value at the boundary
    or if the fit doesn't have a concave-down vertex.
    """
    n = len(Z)
    if idx <= 0 or idx >= n - 1:
        return float(freqs[idx]), float(Z[idx])
    x = np.log(freqs[idx - 1:idx + 2])
    y = np.log(Z[idx - 1:idx + 2])
    # y = c2 x^2 + c1 x + c0 -- vertex at x = -c1/(2 c2)
    coeffs = np.polyfit(x, y, 2)
    if coeffs[0] >= 0 or not np.isfinite(coeffs[0]):
        return float(freqs[idx]), float(Z[idx])
    x_peak = -coeffs[1] / (2.0 * coeffs[0])
    # Clamp vertex to the bracket [x[0], x[2]] -- otherwise the parabola
    # is wagging; the grid value is more trustworthy.
    if not (x[0] <= x_peak <= x[2]):
        return float(freqs[idx]), float(Z[idx])
    y_peak = float(np.polyval(coeffs, x_peak))
    return float(np.exp(x_peak)), float(np.exp(y_peak))


def _safe_remove(path):
    if not path:
        return
    try:
        os.remove(path)
    except OSError:
        pass


# ---------- Driver ----------
def main():
    p = argparse.ArgumentParser(description="SSH circuit Monte Carlo disorder analysis")
    p.add_argument("--m", type=int, default=500, help="realizations per sigma")
    p.add_argument("--sigmas", type=float, nargs="+", default=[0.01, 0.05, 0.10],
                   help="disorder std deviations (fractional)")
    p.add_argument("--seed-base", type=int, default=10_000)
    p.add_argument("--tag", type=str, default="",
                   help="optional suffix on output filenames")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args()

    suffix = f"_{args.tag}" if args.tag else ""

    runner = SimRunner(simulator=LTSPICE, output_folder=RAW_ROOT)

    # Pristine baseline: one deterministic run at sigma=0
    print("Pristine baseline ...")
    t0 = time.perf_counter()
    pristine = run_realization(runner, sigma=0.0, seed=args.seed_base - 1)
    assert pristine is not None, "pristine run failed -- check template"
    print(f"  f_peak = {pristine['f_peak']/1e3:.3f} kHz  "
          f"|Z|_peak = {pristine['z_peak']:.1f}  "
          f"IPR_full = {pristine['IPR_full']:.4f}  "
          f"IPR_bulk = {pristine['IPR_bulk']:.4f}  "
          f"xi = {pristine['xi']:.3f}  "
          f"(elapsed {time.perf_counter()-t0:.2f} s)")

    records = []
    voltages_records = []
    sim_total = 0
    t_start = time.perf_counter()
    for sigma in args.sigmas:
        print(f"\nsigma = {sigma:.3f}  (M = {args.m})")
        t_sigma = time.perf_counter()
        ok = 0
        fail = 0
        for r_idx in range(args.m):
            seed = args.seed_base + int(round(sigma * 1e6)) * 10 + r_idx
            try:
                res = run_realization(runner, sigma=sigma, seed=seed)
            except Exception as e:
                fail += 1
                if fail < 5:
                    print(f"    realization {r_idx} failed: {e}")
                continue
            if res is None:
                fail += 1
                continue
            ok += 1
            sim_total += 1
            records.append(dict(
                sigma=sigma, realization=r_idx, seed=seed,
                f_peak=res["f_peak"], z_peak=res["z_peak"],
                IPR_full=res["IPR_full"], IPR_bulk=res["IPR_bulk"],
                xi=res["xi"],
            ))
            v_row = dict(sigma=sigma, realization=r_idx, seed=seed)
            for k, val in enumerate(res["voltages"], 1):
                v_row[f"V{k}"] = float(val)
            voltages_records.append(v_row)
            if (r_idx + 1) % max(1, args.m // 10) == 0:
                el = time.perf_counter() - t_sigma
                rate = (r_idx + 1) / max(el, 1e-6)
                print(f"  {r_idx + 1:5d}/{args.m}  ({rate:.2f} sim/s)")
        t_dt = time.perf_counter() - t_sigma
        print(f"  done: ok={ok}, fail={fail}, total {t_dt:.1f} s "
              f"(avg {t_dt/max(args.m,1):.2f} s/sim)")

    total_dt = time.perf_counter() - t_start
    print(f"\nFinished {sim_total} simulations in {total_dt:.1f} s "
          f"({sim_total/max(total_dt,1e-6):.2f} sim/s)")

    df = pd.DataFrame(records)
    df_v = pd.DataFrame(voltages_records)
    csv_path = os.path.join(FIG_DIR, f"monte_carlo_results{suffix}.csv")
    csv_v_path = os.path.join(FIG_DIR, f"monte_carlo_voltages{suffix}.csv")
    df.to_csv(csv_path, index=False)
    df_v.to_csv(csv_v_path, index=False)
    print(f"wrote {csv_path}")
    print(f"wrote {csv_v_path}")

    # Summary stats
    summary = df.groupby("sigma").agg(
        f_peak_mean=("f_peak", "mean"),
        f_peak_std=("f_peak", "std"),
        f_peak_min=("f_peak", "min"),
        f_peak_max=("f_peak", "max"),
        z_peak_mean=("z_peak", "mean"),
        z_peak_std=("z_peak", "std"),
        IPR_full_mean=("IPR_full", "mean"),
        IPR_full_std=("IPR_full", "std"),
        IPR_bulk_mean=("IPR_bulk", "mean"),
        IPR_bulk_std=("IPR_bulk", "std"),
        xi_mean=("xi", "mean"),
        xi_std=("xi", "std"),
        n=("f_peak", "count"),
    )
    summary_path = os.path.join(FIG_DIR, f"mc_summary{suffix}.csv")
    summary.to_csv(summary_path)
    print(f"wrote {summary_path}")
    print("\nSummary:")
    with pd.option_context("display.float_format", "{:.4f}".format,
                           "display.width", 160):
        print(summary)

    print(f"\nPristine reference -- "
          f"f_peak = {pristine['f_peak']/1e3:.4f} kHz, "
          f"Z = {pristine['z_peak']:.1f}, "
          f"IPR_full = {pristine['IPR_full']:.4f}, "
          f"IPR_bulk = {pristine['IPR_bulk']:.4f}, "
          f"xi = {pristine['xi']:.3f}")

    if not args.no_plots:
        plot_all(df, pristine, args.sigmas, suffix=suffix)


# ---------- Plotting ----------
def plot_all(df, pristine, sigmas, suffix=""):
    sigmas = sorted(sigmas)

    # f_peak
    fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
    if len(sigmas) == 1:
        axes = [axes]
    for ax, sig in zip(axes, sigmas):
        sub = df[df.sigma == sig]
        ax.hist(sub.f_peak / 1e3, bins=40, alpha=0.75, edgecolor="black")
        ax.axvline(pristine["f_peak"] / 1e3, color="red", ls="--",
                   label=f"pristine = {pristine['f_peak']/1e3:.2f} kHz")
        ax.set_title(f"sigma = {sig*100:.0f}%")
        ax.set_xlabel("f_peak [kHz]")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("count")
    fig.suptitle("Monte Carlo distribution of f_peak under capacitor disorder")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, f"mc_fpeak_histograms{suffix}.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    plt.close(fig)

    # IPR (full + bulk on twin axes)
    fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
    if len(sigmas) == 1:
        axes = [axes]
    for ax, sig in zip(axes, sigmas):
        sub = df[df.sigma == sig]
        ax.hist(sub.IPR_bulk, bins=40, alpha=0.7, edgecolor="black",
                color="C2", label="IPR (N2-N16)")
        ax.axvline(pristine["IPR_bulk"], color="red", ls="--",
                   label=f"pristine = {pristine['IPR_bulk']:.3f}")
        ax.axvline(1 / (N_NODES - 1), color="gray", ls=":",
                   label=f"uniform = {1/(N_NODES-1):.3f}")
        ax.set_title(f"sigma = {sig*100:.0f}%")
        ax.set_xlabel("IPR (bulk)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("count")
    fig.suptitle("Monte Carlo IPR (bulk) -- topological-protection check")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, f"mc_ipr_histograms{suffix}.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    plt.close(fig)

    # Z_peak
    fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
    if len(sigmas) == 1:
        axes = [axes]
    for ax, sig in zip(axes, sigmas):
        sub = df[df.sigma == sig]
        ax.hist(sub.z_peak, bins=40, alpha=0.75, edgecolor="black", color="C1")
        ax.axvline(pristine["z_peak"], color="red", ls="--",
                   label=f"pristine = {pristine['z_peak']:.0f}")
        ax.set_title(f"sigma = {sig*100:.0f}%")
        ax.set_xlabel("|Z|_peak")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("count")
    fig.suptitle("Monte Carlo distribution of peak impedance")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, f"mc_zpeak_histograms{suffix}.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    plt.close(fig)

    # xi
    fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 4), sharey=True)
    if len(sigmas) == 1:
        axes = [axes]
    for ax, sig in zip(axes, sigmas):
        sub = df[df.sigma == sig].dropna(subset=["xi"])
        # cap xi at reasonable upper bound for plotting
        xi_arr = np.clip(sub.xi.values, 0, 30)
        ax.hist(xi_arr, bins=40, alpha=0.75, edgecolor="black", color="C3")
        ax.axvline(pristine["xi"], color="red", ls="--",
                   label=f"pristine = {pristine['xi']:.2f}")
        ax.set_title(f"sigma = {sig*100:.0f}%")
        ax.set_xlabel("xi (localization length)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("count")
    fig.suptitle("Monte Carlo distribution of localization length xi")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, f"mc_xi_histograms{suffix}.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
