"""
plot_repro.py — recreate the two target JPG figures from our ngspice data.

  figures/fig1_AB_BA_pair.png   <- target  b966...JPG  (AB/BA pair impurity)
  figures/fig2_single_AB.png    <- target  edd...JPG   (single A/B impurity)

Each figure is a 2x2 grid:  top row = position sweep, bottom row = cell-3 strength
sweep; normalised |V| on odd nodes (log y) with the pristine 0.667/cell dashed line.
"""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NODES = list(range(1, 16, 2))

# colour palettes matched to the target JPGs
C3 = ["#1f77b4", "#8c564b", "#17becf"]                       # blue, brown, cyan
C4 = ["#1f77b4", "#d62728", "#e377c2", "#17becf"]            # blue, red, pink, cyan


def load():
    data = defaultdict(lambda: defaultdict(dict))   # data[figure][panel][series] = [v...]
    ref = {}
    with open("data/profiles.csv") as fh:
        for r in csv.DictReader(fh):
            if r["figure"] == "ref":
                ref[int(r["node"])] = float(r["voltage"])
                continue
            data[r["figure"]][r["panel"]][r["series"]] = data[r["figure"]][r["panel"]].get(r["series"], {})
            data[r["figure"]][r["panel"]][r["series"]][int(r["node"])] = float(r["voltage"])
    refv = [ref[n] for n in NODES]
    return data, refv


def draw_panel(ax, title, series_dict, refv, palette):
    ax.plot(NODES, refv, "k--", lw=1.6, zorder=1)          # pristine reference
    for i, (label, prof) in enumerate(series_dict.items()):
        ys = [prof[n] for n in NODES]
        ax.plot(NODES, ys, "-o", color=palette[i % len(palette)], ms=6, lw=1.4,
                label=label, zorder=2)
    ax.set_yscale("log")
    ax.set_ylim(3e-2, 1.3)
    ax.set_title(title, fontsize=12)
    ax.set_xticks(NODES)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(fontsize=8, loc="upper right")


def make_figure(data, refv, figkey, panels, palettes, outpath, suptitle):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for ax, (panel, pal) in zip(axes.flat, zip(panels, palettes)):
        draw_panel(ax, panel, data[figkey][panel], refv, pal)
    for ax in axes[:, 0]:
        ax.set_ylabel("Normalised Voltage")
    for ax in axes[1, :]:
        ax.set_xlabel("Node Number")
    fig.suptitle(suptitle, fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outpath, dpi=130)
    print("wrote", outpath)


def main():
    data, refv = load()
    make_figure(data, refv, "fig1",
                ["AB impurity", "BA impurity", "AB cell3", "BA cell3"],
                [C3, C3, C4, C4],
                "figures/fig1_AB_BA_pair.png",
                "SSH edge-state robustness — AB/BA pair impurity (ngspice, 3.3/2.2 nF)")
    make_figure(data, refv, "fig2",
                ["A impurity", "B impurity", "A cell3", "B cell3"],
                [C3, C3, C4, C4],
                "figures/fig2_single_AB.png",
                "SSH edge-state robustness — single A/B impurity (ngspice, 3.3/2.2 nF)")


if __name__ == "__main__":
    main()
