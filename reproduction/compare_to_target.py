"""compare_to_target.py — stack each target JPG above our ngspice reproduction."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

CLAUDE_SSH = "../claude_ssh"
PAIRS = [
    (f"{CLAUDE_SSH}/b966211c-fd3a-414e-ab1c-9bd9cbd23c3f.JPG",
     "figures/fig1_AB_BA_pair.png", "figures/compare_fig1.png",
     "AB / BA pair impurity"),
    (f"{CLAUDE_SSH}/edd1597d-7aa2-4d3d-be13-d60d9e5854be.JPG",
     "figures/fig2_single_AB.png", "figures/compare_fig2.png",
     "single A / B impurity"),
]

for target, ours, out, name in PAIRS:
    fig, ax = plt.subplots(2, 1, figsize=(10, 14))
    ax[0].imshow(mpimg.imread(target)); ax[0].axis("off")
    ax[0].set_title(f"TARGET (reference JPG) — {name}", fontsize=13)
    ax[1].imshow(mpimg.imread(ours)); ax[1].axis("off")
    ax[1].set_title(f"OUR ngspice reproduction — {name}", fontsize=13)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print("wrote", out)
