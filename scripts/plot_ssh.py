import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ── Load data ──────────────────────────────────────────────────────────────────
data = np.loadtxt('/Users/omrishahar/code/ssh_project/netlists/0521.txt', skiprows=1)

time   = data[:, 0]          # seconds
currents = data[:, 1:]       # all 15 diode columns

# Column headers (order as they appear in the file)
labels = [
    "I(D1)",  "I(D2)",  "I(D3)",  "I(D4)",  "I(D5)",  "I(D6)",  "I(D8)",
    "I(D9)",  "I(D10)", "I(D11)", "I(D12)", "I(D13)", "I(D14)", "I(D15)", "I(D16)"
]

# Original column order from the file
original_labels = [
    "I(D3)",  "I(D5)",  "I(D8)",  "I(D10)", "I(D12)", "I(D14)", "I(D16)",
    "I(D1)",  "I(D15)", "I(D13)", "I(D11)", "I(D9)",  "I(D6)",  "I(D4)",  "I(D2)"
]

sorted_indices = [original_labels.index(l) for l in labels]
currents = currents[:, sorted_indices]

# ── Convert units for readability ──────────────────────────────────────────────
time_us   = time * 1e6          # seconds → microseconds
current_mA = currents * 1e3     # amperes → milliamperes

colours = cm.tab20(np.linspace(0, 1, len(labels)))  # ← add this line

# ── Max current per diode ──────────────────────────────────────────────────────
max_currents = np.max(current_mA, axis=0)   # shape: (15,)

fig2, ax2 = plt.subplots(figsize=(10, 5))

ax2.plot(labels, max_currents, color="steelblue", linewidth=1.5,
         marker="o", markersize=6, markerfacecolor="white", markeredgewidth=1.5)

ax2.set_yscale("log")
ax2.set_xlabel("Diode", fontsize=13)
ax2.set_ylabel("Max Current (mA, log scale)", fontsize=13)
ax2.set_title("Maximum Current per Diode", fontsize=14, fontweight="bold")
ax2.tick_params(axis="x", rotation=45)
ax2.grid(True, which="both", linestyle="--", alpha=0.4)

plt.tight_layout()

# ── Plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 6))

for i, (label, colour) in enumerate(zip(labels, colours)):
    ax.plot(time_us, current_mA[:, i], label=label, color=colour, linewidth=1.2)

# ── Labels, title, legend ──────────────────────────────────────────────────────
ax.set_xlabel("Time (µs)", fontsize=13)
ax.set_ylabel("Current (mA)", fontsize=13)
ax.set_title("Diode Currents vs Time (LTSpice Simulation)", fontsize=14, fontweight="bold")

ax.legend(
    loc="upper right",
    ncol=3,
    fontsize=8,
    framealpha=0.9,
    title="Diodes",
    title_fontsize=9
)

ax.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.show()