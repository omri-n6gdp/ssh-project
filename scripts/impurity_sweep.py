import pandas as pd
import matplotlib.pyplot as plt
import os

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, "impurity_sweep.xlsx")

# Read the Excel file without assuming headers
df = pd.read_excel(file_path, sheet_name="Sheet1", header=None)

# Extract cell numbers (first column, rows 1 to 6) -> index 1..6
cell_numbers = df.iloc[1:7, 0].values.astype(int)  # [2,3,4,5,6,7]

# Extract delta values (first row, columns B to F)
delta_values = df.iloc[0, 1:6].values.astype(float)  # [-0.3, -0.1, 0, 0.1, 0.3]

# Extract frequency data: rows 1..6, columns B..F
freq_data = df.iloc[1:7, 1:6].values  # shape (6,5)

# Plot each delta as a separate line
plt.figure(figsize=(8, 5))
for i, delta in enumerate(delta_values):
    plt.plot(cell_numbers, freq_data[:, i], marker='o', label=f"Δ = {delta}")

# Labels and title
plt.xlabel("Cell number")
plt.ylabel("Frequency")
plt.title("Frequency vs. Cell Number for Different Impurity Values (Δ)")
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(title="Impurity parameter")

# Optional: set x-axis ticks to integer cell numbers
plt.xticks(cell_numbers)

# Save the figure
output_path = os.path.join(script_dir, "impurity_sweep_plot.png")
plt.savefig(output_path, dpi=150, bbox_inches="tight")
print(f"Plot saved to: {output_path}")

# Show the plot
plt.show()