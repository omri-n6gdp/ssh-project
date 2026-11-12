import pandas as pd
import matplotlib.pyplot as plt
import os

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

def plot_sweep(prefix, parameter, held_value=0.1):
    """
    Plot frequency vs cell number for a given sweep file.
    prefix: "AB" or "BA"
    parameter: "gamma" or "delta"
    held_value: value of the other parameter (fixed)
    """
    file_name = f"{prefix}_{parameter}_sweep.xlsx"
    file_path = os.path.join(script_dir, file_name)

    # Check if file exists and is non-empty
    if not os.path.exists(file_path):
        print(f"File not found: {file_name} – skipping")
        return
    if os.path.getsize(file_path) == 0:
        print(f"File is empty: {file_name} – skipping")
        return

    try:
        df = pd.read_excel(file_path, sheet_name="Sheet1", header=None)
    except Exception as e:
        print(f"Error reading {file_name}: {e} – skipping")
        return

    # If after reading there are no data rows (only header maybe), skip
    if df.shape[0] < 2 or df.shape[1] < 2:
        print(f"No data in {file_name} – skipping")
        return

    # Extract cell numbers (first column, rows 1 to end)
    cell_numbers = df.iloc[1:, 0].values.astype(int)

    # Extract sweep parameter values (first row, columns B to E)
    param_values = df.iloc[0, 1:5].values.astype(float)

    # Extract frequency data (rows 1 to end, columns B to E)
    freq_data = df.iloc[1:, 1:5].values  # shape (n_cells, n_params)

    # Create plot
    plt.figure(figsize=(8, 5))
    for i, val in enumerate(param_values):
        plt.plot(cell_numbers, freq_data[:, i], marker='o', label=f"{parameter} = {val}")

    plt.xlabel("Cell number")
    plt.ylabel("Frequency")
    title = (f"{prefix.upper()} sweep: frequency vs. cell number – {parameter} varied, "
             f"other fixed at {held_value}")
    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(title="Swept parameter")
    plt.xticks(cell_numbers)

    # Save figure
    output_name = f"{prefix}_{parameter}_sweep_plot.png"
    output_path = os.path.join(script_dir, output_name)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.show()

# Plot all possible combinations (some files may be missing/empty)
for prefix in ["AB", "BA"]:
    for param in ["gamma", "delta"]:
        plot_sweep(prefix, param)