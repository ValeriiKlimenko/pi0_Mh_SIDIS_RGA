import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# --- bin edges (must match your analysis) ---
Z_EDGES  = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
PT2_EDGES= [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.5]

def bin_center(edges, i1b):
    """1-based bin index -> center from edges."""
    if 1 <= i1b < len(edges):
        return 0.5 * (edges[i1b - 1] + edges[i1b])
    return None

# where the CSVs from the ROOT code live
input_dir = "plots_phi_unfold_fit_data_sim/csv"
output_dir = "matplotlib_plots_grid"
os.makedirs(output_dir, exist_ok=True)

# find all xq2 CSV files
csv_files = sorted(glob.glob(os.path.join(input_dir, "xq2_*_data_sim.csv")))
if not csv_files:
    raise RuntimeError(f"No CSV files found in {input_dir}")

for csv_path in csv_files:
    df = pd.read_csv(csv_path)

    # assume xq2 is constant within a file
    xq2 = int(df["xq2"].iloc[0])

    # unique sorted bin indices
    z_bins   = sorted(df["z_bin"].unique())
    pt2_bins = sorted(df["pt2_bin"].unique())
    nrows, ncols = len(z_bins), len(pt2_bins)

    # create subplot grid: share X only (Y different per plot)
    fig, axes = plt.subplots(
        nrows, ncols,
        sharex=True, sharey=False,
        figsize=(ncols * 2.2, nrows * 2.2)
    )

    # make axes 2D for all cases
    if nrows == 1 and ncols == 1:
        axes = [[axes]]
    elif nrows == 1:
        axes = [axes]
    elif ncols == 1:
        axes = [[ax] for ax in axes]

    for i, z in enumerate(z_bins):
        z_ctr = bin_center(Z_EDGES, z)

        for j, pt2 in enumerate(pt2_bins):
            ax = axes[i][j]
            pt2_ctr = bin_center(PT2_EDGES, pt2)

            cell = df[(df["z_bin"] == z) & (df["pt2_bin"] == pt2)]
            if cell.empty:
                ax.set_axis_off()
                continue

            data = cell[cell["source"] == "data"]
            sim  = cell[cell["source"] == "sim_scaled"]

            # --- normalization factors (unit integral per curve) ---
            data_sum = float(data["y"].sum()) if not data.empty else 0.0
            sim_sum  = float(sim["y"].sum())  if not sim.empty  else 0.0

            data_scale = 1.0 / data_sum if data_sum > 0.0 else 1.0
            sim_scale  = 1.0 / sim_sum  if sim_sum  > 0.0 else 1.0

            # data: circles, normalized
            if not data.empty and data_sum > 0.0:
                ax.errorbar(
                    data["phi_deg"],
                    data["y"] * data_scale,
                    yerr=data["yerr"] * data_scale,
                    fmt="o", markersize=3, linewidth=1, label="data"
                )

            # sim: squares, dashed, normalized
            if not sim.empty and sim_sum > 0.0:
                ax.errorbar(
                    sim["phi_deg"],
                    sim["y"] * sim_scale,
                    yerr=sim["yerr"] * sim_scale,
                    fmt="s--", markersize=3, linewidth=1, label="sim (scaled)"
                )

            # x-label only on bottom row (to avoid clutter)
            if i == nrows - 1:
                ax.set_xlabel(r"$\phi_{\mathrm{Trento}}$ [deg]", fontsize=8)

            # Y axis on EVERY subplot
            ax.set_ylabel("Normalized counts", fontsize=8)

            # pT² center above EVERY subplot
            if pt2_ctr is not None:
                ax.set_title(rf"$p_T^2 \approx {pt2_ctr:.2f}$", fontsize=8)

            # row label: z bin center (leftmost column, inside plot)
            if j == 0 and z_ctr is not None:
                ax.text(
                    0.02, 0.95,
                    rf"$z \approx {z_ctr:.2f}$",
                    transform=ax.transAxes,
                    fontsize=8, va="top"
                )

            ax.grid(True, linewidth=0.3)

    # global legend (one per figure)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", linewidth=1, label="data"),
        Line2D([0], [0], marker="s", linestyle="--", linewidth=1, label="sim (scaled)")
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=8)

    fig.suptitle(f"xQ² bin {xq2}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 0.98, 0.95])

    out_png = os.path.join(output_dir, f"xq2_{xq2}_grid.png")
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved {out_png}")
