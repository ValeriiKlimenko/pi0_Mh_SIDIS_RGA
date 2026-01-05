#!/usr/bin/env python3
import csv
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# User options
# ----------------------------------------------------------------------
CSV_FILE = "phi_fit_results_scaled.csv"      # first CSV: fit results
NORM_CSV_FILE = "../unfolding_dis/unfolded_first16.csv"       # second CSV: unfolded spectrum

# z bin edges (must match your analysis!)
Z_EDGES = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
PT2_EDGES = [
    0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50,
    0.65, 0.80, 1.00, 1.50
]

# how many ixq2 bins to show, and layout
N_IX_SHOW = 16
NROWS, NCOLS = 4, 4  # 4x4 grid

# ----------------------------------------------------------------------
# Helper to compute bin width from iz, ipt (for backward compatibility)
# ----------------------------------------------------------------------
def compute_bin_width(iz, ipt):
    """
    iz, ipt: 1-based indices (like in ROOT).
    bin_width = (z_hi - z_lo) * (pt2_hi - pt2_lo)
    """
    if iz < 1 or iz >= len(Z_EDGES):
        return 0.0
    if ipt < 1 or ipt >= len(PT2_EDGES):
        return 0.0
    z_lo = Z_EDGES[iz - 1]
    z_hi = Z_EDGES[iz]
    pt_lo = PT2_EDGES[ipt - 1]
    pt_hi = PT2_EDGES[ipt]
    return (z_hi - z_lo) * (pt_hi - pt_lo)

# ----------------------------------------------------------------------
# Read first CSV into a list of dicts
# Expecting columns:
#   ixq2,iz,ipt,p0,mean_nonzero,no_limit_hit,bin_width
#   (bin_width is new; if missing, we compute it from iz, ipt)
# ----------------------------------------------------------------------
records = []
with open(CSV_FILE, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        ixq2 = int(row["ixq2"])
        iz = int(row["iz"])
        ipt = int(row["ipt"])
        p0 = float(row["p0"])
        mean = float(row["mean_nonzero"])
        no_limit_hit = int(row.get("no_limit_hit", "1"))

        # bin_width from CSV if present, otherwise compute it
        if "bin_width" in row and row["bin_width"] not in (None, "", "nan"):
            bin_width = float(row["bin_width"])
        else:
            bin_width = compute_bin_width(iz, ipt)

        # Guard against zero or negative bin_width (should not happen, but be safe)
        if bin_width > 0.0:
            p0_over_bw = 2 * 3.1415926 * p0 / bin_width
            mean_over_bw =  2 * 3.1415926 * mean / bin_width
        else:
            p0_over_bw = np.nan
            mean_over_bw = np.nan

        rec = {
            "ixq2": ixq2,
            "iz": iz,
            "ipt": ipt,
            "p0": p0,
            "mean": mean,
            "no_limit_hit": no_limit_hit,
            "bin_width": bin_width,
            "p0_over_bw": p0_over_bw,
            "mean_over_bw": mean_over_bw,
        }
        records.append(rec)

if not records:
    raise RuntimeError("No rows read from CSV; check CSV_FILE and header names.")

# ----------------------------------------------------------------------
# Read second CSV (unfolded spectrum) and build normalization map
# We will divide p0_over_bw and mean_over_bw by "content" from this file.
# ixq2 corresponds to (bin + 1) in the second CSV, as requested.
# ----------------------------------------------------------------------
norm_content_by_bin = {}
try:
    with open(NORM_CSV_FILE, newline="") as f_norm:
        norm_reader = csv.DictReader(f_norm)
        for row in norm_reader:
            # expect columns: bin, x_center, content, error
            try:
                b = int(row["bin"])
                c = float(row["content"])
            except (KeyError, ValueError):
                continue
            norm_content_by_bin[b] = c
except FileNotFoundError as e:
    raise RuntimeError(f"Could not open unfolding CSV file '{NORM_CSV_FILE}': {e}")

if not norm_content_by_bin:
    raise RuntimeError(
        "No rows read from unfolding CSV; check NORM_CSV_FILE and header names."
    )

# Map ixq2 -> content (normalization) using: bin = ixq2 + 1
NORM_BY_IXQ2 = {}
all_ixq2_values = {r["ixq2"] for r in records}
for ixq2 in all_ixq2_values:
    bin_idx = ixq2 + 1  # ixq2 corresponds to bin+1 in unfolded CSV
    content = norm_content_by_bin.get(bin_idx)
    if content is None:
        print(
            f"WARNING: no unfolded content for ixq2={ixq2} "
            f"(looking for bin={bin_idx} in {NORM_CSV_FILE}). "
            "This ixq2 bin will NOT be rescaled."
        )
        continue
    if content <= 0.0:
        print(
            f"WARNING: non-positive unfolded content for ixq2={ixq2} "
            f"(bin={bin_idx}, content={content}). "
            "This ixq2 bin will NOT be rescaled."
        )
        continue
    NORM_BY_IXQ2[ixq2] = content

# ----------------------------------------------------------------------
# Global color maps
# Use matplotlib "tab10" discrete colors: blue, orange, green, ...
# ----------------------------------------------------------------------
base_colors = plt.cm.tab10.colors
n_base_colors = len(base_colors)

# --- ipt -> color (for plots vs z) ---
all_ipt_global = sorted({r["ipt"] for r in records})
if len(all_ipt_global) > 2:
    # only bins that will actually be plotted (we drop first and last later)
    plot_ipt_global = all_ipt_global[1:-1]
else:
    plot_ipt_global = all_ipt_global

GLOBAL_COLOR_MAP_IPT = {
    ipt: base_colors[i % n_base_colors] for i, ipt in enumerate(plot_ipt_global)
}

# --- iz -> color (for plots vs pt2) ---
all_iz_global = sorted({r["iz"] for r in records})
if len(all_iz_global) > 2:
    # again, interior z bins, like in the z-plots
    plot_iz_global = all_iz_global[1:-1]
else:
    plot_iz_global = all_iz_global

GLOBAL_COLOR_MAP_Z = {
    iz: base_colors[i % n_base_colors] for i, iz in enumerate(plot_iz_global)
}

# ----------------------------------------------------------------------
# First part: (p0/bin_width)/content, (mean/bin_width)/content vs z (ipt as color)
# ----------------------------------------------------------------------
def make_plot_vs_z(sel_records, tag, title_extra="", norm_by_ixq2=None):
    """
    Make one 4x4 plot:
      (p0/bin_width)/content, (mean/bin_width)/content vs z, ipt as color.

    norm_by_ixq2: dict mapping ixq2 -> content (normalization).
                  If None, no extra normalization is applied.
    """
    if not sel_records:
        print(f"No records for selection '{tag}' (vs z), skipping.")
        return

    # sort records for stable grouping
    sel_records = sorted(sel_records, key=lambda r: (r["ixq2"], r["ipt"], r["iz"]))

    # unique ixq2 values, take the first N_IX_SHOW
    ixq2_values = sorted({r["ixq2"] for r in sel_records})
    ixq2_values = ixq2_values[:N_IX_SHOW]

    # group records by ixq2
    by_ixq2 = defaultdict(list)
    for r in sel_records:
        if r["ixq2"] in ixq2_values:
            by_ixq2[r["ixq2"]].append(r)

    fig, axes = plt.subplots(
        NROWS, NCOLS, figsize=(16, 12), sharex=True, sharey=False
    )
    axes = axes.flatten()

    # Force x-axis tick labels to be shown on every subplot
    for ax in axes:
        ax.tick_params(axis="x", which="both", labelbottom=True)

    max_iz_index = len(Z_EDGES) - 1  # iz runs from 1..max_iz_index

    for idx, ixq2 in enumerate(ixq2_values):
        ax = axes[idx]
        subset = by_ixq2[ixq2]
        if not subset:
            ax.set_visible(False)
            continue

        # Determine normalization factor for this ixq2 bin
        scale = 1.0
        if norm_by_ixq2 is not None and ixq2 in norm_by_ixq2:
            scale = norm_by_ixq2[ixq2]

        # available ipt for this ixq2
        all_ipt_values_local = sorted({r["ipt"] for r in subset})
        # drop very first and very last ipt
        if len(all_ipt_values_local) > 2:
            ipt_values = all_ipt_values_local[1:-1]
        else:
            ipt_values = []  # nothing to plot if only 2 or fewer ipt bins

        for ipt in ipt_values:
            rows_ip = [r for r in subset if r["ipt"] == ipt]
            rows_ip.sort(key=lambda r: r["iz"])

            # z centers from Z_EDGES, drop first and last iz (interior bins)
            z_centers = []
            p0_vals = []
            mean_vals = []

            for r in rows_ip:
                iz = r["iz"]
                # keep only interior z bins: 1 < iz < max_iz_index
                if iz <= 1 or iz >= max_iz_index:
                    continue
                # iz corresponds to interval Z_EDGES[iz-1] .. Z_EDGES[iz]
                zc = 0.5 * (Z_EDGES[iz - 1] + Z_EDGES[iz])
                z_centers.append(zc)
                # divide by unfolding content (scale) if available
                p0_vals.append(r["p0_over_bw"] / scale  )
                mean_vals.append(r["mean_over_bw"] / scale )

            if not z_centers:
                continue

            # Use global ipt color map (consistent across all vs-z plots)
            col = GLOBAL_COLOR_MAP_IPT.get(ipt, base_colors[0])
            pt2_center = 0.5 * (PT2_EDGES[ipt] + PT2_EDGES[ipt - 1])

            # p0/bin_width/content: circles, solid line (these go in the legend)
            ax.plot(
                z_centers,
                p0_vals,
                marker="o",
                linestyle="-",
                color=col,
                label=rf"$p_T^2 = {pt2_center:.3f}$",
            )

            # mean/bin_width/content: squares, dashed line (same color), but hidden from legend
            ax.plot(
                z_centers,
                mean_vals,
                marker="s",
                linestyle="--",
                color=col,
                label=f"_mean_ipt_{ipt}",
            )

        ax.set_title(rf"$x-Q^2 = {ixq2}$")
        ax.set_xlabel("z (bin center)")
        ax.set_ylabel(r"$p_0,\ \langle y\rangle\ /\ (\mathrm{bin\_width}\cdot\mathrm{content})$")
        ax.set_yscale("log")

        # only show legend in the first subplot to avoid clutter
        if idx == 0:
            ax.legend(fontsize=8, ncol=2)

    # hide any unused subplots if ixq2_values < 16
    for j in range(len(ixq2_values), len(axes)):
        axes[j].set_visible(False)

    fig.tight_layout()
    fig.suptitle(
        r"z-dependence of $p_0/(\mathrm{bin\_width}\cdot\mathrm{content})$ and "
        r"$\langle y\rangle/(\mathrm{bin\_width}\cdot\mathrm{content})$ "
        r"for first 16 $x-Q^2$ bins"
        + title_extra,
        y=1.02,
    )

    out_png = f"p0_mean_over_binwidth_over_content_vs_z_4x4_{tag}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Saved figure to {out_png}")
    plt.close(fig)

# ----------------------------------------------------------------------
# Second part: (p0/bin_width)/content, (mean/bin_width)/content vs pT^2 (z as color)
# ----------------------------------------------------------------------
def make_plot_vs_pt2(sel_records, tag, title_extra="", norm_by_ixq2=None):
    """
    Make one 4x4 plot:
      (p0/bin_width)/content, (mean/bin_width)/content vs pT^2, z as color.

    norm_by_ixq2: dict mapping ixq2 -> content (normalization).
                  If None, no extra normalization is applied.
    """
    if not sel_records:
        print(f"No records for selection '{tag}' (vs pT^2), skipping.")
        return

    # sort records for stable grouping
    sel_records = sorted(sel_records, key=lambda r: (r["ixq2"], r["iz"], r["ipt"]))

    # unique ixq2 values, take the first N_IX_SHOW
    ixq2_values = sorted({r["ixq2"] for r in sel_records})
    ixq2_values = ixq2_values[:N_IX_SHOW]

    # group records by ixq2
    by_ixq2 = defaultdict(list)
    for r in sel_records:
        if r["ixq2"] in ixq2_values:
            by_ixq2[r["ixq2"]].append(r)

    fig, axes = plt.subplots(
        NROWS, NCOLS, figsize=(16, 12), sharex=True, sharey=False
    )
    axes = axes.flatten()

    # Force x-axis tick labels to be shown on every subplot
    for ax in axes:
        ax.tick_params(axis="x", which="both", labelbottom=True)

    max_ipt_index = len(PT2_EDGES) - 1  # ipt runs from 1..max_ipt_index

    for idx, ixq2 in enumerate(ixq2_values):
        ax = axes[idx]
        subset = by_ixq2[ixq2]
        if not subset:
            ax.set_visible(False)
            continue

        # Determine normalization factor for this ixq2 bin
        scale = 1.0
        if norm_by_ixq2 is not None and ixq2 in norm_by_ixq2:
            scale = norm_by_ixq2[ixq2]

        # available z for this ixq2
        all_iz_values_local = sorted({r["iz"] for r in subset})
        # use interior z bins (to match above)
        if len(all_iz_values_local) > 2:
            iz_values = all_iz_values_local[1:-1]
        else:
            iz_values = all_iz_values_local

        for iz in iz_values:
            rows_z = [r for r in subset if r["iz"] == iz]
            rows_z.sort(key=lambda r: r["ipt"])

            pt2_centers = []
            p0_vals = []
            mean_vals = []

            for r in rows_z:
                ipt = r["ipt"]
                # keep pt bins including first, excluding last (overflow):
                # 1 <= ipt < max_ipt_index
                if ipt < 1 or ipt >= max_ipt_index:
                    continue

                pt2c = 0.5 * (PT2_EDGES[ipt - 1] + PT2_EDGES[ipt])
                pt2_centers.append(pt2c)
                # divide by unfolding content (scale) if available
                p0_vals.append(r["p0_over_bw"] / scale  )
                mean_vals.append(r["mean_over_bw"] / scale )

            if not pt2_centers:
                continue

            # Use global z color map (consistent across all vs-pt2 plots)
            col = GLOBAL_COLOR_MAP_Z.get(iz, base_colors[0])
            z_label = f"{Z_EDGES[iz - 1]:.2f}–{Z_EDGES[iz]:.2f}"

            # p0/bin_width/content: circles, solid line (these go in the legend)
            ax.plot(
                pt2_centers,
                p0_vals,
                marker="o",
                linestyle="-",
                color=col,
                label=rf"$z = {z_label}$",
            )

            # mean/bin_width/content: squares, dashed line (same color), but hidden from legend
            ax.plot(
                pt2_centers,
                mean_vals,
                marker="s",
                linestyle="--",
                color=col,
                label=f"_mean_iz_{iz}",
            )

        ax.set_title(rf"$x-Q^2 = {ixq2}$")
        ax.set_xlabel(r"$p_T^2$ (bin center)")
        ax.set_ylabel(r"$p_0,\ \langle y\rangle\ /\ (\mathrm{bin\_width}\cdot\mathrm{content})$")
        ax.set_yscale("log")

        # only show legend in the first subplot to avoid clutter
        if idx == 0:
            ax.legend(fontsize=8, ncol=2)

    # hide any unused subplots if ixq2_values < 16
    for j in range(len(ixq2_values), len(axes)):
        axes[j].set_visible(False)

    fig.tight_layout()
    fig.suptitle(
        r"$p_T^2$-dependence of $p_0/(\mathrm{bin\_width}\cdot\mathrm{content})$ and "
        r"$\langle y\rangle/(\mathrm{bin\_width}\cdot\mathrm{content})$ "
        r"for first 16 $x-Q^2$ bins"
        + title_extra,
        y=1.02,
    )

    out_png = f"p0_mean_over_binwidth_over_content_vs_pt2_4x4_{tag}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Saved figure to {out_png}")
    plt.close(fig)

# ----------------------------------------------------------------------
# Build selections and make six plots (3 vs z, 3 vs pT^2)
# ----------------------------------------------------------------------
# 1) Ignore flag: use all records
make_plot_vs_z(records, tag="all", title_extra=r" (all fits)", norm_by_ixq2=NORM_BY_IXQ2)
make_plot_vs_pt2(records, tag="all", title_extra=r" (all fits)", norm_by_ixq2=NORM_BY_IXQ2)

# 2) Only no_limit_hit == 1
records_flag1 = [r for r in records if r["no_limit_hit"] == 1]
make_plot_vs_z(records_flag1, tag="nolimit1", title_extra=r" ($no\_limit\_hit = 1$)", norm_by_ixq2=NORM_BY_IXQ2)
make_plot_vs_pt2(records_flag1, tag="nolimit1", title_extra=r" ($no\_limit\_hit = 1$)", norm_by_ixq2=NORM_BY_IXQ2)

# 3) Only no_limit_hit == 0
records_flag0 = [r for r in records if r["no_limit_hit"] == 0]
make_plot_vs_z(records_flag0, tag="nolimit0", title_extra=r" ($no\_limit\_hit = 0$)", norm_by_ixq2=NORM_BY_IXQ2)
make_plot_vs_pt2(records_flag0, tag="nolimit0", title_extra=r" ($no\_limit\_hit = 0$)", norm_by_ixq2=NORM_BY_IXQ2)
