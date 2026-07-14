#!/usr/bin/env python3
import csv
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------
# User options
# ----------------------------------------------------------------------
CSV_FILE = "phi_fit_twofiles.csv"  # <-- your CSV from the fitter

# z bin edges (must match your analysis!)
Z_EDGES = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
PT2_EDGES = [
    0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50,
    0.65, 0.80, 1.00, 1.50
]

# how many xQ² bins to show, and layout
N_IX_SHOW = 16
NROWS, NCOLS = 4, 4  # 4x4 grid

# ----------------------------------------------------------------------
# Helpers
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
# Read CSV
# Expecting columns from your fitter:
#   xq2_center,xq2_bin,iy,z_bin,pt2_bin,z_label,pt2_label,
#   rootfile,source,npoints,p0,p1,p2,chi2,ndf,status
# ----------------------------------------------------------------------
records = []
with open(CSV_FILE, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        xq2   = int(row["xq2_bin"])
        iz    = int(row["z_bin"])
        ipt   = int(row["pt2_bin"])
        p0    = float(row["p0"])
        src   = row.get("source", "file")

        bw = compute_bin_width(iz, ipt)
        if bw > 0.0:
            p0_over_bw = p0 / bw
        else:
            p0_over_bw = np.nan

        records.append(
            {
                "xq2": xq2,
                "iz": iz,
                "ipt": ipt,
                "p0": p0,
                "p0_over_bw": p0_over_bw,
                "source": src,
            }
        )

if not records:
    raise RuntimeError("No rows read from CSV; check CSV_FILE and header names.")

# unique sources (= your two files)
SOURCES = sorted({r["source"] for r in records})
if len(SOURCES) < 2:
    print("WARNING: found fewer than 2 distinct 'source' labels in CSV.")
print("Sources in CSV:", SOURCES)

# ----------------------------------------------------------------------
# Global color maps (same idea as your old script)
# ----------------------------------------------------------------------
base_colors = plt.cm.tab10.colors
n_base_colors = len(base_colors)

# ipt -> color (for plots vs z)
all_ipt_global = sorted({r["ipt"] for r in records})
if len(all_ipt_global) > 2:
    plot_ipt_global = all_ipt_global[1:-1]  # interior pT² bins
else:
    plot_ipt_global = all_ipt_global

GLOBAL_COLOR_MAP_IPT = {
    ipt: base_colors[i % n_base_colors] for i, ipt in enumerate(plot_ipt_global)
}

# iz -> color (for plots vs pt²)
all_iz_global = sorted({r["iz"] for r in records})
if len(all_iz_global) > 2:
    plot_iz_global = all_iz_global[1:-1]  # interior z bins
else:
    plot_iz_global = all_iz_global

GLOBAL_COLOR_MAP_Z = {
    iz: base_colors[i % n_base_colors] for i, iz in enumerate(plot_iz_global)
}

# source -> line style (solid vs dashed)
LINESTYLES = ["-", "--", ":", "-."]  # more styles in case of >2 sources
SOURCE_LINESTYLE = {
    src: LINESTYLES[i % len(LINESTYLES)] for i, src in enumerate(SOURCES)
}

# ----------------------------------------------------------------------
# Plots vs z: p0/bin_width with pT² as color, files as line style
# ----------------------------------------------------------------------
def make_plot_vs_z(sel_records, tag, title_extra=""):
    if not sel_records:
        print(f"No records for selection '{tag}' (vs z), skipping.")
        return

    sel_records = sorted(sel_records, key=lambda r: (r["xq2"], r["ipt"], r["iz"]))

    # xQ² bins to show
    xq2_values = sorted({r["xq2"] for r in sel_records})[:N_IX_SHOW]

    by_xq2 = defaultdict(list)
    for r in sel_records:
        if r["xq2"] in xq2_values:
            by_xq2[r["xq2"]].append(r)

    fig, axes = plt.subplots(
        NROWS, NCOLS, figsize=(16, 12), sharex=True, sharey=False
    )
    axes = axes.flatten()

    for ax in axes:
        ax.tick_params(axis="x", which="both", labelbottom=True)

    max_iz_index = len(Z_EDGES) - 1  # iz runs from 1..max_iz_index

    for idx, xq2 in enumerate(xq2_values):
        ax = axes[idx]
        subset = by_xq2[xq2]
        if not subset:
            ax.set_visible(False)
            continue

        # local pT² bins available in this xQ²
        all_ipt_vals_local = sorted({r["ipt"] for r in subset})
        if len(all_ipt_vals_local) > 2:
            ipt_values = all_ipt_vals_local[1:-1]  # drop first and last
        else:
            ipt_values = all_ipt_vals_local

        for ipt in ipt_values:
            # For each file/source, plot with different linestyle
            for src in SOURCES:
                rows = [r for r in subset if r["ipt"] == ipt and r["source"] == src]
                rows.sort(key=lambda r: r["iz"])

                z_centers = []
                p0_vals = []

                for r in rows:
                    iz = r["iz"]
                    # keep only interior z bins: 1 < iz < max_iz_index
                    if iz <= 1 or iz >= max_iz_index:
                        continue
                    zc = 0.5 * (Z_EDGES[iz - 1] + Z_EDGES[iz])
                    z_centers.append(zc)
                    p0_vals.append(r["p0_over_bw"])

                if not z_centers:
                    continue

                col = GLOBAL_COLOR_MAP_IPT.get(ipt, base_colors[0])
                pt2_center = 0.5 * (PT2_EDGES[ipt - 1] + PT2_EDGES[ipt])
                ls = SOURCE_LINESTYLE.get(src, "-")

                # Only label ipt once (for the first source) to keep legend small
                if src == SOURCES[0]:
                    label = rf"$p_T^2 = {pt2_center:.3f}$"
                else:
                    # underscore -> hidden from legend
                    label = f"_{src}_ipt_{ipt}"

                ax.plot(
                    z_centers,
                    p0_vals,
                    marker="o",
                    linestyle=ls,
                    color=col,
                    label=label,
                )

        ax.set_title(rf"$x-Q^2 = {xq2}$")
        ax.set_xlabel("z (bin center)")
        ax.set_ylabel(r"$p_0/\mathrm{bin\_width}$")
        ax.set_yscale("log")

        if idx == 0:
            # legend for pT² colors
            leg1 = ax.legend(fontsize=8, ncol=2, loc="upper right")

            # legend for file line styles
            style_handles = [
                Line2D(
                    [0],
                    [0],
                    color="k",
                    linestyle=SOURCE_LINESTYLE.get(src, "-"),
                    label=src,
                )
                for src in SOURCES
            ]
            leg2 = ax.legend(
                handles=style_handles,
                fontsize=8,
                loc="lower left",
                title="source",
            )
            ax.add_artist(leg1)  # keep the first legend too

    # hide unused pads
    for j in range(len(xq2_values), len(axes)):
        axes[j].set_visible(False)

    fig.tight_layout()
    fig.suptitle(
        r"$z$-dependence of $p_0/\mathrm{bin\_width}$ "
        r"(colors = $p_T^2$ bins, solid vs dashed = files)"
        + title_extra,
        y=1.02,
    )

    out_png = f"p0_over_binwidth_vs_z_4x4_{tag}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Saved figure to {out_png}")
    plt.close(fig)


# ----------------------------------------------------------------------
# Plots vs pT²: p0/bin_width with z as color, files as line style
# ----------------------------------------------------------------------
def make_plot_vs_pt2(sel_records, tag, title_extra=""):
    if not sel_records:
        print(f"No records for selection '{tag}' (vs pT^2), skipping.")
        return

    sel_records = sorted(sel_records, key=lambda r: (r["xq2"], r["iz"], r["ipt"]))
    xq2_values = sorted({r["xq2"] for r in sel_records})[:N_IX_SHOW]

    by_xq2 = defaultdict(list)
    for r in sel_records:
        if r["xq2"] in xq2_values:
            by_xq2[r["xq2"]].append(r)

    fig, axes = plt.subplots(
        NROWS, NCOLS, figsize=(16, 12), sharex=True, sharey=False
    )
    axes = axes.flatten()

    for ax in axes:
        ax.tick_params(axis="x", which="both", labelbottom=True)

    max_ipt_index = len(PT2_EDGES) - 1  # ipt 1..max_ipt_index

    for idx, xq2 in enumerate(xq2_values):
        ax = axes[idx]
        subset = by_xq2[xq2]
        if not subset:
            ax.set_visible(False)
            continue

        all_iz_vals_local = sorted({r["iz"] for r in subset})
        # interior z bins (to match vs-z plots)
        if len(all_iz_vals_local) > 2:
            iz_values = all_iz_vals_local[1:-1]
        else:
            iz_values = all_iz_vals_local

        for iz in iz_values:
            for src in SOURCES:
                rows = [
                    r
                    for r in subset
                    if r["iz"] == iz and r["source"] == src
                ]
                rows.sort(key=lambda r: r["ipt"])

                pt2_centers = []
                p0_vals = []

                for r in rows:
                    ipt = r["ipt"]
                    # use 1 <= ipt < max_ipt_index (exclude overflow)
                    if ipt < 1 or ipt >= max_ipt_index:
                        continue
                    pt2c = 0.5 * (PT2_EDGES[ipt - 1] + PT2_EDGES[ipt])
                    pt2_centers.append(pt2c)
                    p0_vals.append(r["p0_over_bw"])

                if not pt2_centers:
                    continue

                col = GLOBAL_COLOR_MAP_Z.get(iz, base_colors[0])
                z_label = f"{Z_EDGES[iz - 1]:.2f}–{Z_EDGES[iz]:.2f}"
                ls = SOURCE_LINESTYLE.get(src, "-")

                # Only label z once (for first source) in the color legend
                if src == SOURCES[0]:
                    label = rf"$z = {z_label}$"
                else:
                    label = f"_{src}_iz_{iz}"

                ax.plot(
                    pt2_centers,
                    p0_vals,
                    marker="o",
                    linestyle=ls,
                    color=col,
                    label=label,
                )

        ax.set_title(rf"$x-Q^2 = {xq2}$")
        ax.set_xlabel(r"$p_T^2$ (bin center)")
        ax.set_ylabel(r"$p_0/\mathrm{bin\_width}$")
        ax.set_yscale("log")

        if idx == 0:
            # legend for z colors
            leg1 = ax.legend(fontsize=8, ncol=2, loc="upper right")

            # legend for file line styles
            style_handles = [
                Line2D(
                    [0],
                    [0],
                    color="k",
                    linestyle=SOURCE_LINESTYLE.get(src, "-"),
                    label=src,
                )
                for src in SOURCES
            ]
            leg2 = ax.legend(
                handles=style_handles,
                fontsize=8,
                loc="lower left",
                title="source",
            )
            ax.add_artist(leg1)

    for j in range(len(xq2_values), len(axes)):
        axes[j].set_visible(False)

    fig.tight_layout()
    fig.suptitle(
        r"$p_T^2$-dependence of $p_0/\mathrm{bin\_width}$ "
        r"(colors = $z$ bins, solid vs dashed = files)"
        + title_extra,
        y=1.02,
    )

    out_png = f"p0_over_binwidth_vs_pt2_4x4_{tag}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Saved figure to {out_png}")
    plt.close(fig)


# ----------------------------------------------------------------------
# Make the two overlays (vs z and vs pT²)
# ----------------------------------------------------------------------
make_plot_vs_z(records, tag="bothfiles", title_extra="")
make_plot_vs_pt2(records, tag="bothfiles", title_extra="")
