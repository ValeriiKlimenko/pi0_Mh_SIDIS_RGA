#!/usr/bin/env python3
import os
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- defaults ----
DEFAULT_PHI_CSV  = "tables_h_meas_data.csv"  # has phi rows
DEFAULT_SPEC_CSV = "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v2/initial_data_spectrum.csv"  # no phi (initial spectrum)
DEFAULT_NOPHI_CSV = "tables_nophi_h_meas_data.csv"  # NEW: already no-phi, no summing needed (override with --nophi_csv)


def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def ensure_int_cols(df: pd.DataFrame, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=list(cols)).copy()
    for c in cols:
        df[c] = df[c].astype(int)
    return df


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Compare 3 sources per (xQ2,z,pt2): "
            "(1) phi-binned CSV summed over phi, "
            "(2) initial spectrum CSV (no phi), "
            "(3) new no-phi CSV (already integrated, no summing). "
            "Plots: one figure per xQ2, subplots per z; top=content vs pt2, bottom=ratios."
        )
    )
    ap.add_argument("--phi_csv", default=DEFAULT_PHI_CSV,
                    help=f"CSV with phi bins (default: {DEFAULT_PHI_CSV})")
    ap.add_argument("--spec_csv", default=DEFAULT_SPEC_CSV,
                    help=f"CSV without phi (default: {DEFAULT_SPEC_CSV})")
    ap.add_argument("--nophi_csv", default=DEFAULT_NOPHI_CSV,
                    help=f"NEW: CSV without phi, already integrated (default: {DEFAULT_NOPHI_CSV})")
    ap.add_argument("--outdir", default="compare_phiSum_vs_spec_vs_nophi",
                    help="Output directory for plots and CSVs")
    ap.add_argument("--phi_content_col", default="content")
    ap.add_argument("--spec_content_col", default="content")
    ap.add_argument("--nophi_content_col", default="content")
    ap.add_argument("--phi_xq2_shift", type=int, default=-1,
                    help="Shift xq2_bin in phi CSV by this integer (default: -1)")
    ap.add_argument("--skip_allzero_rows", action="store_true",
                    help="Drop rows where all contents are zero (after merge)")
    ap.add_argument("--max_xq2", type=int, default=None)
    ap.add_argument("--ncols", type=int, default=4,
                    help="Number of subplot columns for z panels")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------ load ------------------
    df_phi   = load_csv(args.phi_csv)
    df_spec  = load_csv(args.spec_csv)
    df_nophi = load_csv(args.nophi_csv)

    # ------------------ required columns ------------------
    req_phi   = {"xq2_bin", "z_bin", "pt2_bin", "phi_bin", args.phi_content_col}
    req_spec  = {"xq2_bin", "z_bin", "pt2_bin", args.spec_content_col}
    # NEW no-phi CSV: we accept either 0-based (z_bin,pt2_bin) or 1-based (z_bin,pt2_bin) but no phi column
    req_nophi = {"xq2_bin", "z_bin", "pt2_bin", args.nophi_content_col}

    miss_phi   = req_phi - set(df_phi.columns)
    miss_spec  = req_spec - set(df_spec.columns)
    miss_nophi = req_nophi - set(df_nophi.columns)
    if miss_phi:
        raise ValueError(f"phi_csv missing columns: {miss_phi}\nFile: {args.phi_csv}")
    if miss_spec:
        raise ValueError(f"spec_csv missing columns: {miss_spec}\nFile: {args.spec_csv}")
    if miss_nophi:
        raise ValueError(f"nophi_csv missing columns: {miss_nophi}\nFile: {args.nophi_csv}")

    # ------------------ numeric casting ------------------
    for c in ["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]:
        df_phi[c] = pd.to_numeric(df_phi[c], errors="coerce")
    df_phi[args.phi_content_col] = pd.to_numeric(df_phi[args.phi_content_col], errors="coerce")
    df_phi = df_phi.dropna(subset=["xq2_bin", "z_bin", "pt2_bin", "phi_bin", args.phi_content_col]).copy()
    df_phi[["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]] = df_phi[["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]].astype(int)

    for c in ["xq2_bin", "z_bin", "pt2_bin"]:
        df_spec[c] = pd.to_numeric(df_spec[c], errors="coerce")
    df_spec[args.spec_content_col] = pd.to_numeric(df_spec[args.spec_content_col], errors="coerce")
    df_spec = df_spec.dropna(subset=["xq2_bin", "z_bin", "pt2_bin", args.spec_content_col]).copy()
    df_spec[["xq2_bin", "z_bin", "pt2_bin"]] = df_spec[["xq2_bin", "z_bin", "pt2_bin"]].astype(int)

    for c in ["xq2_bin", "z_bin", "pt2_bin"]:
        df_nophi[c] = pd.to_numeric(df_nophi[c], errors="coerce")
    df_nophi[args.nophi_content_col] = pd.to_numeric(df_nophi[args.nophi_content_col], errors="coerce")
    df_nophi = df_nophi.dropna(subset=["xq2_bin", "z_bin", "pt2_bin", args.nophi_content_col]).copy()
    df_nophi[["xq2_bin", "z_bin", "pt2_bin"]] = df_nophi[["xq2_bin", "z_bin", "pt2_bin"]].astype(int)

    # ------------------ shift xQ2 for phi CSV ------------------
    if args.phi_xq2_shift != 0:
        df_phi["xq2_bin"] = df_phi["xq2_bin"] + int(args.phi_xq2_shift)
        df_phi = df_phi[df_phi["xq2_bin"] >= 0].copy()

    # =============================================================================
    # 1) SUM over phi bins for each (xQ2, z, pt2)
    # =============================================================================
    phi_sum = (
        df_phi.groupby(["xq2_bin", "z_bin", "pt2_bin"], as_index=False)
              .agg(phi_content_sum=(args.phi_content_col, "sum"),
                   phi_bin_count=("phi_bin", "size"))
    )

    # Align bin numbering:
    #   phi CSV: z,pt2 start at 1  -> make 0-based for merging
    phi_sum["z0"]   = phi_sum["z_bin"] - 1
    phi_sum["pt20"] = phi_sum["pt2_bin"] - 1
    phi_sum = phi_sum[(phi_sum["z0"] >= 0) & (phi_sum["pt20"] >= 0)].copy()

    out_phi_sum_csv = os.path.join(args.outdir, "phi_summed_by_xq2_z_pt2.csv")
    phi_sum.to_csv(out_phi_sum_csv, index=False)

    # =============================================================================
    # 2) Prepare spec (assumed 0-based) and NEW no-phi CSV (could be 0- or 1-based)
    # =============================================================================
    spec_small = df_spec[["xq2_bin", "z_bin", "pt2_bin", args.spec_content_col]].copy()
    spec_small = spec_small.rename(columns={
        "z_bin": "z0",
        "pt2_bin": "pt20",
        args.spec_content_col: "spec_content"
    })

    # Detect if nophi_csv is 1-based or 0-based by looking at min values
    # If it starts at 1 (and spec starts at 0), shift to 0-based.
    nophi_small = df_nophi[["xq2_bin", "z_bin", "pt2_bin", args.nophi_content_col]].copy()
    zmin = int(nophi_small["z_bin"].min())
    ptmin = int(nophi_small["pt2_bin"].min())
    needs_shift = (zmin >= 1 and ptmin >= 1)

    nophi_small["z0"] = nophi_small["z_bin"] - (1 if needs_shift else 0)
    nophi_small["pt20"] = nophi_small["pt2_bin"] - (1 if needs_shift else 0)
    nophi_small = nophi_small.rename(columns={args.nophi_content_col: "nophi_content"})
    nophi_small = nophi_small.drop(columns=["z_bin", "pt2_bin"])

    # =============================================================================
    # 3) Merge all three on (xq2_bin, z0, pt20)
    # =============================================================================
    phi_small = phi_sum[["xq2_bin", "z0", "pt20", "phi_content_sum", "phi_bin_count"]].copy()

    merged = (
        phi_small
        .merge(spec_small, on=["xq2_bin", "z0", "pt20"], how="outer", indicator=False)
        .merge(nophi_small, on=["xq2_bin", "z0", "pt20"], how="outer", indicator=False)
    )

    # Add ratios (use np.nan where denom is 0)
    merged["ratio(phi/spec)"]   = merged["phi_content_sum"] / merged["spec_content"].replace(0, np.nan)
    merged["ratio(phi/nophi)"]  = merged["phi_content_sum"] / merged["nophi_content"].replace(0, np.nan)
    merged["ratio(nophi/spec)"] = merged["nophi_content"] / merged["spec_content"].replace(0, np.nan)

    if args.skip_allzero_rows:
        z = (merged["phi_content_sum"].fillna(0) == 0) & (merged["spec_content"].fillna(0) == 0) & (merged["nophi_content"].fillna(0) == 0)
        merged = merged[~z].copy()

    out_merged_csv = os.path.join(args.outdir, "phiSum_vs_spec_vs_nophi_merged.csv")
    merged.to_csv(out_merged_csv, index=False)

    # =============================================================================
    # 4) Plots: one figure per xQ2, subplots per z
    #    Each z panel has 2 rows: top content, bottom ratio(s).
    #    Ratio y-range: at least [0.95,1.05] (can expand beyond).
    # =============================================================================
    xq2_values = sorted(set(merged["xq2_bin"].dropna().astype(int).tolist()))
    if args.max_xq2 is not None:
        xq2_values = [x for x in xq2_values if x <= args.max_xq2]

    ncols = max(1, int(args.ncols))

    def ratio_ylim(rvals):
        r = np.asarray(rvals, dtype=float)
        r = r[np.isfinite(r)]
        if r.size == 0:
            return (0.95, 1.05)
        lo = float(np.nanmin(r))
        hi = float(np.nanmax(r))
        lo = min(lo, 0.95)
        hi = max(hi, 1.05)
        # tiny pad
        pad = 0.02 * (hi - lo) if hi > lo else 0.02
        return (lo - pad, hi + pad)

    for xq2 in xq2_values:
        df_x = merged[merged["xq2_bin"] == xq2].copy()
        df_x = df_x.dropna(subset=["z0", "pt20"], how="any").copy()
        if df_x.empty:
            continue

        df_x["z0"] = df_x["z0"].astype(int)
        df_x["pt20"] = df_x["pt20"].astype(int)

        z_values = sorted(df_x["z0"].unique().tolist())
        if not z_values:
            continue

        nrows_grid = int(math.ceil(len(z_values) / ncols))

        # Figure geometry: each z panel is 2 stacked axes -> height grows
        fig_w = max(12, 3.2 * ncols)
        fig_h = max(6,  3.8 * nrows_grid)

        # Create a grid of panels; inside each panel we make 2 axes via gridspec
        fig = plt.figure(figsize=(fig_w, fig_h))
        outer = fig.add_gridspec(nrows_grid, ncols, wspace=0.25, hspace=0.35)

        fig.suptitle(f"xQ2 bin {xq2}: pT² dependence (subplots per z) — phi-sum vs spec vs no-phi", y=0.995)

        # legend handles (global)
        h_phi = h_spec = h_nophi = None
        h_r_phi_spec = h_r_phi_nophi = h_r_nophi_spec = None

        for idx, z in enumerate(z_values):
            r = idx // ncols
            c = idx % ncols
            inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)

            ax_top = fig.add_subplot(inner[0, 0])
            ax_bot = fig.add_subplot(inner[1, 0], sharex=ax_top)

            sub = df_x[df_x["z0"] == z].copy().sort_values("pt20")
            x = sub["pt20"].to_numpy()

            y_phi   = sub["phi_content_sum"].to_numpy()
            y_spec  = sub["spec_content"].to_numpy()
            y_nophi = sub["nophi_content"].to_numpy()

            # ----- top: contents -----
            m = np.isfinite(y_phi)
            if m.any():
                (hp,) = ax_top.plot(x[m], y_phi[m], marker="o", linestyle="-")
                if h_phi is None: h_phi = hp

            m = np.isfinite(y_spec)
            if m.any():
                (hs,) = ax_top.plot(x[m], y_spec[m], marker="x", linestyle="--")
                if h_spec is None: h_spec = hs

            m = np.isfinite(y_nophi)
            if m.any():
                (hn,) = ax_top.plot(x[m], y_nophi[m], marker="s", linestyle=":")
                if h_nophi is None: h_nophi = hn

            ax_top.set_title(f"z bin {z}", fontsize=10)
            ax_top.grid(True, which="both", linewidth=0.5)

            if (idx % ncols) == 0:
                ax_top.set_ylabel("content")

            # ----- bottom: ratios -----
            r_phi_spec   = sub["ratio(phi/spec)"].to_numpy()
            r_phi_nophi  = sub["ratio(phi/nophi)"].to_numpy()
            r_nophi_spec = sub["ratio(nophi/spec)"].to_numpy()

            m = np.isfinite(r_phi_spec)
            if m.any():
                (hr1,) = ax_bot.plot(x[m], r_phi_spec[m], marker="o", linestyle="-")
                if h_r_phi_spec is None: h_r_phi_spec = hr1

            m = np.isfinite(r_phi_nophi)
            if m.any():
                (hr2,) = ax_bot.plot(x[m], r_phi_nophi[m], marker="s", linestyle=":")
                if h_r_phi_nophi is None: h_r_phi_nophi = hr2

            m = np.isfinite(r_nophi_spec)
            if m.any():
                (hr3,) = ax_bot.plot(x[m], r_nophi_spec[m], marker="x", linestyle="--")
                if h_r_nophi_spec is None: h_r_nophi_spec = hr3

            ax_bot.axhline(1.0, linewidth=1.0)
            ax_bot.grid(True, which="both", linewidth=0.5)
            ax_bot.set_ylabel("ratio", fontsize=9)

            # enforce ratio y-range at least [0.95, 1.05]
            ax_bot.set_ylim(*ratio_ylim(np.r_[r_phi_spec, r_phi_nophi, r_nophi_spec]))

            # only show x-labels on bottom row of grid
            if r == (nrows_grid - 1):
                ax_bot.set_xlabel("pT² bin (0-based)")
            else:
                plt.setp(ax_bot.get_xticklabels(), visible=False)

            # hide top x tick labels
            plt.setp(ax_top.get_xticklabels(), visible=False)

        # Turn off unused panels
        total_slots = nrows_grid * ncols
        for idx in range(len(z_values), total_slots):
            r = idx // ncols
            c = idx % ncols
            ax_dummy = fig.add_subplot(outer[r, c])
            ax_dummy.axis("off")

        # Global legend
        handles = []
        labels = []

        if h_phi is not None:
            handles.append(h_phi); labels.append("phi-sum")
        if h_spec is not None:
            handles.append(h_spec); labels.append("spec")
        if h_nophi is not None:
            handles.append(h_nophi); labels.append("no-phi CSV")

        if h_r_phi_spec is not None:
            handles.append(h_r_phi_spec); labels.append("phi/spec")
        if h_r_phi_nophi is not None:
            handles.append(h_r_phi_nophi); labels.append("phi/noPhi")
        if h_r_nophi_spec is not None:
            handles.append(h_r_nophi_spec); labels.append("noPhi/spec")

        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.985), fontsize=9)

        png = os.path.join(args.outdir, f"xq2_{xq2:02d}_pt2_by_z_subplots_with_ratios.png")
        pdf = os.path.join(args.outdir, f"xq2_{xq2:02d}_pt2_by_z_subplots_with_ratios.pdf")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)

    # =============================================================================
    # Summary
    # =============================================================================
    print("Wrote:")
    print("  ", out_phi_sum_csv)
    print("  ", out_merged_csv)
    print("\nNotes:")
    print(f"  nophi_csv interpreted as {'1-based (shifted to 0-based)' if needs_shift else '0-based (no shift)'} for z_bin/pt2_bin")
    print("\nPlots saved in:", args.outdir)


if __name__ == "__main__":
    main()
