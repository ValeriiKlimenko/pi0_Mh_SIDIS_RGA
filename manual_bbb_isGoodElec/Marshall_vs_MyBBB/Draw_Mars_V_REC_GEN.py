#!/usr/bin/env python3
import os
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULTS = {
    "rec": {
        "phi_csv": "tables_h_meas_mc.csv",
        "spec_csv": "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v2/initial_rec_spectrum.csv",
        "outdir": "compare_REC_phiSum_vs_initialREC",
    },
    "gen": {
        "phi_csv": "tables_h_true_mc.csv",
        "spec_csv": "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v2/initial_gen_spectrum.csv",
        "outdir": "compare_GEN_phiSum_vs_initialGEN",
    },
    # keep data here too if you ever want it
    "data": {
        "phi_csv": "tables_h_meas_data.csv",
        "spec_csv": "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v2/initial_data_spectrum.csv",
        "outdir": "compare_DATA_phiSum_vs_initialDATA",
    },
}


def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def compare_pair(phi_csv: str,
                 spec_csv: str,
                 outdir: str,
                 phi_content_col: str = "content",
                 spec_content_col: str = "content",
                 phi_error_col: str = "error",
                 spec_error_col: str = "error",
                 phi_xq2_shift: int = -1,
                 skip_allzero_rows: bool = False,
                 max_xq2: int | None = None,
                 ncols: int = 4):
    os.makedirs(outdir, exist_ok=True)

    # ------------------ load ------------------
    df_phi = load_csv(phi_csv)
    df_spec = load_csv(spec_csv)

    # ------------------ required columns ------------------
    req_phi = {"xq2_bin", "z_bin", "pt2_bin", "phi_bin", phi_content_col}
    req_spec = {"xq2_bin", "z_bin", "pt2_bin", spec_content_col}
    miss_phi = req_phi - set(df_phi.columns)
    miss_spec = req_spec - set(df_spec.columns)
    if miss_phi:
        raise ValueError(f"phi_csv missing columns: {miss_phi}\nFile: {phi_csv}")
    if miss_spec:
        raise ValueError(f"spec_csv missing columns: {miss_spec}\nFile: {spec_csv}")

    # ------------------ numeric casting ------------------
    for c in ["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]:
        df_phi[c] = pd.to_numeric(df_phi[c], errors="coerce")
    df_phi[phi_content_col] = pd.to_numeric(df_phi[phi_content_col], errors="coerce")
    if phi_error_col in df_phi.columns:
        df_phi[phi_error_col] = pd.to_numeric(df_phi[phi_error_col], errors="coerce")

    for c in ["xq2_bin", "z_bin", "pt2_bin"]:
        df_spec[c] = pd.to_numeric(df_spec[c], errors="coerce")
    df_spec[spec_content_col] = pd.to_numeric(df_spec[spec_content_col], errors="coerce")
    if spec_error_col in df_spec.columns:
        df_spec[spec_error_col] = pd.to_numeric(df_spec[spec_error_col], errors="coerce")

    drop_phi_cols = ["xq2_bin", "z_bin", "pt2_bin", "phi_bin", phi_content_col]
    drop_spec_cols = ["xq2_bin", "z_bin", "pt2_bin", spec_content_col]
    df_phi = df_phi.dropna(subset=drop_phi_cols).copy()
    df_spec = df_spec.dropna(subset=drop_spec_cols).copy()

    # cast bin columns to int
    for c in ["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]:
        df_phi[c] = df_phi[c].astype(int)
    for c in ["xq2_bin", "z_bin", "pt2_bin"]:
        df_spec[c] = df_spec[c].astype(int)

    # ------------------ shift xQ2 for phi CSV ------------------
    if phi_xq2_shift != 0:
        df_phi["xq2_bin"] = df_phi["xq2_bin"] + int(phi_xq2_shift)
        df_phi = df_phi[df_phi["xq2_bin"] >= 0].copy()

    # =============================================================================
    # 1) SUM over phi bins for each (xQ2, z, pt2)
    # =============================================================================
    # Optional: propagate phi errors in quadrature (useful sanity check)
    if phi_error_col in df_phi.columns:
        df_phi["_err2"] = df_phi[phi_error_col] ** 2
        phi_sum = (
            df_phi.groupby(["xq2_bin", "z_bin", "pt2_bin"], as_index=False)
                  .agg(phi_content_sum=(phi_content_col, "sum"),
                       phi_err2_sum=("_err2", "sum"),
                       phi_bin_count=("phi_bin", "size"))
        )
        phi_sum["phi_error_sum"] = np.sqrt(phi_sum["phi_err2_sum"])
        phi_sum = phi_sum.drop(columns=["phi_err2_sum"])
    else:
        phi_sum = (
            df_phi.groupby(["xq2_bin", "z_bin", "pt2_bin"], as_index=False)
                  .agg(phi_content_sum=(phi_content_col, "sum"),
                       phi_bin_count=("phi_bin", "size"))
        )
        phi_sum["phi_error_sum"] = np.nan

    # Align bin numbering:
    #   phi CSV: z,pt2 start at 1
    #   spec CSV: z,pt2 start at 0
    phi_sum["z_bin0"]   = phi_sum["z_bin"] - 1
    phi_sum["pt2_bin0"] = phi_sum["pt2_bin"] - 1
    phi_sum = phi_sum[(phi_sum["z_bin0"] >= 0) & (phi_sum["pt2_bin0"] >= 0)].copy()

    out_phi_sum_csv = os.path.join(outdir, "phi_summed_by_xq2_z_pt2.csv")
    phi_sum.to_csv(out_phi_sum_csv, index=False)

    # =============================================================================
    # 2) MERGE with spectrum CSV (avoid suffix mess by using z_bin0/pt2_bin0)
    # =============================================================================
    df_spec_small = df_spec[["xq2_bin", "z_bin", "pt2_bin", spec_content_col]].copy()
    df_spec_small = df_spec_small.rename(columns={spec_content_col: "spec_content"})

    if spec_error_col in df_spec.columns:
        df_spec_small["spec_error"] = df_spec[spec_error_col].values
    else:
        df_spec_small["spec_error"] = np.nan

    phi_sum_small = phi_sum[["xq2_bin", "z_bin0", "pt2_bin0", "phi_content_sum", "phi_error_sum", "phi_bin_count"]].copy()

    merged = phi_sum_small.merge(
        df_spec_small,
        left_on=["xq2_bin", "z_bin0", "pt2_bin0"],
        right_on=["xq2_bin", "z_bin", "pt2_bin"],
        how="outer",
        indicator=True
    )

    merged["diff(phiSum - spec)"] = merged["phi_content_sum"] - merged["spec_content"]
    merged["ratio(phiSum/spec)"] = merged["phi_content_sum"] / merged["spec_content"].replace(0, np.nan)

    if skip_allzero_rows:
        merged = merged[~((merged["phi_content_sum"].fillna(0) == 0) &
                          (merged["spec_content"].fillna(0) == 0))].copy()

    out_merged_csv = os.path.join(outdir, "phiSum_vs_spec_merged.csv")
    merged.to_csv(out_merged_csv, index=False)

    # =============================================================================
    # 3) Unified 0-based z/pt2 for plotting
    # =============================================================================
    merged["z0"] = merged["z_bin"].fillna(merged["z_bin0"])
    merged["pt20"] = merged["pt2_bin"].fillna(merged["pt2_bin0"])
    merged["z0"] = pd.to_numeric(merged["z0"], errors="coerce")
    merged["pt20"] = pd.to_numeric(merged["pt20"], errors="coerce")

    # =============================================================================
    # 4) PLOTS: one figure per xQ2, subplots per z, pT2 dependence overlay
    # =============================================================================
    xq2_values = sorted(set(merged["xq2_bin"].dropna().astype(int).tolist()))
    if max_xq2 is not None:
        xq2_values = [x for x in xq2_values if x <= max_xq2]

    ncols = max(1, int(ncols))

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

        nrows = int(math.ceil(len(z_values) / ncols))
        fig_w = max(12, 3.2 * ncols)
        fig_h = max(6,  2.8 * nrows)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_w, fig_h), sharex=True, sharey=False)
        axes = np.array(axes).reshape(-1)

        fig.suptitle(f"xQ2 bin {xq2}: pT² dependence (subplots per z) — phi-sum vs spec", y=0.995)

        handle_phi = None
        handle_spec = None

        for idx, z in enumerate(z_values):
            ax = axes[idx]
            sub = df_x[df_x["z0"] == z].copy().sort_values("pt20")

            x = sub["pt20"].to_numpy()
            y_phi = sub["phi_content_sum"].to_numpy()
            y_spec = sub["spec_content"].to_numpy()

            mphi = np.isfinite(y_phi)
            if mphi.any():
                (hphi,) = ax.plot(x[mphi], y_phi[mphi], marker="o", linestyle="-")
                if handle_phi is None:
                    handle_phi = hphi

            mspec = np.isfinite(y_spec)
            if mspec.any():
                (hspec,) = ax.plot(x[mspec], y_spec[mspec], marker="x", linestyle="--")
                if handle_spec is None:
                    handle_spec = hspec

            ax.set_title(f"z bin {z}", fontsize=10)
            ax.grid(True, which="both", linewidth=0.5)

            if (idx % ncols) == 0:
                ax.set_ylabel("content")

        for j in range(len(z_values), len(axes)):
            axes[j].axis("off")

        for ax in axes[max(0, (nrows - 1) * ncols): (nrows * ncols)]:
            if ax.has_data():
                ax.set_xlabel("pT² bin (0-based)")

        handles, labels = [], []
        if handle_phi is not None:
            handles.append(handle_phi); labels.append("phi-sum")
        if handle_spec is not None:
            handles.append(handle_spec); labels.append("spec")
        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.98), fontsize=10)

        fig.tight_layout(rect=[0, 0, 0.98, 0.97])

        png = os.path.join(outdir, f"xq2_{xq2:02d}_pt2_by_z_subplots.png")
        pdf = os.path.join(outdir, f"xq2_{xq2:02d}_pt2_by_z_subplots.pdf")
        fig.savefig(png, dpi=150)
        fig.savefig(pdf)
        plt.close(fig)

    print("\nWrote:")
    print("  ", out_phi_sum_csv)
    print("  ", out_merged_csv)
    print("\nMerge status counts:")
    print(merged["_merge"].value_counts(dropna=False))
    print("\nPlots saved in:", outdir)


def main():
    ap = argparse.ArgumentParser(
        description="Compare (phi-binned hist CSV) summed over phi vs (no-phi initial spectrum CSV). "
                    "Makes one figure per xQ2 with subplots per z showing pT2 dependence."
    )
    ap.add_argument("--mode", choices=sorted(DEFAULTS.keys()), default="rec",
                    help="Which default pair to compare (rec/gen/data).")
    ap.add_argument("--phi_csv", default=None, help="Override phi CSV path.")
    ap.add_argument("--spec_csv", default=None, help="Override spec CSV path.")
    ap.add_argument("--outdir", default=None, help="Override output directory.")
    ap.add_argument("--phi_content_col", default="content")
    ap.add_argument("--spec_content_col", default="content")
    ap.add_argument("--phi_error_col", default="error")
    ap.add_argument("--spec_error_col", default="error")
    ap.add_argument("--phi_xq2_shift", type=int, default=-1,
                    help="Shift xq2_bin in phi CSV by this integer (default: -1).")
    ap.add_argument("--skip_allzero_rows", action="store_true")
    ap.add_argument("--max_xq2", type=int, default=None)
    ap.add_argument("--ncols", type=int, default=4)
    args = ap.parse_args()

    d = DEFAULTS[args.mode]
    phi_csv = args.phi_csv if args.phi_csv is not None else d["phi_csv"]
    spec_csv = args.spec_csv if args.spec_csv is not None else d["spec_csv"]
    outdir = args.outdir if args.outdir is not None else d["outdir"]

    compare_pair(
        phi_csv=phi_csv,
        spec_csv=spec_csv,
        outdir=outdir,
        phi_content_col=args.phi_content_col,
        spec_content_col=args.spec_content_col,
        phi_error_col=args.phi_error_col,
        spec_error_col=args.spec_error_col,
        phi_xq2_shift=args.phi_xq2_shift,
        skip_allzero_rows=args.skip_allzero_rows,
        max_xq2=args.max_xq2,
        ncols=args.ncols,
    )


if __name__ == "__main__":
    main()
