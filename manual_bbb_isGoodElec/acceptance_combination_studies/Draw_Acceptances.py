#!/usr/bin/env python3
import os
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SCALE_CSV_FILE_DEFAULT = "binbybin_unfolding.csv"

def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def load_scale_map(scale_csv_file: str) -> dict[int, tuple[float, float]]:
    """
    Supports two formats:

    (A) Old:
        columns: bin, content, acceptance
        returns: bin -> (content_scale, acceptance_scale)

    (B) New (binbybin_unfolding.csv):
        columns include: 'rec/gen' (used as acceptance_scale)
        bin is inferred from row number: row i -> bin=i+1
        content_scale is taken from 'data*gen/rec' if present, else NaN.
    """
    # read CSV (comma) first; if it looks wrong, fall back to whitespace
    df = pd.read_csv(scale_csv_file)
    df.columns = [c.strip() for c in df.columns]
    if df.shape[1] == 1:
        df = pd.read_csv(scale_csv_file, sep=r"\s+", engine="python")
        df.columns = [c.strip() for c in df.columns]

    # ---- Old format ----
    if {"bin", "content", "acceptance"}.issubset(df.columns):
        m: dict[int, tuple[float, float]] = {}
        for _, r in df.iterrows():
            b = int(r["bin"])
            c = float(r["content"])
            a = float(r["acceptance"])
            m[b] = (c, a)
        return m

    # ---- New bin-by-bin format ----
    if "rec/gen" in df.columns:
        # infer 1-based "bin" from row number
        df = df.copy()
        df["bin"] = np.arange(1, len(df) + 1, dtype=int)

        # acceptance_scale := rec/gen
        df["acceptance"] = pd.to_numeric(df["rec/gen"], errors="coerce")

        # optional "content" (not used by your script, but keep tuple structure)
        if "data*gen/rec" in df.columns:
            df["content"] = pd.to_numeric(df["data*gen/rec"], errors="coerce")
        else:
            df["content"] = np.nan

        m: dict[int, tuple[float, float]] = {}
        for _, r in df.iterrows():
            b = int(r["bin"])
            c = float(r["content"]) if np.isfinite(r["content"]) else float("nan")
            a = float(r["acceptance"]) if np.isfinite(r["acceptance"]) else float("nan")
            m[b] = (c, a)
        return m

    raise ValueError(
        f"SCALE: unsupported columns in {scale_csv_file}. "
        f"Expected either (bin,content,acceptance) or a column named 'rec/gen'. "
        f"Found: {list(df.columns)}"
    )


def phi_sum(df: pd.DataFrame,
            content_col: str = "content",
            error_col: str = "error",
            xq2_shift: int = 0,
            max_xq2: int | None = None,
            skip_zeros: bool = False) -> pd.DataFrame:
    req = {"xq2_bin", "z_bin", "pt2_bin", "phi_bin", content_col}
    missing = req - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    for c in ["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df[content_col] = pd.to_numeric(df[content_col], errors="coerce")
    if error_col in df.columns:
        df[error_col] = pd.to_numeric(df[error_col], errors="coerce")
    else:
        df[error_col] = np.nan

    df = df.dropna(subset=["xq2_bin", "z_bin", "pt2_bin", "phi_bin", content_col]).copy()

    for c in ["xq2_bin", "z_bin", "pt2_bin", "phi_bin"]:
        df[c] = df[c].astype(int)

    if xq2_shift != 0:
        df["xq2_bin"] = df["xq2_bin"] + int(xq2_shift)
        df = df[df["xq2_bin"] >= 0].copy()

    if max_xq2 is not None:
        df = df[df["xq2_bin"] <= int(max_xq2)].copy()

    if skip_zeros:
        e = df[error_col].fillna(0.0)
        df = df[~((df[content_col] == 0.0) & (e == 0.0))].copy()

    df["_err2"] = df[error_col].fillna(0.0) ** 2
    out = (df.groupby(["xq2_bin", "z_bin", "pt2_bin"], as_index=False)
             .agg(content_sum=(content_col, "sum"),
                  err2_sum=("_err2", "sum"),
                  phi_bin_count=("phi_bin", "size")))
    out["error_sum"] = np.sqrt(out["err2_sum"])
    out = out.drop(columns=["err2_sum"])
    return out


def compute_eff(rec_sum: pd.DataFrame, gen_sum: pd.DataFrame, tag: str) -> pd.DataFrame:
    a = rec_sum.rename(columns={"content_sum": f"rec_{tag}", "error_sum": f"recerr_{tag}"})
    b = gen_sum.rename(columns={"content_sum": f"gen_{tag}", "error_sum": f"generr_{tag}"})

    m = a.merge(b, on=["xq2_bin", "z_bin", "pt2_bin"], how="outer", indicator=False)

    rec = m[f"rec_{tag}"].to_numpy(dtype=float)
    gen = m[f"gen_{tag}"].to_numpy(dtype=float)
    recerr = m[f"recerr_{tag}"].to_numpy(dtype=float)
    generr = m[f"generr_{tag}"].to_numpy(dtype=float)

    eff = np.divide(rec, gen, out=np.full_like(rec, np.nan), where=(gen != 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        term_rec = np.where(rec != 0.0, (recerr / rec) ** 2, np.nan)
        term_gen = np.where(gen != 0.0, (generr / gen) ** 2, np.nan)
        efferr = np.abs(eff) * np.sqrt(term_rec + term_gen)

    m[f"eff_{tag}"] = eff
    m[f"efferr_{tag}"] = efferr
    return m


def main():
    ap = argparse.ArgumentParser(
        description="Compare REC/GEN between two sets of ROOT-exported phi-binned CSVs "
                    "(h_meas_mc vs h_true_mc), with plots per xQ2 and subplots per z."
    )

    ap.add_argument("--A_rec", default="noGoodElec_h_meas_mc.csv")
    ap.add_argument("--A_gen", default="noGoodElec_h_true_mc.csv")
    ap.add_argument("--B_rec", default="isGoodElec_h_meas_mc.csv")
    ap.add_argument("--B_gen", default="isGoodElec_h_true_mc.csv")

    ap.add_argument("--outdir_noacc", default=None,
                  help="Output dir for plots WITHOUT acceptance correction to A. "
                       "Default: <outdir>_noAcc")

    ap.add_argument("--outdir", default="compare_RECGEN_A_divAcc_vs_B")
    ap.add_argument("--xq2_shift", type=int, default=0)
    ap.add_argument("--max_xq2", type=int, default=None)
    ap.add_argument("--skip_zeros", action="store_true")

    ap.add_argument("--z0", action="store_true")
    ap.add_argument("--pt20", action="store_true")

    ap.add_argument("--ncols", type=int, default=4)
    ap.add_argument("--ratio_min", type=float, default=None)
    ap.add_argument("--ratio_max", type=float, default=None)

    # acceptance scale map (bin == xq2_bin+1)
    ap.add_argument("--scale_csv", default=SCALE_CSV_FILE_DEFAULT,
                    help="CSV providing correction factors. If using binbybin_unfolding.csv, "
                         "it must contain column 'rec/gen'. Bin is inferred from row number "
                         "(row i -> bin=i+1), and the script maps bin == xq2_bin+1.")
    args = ap.parse_args()

    A_rec_csv, A_gen_csv = args.A_rec, args.A_gen
    B_rec_csv, B_gen_csv = args.B_rec, args.B_gen
    labelA = os.path.basename(A_rec_csv).replace("_h_meas_mc.csv", "")
    labelB = os.path.basename(B_rec_csv).replace("_h_meas_mc.csv", "")

    os.makedirs(args.outdir, exist_ok=True)

    outdir_noacc = args.outdir_noacc or (args.outdir + "_noAcc")
    os.makedirs(outdir_noacc, exist_ok=True)

    dfA_rec = load_csv(A_rec_csv)
    dfA_gen = load_csv(A_gen_csv)
    dfB_rec = load_csv(B_rec_csv)
    dfB_gen = load_csv(B_gen_csv)

    A_rec_sum = phi_sum(dfA_rec, xq2_shift=args.xq2_shift, max_xq2=args.max_xq2, skip_zeros=args.skip_zeros)
    A_gen_sum = phi_sum(dfA_gen, xq2_shift=args.xq2_shift, max_xq2=args.max_xq2, skip_zeros=args.skip_zeros)
    B_rec_sum = phi_sum(dfB_rec, xq2_shift=args.xq2_shift, max_xq2=args.max_xq2, skip_zeros=args.skip_zeros)
    B_gen_sum = phi_sum(dfB_gen, xq2_shift=args.xq2_shift, max_xq2=args.max_xq2, skip_zeros=args.skip_zeros)

    effA = compute_eff(A_rec_sum, A_gen_sum, "A")
    effB = compute_eff(B_rec_sum, B_gen_sum, "B")

    merged = effA.merge(
        effB,
        on=["xq2_bin", "z_bin", "pt2_bin"],
        how="outer",
        suffixes=("_A", "_B"),
    )

    # ---------------- apply acceptance_scale to A by DIVISION ----------------
    scale_map = load_scale_map(args.scale_csv)

    merged["acceptance_scale"] = merged["xq2_bin"].map(
        lambda ix: scale_map.get(int(ix + 1), (np.nan, np.nan))[1]
    ).astype(float)

    acc = merged["acceptance_scale"].to_numpy(dtype=float)

    merged["eff_A_divAcc"] = merged["eff_A"] / np.where(acc != 0.0, acc, np.nan)
    merged["efferr_A_divAcc"] = merged["efferr_A"] / np.where(acc != 0.0, np.abs(acc), np.nan)

    # Ratio (A/acc) / B
    merged["eff_ratio_AdivAcc_over_B"] = merged["eff_A_divAcc"] / merged["eff_B"].replace(0, np.nan)


    # ---------------- ALSO: ratio WITHOUT acceptance correction to A ----------------
    merged["eff_ratio_A_over_B"] = merged["eff_A"] / merged["eff_B"].replace(0, np.nan)
    
    with np.errstate(divide="ignore", invalid="ignore"):
        ra = merged["eff_A"].to_numpy(dtype=float)
        rb = merged["eff_B"].to_numpy(dtype=float)
        ea = merged["efferr_A"].to_numpy(dtype=float)
        eb = merged["efferr_B"].to_numpy(dtype=float)
        r  = merged["eff_ratio_A_over_B"].to_numpy(dtype=float)
    
        term_a = np.where(ra != 0.0, (ea / ra) ** 2, np.nan)
        term_b = np.where(rb != 0.0, (eb / rb) ** 2, np.nan)
        merged["eff_ratio_err_noacc"] = np.abs(r) * np.sqrt(term_a + term_b)


  
    # ratio error propagation
    with np.errstate(divide="ignore", invalid="ignore"):
        ra = merged["eff_A_divAcc"].to_numpy(dtype=float)
        rb = merged["eff_B"].to_numpy(dtype=float)
        ea = merged["efferr_A_divAcc"].to_numpy(dtype=float)
        eb = merged["efferr_B"].to_numpy(dtype=float)
        r = merged["eff_ratio_AdivAcc_over_B"].to_numpy(dtype=float)

        term_a = np.where(ra != 0.0, (ea / ra) ** 2, np.nan)
        term_b = np.where(rb != 0.0, (eb / rb) ** 2, np.nan)
        merged["eff_ratio_err"] = np.abs(r) * np.sqrt(term_a + term_b)

    merged["z_plot"] = merged["z_bin"] - 1 if args.z0 else merged["z_bin"]
    merged["pt2_plot"] = merged["pt2_bin"] - 1 if args.pt20 else merged["pt2_bin"]

    merged_csv = os.path.join(args.outdir, "rec_gen_eff_compare_merged.csv")
    merged.to_csv(merged_csv, index=False)

    xq2_values = sorted(pd.to_numeric(merged["xq2_bin"], errors="coerce").dropna().astype(int).unique().tolist())
    if args.max_xq2 is not None:
        xq2_values = [x for x in xq2_values if x <= args.max_xq2]

    ncols = max(1, int(args.ncols))

    for xq2 in xq2_values:
        df_x = merged[merged["xq2_bin"] == xq2].copy()
        df_x = df_x.dropna(subset=["z_plot", "pt2_plot"], how="any").copy()
        if df_x.empty:
            continue

        df_x["z_plot"] = pd.to_numeric(df_x["z_plot"], errors="coerce").astype(int)
        df_x["pt2_plot"] = pd.to_numeric(df_x["pt2_plot"], errors="coerce").astype(int)

        z_values = sorted(df_x["z_plot"].unique().tolist())
        if not z_values:
            continue

        nrows = int(math.ceil(len(z_values) / ncols))
        fig_w = max(12, 3.2 * ncols)
        fig_h = max(6,  2.8 * nrows)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_w, fig_h), sharex=True, sharey=False)
        axes = np.array(axes).reshape(-1)

        acc_val = scale_map.get(int(xq2 + 1), (np.nan, np.nan))[1]
        fig.suptitle(
            f"xQ2 bin {xq2}: pT² dependence — ({labelA} / acc_scale) vs {labelB}  [acc_scale={acc_val:.4g}] + ratio((A/acc)/B)",
            y=0.995
        )

        handle_A = None
        handle_B = None
        handle_R = None

        for idx, z in enumerate(z_values):
            ax = axes[idx]
            sub = df_x[df_x["z_plot"] == z].copy().sort_values("pt2_plot")

            x = sub["pt2_plot"].to_numpy()
            yA = sub["eff_A_divAcc"].to_numpy()
            yB = sub["eff_B"].to_numpy()
            r  = sub["eff_ratio_AdivAcc_over_B"].to_numpy()

            mA = np.isfinite(yA)
            if mA.any():
                (hA,) = ax.plot(x[mA], yA[mA], marker="o", linestyle="-", label=f"{labelA} Separate")
                if handle_A is None:
                    handle_A = hA

            mB = np.isfinite(yB)
            if mB.any():
                (hB,) = ax.plot(x[mB], yB[mB], marker="x", linestyle="--", label=f"{labelB} Combined")
                if handle_B is None:
                    handle_B = hB

            ax.set_title(f"z bin {z}", fontsize=10)
            ax.grid(True, which="both", linewidth=0.5)

            if (idx % ncols) == 0:
                ax.set_ylabel("REC/GEN")

            ax2 = ax.twinx()
            mR = np.isfinite(r)
            if mR.any():
                (hR,) = ax2.plot(
                    x[mR], r[mR],
                    marker="s", linestyle=":",
                    color="green",          # or "tab:green"
                    label="Separate/Combined"
                )
                if handle_R is None:
                    handle_R = hR

                if args.ratio_min is not None or args.ratio_max is not None:
                    ylo = args.ratio_min if args.ratio_min is not None else float(np.nanmin(r[mR]))
                    yhi = args.ratio_max if args.ratio_max is not None else float(np.nanmax(r[mR]))
                else:
                    rmin = float(np.nanmin(r[mR]))
                    rmax = float(np.nanmax(r[mR]))
                    ylo = min(0.9, rmin)
                    yhi = max(1.07, rmax)

                if not np.isfinite(ylo) or not np.isfinite(yhi) or yhi <= ylo:
                    ylo, yhi = 0.75, 1.0
                ax2.set_ylim(ylo, yhi)

            ax2.axhline(1.0, linestyle=":", linewidth=1)

            if (idx % ncols) == (ncols - 1):
                ax2.set_ylabel("ratio")
            else:
                ax2.set_ylabel("")

        for j in range(len(z_values), len(axes)):
            axes[j].axis("off")

        for ax in axes[max(0, (nrows - 1) * ncols): (nrows * ncols)]:
            if ax.has_data():
                ax.set_xlabel("pT² bin" + (" (0-based)" if args.pt20 else ""))

        handles, labels = [], []
        if handle_A is not None:
            handles.append(handle_A); labels.append(f"Separate")
        if handle_B is not None:
            handles.append(handle_B); labels.append(f"Combined")
        if handle_R is not None:
            handles.append(handle_R); labels.append("Separate/Combined")
        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.98), fontsize=10)

        fig.tight_layout(rect=[0, 0, 0.98, 0.97])

        png = os.path.join(args.outdir, f"xq2_{xq2:02d}_pt2_by_z_recgen_compare.png")
        pdf = os.path.join(args.outdir, f"xq2_{xq2:02d}_pt2_by_z_recgen_compare.pdf")
        fig.savefig(png, dpi=150)
        fig.savefig(pdf)
        plt.close(fig)


        # ============================================================
        # SECOND FIGURE: WITHOUT acceptance correction to A
        # ============================================================
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_w, fig_h), sharex=True, sharey=False)
        axes = np.array(axes).reshape(-1)
        
        fig.suptitle(
            f"xQ2 bin {xq2}: pT² dependence — {labelA} vs {labelB} (NO acc correction) + ratio(A/B)",
            y=0.995
        )
        
        handle_A = None
        handle_B = None
        handle_R = None
        
        for idx, z in enumerate(z_values):
            ax = axes[idx]
            sub = df_x[df_x["z_plot"] == z].copy().sort_values("pt2_plot")
        
            x = sub["pt2_plot"].to_numpy()
            yA = sub["eff_A"].to_numpy()                    # <-- NO acc correction
            yB = sub["eff_B"].to_numpy()
            r  = sub["eff_ratio_A_over_B"].to_numpy()       # <-- NO acc correction ratio
        
            mA = np.isfinite(yA)
            if mA.any():
                (hA,) = ax.plot(x[mA], yA[mA], marker="o", linestyle="-", label=f"{labelA} Separate")
                if handle_A is None:
                    handle_A = hA
        
            mB = np.isfinite(yB)
            if mB.any():
                (hB,) = ax.plot(x[mB], yB[mB], marker="x", linestyle="--", label=f"{labelB} Combined")
                if handle_B is None:
                    handle_B = hB
        
            ax.set_title(f"z bin {z}", fontsize=10)
            ax.grid(True, which="both", linewidth=0.5)
        
            if (idx % ncols) == 0:
                ax.set_ylabel("REC/GEN")
        
            ax2 = ax.twinx()
            mR = np.isfinite(r)
            if mR.any():
                (hR,) = ax2.plot(
                    x[mR], r[mR],
                    marker="s", linestyle=":",
                    color="green",
                    label="Separate/Combined"
                )
                if handle_R is None:
                    handle_R = hR
        
                # same ratio y-limits logic as your acc-corrected plot
                if args.ratio_min is not None or args.ratio_max is not None:
                    ylo = args.ratio_min if args.ratio_min is not None else float(np.nanmin(r[mR]))
                    yhi = args.ratio_max if args.ratio_max is not None else float(np.nanmax(r[mR]))
                else:
                    rmin = float(np.nanmin(r[mR]))
                    rmax = float(np.nanmax(r[mR]))
                    ylo = min(0.5, rmin)
                    yhi = max(0.8, rmax)
        
                if not np.isfinite(ylo) or not np.isfinite(yhi) or yhi <= ylo:
                    ylo, yhi = 0.75, 1.0
                ax2.set_ylim(ylo, yhi)
        
            ax2.axhline(1.0, linestyle=":", linewidth=1)
            if (idx % ncols) == (ncols - 1):
                ax2.set_ylabel("ratio")
            else:
                ax2.set_ylabel("")
        
        # Hide unused axes
        for j in range(len(z_values), len(axes)):
            axes[j].axis("off")
        
        # x-label only on bottom row
        for ax in axes[max(0, (nrows - 1) * ncols): (nrows * ncols)]:
            if ax.has_data():
                ax.set_xlabel("pT² bin" + (" (0-based)" if args.pt20 else ""))
        
        handles, labels = [], []
        if handle_A is not None:
            handles.append(handle_A); labels.append("Separate")
        if handle_B is not None:
            handles.append(handle_B); labels.append("Combined")
        if handle_R is not None:
            handles.append(handle_R); labels.append("Separate/Combined")
        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.98), fontsize=10)
        
        fig.tight_layout(rect=[0, 0, 0.98, 0.97])
        
        png = os.path.join(outdir_noacc, f"xq2_{xq2:02d}_pt2_by_z_recgen_compare.png")
        pdf = os.path.join(outdir_noacc, f"xq2_{xq2:02d}_pt2_by_z_recgen_compare.pdf")
        fig.savefig(png, dpi=150)
        fig.savefig(pdf)
        plt.close(fig)

  

    print("Wrote:")
    print("  ", merged_csv)
    print("Plots saved in:", args.outdir)


if __name__ == "__main__":
    main()