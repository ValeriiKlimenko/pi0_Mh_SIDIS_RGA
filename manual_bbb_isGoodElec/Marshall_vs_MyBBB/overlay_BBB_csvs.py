#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SCALE_CSV_FILE_DEFAULT = "unfolded_first16.csv"

PT_EDGES = np.array(
    [0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00, 1.50],
    dtype=float
)

Z_EDGES = np.array(
    [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0],
    dtype=float
)

SKIP_Z_BINS = {1, 8}
SKIP_PT2_BINS = {1, 11}

TWO_PI = 2.0 * np.pi


def _sym_ref_err(eyl: pd.Series, eyh: pd.Series) -> np.ndarray:
    return 0.5 * (eyl.to_numpy(dtype=float) + eyh.to_numpy(dtype=float))


def _ratio_with_err(num: np.ndarray, enum: np.ndarray,
                    den: np.ndarray, eden: np.ndarray):
    """
    r = num/den
    σr = |r| * sqrt( (σnum/num)^2 + (σden/den)^2 )
    Returns (r, σr) with NaN where invalid.
    """
    num = np.asarray(num, dtype=float)
    enum = np.asarray(enum, dtype=float)
    den = np.asarray(den, dtype=float)
    eden = np.asarray(eden, dtype=float)

    r = np.full_like(num, np.nan, dtype=float)
    er = np.full_like(num, np.nan, dtype=float)

    good = np.isfinite(num) & np.isfinite(enum) & np.isfinite(den) & np.isfinite(eden) & (den != 0.0)
    if not np.any(good):
        return r, er

    r_good = num[good] / den[good]

    # protect against num==0 in the (σnum/num) term
    rel_num2 = np.zeros_like(r_good)
    nnz = (num[good] != 0.0)
    rel_num2[nnz] = (enum[good][nnz] / np.abs(num[good][nnz])) ** 2

    rel_den2 = (eden[good] / np.abs(den[good])) ** 2
    er_good = np.abs(r_good) * np.sqrt(rel_num2 + rel_den2)

    r[good] = r_good
    er[good] = er_good
    return r, er


def pt2_value_to_bin(x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    b = np.digitize(x, edges, right=False)              # left-inclusive
    b = np.where(b == len(edges), len(edges) - 1, b)    # clamp x==last_edge
    b = np.where((b >= 1) & (b <= len(edges) - 1), b, 0)
    return b.astype(int)


def pt2_bin_center(pt2_bin: np.ndarray, edges: np.ndarray) -> np.ndarray:
    pt2_bin = np.asarray(pt2_bin, dtype=int)
    centers = np.full(pt2_bin.shape, np.nan, dtype=float)
    ok = (pt2_bin >= 1) & (pt2_bin <= (len(edges) - 1))
    i = pt2_bin[ok] - 1
    centers[ok] = 0.5 * (edges[i] + edges[i + 1])
    return centers


def z_bin_center(z_bin: np.ndarray, edges: np.ndarray) -> np.ndarray:
    z_bin = np.asarray(z_bin, dtype=int)
    centers = np.full(z_bin.shape, np.nan, dtype=float)
    ok = (z_bin >= 1) & (z_bin <= (len(edges) - 1))
    i = z_bin[ok] - 1
    centers[ok] = 0.5 * (edges[i] + edges[i + 1])
    return centers


def read_method_csv(path: str, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"xq2_bin", "z_bin", "pt2_bin", "content", "error"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{label}: missing columns {sorted(missing)} in {path}")

    out = df[["xq2_bin", "z_bin", "pt2_bin", "content", "error"]].copy()
    out["xq2_bin"] = out["xq2_bin"].astype(int)
    out["z_bin"]   = out["z_bin"].astype(int)
    out["pt2_bin"] = out["pt2_bin"].astype(int)
    out = out.rename(columns={"content": f"{label}_content", "error": f"{label}_error"})
    return out


def read_reference_table(path: str, z_offset: int = 2) -> pd.DataFrame:
    """
    Reads table with columns: ixq2 iz point x y exl exh eyl eyh
    Accepts comma-CSV or whitespace/tab.
    """
    needed = {"ixq2", "iz", "x", "y", "eyl", "eyh"}

    df = pd.read_csv(path)
    if not needed.issubset(df.columns):
        df = pd.read_csv(path, sep=r"\s+", engine="python")

    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"REF: missing columns {sorted(missing)} in {path}. Found: {list(df.columns)}")

    out = df.copy()
    out["xq2_bin"] = out["ixq2"].astype(int)
    out["z_bin"]   = out["iz"].astype(int) + int(z_offset)

    out["pt2_bin"] = pt2_value_to_bin(out["x"].to_numpy(), PT_EDGES)
    out = out[out["pt2_bin"] != 0].copy()

    out = out.rename(columns={"x": "ref_x", "y": "ref_y", "eyl": "ref_eyl", "eyh": "ref_eyh"})
    return out[["xq2_bin", "z_bin", "pt2_bin", "ref_x", "ref_y", "ref_eyl", "ref_eyh"]]


def load_scale_map(scale_csv_file: str) -> dict[int, tuple[float, float]]:
    needed = {"bin", "content", "acceptance"}
    df = pd.read_csv(scale_csv_file)
    if not needed.issubset(df.columns):
        df = pd.read_csv(scale_csv_file, sep=r"\s+", engine="python")

    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"SCALE: missing columns {sorted(missing)} in {scale_csv_file}. Found: {list(df.columns)}")

    m: dict[int, tuple[float, float]] = {}
    for _, r in df.iterrows():
        b = int(r["bin"])
        c = float(r["content"])
        a = float(r["acceptance"])
        m[b] = (c, a)
    return m


def apply_scaling_to_methods(
    merged: pd.DataFrame,
    scale_map: dict[int, tuple[float, float]],
    bin_width: np.ndarray,
) -> pd.DataFrame:
    """
    Scale m1 and m2 only:

      scale_factor = 1 / (bin_width[pt2_bin-1] * content_scale * acceptance_scale)

    where content_scale and acceptance_scale come from SCALE_CSV (row with bin = xq2_bin+1).
    """
    out = merged.copy()

    content_scale = out["xq2_bin"].map(lambda ix: scale_map.get(int(ix + 1), (np.nan, np.nan))[0]).astype(float)
    acc_scale     = out["xq2_bin"].map(lambda ix: scale_map.get(int(ix + 1), (np.nan, np.nan))[1]).astype(float)

    out["scale_content"] = content_scale
    out["scale_acc"]     = acc_scale

    pt = out["pt2_bin"].astype(int).to_numpy()
    bw = np.full(len(out), np.nan, dtype=float)
    ok_pt = (pt >= 1) & (pt <= len(bin_width))
    bw[ok_pt] = bin_width[pt[ok_pt] - 1]
    out["bin_width"] = bw

    sf = np.full(len(out), np.nan, dtype=float)
    ok = (
        np.isfinite(bw) & (bw != 0.0) &
        np.isfinite(content_scale) & (content_scale != 0.0) &
        np.isfinite(acc_scale) & (acc_scale != 0.0)
    )
    sf[ok] = 1.0 / (bw[ok] * content_scale.to_numpy()[ok] * acc_scale.to_numpy()[ok])

    out["scale_factor"] = sf

    # Apply to m1 and m2 only (NOT m3)
    for lbl in ["m1", "m2", "m3"]:
        yc = f"{lbl}_content"
        ye = f"{lbl}_error"
        if yc in out.columns:
            out[yc] = out[yc] * out["scale_factor"]
        if ye in out.columns:
            out[ye] = out[ye] * np.abs(out["scale_factor"])

    return out


def choose_subplot_grid(n: int) -> tuple[int, int]:
    ncols = min(4, max(1, n))
    nrows = int(np.ceil(n / ncols))
    return nrows, ncols


def main():
    ap = argparse.ArgumentParser(description="Overlay m1/m2/m3 with reference; ratio plots m1/ref, m2/ref, m3/ref.")
    ap.add_argument("--m1", default="dataGenOverRec_method1_sumPhi_then_ratio.csv", help="method1 CSV")
    ap.add_argument("--m2", default="dataGenOverRec_method2_ratio_then_avgPhi_unweighted.csv", help="method2 CSV")
    ap.add_argument("--m3", default="dataGenOverRec_method3_sumPhi_all_then_ratio.csv", help="method3 CSV")
    ap.add_argument("--ref", default="tgp_2026_fullData_v3.csv", help="reference table (comma or whitespace/tab)")
    ap.add_argument("--scale_csv", default=SCALE_CSV_FILE_DEFAULT,
                    help="scale CSV with columns bin, content, acceptance (bin == xq2_bin+1)")
    ap.add_argument("--z_offset", type=int, default=2, help="z_bin = iz + z_offset for reference (default 2)")
    ap.add_argument("--out_csv", default="overlay_merged.csv", help="output merged CSV")
    ap.add_argument("--plots_dir", default="Mh", help="if set: write PNG plots into this dir (one per xq2)")
    args = ap.parse_args()

    bin_width = np.asarray([
        0.005, 0.005, 0.005, 0.005,
        0.01,  0.01,  0.01,
        0.015, 0.015,
        0.02,
        0.05,
    ], dtype=float)

    df1 = read_method_csv(args.m1, "m1")
    df2 = read_method_csv(args.m2, "m2")
    df3 = read_method_csv(args.m3, "m3")
    dfr = read_reference_table(args.ref, z_offset=args.z_offset)

    # Average reference if multiple rows per (xq2,z,pt2)
    g = dfr.groupby(["xq2_bin", "z_bin", "pt2_bin"], as_index=False)
    dfr_b = g.agg(
        ref_x=("ref_x", "mean"),
        ref_y=("ref_y", "mean"),
        ref_eyl=("ref_eyl", lambda s: float(np.sqrt(np.mean(np.square(s))))),
        ref_eyh=("ref_eyh", lambda s: float(np.sqrt(np.mean(np.square(s))))),
    )

    key = ["xq2_bin", "z_bin", "pt2_bin"]
    merged = df1.merge(df2, on=key, how="outer")
    merged = merged.merge(df3, on=key, how="outer")
    merged = merged.merge(dfr_b, on=key, how="outer")

    # Skip requested bins globally
    merged = merged[
        (~merged["z_bin"].isin(SKIP_Z_BINS)) &
        (~merged["pt2_bin"].isin(SKIP_PT2_BINS))
    ].copy()

    merged["pt2_center"] = pt2_bin_center(merged["pt2_bin"].to_numpy(), PT_EDGES)
    merged["z_center"] = z_bin_center(merged["z_bin"].to_numpy(), Z_EDGES)

    # scaling (m1/m2 only)
    scale_map = load_scale_map(args.scale_csv)
    merged = apply_scaling_to_methods(merged, scale_map=scale_map, bin_width=bin_width)

    merged["m2_minus_m1"] = merged["m2_content"] - merged["m1_content"]
    merged["ref_minus_m1"] = merged["ref_y"] - merged["m1_content"]
    merged["ref_minus_m2"] = merged["ref_y"] - merged["m2_content"]

    merged = merged.sort_values(["xq2_bin", "z_bin", "pt2_bin"]).reset_index(drop=True)
    merged.to_csv(args.out_csv, index=False)
    print(f"Wrote merged overlay CSV: {args.out_csv}")

    if args.plots_dir:
        pdir = Path(args.plots_dir)
        pdir.mkdir(parents=True, exist_ok=True)

        ratio_dir = pdir / "ratios"
        ratio_dir.mkdir(parents=True, exist_ok=True)

        DRAW_COEF = TWO_PI

        for xq2, df_x in merged.groupby("xq2_bin"):
            z_bins = sorted([int(z) for z in df_x["z_bin"].dropna().unique()])
            if not z_bins:
                continue

            cs, ca = scale_map.get(int(xq2 + 1), (np.nan, np.nan))
            nrows, ncols = choose_subplot_grid(len(z_bins))

            # ----------- Overlay figure -----------
            W_PER_COL = 5.5
            H_PER_ROW = 3.2

            fig, axes = plt.subplots(
                nrows=nrows, ncols=ncols, sharex=True,
                figsize=(W_PER_COL * ncols, H_PER_ROW * nrows)
            )
            axes = np.array(axes).reshape(-1)

            # ----------- Ratio figure -----------
            figR, axesR = plt.subplots(
                nrows=nrows, ncols=ncols, sharex=True,
                figsize=(W_PER_COL * ncols, H_PER_ROW * nrows)
            )
            axesR = np.array(axesR).reshape(-1)

            for i, z in enumerate(z_bins):
                ax = axes[i]
                axR = axesR[i]

                sub = df_x[df_x["z_bin"] == z].sort_values("pt2_center")
                zc = z_bin_center(np.array([z]), Z_EDGES)[0]

                # ---------- overlay: m3 ----------
                m = sub["m3_content"].notna()
                if m.any():
                    ax.errorbar(
                        sub.loc[m, "pt2_center"],
                        sub.loc[m, "m3_content"],
                        yerr=sub.loc[m, "m3_error"],
                        fmt="D",
                        label="Primary workflow"
                    )

                # ---------- overlay: reference ----------
                mref = sub["ref_y"].notna()
                if mref.any():
                    xref = sub.loc[mref, "ref_x"].copy()
                    xref = xref.fillna(sub.loc[mref, "pt2_center"])
                    yerr_ref = np.vstack([sub.loc[mref, "ref_eyl"], sub.loc[mref, "ref_eyh"]])
                    ax.errorbar(
                        xref,
                        sub.loc[mref, "ref_y"],
                        yerr=yerr_ref,
                        fmt="^",
                        label="Legacy"
                    )

                ax.set_title(f"z = {zc:.2f}")
                ax.grid(True, alpha=0.3)

                # ---------------- Ratios m3/ref ----------------
                x = sub["pt2_center"].to_numpy(dtype=float)
                den = sub["ref_y"].to_numpy(dtype=float)
                eden = _sym_ref_err(sub["ref_eyl"], sub["ref_eyh"])

                num3 = sub["m3_content"].to_numpy(dtype=float)
                en3  = sub["m3_error"].to_numpy(dtype=float)
                r3, er3 = _ratio_with_err(num3, en3, den, eden)
                ok3 = np.isfinite(x) & np.isfinite(r3) & np.isfinite(er3)
                if np.any(ok3):
                    axR.errorbar(x[ok3], r3[ok3], yerr=er3[ok3], fmt="D", label="Primary / Legacy workflows")

                axR.axhline(1.0, linestyle="--", linewidth=1)
                axR.set_title(f"z = {zc:.2f}")
                axR.grid(True, alpha=0.3)

                if ax.get_legend_handles_labels()[0]:
                    ax.legend(fontsize=8)
                if axR.get_legend_handles_labels()[0]:
                    axR.legend(fontsize=8)

            # Hide unused axes
            for j in range(len(z_bins), len(axes)):
                axes[j].axis("off")
                axesR[j].axis("off")

            fig.suptitle(
                f"Overlay vs pT2 (xq2_bin={xq2}, scale_content={cs:.3g}, scale_acc={ca:.3g}, DRAW_COEF={DRAW_COEF:.6g}; m3 no coef)"
            )
            fig.supxlabel("pT2")
            fig.supylabel("Scaled value (as drawn)")
            fig.tight_layout(rect=[0, 0, 1, 0.95])

            out_png = pdir / f"overlay_xq2_{xq2}.png"
            fig.savefig(out_png, dpi=150)
            plt.close(fig)

            figR.suptitle(
                f"Ratios vs pT2 (xq2_bin={xq2}, scale_content={cs:.3g}, scale_acc={ca:.3g}, DRAW_COEF={DRAW_COEF:.6g}; m3 no coef)"
            )
            figR.supxlabel("pT2")
            figR.supylabel("ratio")
            figR.tight_layout(rect=[0, 0, 1, 0.95])

            out_png_R = ratio_dir / f"ratio_xq2_{xq2}.png"
            figR.savefig(out_png_R, dpi=150)
            plt.close(figR)

        print(f"Wrote PNG plots into: {args.plots_dir}")
        print(f"Wrote ratio PNG plots into: {ratio_dir}")


if __name__ == "__main__":
    main()