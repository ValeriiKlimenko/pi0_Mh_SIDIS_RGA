#!/usr/bin/env python3
from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_scale_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Drop accidental unnamed index col if present
    if "bin" not in df.columns and df.columns[0].lower().startswith("unnamed"):
        df = df.drop(columns=[df.columns[0]])

    needed = {"bin", "content", "acceptance"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Scale CSV missing columns {sorted(missing)}. Found: {list(df.columns)}")

    df = df.copy()
    df["bin"] = pd.to_numeric(df["bin"], errors="coerce")
    df["content"] = pd.to_numeric(df["content"], errors="coerce")
    df["acceptance"] = pd.to_numeric(df["acceptance"], errors="coerce")
    df = df.dropna(subset=["bin"]).copy()
    df["bin"] = df["bin"].astype(int)
    return df


def load_parvec_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Drop accidental unnamed index col if present
    if "xq2_bin" not in df.columns and df.columns[0].lower().startswith("unnamed"):
        df = df.drop(columns=[df.columns[0]])

    needed = {"xq2_bin", "dis"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Parvec CSV missing columns {sorted(missing)}. Found: {list(df.columns)}")

    df = df.copy()
    df["xq2_bin"] = pd.to_numeric(df["xq2_bin"], errors="coerce")
    df["dis"] = pd.to_numeric(df["dis"], errors="coerce")
    df = df.dropna(subset=["xq2_bin"]).copy()
    df["xq2_bin"] = df["xq2_bin"].astype(int)
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot (content*acceptance), DIS, and ratio DIS/(content*acceptance) vs xq2_bin.")
    ap.add_argument("--scale", default="unfolded_first16.csv", help="Scale CSV (bin, content, acceptance)")
    ap.add_argument("--parvec", default="dis_parvec.csv", help="Parvec CSV (xq2_bin, dis)")
    ap.add_argument("--bin_offset", type=int, default=2, help="Mapping: xq2_bin = bin - bin_offset (default 2)")
    ap.add_argument("--out_png", default="content_times_acceptance_vs_dis_with_ratio.png", help="Output PNG filename")
    ap.add_argument("--no_show", action="store_true", help="Do not open an interactive window")
    args = ap.parse_args()

    scale = load_scale_csv(args.scale)
    parv = load_parvec_csv(args.parvec)

    # Map: xq2_bin = bin - 2
    scale = scale.copy()
    scale["xq2_bin"] = scale["bin"] - int(args.bin_offset)

    # Compute content*acceptance
    scale["content_x_acceptance"] = scale["content"] * scale["acceptance"]

    # Merge for plotting on same x
    merged = scale[["xq2_bin", "content_x_acceptance"]].merge(
        parv[["xq2_bin", "dis"]],
        on="xq2_bin",
        how="outer",
    ).sort_values("xq2_bin")

    # Ratio: dis / (content*acceptance)
    num = merged["dis"].to_numpy(dtype=float)
    den = merged["content_x_acceptance"].to_numpy(dtype=float)
    ratio = np.full_like(num, np.nan, dtype=float)
    good = np.isfinite(num) & np.isfinite(den) & (den != 0.0)
    ratio[good] = num[good] / den[good]
    merged["ratio_dis_over_scale"] = ratio

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    axr = ax.twinx()  # second y-axis for ratio

    # content*acceptance
    m1 = np.isfinite(merged["content_x_acceptance"])
    if m1.any():
        ax.plot(
            merged.loc[m1, "xq2_bin"],
            merged.loc[m1, "content_x_acceptance"],
            marker="o",
            linestyle="-",
            label="DIS Valerii (acc*unf)",
        )

    # dis
    m2 = np.isfinite(merged["dis"])
    if m2.any():
        ax.plot(
            merged.loc[m2, "xq2_bin"],
            merged.loc[m2, "dis"],
            marker="s",
            linestyle="-",
            label="DIS Marshall",
        )

    # ratio on right axis
    mr = np.isfinite(merged["ratio_dis_over_scale"])
    if mr.any():
        axr.plot(
            merged.loc[mr, "xq2_bin"],
            merged.loc[mr, "ratio_dis_over_scale"],
            marker="^",
            linestyle="--",
            color = "#2ca02c",
            label="ratio: DIS / (acc*unf)",
        )
        axr.axhline(1.0, linestyle=":", linewidth=1)

    ax.set_xlabel("xq2_bin")
    ax.set_ylabel("value")
    axr.set_ylabel("ratio")

    ax.set_title("content*acceptance vs DIS (+ ratio DIS/(content*acceptance))")
    ax.grid(True, alpha=0.3)

    # Combined legend (left + right axis)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    if h1 or h2:
        ax.legend(h1 + h2, l1 + l2, frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(args.out_png, dpi=200)
    print(f"Wrote: {args.out_png}")

    if not args.no_show:
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
