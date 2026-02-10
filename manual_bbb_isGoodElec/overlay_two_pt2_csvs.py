#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_pt_dependence_csv(path: str) -> pd.DataFrame:
    """
    Reads pt_dependence.csv-like files with columns:
      ix, iz, iPt, ratio_corr

    Accepts comma/semicolon/whitespace separators and ignores lines starting with '#'.
    """
    df = pd.read_csv(
        path,
        comment="#",
        header=None,
        sep=r"[\s,;]+",
        engine="python",
        names=["ix", "iz", "iPt", "ratio_corr"],
        usecols=[0, 1, 2, 3],
    )

    # Coerce numeric and drop bad lines
    for c in ["ix", "iz", "iPt", "ratio_corr"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["ix", "iz", "iPt", "ratio_corr"])

    df["ix"] = df["ix"].astype(int)
    df["iz"] = df["iz"].astype(int)
    df["iPt"] = df["iPt"].astype(int)
    return df


def _prep_axes(nplots: int, ncols: int):
    ncols_eff = max(1, min(ncols, nplots))
    nrows_eff = int(np.ceil(nplots / ncols_eff))

    fig_w = max(10, 4.0 * ncols_eff)
    fig_h = max(6, 3.0 * nrows_eff)
    fig, axes = plt.subplots(nrows_eff, ncols_eff, figsize=(fig_w, fig_h), sharex=True)

    axes = np.atleast_1d(axes).reshape(nrows_eff, ncols_eff)
    return fig, axes, nrows_eff, ncols_eff


def drop_first_last_ipt(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drop the smallest and largest iPt rows in a (ix, iz) series.
    If fewer than 3 points exist, returns empty (or keeps none) to avoid misleading plots.
    """
    if df.empty:
        return df
    df = df.sort_values("iPt")
    if len(df) < 3:
        # not enough points to drop first+last safely
        return df.iloc[0:0].copy()
    return df.iloc[1:-1].copy()


def plot_overlay_per_ix(df1: pd.DataFrame, df2: pd.DataFrame, outdir: str,
                        label1: str, label2: str, ncols: int = 3):
    os.makedirs(outdir, exist_ok=True)

    # Drop iPt == 11
    df1 = df1[df1["iPt"] != 11].copy()
    df2 = df2[df2["iPt"] != 11].copy()

    # Which ix values to plot (union)
    ix_vals = sorted(set(df1["ix"].unique()).union(set(df2["ix"].unique())))

    for ix in ix_vals:
        d1_ix = df1[df1["ix"] == ix]
        d2_ix = df2[df2["ix"] == ix]

        # iz values for this ix (union), skipping iz=1 and iz=8
        iz_vals = sorted(set(d1_ix["iz"].unique()).union(set(d2_ix["iz"].unique())))
        iz_vals = [iz for iz in iz_vals if iz not in (1, 8)]
        if not iz_vals:
            continue

        nplots = len(iz_vals)
        fig, axes, nrows_eff, ncols_eff = _prep_axes(nplots, ncols)

        for k, iz in enumerate(iz_vals):
            r, c = divmod(k, ncols_eff)
            ax = axes[r, c]

            s1 = d1_ix[d1_ix["iz"] == iz][["iPt", "ratio_corr"]].copy()
            s2 = d2_ix[d2_ix["iz"] == iz][["iPt", "ratio_corr"]].copy()

            s1 = drop_first_last_ipt(s1)
            s2 = drop_first_last_ipt(s2)

            if not s1.empty:
                ax.plot(s1["iPt"], s1["ratio_corr"], marker="o", linewidth=1.5, label=label1)
            if not s2.empty:
                ax.plot(s2["iPt"], s2["ratio_corr"], marker="s", linewidth=1.5, label=label2)

            ax.set_title(f"iz = {iz}")
            ax.grid(True, alpha=0.3)
            if c == 0:
                ax.set_ylabel("Acceptance")

        # Turn off any unused subplots
        for k in range(nplots, nrows_eff * ncols_eff):
            r, c = divmod(k, ncols_eff)
            axes[r, c].axis("off")

        # Common x-label
        for c in range(ncols_eff):
            axes[-1, c].set_xlabel("iPt")

        # One legend for the whole figure
        handles, labels = [], []
        for ax in fig.axes:
            h, l = ax.get_legend_handles_labels()
            for hh, ll in zip(h, l):
                if ll not in labels:
                    handles.append(hh)
                    labels.append(ll)
        if handles:
            fig.legend(handles, labels, loc="upper right")

        fig.suptitle(f"Acceptance vs iPt (iPt!=11; drop first/last iPt), ix = {ix}", fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        outpath = os.path.join(outdir, f"overlay_ix{ix:02d}.png")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)


def plot_ratio_per_ix(df1: pd.DataFrame, df2: pd.DataFrame, outdir: str,
                      label1: str, label2: str, ncols: int = 3):
    """
    For each ix, create a figure with iz panels showing:
      (ratio_corr from csv1) / (ratio_corr from csv2)
    computed point-by-point after matching on iPt.

    The first and last iPt points are dropped AFTER matching (common iPt range).
    """
    os.makedirs(outdir, exist_ok=True)

    # Drop iPt == 11
    df1 = df1[df1["iPt"] != 11].copy()
    df2 = df2[df2["iPt"] != 11].copy()

    ix_vals = sorted(set(df1["ix"].unique()).union(set(df2["ix"].unique())))

    for ix in ix_vals:
        d1_ix = df1[df1["ix"] == ix]
        d2_ix = df2[df2["ix"] == ix]

        iz_vals = sorted(set(d1_ix["iz"].unique()).union(set(d2_ix["iz"].unique())))
        iz_vals = [iz for iz in iz_vals if iz not in (1, 8)]
        if not iz_vals:
            continue

        nplots = len(iz_vals)
        fig, axes, nrows_eff, ncols_eff = _prep_axes(nplots, ncols)

        for k, iz in enumerate(iz_vals):
            r, c = divmod(k, ncols_eff)
            ax = axes[r, c]

            s1 = d1_ix[d1_ix["iz"] == iz][["iPt", "ratio_corr"]].sort_values("iPt")
            s2 = d2_ix[d2_ix["iz"] == iz][["iPt", "ratio_corr"]].sort_values("iPt")

            # Enforce fixed ratio range on every subplot
            ax.set_ylim(0.9, 1.1)

            if s1.empty or s2.empty:
                ax.text(0.5, 0.5, "missing data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"iz = {iz}")
                ax.grid(True, alpha=0.3)
                if c == 0:
                    ax.set_ylabel(f"{label1}/{label2}")
                continue

            m = pd.merge(s1, s2, on="iPt", how="inner", suffixes=("_1", "_2"))
            m = m[m["ratio_corr_2"] != 0].sort_values("iPt")

            # Drop first/last iPt on the matched set
            m = drop_first_last_ipt(m)

            if m.empty:
                ax.text(0.5, 0.5, "no overlap", ha="center", va="center", transform=ax.transAxes)
            else:
                ratio = m["ratio_corr_1"] / m["ratio_corr_2"]
                ax.plot(m["iPt"], ratio, marker="o", linewidth=1.5)
                ax.axhline(1.0, linewidth=1.0, linestyle="--")

            ax.set_title(f"iz = {iz}")
            ax.grid(True, alpha=0.3)
            if c == 0:
                ax.set_ylabel(f"{label1}/{label2}")

        for k in range(nplots, nrows_eff * ncols_eff):
            r, c = divmod(k, ncols_eff)
            axes[r, c].axis("off")

        for c in range(ncols_eff):
            axes[-1, c].set_xlabel("iPt")

        fig.suptitle(f"({label1})/({label2}) vs iPt (iPt!=11; drop first/last iPt), ix = {ix}", fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        outpath = os.path.join(outdir, f"ratio_ix{ix:02d}.png")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv1", nargs="?", default="pt_dependence_match.csv",
                    help="First pt_dependence.csv-like file")
    ap.add_argument("csv2", nargs="?", default="pt_dependence_noGoodElec.csv",
                    help="Second pt_dependence.csv-like file")

    ap.add_argument("-o", "--outdir", default="plots", help="Base output directory for PNGs")
    ap.add_argument("--label1", default="Combined", help="Legend label for first CSV")
    ap.add_argument("--label2", default="Separate", help="Legend label for second CSV")
    ap.add_argument("--ncols", type=int, default=3, help="Number of subplot columns per ix figure")
    args = ap.parse_args()

    df1 = read_pt_dependence_csv(args.csv1)
    df2 = read_pt_dependence_csv(args.csv2)

    overlay_dir = os.path.join(args.outdir, "overlay")
    ratio_dir   = os.path.join(args.outdir, "ratio_csv1_over_csv2")

    plot_overlay_per_ix(df1, df2, overlay_dir, args.label1, args.label2, ncols=args.ncols)
    plot_ratio_per_ix(df1, df2, ratio_dir, args.label1, args.label2, ncols=args.ncols)

    print(f"Saved overlay figures into: {overlay_dir}")
    print(f"Saved ratio (csv1/csv2) figures into: {ratio_dir}")


if __name__ == "__main__":
    main()
