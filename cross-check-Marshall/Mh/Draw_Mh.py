#!/usr/bin/env python3
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import uproot


SCALE_CSV_FILE = "unfolded_first16.csv"
OVERLAY_CSV_FILE = "tgp_2026_fullData.csv"

PT_EDGES = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00, 1.50], dtype=float)


def _read_table(csv_file: str):
    """Try comma-delimited then whitespace-delimited. Returns structured array with names."""
    try:
        arr = np.genfromtxt(csv_file, delimiter=",", names=True, dtype=None, encoding=None)
        if arr is None or (hasattr(arr, "size") and arr.size == 0):
            raise ValueError("empty")
        return arr
    except Exception:
        arr = np.genfromtxt(csv_file, delimiter=None, names=True, dtype=None, encoding=None)
        if arr is None or (hasattr(arr, "size") and arr.size == 0):
            raise ValueError(f"Could not parse '{csv_file}'")
        return arr


def _load_ix_scale_map(scale_csv_file: str) -> dict[int, tuple[float, float]]:
    """
    Returns: bin -> (content, acceptance)
    """
    arr = _read_table(scale_csv_file)
    required = {"bin", "content", "acceptance"}
    if not required.issubset(set(arr.dtype.names or [])):
        raise ValueError(f"'{scale_csv_file}' must contain columns: {sorted(required)}")

    rows = arr if np.ndim(arr) > 0 else np.array([arr], dtype=arr.dtype)
    m = {}
    for r in rows:
        b = int(r["bin"])
        c = float(r["content"])
        a = float(r["acceptance"])
        m[b] = (c, a)
    return m


def _load_overlay_points(overlay_csv_file: str):
    arr = _read_table(overlay_csv_file)

    needed = {"ixq2", "iz", "x", "y", "eyl", "eyh"}
    if not needed.issubset(set(arr.dtype.names or [])):
        raise ValueError(
            f"'{overlay_csv_file}' must contain columns at least: {sorted(needed)}; "
            f"found: {arr.dtype.names}"
        )
    return arr if np.ndim(arr) > 0 else np.array([arr], dtype=arr.dtype)


def _x_to_iPt(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    idx = np.searchsorted(PT_EDGES, x, side="right") - 1
    nbins = len(PT_EDGES) - 1
    valid = (idx >= 0) & (idx < nbins) & np.isfinite(x)
    iPt = idx + 1
    return iPt.astype(int), valid


def _aggregate_overlay_by_iPt(iPt_ov: np.ndarray, y_ov: np.ndarray, eyl: np.ndarray, eyh: np.ndarray):
    """
    Aggregate overlay points that land in the same iPt bin.
    Returns dict: iPt -> (mean_y, mean_err)
    Uses weighted mean with sigma = (eyl+eyh)/2 when sigma>0, else falls back to simple mean.
    """
    out = {}
    iPt_ov = np.asarray(iPt_ov, dtype=int)
    y_ov = np.asarray(y_ov, dtype=float)
    sigma = 0.5 * (np.asarray(eyl, dtype=float) + np.asarray(eyh, dtype=float))

    for ip in np.unique(iPt_ov):
        m = (iPt_ov == ip) & np.isfinite(y_ov)
        if not np.any(m):
            continue

        yy = y_ov[m]
        ss = sigma[m]

        # weights where sigma>0 and finite
        ok_w = np.isfinite(ss) & (ss > 0)
        if np.any(ok_w):
            w = 1.0 / (ss[ok_w] ** 2)
            mean = np.sum(w * yy[ok_w]) / np.sum(w)
            err = np.sqrt(1.0 / np.sum(w))
        else:
            mean = float(np.mean(yy))
            # crude fallback uncertainty: std/sqrt(N) if possible, else 0
            err = float(np.std(yy) / np.sqrt(len(yy))) if len(yy) > 1 else 0.0

        out[int(ip)] = (float(mean), float(err))

    return out


def plot_p0_vs_iPt_by_ix(
    root_file="phi_fit_results_1D.root",
    tree_name="phi_fit_p0",
    outdir="p0_vs_iPt_by_ix",
    ncols=4,
    iz_mode="all",
    bin_width=None,
    scale_csv_file=SCALE_CSV_FILE,
    overlay_csv_file=OVERLAY_CSV_FILE,
):
    if bin_width is None:
        bin_width = [
            0.005, 0.005, 0.005, 0.005,
            0.01,  0.01,  0.01,
            0.015, 0.015,
            0.02,
            0.05,
        ]
    bin_width = np.asarray(bin_width, dtype=float)

    ix_to_scale = _load_ix_scale_map(scale_csv_file)  # bin -> (content, acceptance)
    overlay = _load_overlay_points(overlay_csv_file)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ratio_dir = outdir / "ratios"
    ratio_dir.mkdir(parents=True, exist_ok=True)

    with uproot.open(root_file) as f:
        t = f[tree_name]
        a = t.arrays(["ix", "iz", "iPt", "p0"], library="np")

    ix_all = a["ix"].astype(int)
    iz_all = a["iz"].astype(int)
    ipt_all = a["iPt"].astype(int)
    p0_all = a["p0"].astype(float)

    min_ipt = int(np.min(ipt_all))
    max_ipt = int(np.max(ipt_all))
    if min_ipt < 1:
        raise ValueError(f"Found iPt={min_ipt}, but iPt is expected to start at 1.")
    if max_ipt > len(bin_width):
        raise ValueError(
            f"Found iPt up to {max_ipt}, but bin_width has only {len(bin_width)} entries."
        )

    ix_vals = np.unique(ix_all)
    ix_vals.sort()

    for ix in ix_vals:
        m_ix = (ix_all == ix)
        if not np.any(m_ix):
            continue

        if ix not in ix_to_scale:
            print(f"WARNING: skipping ix={ix}: no row with bin={ix} in {scale_csv_file}")
            continue

        content, acceptance = ix_to_scale[ix]

        if (not np.isfinite(content)) or content == 0.0:
            print(f"WARNING: skipping ix={ix}: invalid content for bin={ix} -> {content}")
            continue
        if (not np.isfinite(acceptance)):
            print(f"WARNING: skipping ix={ix}: invalid acceptance for bin={ix} -> {acceptance}")
            continue

        if iz_mode == "present":
            iz_vals = np.unique(iz_all[m_ix])
            iz_vals.sort()
        else:
            iz_max = int(np.max(iz_all[m_ix]))
            iz_vals = np.arange(1, iz_max + 1, dtype=int)

        nsub = len(iz_vals)
        if nsub == 0:
            continue

        ncols_eff = min(ncols, nsub)
        nrows_eff = math.ceil(nsub / ncols_eff)

        # --- Main figure (overlay + ROOT) ---
        fig, axes = plt.subplots(
            nrows_eff,
            ncols_eff,
            figsize=(4.2 * ncols_eff, 3.2 * nrows_eff),
            sharex=True,
        )
        axes = np.atleast_1d(axes).ravel()

        # --- Ratio figure (ROOT/CSV) ---
        figR, axesR = plt.subplots(
            nrows_eff,
            ncols_eff,
            figsize=(4.2 * ncols_eff, 3.2 * nrows_eff),
            sharex=True,
        )
        axesR = np.atleast_1d(axesR).ravel()

        for k, iz in enumerate(iz_vals):
            ax = axes[k]
            axR = axesR[k]

            m = m_ix & (iz_all == iz)

            # labels
            ax.set_title(f"iz = {iz}")
            ax.set_xlabel("iPt")
            ax.set_ylabel("p0 / (bin_width * content)")

            axR.set_title(f"iz = {iz}")
            axR.set_xlabel("iPt")
            axR.set_ylabel("ROOT / CSV")
            axR.axhline(1.0, linestyle="--", linewidth=1)

            # --- ROOT curve (scaled) ---
            root_map = {}  # iPt -> y
            if np.any(m):
                iPt = ipt_all[m]
                y = p0_all[m]

                order = np.argsort(iPt)
                iPt = iPt[order]
                y = y[order]

                y = y / bin_width[iPt - 1]
                y = y / content
                y = y*2*3.1415926
                # if you later want acceptance: y = y * acceptance

                ax.plot(iPt, y, marker="o", linestyle="-", label="5D unfolding * pi")
                root_map = {int(ip): float(yy) for ip, yy in zip(iPt, y)}
            else:
                ax.text(0.5, 0.6, "no ROOT data", ha="center", va="center", transform=ax.transAxes)
                axR.text(0.5, 0.6, "no ROOT data", ha="center", va="center", transform=axR.transAxes)

            # --- Overlay points (x -> iPt) ---
            ixq2_target = ix - 1
            iz_csv_target = iz -2
            mo = (overlay["ixq2"].astype(int) == ixq2_target) & (overlay["iz"].astype(int) == iz_csv_target)

            overlay_map = {}  # iPt -> (mean_y, mean_err)
            if np.any(mo):
                o = overlay[mo]
                x_raw = o["x"].astype(float)
                iPt_ov, valid = _x_to_iPt(x_raw)

                if np.any(~valid):
                    o = o[valid]
                    iPt_ov = iPt_ov[valid]

                if o.size > 0:
                    y_ov = o["y"].astype(float)
                    eyl = o["eyl"].astype(float)
                    eyh = o["eyh"].astype(float)

                    # plot raw overlay points
                    yerr = np.vstack([eyl, eyh]) if (np.any(eyl) or np.any(eyh)) else None
                    ax.errorbar(
                        iPt_ov, y_ov,
                        yerr=yerr,
                        fmt="s", linestyle="none",
                        label="Marshall",
                        capsize=2,
                    )

                    # build aggregated overlay per iPt for ratio
                    overlay_map = _aggregate_overlay_by_iPt(iPt_ov, y_ov, eyl, eyh)

            # --- Ratio ROOT/CSV at common iPt ---
            if root_map and overlay_map:
                common = sorted(set(root_map.keys()) & set(overlay_map.keys()))
                if len(common) > 0:
                    xr = np.array(common, dtype=int)
                    yr = np.array([root_map[ip] for ip in common], dtype=float)
                    yo = np.array([overlay_map[ip][0] for ip in common], dtype=float)
                    eo = np.array([overlay_map[ip][1] for ip in common], dtype=float)

                    good = np.isfinite(yr) & np.isfinite(yo) & (yo != 0.0)
                    xr, yr, yo, eo = xr[good], yr[good], yo[good], eo[good]

                    if xr.size > 0:
                        ratio = yr / yo
                        # propagate only overlay uncertainty (ROOT uncertainty not provided)
                        ratio_err = np.abs(ratio) * (eo / np.abs(yo))

                        axR.errorbar(
                            xr, ratio,
                            yerr=ratio_err,
                            fmt="o", linestyle="none",
                            capsize=2,
                            label="(ROOT/CSV)",
                        )
                        axR.legend(fontsize=8)
                    else:
                        axR.text(0.5, 0.6, "no valid overlap", ha="center", va="center", transform=axR.transAxes)
                else:
                    axR.text(0.5, 0.6, "no iPt overlap", ha="center", va="center", transform=axR.transAxes)
            else:
                axR.text(0.5, 0.6, "missing series", ha="center", va="center", transform=axR.transAxes)

            # legend on main plot
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8)

        # turn off unused subplots
        for k in range(nsub, len(axes)):
            axes[k].axis("off")
            axesR[k].axis("off")

        # main fig save
        fig.suptitle(f"Scaled p0 vs iPt (ix={ix}, content={content:.3g}, acc={acceptance:.3g}) with overlay")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        outpath = outdir / f"p0_vs_iPt_ix{ix:02d}.png"
        fig.savefig(outpath, dpi=200)
        plt.close(fig)
        print(f"wrote: {outpath}")

        # ratio fig save
        figR.suptitle(f"Ratio ROOT/CSV vs iPt (ix={ix})")
        figR.tight_layout(rect=(0, 0, 1, 0.95))
        outpathR = ratio_dir / f"ratio_root_over_csv_ix{ix:02d}.png"
        figR.savefig(outpathR, dpi=200)
        plt.close(figR)
        print(f"wrote: {outpathR}")


if __name__ == "__main__":
    plot_p0_vs_iPt_by_ix(
        root_file="phi_fit_results_1D.root",
        tree_name="phi_fit_p0",
        outdir="p0_vs_iPt_by_ix",
        ncols=4,
        iz_mode="all",
        scale_csv_file=SCALE_CSV_FILE,
        overlay_csv_file=OVERLAY_CSV_FILE,
    )
