#!/usr/bin/env python3
"""
For EACH (xQ2 bin, z bin) make ONE figure with subplots for pT^2 bins,
and in each subplot draw closure ratio (unfolded/truth) as a function of phi.

Binning decode (matches your flattening):
  global_bin = xq2bin * nZ + z_pt2_phi_bin
  z_pt2_phi_bin = ( (zbin*nPt2 + pt2bin) * nPhi ) + phi   , with phi = 1..nPhi

Cuts / drops:
  - drop first and last z bins:    zbin = 1..(N_Zbins-2)
  - drop first and last pT^2 bins: pt2bin = 1..(N_pT2-2)
  - require U>0 and T>0 for a point to be drawn
  - optional: skip points with ratio < min_ratio (default disabled unless set)

Usage:
  python3 draw_ratio_phi_by_pt2.py --rootfile closure_test_bayes.root -o ratio_phi
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ---------------- defaults from your zpt2phit_8x8x9_binning.cxx ----------------
DEFAULT_Z_EDGES   = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0], dtype=float)  # 8 bins
DEFAULT_PT2_EDGES = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00, 1.50], dtype=float)  # 11 bins
DEFAULT_NPHI = 8
DEFAULT_NX   = 21  # xq2bin: 0..20 (0 underflow; 1..20 physical)
# -----------------------------------------------------------------------------


def make_ratio_hist_per_x(
    true_vals: np.ndarray,
    unf_vals: np.ndarray,
    x_bins: list[int],
    z_bins_to_plot: list[int],
    pt2_bins_to_plot: list[int],
    nPhi: int,
    nPt2: int,
    nZ: int,
    nTot: int,
    outdir: Path,
    min_ratio: float | None = None,
    ratio_band: float = 0.10,
):
    """
    For each xQ2 bin:
      - compute ratio U/T separately for every allowed (z, pt2, phi) bin
      - fill all ratio values into one 1D histogram
    """

    def zphi_bin(iz: int, ipt: int, iphi_1based: int) -> int:
        return (iz * nPt2 + ipt) * nPhi + iphi_1based

    for x in x_bins:
        ratios = []

        for iz in z_bins_to_plot:
            for ipt in pt2_bins_to_plot:
                for iphi in range(1, nPhi + 1):
                    zphi = zphi_bin(iz, ipt, iphi)
                    g = global_index(x, zphi, nZ)

                    if not (0 <= g < nTot):
                        continue

                    U = unf_vals[g]
                    T = true_vals[g]

                    if (U > 0) and (T > 0):
                        r = U / T
                        if np.isfinite(r):
                            if min_ratio is None or r >= min_ratio:
                                ratios.append(r)

        ratios = np.asarray(ratios, dtype=float)

        fig, ax = plt.subplots(figsize=(8, 5.5))
        ax.grid(True, alpha=0.3)

        if ratios.size == 0:
            ax.text(0.5, 0.5, "No valid ratio entries", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"xQ2 bin {x}: ratio distribution")
            ax.set_xlabel("Unfolded / Truth")
            ax.set_ylabel("Count")
        else:
            # Choose a reasonable binning around the observed spread
            rmin = np.nanmin(ratios)
            rmax = np.nanmax(ratios)

            # If all ratios are the same (or almost the same), force a wide central bin
            if np.allclose(ratios, ratios[0], rtol=0.0, atol=1e-12):
                center = ratios[0]
                half_width = 0.01   # gives one visible pillar from center-0.05 to center+0.05
                bins = np.array([center - half_width, center + half_width])

            elif np.isfinite(rmin) and np.isfinite(rmax) and (rmax - rmin) < 1e-3:
                center = 0.5 * (rmin + rmax)
                half_width = 0.01
                bins = np.linspace(center - half_width, center + half_width, 9)

            else:
                pad = 0.08 * (rmax - rmin)
                nbins = min(40, max(15, int(np.sqrt(ratios.size) * 2.5)))
                bins = np.linspace(rmin - pad, rmax + pad, nbins + 1)

            ax.hist(ratios, bins=bins, rwidth=0.9)

            ax.axvline(1.0, linestyle="--", linewidth=1.5)
            if ratio_band and ratio_band > 0:
                ax.axvline(1.0 + ratio_band, linestyle=":", linewidth=1.0)
                ax.axvline(1.0 - ratio_band, linestyle=":", linewidth=1.0)

            text = "\n".join([
                f"N = {ratios.size}",
                f"mean = {np.mean(ratios):.4f}",
                f"std = {np.std(ratios, ddof=1):.4f}" if ratios.size > 1 else "std = n/a",
                f"RMS = {np.sqrt(np.mean(ratios**2)):.4f}",
            ])
            ax.text(
                0.98,
                0.97,
                text,
                ha="right",
                va="top",
                transform=ax.transAxes,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
            )

            ax.set_title(f"xQ2 bin {x}: distribution of per-bin closure ratios")
            ax.set_xlabel("Unfolded / Truth")
            ax.set_ylabel("Count")

        fig.tight_layout()
        fig.savefig(outdir / f"closure_ratio_hist_xQ2bin{x:02d}.png", dpi=200)
        fig.savefig(outdir / f"closure_ratio_hist_xQ2bin{x:02d}.pdf")
        plt.close(fig)

def make_ratio_hist_all_x_1to16(
    true_vals: np.ndarray,
    unf_vals: np.ndarray,
    z_bins_to_plot: list[int],
    pt2_bins_to_plot: list[int],
    nPhi: int,
    nPt2: int,
    nZ: int,
    nTot: int,
    outdir: Path,
    min_ratio: float | None = None,
    ratio_band: float = 0.10,
):
    """
    Build one 1D histogram of per-bin ratios U/T, combining all allowed
    (z, pt2, phi) bins from xQ2 bins 1..16 inclusive.
    """

    def zphi_bin(iz: int, ipt: int, iphi_1based: int) -> int:
        return (iz * nPt2 + ipt) * nPhi + iphi_1based

    ratios = []

    for x in range(1, 17):   # include both 1 and 16
        for iz in z_bins_to_plot:
            for ipt in pt2_bins_to_plot:
                for iphi in range(1, nPhi + 1):
                    zphi = zphi_bin(iz, ipt, iphi)
                    g = global_index(x, zphi, nZ)

                    if not (0 <= g < nTot):
                        continue

                    U = unf_vals[g]
                    T = true_vals[g]

                    if (U > 0) and (T > 0):
                        r = U / T
                        if np.isfinite(r):
                            if min_ratio is None or r >= min_ratio:
                                ratios.append(r)

    ratios = np.asarray(ratios, dtype=float)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.grid(True, alpha=0.3)

    if ratios.size == 0:
        ax.text(0.5, 0.5, "No valid ratio entries", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Combined xQ2 bins 1-16: ratio distribution")
        ax.set_xlabel("Unfolded / Truth")
        ax.set_ylabel("Count")
    else:
        rmin = np.nanmin(ratios)
        rmax = np.nanmax(ratios)

        # Special handling when everything is ~1 so the bar is visible
        if np.allclose(ratios, ratios[0], rtol=0.0, atol=1e-12):
            center = ratios[0]
            half_width = 0.01
            bins = np.array([center - half_width, center + half_width])

        elif np.isfinite(rmin) and np.isfinite(rmax) and (rmax - rmin) < 1e-3:
            center = 0.5 * (rmin + rmax)
            half_width = 0.01
            bins = np.linspace(center - half_width, center + half_width, 9)

        else:
            pad = 0.08 * (rmax - rmin)
            nbins = min(40, max(15, int(np.sqrt(ratios.size) * 2.5)))
            bins = np.linspace(rmin - pad, rmax + pad, nbins + 1)

        ax.hist(ratios, bins=bins, rwidth=0.9)

        ax.axvline(1.0, linestyle="--", linewidth=1.5)
        if ratio_band and ratio_band > 0:
            ax.axvline(1.0 + ratio_band, linestyle=":", linewidth=1.0)
            ax.axvline(1.0 - ratio_band, linestyle=":", linewidth=1.0)

        text = "\n".join([
            f"N = {ratios.size}",
            f"mean = {np.mean(ratios):.4f}",
            f"std = {np.std(ratios, ddof=1):.4f}" if ratios.size > 1 else "std = n/a",
            f"RMS = {np.sqrt(np.mean(ratios**2)):.4f}",
        ])
        ax.text(
            0.98,
            0.97,
            text,
            ha="right",
            va="top",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
        )

        ax.set_title("Distribution of per-bin closure ratios for xQ2 bins 1-16")
        ax.set_xlabel("Unfolded / Truth")
        ax.set_ylabel("Count")

    fig.tight_layout()
    fig.savefig(outdir / "closure_ratio_hist_xQ2bin01to16.png", dpi=200)
    fig.savefig(outdir / "closure_ratio_hist_xQ2bin01to16.pdf")
    plt.close(fig)



def _find_hist(upfile, names):
    for n in names:
        if n in upfile:
            return upfile[n], n
    raise KeyError(f"None of these histograms found: {names}")

def _read_th1_as_arrays(upfile, names):
    h, used = _find_hist(upfile, names)
    vals = h.values(flow=False)  # excludes under/overflow
    vars_ = h.variances(flow=False)
    errs = np.zeros_like(vals, dtype=float) if vars_ is None else np.sqrt(np.maximum(vars_, 0.0))
    return vals.astype(float), errs.astype(float), used

def global_index(xbin: int, z_pt2_phi_bin: int, nZ: int) -> int:
    return xbin * nZ + z_pt2_phi_bin

def subplot_grid(n: int):
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))
    return nrows, ncols

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rootfile", default="closure_test_bayes.root",
                    help="closure_test_bayes.root or qa_response.root")
    ap.add_argument("-o", "--outdir", default="ratio_phi_plots", help="output directory")
    ap.add_argument("--ratio-band", type=float, default=0.10,
                    help="draw +/- band around 1 (default 0.10)")
    ap.add_argument("--min-ratio", type=float, default=None,
                    help="if set, skip points with ratio < min-ratio (e.g. 0.6). Default: no cut.")
    ap.add_argument("--include-underflow-x", action="store_true",
                    help="include xQ2 bin 0 (underflow)")
    ap.add_argument("--nx", type=int, default=DEFAULT_NX,
                    help=f"number of xQ2 bins including underflow (default {DEFAULT_NX})")
    ap.add_argument("--nphi", type=int, default=DEFAULT_NPHI,
                    help=f"number of phi bins (default {DEFAULT_NPHI})")
    ap.add_argument("--x-bins", default="all",
                    help="x bins to plot: 'all' or comma list like '1,2,3' (0/1-based as you use in code)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    z_edges   = DEFAULT_Z_EDGES
    pt2_edges = DEFAULT_PT2_EDGES
    nZbins = len(z_edges) - 1
    nPt2   = len(pt2_edges) - 1
    nPhi   = args.nphi

    # drop first and last z, pt2 bins
    z_bins_to_plot   = list(range(1, nZbins - 1))   # 1..N-2
    pt2_bins_to_plot = list(range(1, nPt2   - 1))   # 1..N-2

    # phi geometry (0..360)
    phi_edges   = np.linspace(0.0, 360.0, nPhi + 1)
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])

    # nZ in unfolding code: N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1
    nZ = nZbins * nPt2 * nPhi + 1

    try:
        import uproot
    except ImportError as e:
        raise SystemExit("ERROR: need uproot (pip install uproot) or rewrite using PyROOT.") from e

    with uproot.open(args.rootfile) as f:
        true_vals, true_errs, true_name = _read_th1_as_arrays(f, ["h_true_mc", "h_truth_trained"])
        unf_vals,  unf_errs,  unf_name  = _read_th1_as_arrays(f, ["h_unfold_closure"])

    nTot = len(true_vals)
    if len(unf_vals) != nTot:
        raise SystemExit(f"ERROR: histogram lengths differ: truth={len(true_vals)} unfolded={len(unf_vals)}")

    expected = args.nx * nZ
    if nTot != expected:
        if nTot % nZ == 0:
            nx = nTot // nZ
            print(f"WARNING: nbins={nTot} != nx*nZ={expected}. Inferring nx={nx} from nbins/nZ.")
        else:
            raise SystemExit(f"ERROR: nbins={nTot} not compatible with nZ={nZ}. Check nx/nphi/pt2/z settings.")
    else:
        nx = args.nx

    # x bins selection
    if args.x_bins.strip().lower() == "all":
        x_start = 0 if args.include_underflow_x else 1
        x_bins = list(range(x_start, nx))
    else:
        x_bins = [int(x.strip()) for x in args.x_bins.split(",") if x.strip()]
        if not args.include_underflow_x:
            x_bins = [x for x in x_bins if x != 0]

    print(f"Using histograms: truth='{true_name}', unfolded='{unf_name}'")
    print(f"Decoded geometry: nx={nx}, nZ={nZ} (z={nZbins}, pt2={nPt2}, phi={nPhi})")
    print(f"Plot z bins (0-based): {z_bins_to_plot}")
    print(f"Plot pt2 bins (0-based): {pt2_bins_to_plot}")
    print(f"Plot x bins: {x_bins}")
    print(f"Writing plots to: {outdir}")

    nrows, ncols = subplot_grid(len(pt2_bins_to_plot))

    # helper: compute z_pt2_phi_bin (1-based z_pt2_phi_bin, phi is 1..nPhi)
    # z_pt2_phi_bin = ( (zbin*nPt2 + pt2bin) * nPhi ) + phi
    def zphi_bin(iz: int, ipt: int, iphi_1based: int) -> int:
        return (iz * nPt2 + ipt) * nPhi + iphi_1based

    for x in x_bins:
        for iz in z_bins_to_plot:
            fig, axs = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.8 * nrows),
                                    sharex=True, sharey=True)
            axs = np.atleast_1d(axs).ravel()

            zlo, zhi = z_edges[iz], z_edges[iz + 1]

            for k, ipt in enumerate(pt2_bins_to_plot):
                ax = axs[k]
                ax.grid(True, alpha=0.3)
                ax.axhline(1.0, linestyle="--", linewidth=1)

                band = args.ratio_band
                if band and band > 0:
                    ax.axhline(1.0 + band, linestyle=":", linewidth=1)
                    ax.axhline(1.0 - band, linestyle=":", linewidth=1)

                # build arrays over phi
                U = np.zeros(nPhi, dtype=float)
                T = np.zeros(nPhi, dtype=float)
                eU = np.zeros(nPhi, dtype=float)
                eT = np.zeros(nPhi, dtype=float)

                for iphi in range(1, nPhi + 1):
                    zphi = zphi_bin(iz, ipt, iphi)
                    g = global_index(x, zphi, nZ)
                    if 0 <= g < nTot:
                        U[iphi - 1] = unf_vals[g]
                        T[iphi - 1] = true_vals[g]
                        eU[iphi - 1] = unf_errs[g]
                        eT[iphi - 1] = true_errs[g]

                good = (U > 0) & (T > 0)
                ratio = np.full(nPhi, np.nan, dtype=float)
                eratio = np.full(nPhi, np.nan, dtype=float)

                ratio[good] = U[good] / T[good]
                eratio[good] = ratio[good] * np.sqrt((eU[good] / U[good]) ** 2 + (eT[good] / T[good]) ** 2)

                keep = good & np.isfinite(ratio)
                if args.min_ratio is not None:
                    keep &= (ratio >= args.min_ratio)

                if np.any(keep):
                    ax.errorbar(phi_centers[keep], ratio[keep], yerr=eratio[keep],
                                fmt="o", markersize=3, capsize=2)

                ptlo, pthi = pt2_edges[ipt], pt2_edges[ipt + 1]
                ax.set_title(rf"$p_T^2 \in [{ptlo:g}, {pthi:g}]$")

                if (k % ncols) == 0:
                    ax.set_ylabel("Unfolded / Truth")
                if k >= (len(pt2_bins_to_plot) - ncols):
                    ax.set_xlabel(r"$\phi$ [deg]")

            # hide unused axes
            for k in range(len(pt2_bins_to_plot), len(axs)):
                axs[k].axis("off")

            fig.suptitle(f"Closure ratio vs $\\phi$  (xQ2 bin {x}, z ∈ [{zlo:g}, {zhi:g}])", y=0.995)
            fig.tight_layout(rect=[0, 0, 1, 0.985])

            base = f"closure_ratio_phi_xQ2bin{x:02d}_z{iz:02d}"
            fig.savefig(outdir / f"{base}.png", dpi=200)
            fig.savefig(outdir / f"{base}.pdf")
            plt.close(fig)

    make_ratio_hist_per_x(
        true_vals=true_vals,
        unf_vals=unf_vals,
        x_bins=x_bins,
        z_bins_to_plot=z_bins_to_plot,
        pt2_bins_to_plot=pt2_bins_to_plot,
        nPhi=nPhi,
        nPt2=nPt2,
        nZ=nZ,
        nTot=nTot,
        outdir=outdir,
        min_ratio=args.min_ratio,
        ratio_band=args.ratio_band,
    ) 

    make_ratio_hist_all_x_1to16(
        true_vals=true_vals,
        unf_vals=unf_vals,
        z_bins_to_plot=z_bins_to_plot,
        pt2_bins_to_plot=pt2_bins_to_plot,
        nPhi=nPhi,
        nPt2=nPt2,
        nZ=nZ,
        nTot=nTot,
        outdir=outdir,
        min_ratio=args.min_ratio,
        ratio_band=args.ratio_band,
    )
  

    print("Done.")

if __name__ == "__main__":
    main()
