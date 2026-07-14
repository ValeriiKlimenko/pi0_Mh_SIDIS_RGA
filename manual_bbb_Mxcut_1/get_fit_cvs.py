#!/usr/bin/env python3
import os
import argparse
import math
import ROOT

# ---------------- binning (must match your analysis!) ----------------
N_ZBINS = 8
N_PT2_BINS_WITH_OVERFLOW = 11  # 10 + overflow
N_PHI = 8

Z_EDGES  = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
PT2_EDGES= [0, 0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1.0,1.5]


def phi_center_from_bin(phi_bin_1b: int) -> float:
    """Return the center of a 1-based phi bin in degrees."""
    w = 360.0 / N_PHI
    return (phi_bin_1b - 0.5) * w


def bin_label(edges, i1b: int) -> str:
    """Pretty label 'a–b' for 1-based bin index."""
    if i1b < 1 or i1b >= len(edges):
        return "out"
    a, b = edges[i1b-1], edges[i1b]
    return f"{a:.2f}–{b:.2f}"


def get_phi_points_for_cell(h2: ROOT.TH2,
                            ybin: int,
                            z_bin: int,
                            pt2_bin: int):
    """
    Extract phi points for a given (z_bin, pt2_bin, ybin) from a TH2 with:
      X = composite index z_pt2_phi_bin (1..N_ZBINS*N_PT2*N_PHI, plus +1)
      Y = xq2bin (0..20)
    Returns list of (phi_deg, y, yerr).
    """
    pts = []
    for phi_bin in range(1, N_PHI+1):
        # zpt2_bin: 1..(N_ZBINS*N_PT2_BINS_WITH_OVERFLOW)
        zpt2_bin = (z_bin - 1) * N_PT2_BINS_WITH_OVERFLOW + pt2_bin
        comp = (zpt2_bin - 1) * N_PHI + phi_bin  # 1..N_ZBINS*N_PT2*N_PHI

        xcoord = float(comp)
        bx = h2.GetXaxis().FindBin(xcoord)
        if bx < 1 or bx > h2.GetNbinsX():
            continue

        y = h2.GetBinContent(bx, ybin)
        e = h2.GetBinError(bx, ybin)
        phi_deg = phi_center_from_bin(phi_bin)
        pts.append((phi_deg, y, e))
    return pts


def fit_phi_points(pts,
                   empty_eps: float = 0.0,
                   min_points: int = 6):
    """
    Fit a list of (phi_deg, y, yerr) with:
      f(phi) = p0 + p1*cos(phi) + p2*cos(2*phi)
    phi in degrees.

    Returns dict:
      {
        "success": bool,
        "npts": int,
        "p0": float,
        "p1": float,
        "p2": float,
        "chi2": float,
        "ndf": int,
        "status": int
      }
    """
    result = {
        "success": False,
        "npts": 0,
        "p0": float("nan"),
        "p1": float("nan"),
        "p2": float("nan"),
        "chi2": float("nan"),
        "ndf": -1,
        "status": -1,
    }

    # keep only non-empty points
    filtered = [(phi, y, e) for (phi, y, e) in pts
                if not (abs(y) <= empty_eps and abs(e) <= empty_eps)]
    if len(filtered) < min_points:
        return result

    # histogram average over non-empty bins
    vals = [y for (_, y, _) in filtered]
    if not vals:
        return result
    histAvg = sum(vals) / len(vals)

    # Simple Fourier-like initial guesses for p1,p2
    sumC1 = sumC2 = 0.0
    sumC1C1 = sumC2C2 = 0.0

    for (phi_deg, v, e) in filtered:
        phi = phi_deg * math.pi / 180.0
        w = 1.0 / (e * e) if (abs(e) > empty_eps) else 1.0

        c1 = math.cos(phi)
        c2 = math.cos(2.0 * phi)
        y  = v - histAvg

        sumC1   += w * y * c1
        sumC2   += w * y * c2
        sumC1C1 += w * c1 * c1
        sumC2C2 += w * c2 * c2

    p0_init = max(1e-12, max(0.0, histAvg))
    p1_init = sumC1 / sumC1C1 if sumC1C1 > 0 else 0.0
    p2_init = sumC2 / sumC2C2 if sumC2C2 > 0 else 0.0

    # Build TGraphErrors
    npts = len(filtered)
    gr = ROOT.TGraphErrors(npts)
    for i, (phi_deg, v, e) in enumerate(filtered):
        gr.SetPoint(i, phi_deg, v)
        gr.SetPointError(i, 0.0, e if e > empty_eps else 1.0)

    # Define fit function: phi in degrees
    # Use numeric pi to avoid TMath in formula string
    pi_val = math.pi
    fphi = ROOT.TF1("fphi_tmp",
                    f"[0] + [1]*cos(x*{pi_val}/180.) + "
                    f"[2]*cos(2.*x*{pi_val}/180.)",
                    0.0, 360.0)

    fphi.SetParameters(p0_init, p1_init, p2_init)
    scale = max(histAvg, 1.0)
    fphi.SetParLimits(0, 0.0,        20.0 * scale)
    fphi.SetParLimits(1, -10.0*scale, 10.0 * scale)
    fphi.SetParLimits(2, -10.0*scale, 10.0 * scale)

    # Fit quietly and store result ("QS": Quiet, Store)
    fitres = gr.Fit(fphi, "QS")
    # fitres is a TFitResultPtr
    try:
        status = int(fitres)
        chi2   = float(fitres.Chi2())
        ndf    = int(fitres.Ndf())
    except Exception:
        # fallback if PyROOT semantics differ
        status, chi2, ndf = -1, float("nan"), -1

    p0 = fphi.GetParameter(0)
    p1 = fphi.GetParameter(1)
    p2 = fphi.GetParameter(2)

    result.update({
        "success": True,
        "npts": npts,
        "p0": p0,
        "p1": p1,
        "p2": p2,
        "chi2": chi2,
        "ndf": ndf,
        "status": status,
    })
    return result


def main():
    ap = argparse.ArgumentParser(
        description="Fit φ-dependence in two TH2 histograms "
                    "and save p0,p1,p2 for each (xQ²,z,pT²) cell to CSV."
    )
    ap.add_argument("--root1", required=False,
                    default="../manual_bbb_Mxcut_1/measData_times_truth_over_meas.root",
                    help="First ROOT file")
    ap.add_argument("--root2", required=False,
                    default="../measData_times_truth_over_meas.root",
                    help="Second ROOT file")
    ap.add_argument("--hname", default="h_measData_times_truth_over_meas",
                    help="Name of TH2 in both ROOT files")
    ap.add_argument("--outcsv", default="phi_fit_twofiles.csv",
                    help="Output CSV file")
    ap.add_argument("--max_ix", type=int, default=None,
                    help="Optional maximum xQ² bin (center value) to fit")
    ap.add_argument("--label1", default="file1",
                    help="Label for first file (goes in 'source' column)")
    ap.add_argument("--label2", default="file2",
                    help="Label for second file (goes in 'source' column)")
    ap.add_argument("--min_points", type=int, default=6,
                    help="Minimum number of non-empty φ points to attempt a fit")
    ap.add_argument("--empty_eps", type=float, default=0.0,
                    help="Threshold below which y and err are treated as empty")
    args = ap.parse_args()

    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    # Open files and histograms
    f1 = ROOT.TFile.Open(args.root1, "READ")
    if not f1 or f1.IsZombie():
        raise RuntimeError(f"Cannot open file {args.root1}")

    f2 = ROOT.TFile.Open(args.root2, "READ")
    if not f2 or f2.IsZombie():
        raise RuntimeError(f"Cannot open file {args.root2}")

    h1 = f1.Get(args.hname)
    if not h1 or not isinstance(h1, ROOT.TH2):
        raise RuntimeError(f"{args.hname} not found as TH2 in {args.root1}")

    h2 = f2.Get(args.hname)
    if not h2 or not isinstance(h2, ROOT.TH2):
        raise RuntimeError(f"{args.hname} not found as TH2 in {args.root2}")

    # Basic sanity check: same binning
    if h1.GetNbinsX() != h2.GetNbinsX() or h1.GetNbinsY() != h2.GetNbinsY():
        raise RuntimeError("Histograms have different binning in X or Y")

    ny = h1.GetNbinsY()
    yaxis = h1.GetYaxis()

    # Prepare CSV
    outdir = os.path.dirname(os.path.abspath(args.outcsv))
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    with open(args.outcsv, "w") as csv:
        csv.write(
            "xq2_center,xq2_bin,iy,"
            "z_bin,pt2_bin,z_label,pt2_label,"
            "rootfile,source,"
            "npoints,p0,p1,p2,chi2,ndf,status\n"
        )

        # Loop over all xQ² (Y) bins
        for iy in range(1, ny + 1):
            xq2_center = yaxis.GetBinCenter(iy)
            xq2_int = int(round(xq2_center))

            if args.max_ix is not None and xq2_int > args.max_ix:
                continue

            # Optionally skip bins with no content at all in both hists
            proj1 = h1.ProjectionX("_px_tmp1", iy, iy)
            proj2 = h2.ProjectionX("_px_tmp2", iy, iy)
            if proj1.Integral() == 0.0 and proj2.Integral() == 0.0:
                continue

            for iz in range(1, N_ZBINS + 1):
                for ipt in range(1, N_PT2_BINS_WITH_OVERFLOW + 1):
                    # Extract φ distributions for this (z, pT², xQ²)
                    pts1 = get_phi_points_for_cell(h1, iy, iz, ipt)
                    pts2 = get_phi_points_for_cell(h2, iy, iz, ipt)

                    # Fit each file separately
                    fit1 = fit_phi_points(
                        pts1,
                        empty_eps=args.empty_eps,
                        min_points=args.min_points
                    ) if pts1 else {"success": False}

                    fit2 = fit_phi_points(
                        pts2,
                        empty_eps=args.empty_eps,
                        min_points=args.min_points
                    ) if pts2 else {"success": False}

                    zlab  = bin_label(Z_EDGES,  iz)
                    ptlab = bin_label(PT2_EDGES, ipt)

                    # Write CSV rows (only for successful fits)
                    if fit1.get("success"):
                        csv.write(
                            "{:.6f},{:d},{:d},"
                            "{:d},{:d},{},{},"
                            "{},{},"
                            "{:d},{:.10g},{:.10g},{:.10g},{:.10g},{:d},{:d}\n".format(
                                xq2_center, xq2_int, iy,
                                iz, ipt, zlab, ptlab,
                                os.path.basename(args.root1), args.label1,
                                fit1["npts"],
                                fit1["p0"], fit1["p1"], fit1["p2"],
                                fit1["chi2"], fit1["ndf"], fit1["status"]
                            )
                        )

                    if fit2.get("success"):
                        csv.write(
                            "{:.6f},{:d},{:d},"
                            "{:d},{:d},{},{},"
                            "{},{},"
                            "{:d},{:.10g},{:.10g},{:.10g},{:.10g},{:d},{:d}\n".format(
                                xq2_center, xq2_int, iy,
                                iz, ipt, zlab, ptlab,
                                os.path.basename(args.root2), args.label2,
                                fit2["npts"],
                                fit2["p0"], fit2["p1"], fit2["p2"],
                                fit2["chi2"], fit2["ndf"], fit2["status"]
                            )
                        )

            print(f"Finished fits for xQ^2 bin (center) = {xq2_int}")

    f1.Close()
    f2.Close()
    print(f"All fit results written to: {args.outcsv}")


if __name__ == "__main__":
    main()
