#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sum_slice_fit.py
----------------
Sum TH2 histograms across multiple ROOT files, then slice along X and fit each
Y-projection (pi0 mass) with a Gaussian signal + polynomial background.

Defaults tailored for pi0:
  - expected mean ~ 0.133
  - sigma constrained to [0.003, 0.03]
  - mean constrained to [mean - mean_window, mean + mean_window] (default window = 0.010)

Outputs:
  - ROOT file with summed TH2s and per-slice projection hists + TF1 fits
  - CSV file with per-slice fit parameters (amplitude, mean, sigma, yield, chi2/ndf, background coeffs)
  - ROOT TTree ("fit_results") mirroring the CSV
  - PNGs: each fitted slice saved under --png-dir/<hist_name>/<hist_name>_py_<xbin>_fit.png

Usage (example):
  python3 sum_slice_fit.py \
    --inputs pi0_data_v3_unfold_0.root pi0_data_v3_unfold_50.root pi0_data_v3_unfold_100.root pi0_data_v3_unfold_150.root \
    --poly 1 \
    --mean 0.133 --mean-window 0.010 --sigma-min 0.003 --sigma-max 0.03 \
    --out pi0_unfold_sum_fits.root \
    --csv pi0_unfold_sum_fits.csv \
    --png-dir fit_png
"""
import os
import sys
import math
import csv
import argparse

try:
    import ROOT
except Exception as e:
    sys.stderr.write("ERROR: This script requires PyROOT (ROOT with Python bindings).\n")
    raise

# Batch mode for headless plotting/saving
ROOT.gROOT.SetBatch(True)

# Make fits stable
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2", "Migrad")
ROOT.Math.MinimizerOptions.SetDefaultStrategy(1)

def is_th2(obj):
    try:
        return bool(obj.InheritsFrom("TH2"))
    except Exception:
        return False

def list_target_hists(root_file, prefix="d9_"):
    names = []
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if obj and is_th2(obj) and obj.GetName().startswith(prefix):
            names.append(obj.GetName())
    names.sort(key=lambda n: int(n.split("_")[1]) if "_" in n and n.split("_")[1].isdigit() else n)
    return names

def clone_with_dir(h, newname):
    clone = h.Clone(newname)
    clone.SetDirectory(0)
    return clone

def sum_histograms(files, hist_names):
    sums = {}
    for name in hist_names:
        summed = None
        for f in files:
            h = f.Get(name)
            if not h or not is_th2(h):
                continue
            if summed is None:
                summed = clone_with_dir(h, name)
            else:
                summed.Add(h)
        if summed is not None:
            sums[name] = summed
    return sums

def initial_guesses(h1, expected_mean, sigma_min, sigma_max):
    """Return (amp, mean, sigma) initial guesses from a TH1 but anchored to expected mean/sigma bounds."""
    if h1.GetEntries() <= 0:
        return (1.0, expected_mean, max(0.5*(sigma_min+sigma_max), 1e-6))
    max_bin = h1.GetMaximumBin()
    amp0 = max(h1.GetMaximum(), 1.0)
    # Use expected mean, but keep within axis range
    mean0 = expected_mean
    ax = h1.GetXaxis()
    mean0 = min(max(mean0, ax.GetXmin()+1e-6), ax.GetXmax()-1e-6)
    # Start sigma in the middle of the allowed band, or from RMS but clipped to [sigma_min, sigma_max]
    rms = h1.GetRMS()
    if rms <= 0:
        sig0 = 0.5*(sigma_min + sigma_max)
    else:
        sig0 = max(min(rms*0.6, sigma_max), sigma_min)
    return (amp0, mean0, sig0)

def build_model(pol_order, ylow, yhigh, expected_mean, mean_window, sigma_min, sigma_max):
    """
    Create TF1: gaus(0) + polN(3) and impose constraints:
      mean in [expected_mean - mean_window, expected_mean + mean_window]
      sigma in [sigma_min, sigma_max]
    """
    if pol_order < 0 or pol_order > 5:
        raise ValueError("pol_order must be between 0 and 5")
    expr = f"gaus(0)+pol{pol_order}(3)" if pol_order >= 0 else "gaus(0)"
    f = ROOT.TF1("fitfunc", expr, ylow, yhigh)
    pnames = ["A", "mean", "sigma"] + [f"b{i}" for i in range(pol_order + 1)]

    # Constrain mean and sigma
    f.SetParLimits(1, expected_mean - mean_window, expected_mean + mean_window)  # mean
    f.SetParLimits(2, sigma_min, sigma_max)                                      # sigma

    return f, pnames

def fit_projection(hproj, pol_order, expected_mean, mean_window, sigma_min, sigma_max, opts="RSQ"):
    """Fit a Y-projection with gaus + polN background under the requested constraints."""
    ylow = hproj.GetXaxis().GetXmin()
    yhigh = hproj.GetXaxis().GetXmax()
    func, pnames = build_model(pol_order, ylow, yhigh, expected_mean, mean_window, sigma_min, sigma_max)

    # Initial guesses
    A0, mu0, sig0 = initial_guesses(hproj, expected_mean, sigma_min, sigma_max)
    func.SetParameter(0, A0)
    func.SetParameter(1, mu0)
    func.SetParameter(2, sig0)
    for i in range(pol_order + 1):
        func.SetParameter(3 + i, 0.0)

    # Fit (R=use range, S=store, Q=quiet)
    res = hproj.Fit(func, opts)
    status = int(res)
    chi2 = res.Chi2() if hasattr(res, "Chi2") else func.GetChisquare()
    ndf  = res.Ndf()  if hasattr(res, "Ndf")  else func.GetNDF()

    pars = [func.GetParameter(i) for i in range(func.GetNpar())]
    errs = [func.GetParError(i) for i in range(func.GetNpar())]

    A, mean, sigma = pars[0], pars[1], max(pars[2], 1e-12)
    Aerr, meanerr, sigmaerr = errs[0], errs[1], errs[2]

    # Yield = A * sigma * sqrt(2*pi)
    yield_val = A * sigma * math.sqrt(2.0 * math.pi)
    try:
        cov = res.GetCovarianceMatrix()
        try:
            cov_A_sigma = cov[0][2]
        except Exception:
            cov_A_sigma = cov(0, 2)
        var_prod = (sigma**2) * (Aerr**2) + (A**2) * (sigmaerr**2) + 2.0 * A * sigma * cov_A_sigma
        var_prod = max(var_prod, 0.0)
        yield_err = math.sqrt(2.0 * math.pi) * math.sqrt(var_prod)
    except Exception:
        yield_err = math.sqrt((sigma * Aerr)**2 + (A * sigmaerr)**2) * math.sqrt(2.0 * math.pi)

    return {
        "status": status,
        "chi2": chi2,
        "ndf": ndf,
        "A": A, "A_err": Aerr,
        "mean": mean, "mean_err": meanerr,
        "sigma": sigma, "sigma_err": sigmaerr,
        "yield": yield_val, "yield_err": yield_err,
        "pars": pars, "errs": errs,
        "pnames": pnames,
        "func": func,
        "result": res,
    }

def parse_d9_index(name):
    try:
        return int(name.split("_")[1])
    except Exception:
        return -1

def sanitize_filename(s: str) -> str:
    return "".join(c if (c.isalnum() or c in ("-", "_", ".")) else "_" for c in s)

def main():
    parser = argparse.ArgumentParser(description="Sum TH2s across ROOT files and fit Y-slices with gaus+pol background.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input ROOT files to sum (same hist names).")
    parser.add_argument("--out", default="pi0_unfold_sum_fits.root", help="Output ROOT file name.")
    parser.add_argument("--csv", default="pi0_unfold_sum_fits.csv", help="Output CSV summary filename.")
    parser.add_argument("--prefix", default="d9_", help="Prefix for target histograms (default: d9_).")
    parser.add_argument("--poly", type=int, default=1, help="Polynomial background order (0..5). Default: 1 (pol1).")
    parser.add_argument("--min-entries", type=int, default=30, help="Minimum entries in slice to attempt a fit. Default: 30.")
    parser.add_argument("--png-dir", default="fit_png", help="Directory for per-slice PNGs (subdirs per histogram).")
    # π0 constraints
    parser.add_argument("--mean", type=float, default=0.133, help="Expected Gaussian mean (default 0.133).")
    parser.add_argument("--mean-window", type=float, default=0.010, help="Allowed ±window around mean (default 0.010).")
    parser.add_argument("--sigma-min", type=float, default=0.003, help="Minimum sigma (default 0.003).")
    parser.add_argument("--sigma-max", type=float, default=0.030, help="Maximum sigma (default 0.030).")
    args = parser.parse_args()

    # Ensure PNG base directory exists
    os.makedirs(args.png_dir, exist_ok=True)

    # Open input files
    files = []
    for path in args.inputs:
        if not os.path.exists(path):
            sys.stderr.write(f"WARNING: Input file not found: {path}\n")
            continue
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            sys.stderr.write(f"WARNING: Could not open: {path}\n")
            continue
        files.append(f)
    if not files:
        sys.stderr.write("ERROR: No input files could be opened.\n")
        sys.exit(1)

    # Determine hist names from the first file
    hist_names = list_target_hists(files[0], prefix=args.prefix)
    if not hist_names:
        sys.stderr.write(f"ERROR: No target TH2 histograms with prefix '{args.prefix}' in {files[0].GetName()}.\n")
        sys.exit(1)

    print(f"Found {len(hist_names)} target TH2 histograms: {', '.join(hist_names)}")

    # Sum
    sums = sum_histograms(files, hist_names)
    if not sums:
        sys.stderr.write("ERROR: No histograms could be summed.\n")
        sys.exit(1)

    # Output ROOT
    outf = ROOT.TFile(args.out, "RECREATE")
    outf.mkdir("sums")
    outf.cd("sums")
    for name, h in sums.items():
        h.Write(name)
    outf.cd()

    # CSV + TTree setup
    csv_fields = [
        "hist_name","binset","xbin","xbin_low","xbin_high","n_entries",
        "status","chi2","ndf",
        "A","A_err","mean","mean_err","sigma","sigma_err","yield","yield_err"
    ]
    max_pol = 5
    for i in range(max_pol + 1):
        csv_fields += [f"b{i}", f"b{i}_err"]

    csv_file = open(args.csv, "w", newline="")
    writer = csv.DictWriter(csv_file, fieldnames=csv_fields)
    writer.writeheader()

    tree = ROOT.TTree("fit_results", "Per-slice fit results")
    from array import array
    i_binset  = array('i', [0])
    i_xbin    = array('i', [0])
    d_xlow    = array('d', [0.0])
    d_xhigh   = array('d', [0.0])
    i_entries = array('i', [0])
    i_status  = array('i', [0])
    d_chi2    = array('d', [0.0])
    i_ndf     = array('i', [0])
    d_A       = array('d', [0.0])
    d_Ae      = array('d', [0.0])
    d_mean    = array('d', [0.0])
    d_meane   = array('d', [0.0])
    d_sigma   = array('d', [0.0])
    d_sigmae  = array('d', [0.0])
    d_yield   = array('d', [0.0])
    d_yielde  = array('d', [0.0])
    b = [array('d', [0.0]) for _ in range(max_pol + 1)]
    be = [array('d', [0.0]) for _ in range(max_pol + 1)]

    tree.Branch("binset",  i_binset,  "binset/I")
    tree.Branch("xbin",    i_xbin,    "xbin/I")
    tree.Branch("xlow",    d_xlow,    "xlow/D")
    tree.Branch("xhigh",   d_xhigh,   "xhigh/D")
    tree.Branch("n_entries", i_entries, "n_entries/I")
    tree.Branch("status",  i_status,  "status/I")
    tree.Branch("chi2",    d_chi2,    "chi2/D")
    tree.Branch("ndf",     i_ndf,     "ndf/I")
    tree.Branch("A",       d_A,       "A/D")
    tree.Branch("A_err",   d_Ae,      "A_err/D")
    tree.Branch("mean",    d_mean,    "mean/D")
    tree.Branch("mean_err",d_meane,   "mean_err/D")
    tree.Branch("sigma",   d_sigma,   "sigma/D")
    tree.Branch("sigma_err",d_sigmae, "sigma_err/D")
    tree.Branch("yield",   d_yield,   "yield/D")
    tree.Branch("yield_err",d_yielde, "yield_err/D")
    for i in range(max_pol + 1):
        tree.Branch(f"b{i}",     b[i],  f"b{i}/D")
        tree.Branch(f"b{i}_err", be[i], f"b{i}_err/D")

    # Canvas for PNG saving (reused)
    canvas = ROOT.TCanvas("c_fit", "c_fit", 900, 700)

    # Per-hist loop
    outf.mkdir("fits")
    for name, h2 in sums.items():
        binset_idx = parse_d9_index(name)
        print(f"Processing {name} (binset {binset_idx}) ...")
        outf.cd("fits")
        if not ROOT.gDirectory.GetDirectory(name):
            ROOT.gDirectory.mkdir(name)
        ROOT.gDirectory.cd(name)

        # Ensure per-hist PNG directory exists
        hist_png_dir = os.path.join(args.png_dir, sanitize_filename(name))
        os.makedirs(hist_png_dir, exist_ok=True)

        nx = h2.GetNbinsX()
        xax = h2.GetXaxis()
        for ix in range(1, nx + 1):
            proj_name = f"{name}_py_{ix:03d}"
            hproj = h2.ProjectionY(proj_name, ix, ix, "e")
            hproj.SetTitle(f"{name} slice xbin={ix}")
            hproj.Write()

            i_binset[0]  = binset_idx
            i_xbin[0]    = ix
            d_xlow[0]    = xax.GetBinLowEdge(ix)
            d_xhigh[0]   = xax.GetBinUpEdge(ix)
            i_entries[0] = int(hproj.GetEntries())

            # If not enough entries: record blanks, skip fit & PNG
            if hproj.GetEntries() < args.min_entries:
                i_status[0] = -1
                d_chi2[0] = 0.0
                i_ndf[0] = 0
                d_A[0]=d_Ae[0]=d_mean[0]=d_meane[0]=d_sigma[0]=d_sigmae[0]=d_yield[0]=d_yielde[0]=0.0
                for k in range(6):
                    b[k][0]=0.0; be[k][0]=0.0
                tree.Fill()
                continue

            # Fit
            fitres = fit_projection(
                hproj, args.poly,
                expected_mean=args.mean,
                mean_window=args.mean_window,
                sigma_min=args.sigma_min,
                sigma_max=args.sigma_max,
                opts="RSQ"
            )
            func = fitres["func"]
            func.SetName(f"{proj_name}_fit")
            func.SetLineColor(ROOT.kRed)
            func.Write()

            # Fill outputs
            i_status[0] = int(fitres["status"])
            d_chi2[0]   = float(fitres["chi2"])
            i_ndf[0]    = int(fitres["ndf"])
            d_A[0]      = float(fitres["A"])
            d_Ae[0]     = float(fitres["A_err"])
            d_mean[0]   = float(fitres["mean"])
            d_meane[0]  = float(fitres["mean_err"])
            d_sigma[0]  = float(fitres["sigma"])
            d_sigmae[0] = float(fitres["sigma_err"])
            d_yield[0]  = float(fitres["yield"])
            d_yielde[0] = float(fitres["yield_err"])

            pars = fitres["pars"]
            errs = fitres["errs"]
            for k in range(6):
                b[k][0]=0.0; be[k][0]=0.0
            for k in range(args.poly + 1):
                b[k][0]  = float(pars[3 + k])
                be[k][0] = float(errs[3 + k])

            row = {
                "hist_name": name,
                "binset": binset_idx,
                "xbin": ix,
                "xbin_low": d_xlow[0],
                "xbin_high": d_xhigh[0],
                "n_entries": i_entries[0],
                "status": i_status[0],
                "chi2": d_chi2[0],
                "ndf": i_ndf[0],
                "A": d_A[0], "A_err": d_Ae[0],
                "mean": d_mean[0], "mean_err": d_meane[0],
                "sigma": d_sigma[0], "sigma_err": d_sigmae[0],
                "yield": d_yield[0], "yield_err": d_yielde[0],
            }
            for kk in range(6):
                row[f"b{kk}"] = b[kk][0]
                row[f"b{kk}_err"] = be[kk][0]
            writer.writerow(row)

            tree.Fill()

            # ---- Save PNG for this fitted slice ----
            canvas.cd()
            hproj.SetLineColor(ROOT.kBlack)
            hproj.Draw("E1")
            func.Draw("SAME")
            canvas.Update()

            png_name = f"{sanitize_filename(proj_name)}_fit.png"
            png_path = os.path.join(hist_png_dir, png_name)
            canvas.Print(png_path)
            # ---------------------------------------

        outf.cd()

    outf.cd()
    log = ROOT.TNamed("command", " ".join(sys.argv))
    log.Write()
    tree.Write()
    csv_file.close()
    outf.Close()

    for f in files:
        f.Close()

    print("Done.")
    print(f"ROOT output: {args.out}")
    print(f"CSV output : {args.csv}")
    print(f"PNGs saved under: {args.png_dir}")

if __name__ == "__main__":
    main()
