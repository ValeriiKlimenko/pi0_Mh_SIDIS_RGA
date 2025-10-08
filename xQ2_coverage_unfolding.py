#!/usr/bin/env python3
import os, glob, argparse
import ROOT

def main():
    ap = argparse.ArgumentParser(description="Plot xB vs Q2 for -5 bins")
    ap.add_argument("--input_dir", default="/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/rec/", help="Folder with ROOT files")
    ap.add_argument("--tree", default="h22", help="TTree name (default: h22)")
    ap.add_argument("--out", default="plots", help="Output directory")
    ap.add_argument("--threads", type=int, default=6,
                    help="EnableImplicitMT with N threads (0=single-thread)")
    ap.add_argument("--nx", type=int, default=120, help="xB bins (default: 100)")
    ap.add_argument("--ny", type=int, default=120, help="Q2 bins (default: 100)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "*.root")))
    if not files:
        raise SystemExit(f"No ROOT files found under: {args.input_dir}")

    # (Optional) multithreading
    if args.threads > 0:
        ROOT.ROOT.EnableImplicitMT(args.threads)

    # Build a single dataframe over all files
    df = ROOT.RDataFrame(args.tree, files)

    # Sanity check columns up front
    needed = {"xB", "Q2", "bin_xBQ2_Valerii", "bin_xBQ2_Valeriim"}
    have = set(df.GetColumnNames())
    missing = needed - have
    if missing:
        raise SystemExit(f"Missing columns in tree '{args.tree}': {sorted(missing)}")

    # Discover ranges from the data (one pass for all 4)
    xmin_r = df.Min("xB"); xmax_r = df.Max("xB")
    q2min_r = df.Min("Q2"); q2max_r = df.Max("Q2")
    ROOT.RDF.RunGraphs([xmin_r, xmax_r, q2min_r, q2max_r])

    xMin = float(xmin_r.GetValue()); xMax = float(xmax_r.GetValue())
    yMin = float(q2min_r.GetValue()); yMax = float(q2max_r.GetValue())
    if xMax <= xMin: xMax = xMin + 1e-6
    if yMax <= yMin: yMax = yMin + 1e-6

    # Book both histograms with same binning/model
    title = "xB vs Q^{2};xB;Q^{2} [GeV^{2}]"
    model = ROOT.ROOT.RDF.TH2DModel("hxbq2", title, args.nx, xMin, xMax, args.ny, yMin, yMax)

    h_valerii  = df.Filter("bin_xBQ2_Valerii == -5").Histo2D(model, "xB", "Q2")
    h_valerii_gen  = df.Filter("bin_xBQ2_Valerii == -5").Histo2D(model, "xBm", "Q2m")
  
    h_valeriim = df.Filter("bin_xBQ2_Valeriim == -5").Histo2D(model, "xB", "Q2")
    h_valeriim_gen = df.Filter("bin_xBQ2_Valeriim == -5").Histo2D(model, "xBm", "Q2m")

    # Fill both in one event loop
    ROOT.RDF.RunGraphs([h_valerii, h_valeriim, h_valerii_gen, h_valeriim_gen])

    # Make output dir
    os.makedirs(args.out, exist_ok=True)

    # Style & draw
    ROOT.gStyle.SetOptStat(0)

    c1 = ROOT.TCanvas("c_valerii", "xB vs Q2 (bin_xBQ2_Valerii == -5)", 900, 800)
    h_valerii.Draw("COLZ")
    c1.SetRightMargin(0.15); c1.SetLogz()
    c1.SaveAs(os.path.join(args.out, "xb_q2_bin_xBQ2_Valerii_eq_minus5_withMigr.png"))

    c2 = ROOT.TCanvas("c_valeriim", "xB vs Q2 (bin_xBQ2_Valeriim == -5)", 900, 800)
    h_valeriim.Draw("COLZ")
    c2.SetRightMargin(0.15); c2.SetLogz()
    c2.SaveAs(os.path.join(args.out, "xb_q2_bin_xBQ2_Valeriim_eq_minus5_withMigr.png"))


    c3 = ROOT.TCanvas("c_valerii_gen", "xB vs Q2 (bin_xBQ2_Valerii == -5)", 900, 800)
    h_valerii.Draw("COLZ")
    c3.SetRightMargin(0.15); c3.SetLogz()
    c3.SaveAs(os.path.join(args.out, "Cut_REC_GeN_distrib.png"))

    c4 = ROOT.TCanvas("c_valeriim_gen", "xB vs Q2 (bin_xBQ2_Valeriim == -5)", 900, 800)
    h_valeriim.Draw("COLZ")
    c4.SetRightMargin(0.15); c4.SetLogz()
    c4.SaveAs(os.path.join(args.out, "Cut_GEN_GEN_distrib.png"))

    # Save histograms to a ROOT file too
    fout = ROOT.TFile(os.path.join(args.out, "xb_q2_minus5_hists.root"), "RECREATE")
    h_valerii.Write("h_xb_q2_Valerii_m5")
    h_valeriim.Write("h_xb_q2_Valeriim_m5")
    fout.Close()

    print("Done. Outputs in:", args.out)
    print("Entries:", "Valerii", int(h_valerii.GetEntries()),
          "| Valeriim", int(h_valeriim.GetEntries()))

if __name__ == "__main__":
    main()
