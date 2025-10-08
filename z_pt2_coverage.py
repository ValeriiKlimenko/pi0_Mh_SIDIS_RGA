#!/usr/bin/env python3
import os, glob, argparse
import ROOT

def main():
    ap = argparse.ArgumentParser(description="Plot z vs pi0_sidis_PT2 for negative zpt2phit flags")
    ap.add_argument("--input_dir", default="/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/rec/",
                    help="Folder with ROOT files")
    ap.add_argument("--tree", default="h22", help="TTree name (default: h22)")
    ap.add_argument("--out", default="plots_z_pt2", help="Output directory")
    ap.add_argument("--threads", type=int, default=6,
                    help="EnableImplicitMT with N threads (0=single-thread)")
    ap.add_argument("--nx", type=int, default=120, help="z bins (default: 120)")
    ap.add_argument("--ny", type=int, default=120, help="pi0_sidis_PT2 bins (default: 120)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "*.root")))
    if not files:
        raise SystemExit(f"No ROOT files found under: {args.input_dir}")

    if args.threads > 0:
        ROOT.ROOT.EnableImplicitMT(args.threads)

    # Build a single dataframe over all files
    df = ROOT.RDataFrame(args.tree, files)

    # >>> Print all columns (one line) <<<
    print("All columns:", ", ".join(str(c) for c in df.GetColumnNames()))

    # Sanity check columns up front
    needed = {"z", "pi0_sidis_PT2", "zpt2phit_8x8x9", "zpt2phit_8x8x9m"}
    have = set(df.GetColumnNames())
    missing = needed - have
    if missing:
        raise SystemExit(f"Missing columns in tree '{args.tree}': {sorted(missing)}")

    # Discover axes ranges from the data (single pass)
    zmin_r  = df.Min("z");  zmax_r  = df.Max("z")
    pt2min_r = df.Min("pi0_sidis_PT2"); pt2max_r = df.Max("pi0_sidis_PT2")
    ROOT.RDF.RunGraphs([zmin_r, zmax_r, pt2min_r, pt2max_r])

    zMin = float(zmin_r.GetValue());  zMax = float(zmax_r.GetValue())
    yMin = float(pt2min_r.GetValue()); yMax = float(pt2max_r.GetValue())
    if zMax <= zMin: zMax = zMin + 1e-6
    if yMax <= yMin: yMax = yMin + 1e-6

    # Book the 2D model once and reuse
    title = "z vs p_{T}^{2};z;#pi^{0} p_{T}^{2} [GeV^{2}]"
    model = ROOT.RDF.TH2DModel("hz_pt2", title, args.nx, zMin, zMax, args.ny, yMin, yMax)

    # Two histograms with the requested filters
    h_rec = df.Filter("zpt2phit_8x8x9 < 0").Histo2D(model, "z", "pi0_sidis_PT2")
    h_gen = df.Filter("zpt2phit_8x8x9m < 0").Histo2D(model, "z", "pi0_sidis_PT2")

    # Fill both in one event loop
    ROOT.RDF.RunGraphs([h_rec, h_gen])

    # Output dir
    os.makedirs(args.out, exist_ok=True)

    ROOT.gStyle.SetOptStat(0)

    c1 = ROOT.TCanvas("c_rec", "z vs pi0_sidis_PT2 (zpt2phit_8x8x9 < 0)", 900, 800)
    h_rec.Draw("COLZ")
    c1.SetRightMargin(0.15); c1.SetLogz()
    c1.SaveAs(os.path.join(args.out, "z_pt2_zpt2phit_8x8x9_lt0.png"))

    c2 = ROOT.TCanvas("c_gen", "z vs pi0_sidis_PT2 (zpt2phit_8x8x9m < 0)", 900, 800)
    h_gen.Draw("COLZ")
    c2.SetRightMargin(0.15); c2.SetLogz()
    c2.SaveAs(os.path.join(args.out, "z_pt2_zpt2phit_8x8x9m_lt0.png"))

    # Save histograms to a ROOT file too
    fout = ROOT.TFile(os.path.join(args.out, "z_pt2_negative_flag_hists.root"), "RECREATE")
    h_rec.Write("h_z_pt2_rec_flag_neg")
    h_gen.Write("h_z_pt2_gen_flag_neg")
    fout.Close()

    print("Done. Outputs in:", args.out)
    print("Entries:",
          "rec(zpt2phit_8x8x9<0) =", int(h_rec.GetEntries()),
          "| gen(zpt2phit_8x8x9m<0) =", int(h_gen.GetEntries()))

if __name__ == "__main__":
    main()
