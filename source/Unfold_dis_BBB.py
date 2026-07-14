#!/usr/bin/env python3
import csv
import math
import argparse
import ROOT

def get_hist(tfile, name):
    h = tfile.Get(name)
    if not h:
        raise RuntimeError(f"Histogram '{name}' not found in file '{tfile.GetName()}'")
    return h

def safe_div(num, den):
    if den == 0:
        return float("nan")
    return num / den

def main():
    ap = argparse.ArgumentParser(description="Bin-by-bin unfolding CSV: data*gen/rec")
    ap.add_argument("--data", default="unfolding_dis/dis_data_hists_all.root")
    ap.add_argument("--rec",  default="unfolding_dis/dis_rec_hists_all.root")
    ap.add_argument("--gen",  default="unfolding_dis/dis_gen_hists_all.root")
    ap.add_argument("--out",  default="unfolding_dis/binbybin_unfolding.csv")
    args = ap.parse_args()

    f_data = ROOT.TFile.Open(args.data, "READ")
    f_rec  = ROOT.TFile.Open(args.rec,  "READ")
    f_gen  = ROOT.TFile.Open(args.gen,  "READ")
    if not f_data or f_data.IsZombie(): raise RuntimeError(f"Cannot open {args.data}")
    if not f_rec  or f_rec.IsZombie():  raise RuntimeError(f"Cannot open {args.rec}")
    if not f_gen  or f_gen.IsZombie():  raise RuntimeError(f"Cannot open {args.gen}")

    h_data = get_hist(f_data, "h_xQ2_data")
    h_rec  = get_hist(f_rec,  "h_xQ2_rec")
    h_gen  = get_hist(f_gen,  "h_xQ2_gen")

    nb = h_data.GetNbinsX()
    if h_rec.GetNbinsX() != nb or h_gen.GetNbinsX() != nb:
        raise RuntimeError("Histogram binnings do not match (different number of bins).")

    # Write CSV with requested columns
    with open(args.out, "w", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(["data", "rec", "gen", "rec/gen", "data*gen/rec"])

        for i in range(1, nb + 1):  # exclude under/overflow (0 and nb+1)
            data = float(h_data.GetBinContent(i))
            rec  = float(h_rec.GetBinContent(i))
            gen  = float(h_gen.GetBinContent(i))

            rec_over_gen = safe_div(rec, gen)
            unfolded     = safe_div(data * gen, rec)

            w.writerow([data, rec, gen, rec_over_gen, unfolded])

    print(f"Wrote: {args.out}")

    f_data.Close()
    f_rec.Close()
    f_gen.Close()

if __name__ == "__main__":
    main()