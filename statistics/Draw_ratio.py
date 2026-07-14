#!/usr/bin/env python3
import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DEFAULT_FILE = "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/statistics/pi0_nPions_ratio.root"

# keep TF1s alive until canvases are saved
_GLOBAL_FITS_INDIV = []
_GLOBAL_FITS_GRID  = []

def is_hist(obj) -> bool:
    return bool(obj.InheritsFrom("TH1"))

def style_hist(h):
    h.SetLineWidth(2)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.7)
    # ratio-like plots: helpful y-range
    nm = h.GetName().lower()
    if "ratio" in nm or "mc_over_data" in nm:
        h.SetMinimum(0)
        if h.GetMaximum() < 2.0:
            h.SetMaximum(2.0)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)

def fit_horizontal(h, fit_color=ROOT.kRed+1):
    """
    Fit histogram with a horizontal line (pol0) over full X range.
    Returns (tf1, value, err, chi2, ndf, pval).
    """
    xmin = h.GetXaxis().GetXmin()
    xmax = h.GetXaxis().GetXmax()
    f = ROOT.TF1(f"fconst_{h.GetName()}", "pol0", xmin, xmax)
    f.SetLineColor(fit_color)
    f.SetLineWidth(2)

    # Quiet + Store result; don't suppress drawing so the function appears
    fitres = h.Fit(f, "QS")
    # In some ROOT versions fitres may be None if fit failed silently
    chi2 = fitres.Chi2() if fitres else f.GetChisquare()
    ndf  = fitres.Ndf()  if fitres else max(1, f.GetNDF())
    val  = f.GetParameter(0)
    err  = f.GetParError(0)
    pval = ROOT.TMath.Prob(chi2, ndf) if ndf > 0 else 0.0
    return f, val, err, chi2, ndf, pval

def draw_and_save(h, out_dir: str):
    """
    Draw a single histogram with a horizontal fit.
    Legend: 'MC/Data = <fit result>' (no error shown).
    """
    style_hist(h)
    c = ROOT.TCanvas("c_" + h.GetName(), "", 900, 700)
    c.SetGrid()
    h.Draw("E1")
    f, val, err, chi2, ndf, pval = fit_horizontal(h)
    _GLOBAL_FITS_INDIV.append(f)  # keep alive
    f.Draw("same")

    # Legend with MC/Data = fit result (no error)
    leg = ROOT.TLegend(0.60, 0.80, 0.88, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(0, f"MC/Data = {val:.3g}", "")
    leg.Draw()

    out_path = os.path.join(out_dir, f"{h.GetName()}.png")
    c.SaveAs(out_path)
    c.Close()

def visit_dir(tdir, out_dir: str, firstN_list=None, Nwant=16):
    """Save each hist; optionally collect first N hists into firstN_list."""
    keys = tdir.GetListOfKeys()
    if not keys:
        return
    for key in keys:
        obj = key.ReadObj()
        if obj.InheritsFrom("TDirectory"):
            sub_out = os.path.join(out_dir, obj.GetName())
            os.makedirs(sub_out, exist_ok=True)
            visit_dir(obj, sub_out, firstN_list, Nwant)
        elif is_hist(obj):
            draw_and_save(obj, out_dir)
            if firstN_list is not None and len(firstN_list) < Nwant:
                clone = obj.Clone(obj.GetName() + "_clone_for_grid")
                clone.SetDirectory(0)
                firstN_list.append(clone)

def draw_grid_4x4(hists, out_png):
    """
    Draw up to 16 histograms on a 4x4 canvas, each fitted with a horizontal line.
    - Title in each panel: 'x-Q2 = panel number [1:16]'
    - Remove histogram name from the panel.
    - Legend text: 'MC/Data = <fit result>' (no error).
    """
    if not hists:
        print("[WARN] No histograms to draw on grid.")
        return

    c = ROOT.TCanvas("c_grid_first16", "", 1600, 1200)
    c.Divide(4, 4, 0.001, 0.001)

    for i, h in enumerate(hists[:16]):
        pad = c.cd(i+1)
        pad.SetGrid()
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.12)
        pad.SetTopMargin(0.08)

        style_hist(h)
        h.Draw("E1")

        f, val, err, chi2, ndf, pval = fit_horizontal(h)
        _GLOBAL_FITS_GRID.append(f)  # keep alive
        f.Draw("same")

        # Panel title: x-Q2 = panel number [1:16]
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextSize(0.040)
        title.DrawLatex(0.15, 0.92, f"x-Q2 = {i+1}")

        # Legend with MC/Data = fit result (no error)
        leg = ROOT.TLegend(0.55, 0.78, 0.88, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(0, f"MC/Data = {val:.3g}", "")
        leg.Draw()

    c.SaveAs(out_png)
    c.Close()

def main():
    in_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_FILE
    if not os.path.exists(in_file):
        sys.exit(f"[ERROR] File not found: {in_file}")

    out_dir = os.path.join(os.path.dirname(in_file), "plots_ratio")
    os.makedirs(out_dir, exist_ok=True)

    f = ROOT.TFile.Open(in_file, "READ")
    if not f or f.IsZombie():
        sys.exit(f"[ERROR] Failed to open ROOT file: {in_file}")

    first16 = []
    visit_dir(f, out_dir, firstN_list=first16, Nwant=16)
    f.Close()

    grid_png = os.path.join(out_dir, "grid_first16_with_const_fits.png")
    draw_grid_4x4(first16, grid_png)

    print(f"[OK] Saved individual histograms (with const fit) to: {out_dir}")
    print(f"[OK] Saved 4x4 grid of first 16 histograms (with const fit) to: {grid_png}")

if __name__ == "__main__":
    main()
