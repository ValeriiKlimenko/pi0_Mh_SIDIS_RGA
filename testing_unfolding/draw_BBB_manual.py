#!/usr/bin/env python3
import os
import argparse
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


def draw_one_xq2_overlay(iy: int,
                         h1: ROOT.TH2,
                         h2: ROOT.TH2,
                         outdir: str,
                         label1: str = "file1",
                         label2: str = "file2"):
    """
    Draw phi vs counts for all (z, pt2) bins for a single Y-bin index 'iy'
    (one xQ^2 bin), overlaying distributions from h1 and h2 in different colors.

    One canvas per xQ^2 bin:
      - rows    = z bins
      - columns = pt2 bins (including overflow)
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetEndErrorSize(4)

    yaxis = h1.GetYaxis()
    xq2_center = yaxis.GetBinCenter(iy)
    xq2_bin_int = int(round(xq2_center))

    width  = max(1400, int(N_PT2_BINS_WITH_OVERFLOW * 150))
    height = max(900,  int(N_ZBINS * 150))
    c = ROOT.TCanvas(f"c_overlay_xq2_{xq2_bin_int}",
                     f"overlay: xQ^2 bin {xq2_bin_int}",
                     width, height)
    c.Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS)

    latex_objs = []
    graphs = []

    phi_min, phi_max = 0.0, 360.0

    for iz in range(1, N_ZBINS+1):
        for ipt in range(1, N_PT2_BINS_WITH_OVERFLOW+1):
            padnum = (iz-1)*N_PT2_BINS_WITH_OVERFLOW + ipt
            c.cd(padnum)
            pad = ROOT.gPad

            pad.SetFillColor(0)
            pad.SetFrameFillColor(0)
            pad.SetTicks(1,1)
            pad.SetLeftMargin(0.18)
            pad.SetBottomMargin(0.22 if iz==N_ZBINS else 0.08)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBit(ROOT.TPad.kClipFrame, True)

            # Extract phi distributions for this (z, pt2, xq2) for both files
            pts1 = get_phi_points_for_cell(h1, iy, iz, ipt)
            pts2 = get_phi_points_for_cell(h2, iy, iz, ipt)

            # Per-pad Y max (shared axis for both)
            ymax = 1.0
            all_y = [y for (_, y, _) in pts1] + [y for (_, y, _) in pts2]
            if all_y:
                ymax = max(all_y)
            if ymax <= 0:
                ymax = 1.0
            ymax *= 1.2

            frame = pad.DrawFrame(phi_min, 0.0, phi_max, ymax)
            frame.GetXaxis().SetTitle("#phi_{Trento} [deg]")
            frame.GetYaxis().SetTitle("")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)

            # X axis styling
            if iz == N_ZBINS:
                frame.GetXaxis().SetLabelSize(0.08)
                frame.GetXaxis().SetTitleSize(0.08)
            else:
                frame.GetXaxis().SetLabelSize(0.0)
                frame.GetXaxis().SetTitleSize(0.0)
            frame.GetXaxis().SetTitleOffset(1.0)

            # Y axis: keep labels on every pad, larger font, no title
            frame.GetYaxis().SetLabelSize(0.08)
            frame.GetYaxis().SetTitleSize(0.01)
            frame.GetYaxis().SetTitleOffset(0.6)

            # First file (blue)
            if pts1:
                g1 = ROOT.TGraphErrors(len(pts1))
                for i,(x,y,ye) in enumerate(pts1):
                    g1.SetPoint(i, x, y)
                    g1.SetPointError(i, 0.0, ye)
                g1.SetMarkerStyle(20)
                g1.SetMarkerSize(0.7)
                g1.SetLineWidth(1)
                g1.SetLineColor(ROOT.kBlue+1)
                g1.SetMarkerColor(ROOT.kBlue+1)
                g1.Draw("PC SAME")
                graphs.append(g1)

            # Second file (red)
            if pts2:
                g2 = ROOT.TGraphErrors(len(pts2))
                for i,(x,y,ye) in enumerate(pts2):
                    g2.SetPoint(i, x, y)
                    g2.SetPointError(i, 0.0, ye)
                g2.SetMarkerStyle(21)
                g2.SetMarkerSize(0.7)
                g2.SetLineWidth(1)
                g2.SetLineColor(ROOT.kRed+1)
                g2.SetMarkerColor(ROOT.kRed+1)
                g2.Draw("PC SAME")
                graphs.append(g2)

            # top row pT^2 label
            if iz == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.05)
                t.DrawLatex(0.18, 0.92,
                            f"p_{{T}}^{{2}} {bin_label(PT2_EDGES, ipt)}")
                latex_objs.append(t)

            # first column z label
            if ipt == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.05)
                t.DrawLatex(0.20, 0.85,
                            f"z {bin_label(Z_EDGES, iz)}")
                latex_objs.append(t)

    # global title
    c.cd(0)
    st = ROOT.TLatex()
    st.SetNDC(); st.SetTextSize(0.045)
    st.DrawLatex(0.02, 0.98, f"xQ^2 bin (center): {xq2_bin_int}")

    # simple legend in pad 1
    c.cd(1)
    leg = ROOT.TLegend(0.10, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    dummy1 = ROOT.TH1F("dummy1", "", 1, 0, 1)
    dummy1.SetLineColor(ROOT.kBlue+1)
    dummy1.SetMarkerColor(ROOT.kBlue+1)
    dummy1.SetMarkerStyle(20)

    dummy2 = ROOT.TH1F("dummy2", "", 1, 0, 1)
    dummy2.SetLineColor(ROOT.kRed+1)
    dummy2.SetMarkerColor(ROOT.kRed+1)
    dummy2.SetMarkerStyle(21)

    leg.AddEntry(dummy1, label1, "lep")
    leg.AddEntry(dummy2, label2, "lep")
    leg.Draw()

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"xq2_{xq2_bin_int:02d}_overlay.png")
    pdf = os.path.join(outdir, f"xq2_{xq2_bin_int:02d}_overlay.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"  saved {png}, {pdf} for this xQ^2 bin")


def main():
    ap = argparse.ArgumentParser(
        description="Overlay φ-dependence for all z–pt2 bins per xQ^2 bin "
                    "from two TH2 histograms (same structure) in two ROOT files."
    )
    ap.add_argument("--root1", required=False, default = "../manual_bbb_Mxcut_1/measData_times_truth_over_meas.root",
                    help="First ROOT file")
    ap.add_argument("--root2", required=False,default = "../measData_times_truth_over_meas.root",
                    help="Second ROOT file")
    ap.add_argument("--hname", default="h_measData_times_truth_over_meas",
                    help="Name of TH2 in both ROOT files")
    ap.add_argument("--outdir", default="plots_phi_overlay",
                    help="Output directory for PNG/PDF grids")
    ap.add_argument("--max_ix", type=int, default=None,
                    help="Optional maximum xQ2 bin (center value) to draw")
    ap.add_argument("--label1", default="file1",
                    help="Legend label for first file")
    ap.add_argument("--label2", default="file2",
                    help="Legend label for second file")
    args = ap.parse_args()

    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    f1 = ROOT.TFile.Open(args.root1, "READ")
    if not f1 or f1.IsZombie():
        raise RuntimeError(f"Cannot open file {args.root1}")

    f2 = ROOT.TFile.Open(args.root2, "READ")
    if not f2 or f2.IsZombie():
        raise RuntimeError(f"Cannot open file {args.root2}")

    h1 = f1.Get(args.hname)
    if not h1:
        raise RuntimeError(f"{args.hname} not found in {args.root1}")

    h2 = f2.Get(args.hname)
    if not h2:
        raise RuntimeError(f"{args.hname} not found in {args.root2}")

    # Basic sanity check: same binning
    if h1.GetNbinsX() != h2.GetNbinsX() or h1.GetNbinsY() != h2.GetNbinsY():
        raise RuntimeError("Histograms have different binning in X or Y")

    ny = h1.GetNbinsY()
    yaxis = h1.GetYaxis()

    # Loop over all xQ2 Y bins
    for iy in range(1, ny+1):
        xq2_center = yaxis.GetBinCenter(iy)
        xq2_int = int(round(xq2_center))
        if args.max_ix is not None and xq2_int > args.max_ix:
            continue

        # Optionally skip bins with no content at all in both hists
        proj1 = h1.ProjectionX("_px_tmp1", iy, iy)
        proj2 = h2.ProjectionX("_px_tmp2", iy, iy)
        if proj1.Integral() == 0 and proj2.Integral() == 0:
            continue

        draw_one_xq2_overlay(iy, h1, h2, args.outdir,
                             label1=args.label1,
                             label2=args.label2)

    f1.Close()
    f2.Close()
    print(f"Finished overlay plots in directory: {args.outdir}")


if __name__ == "__main__":
    main()
