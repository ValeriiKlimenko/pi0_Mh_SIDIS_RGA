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
        # Encode composite index: same convention as your C++/Python
        # zpt2_bin: 1..(N_ZBINS*N_PT2_BINS_WITH_OVERFLOW)
        zpt2_bin = (z_bin - 1) * N_PT2_BINS_WITH_OVERFLOW + pt2_bin
        comp = (zpt2_bin - 1) * N_PHI + phi_bin  # 1..N_ZBINS*N_PT2*N_PHI

        # Coordinate is the integer comp; TH2 axis is uniform -0.5..nZ-0.5
        xcoord = float(comp)
        bx = h2.GetXaxis().FindBin(xcoord)
        if bx < 1 or bx > h2.GetNbinsX():
            continue

        y = h2.GetBinContent(bx, ybin)
        e = h2.GetBinError(bx, ybin)
        phi_deg = phi_center_from_bin(phi_bin)
        pts.append((phi_deg, y, e))
    return pts


def scale_meas_to_truth(pts_meas, pts_truth):
    """
    Rescale measured points so that:
        sum_y(meas_scaled) = sum_y(truth)
    Truth points are NOT modified.

    pts_* are lists of (phi, y, yerr).
    """
    if not pts_meas or not pts_truth:
        return pts_meas

    sum_meas  = sum(y for (_, y, _) in pts_meas)
    sum_truth = sum(y for (_, y, _) in pts_truth)

    if sum_meas <= 0.0 or sum_truth <= 0.0:
        # Nothing sensible to do; leave measured as-is
        return pts_meas

    factor = sum_truth / sum_meas
    return [(phi, y * factor, err * factor) for (phi, y, err) in pts_meas]


def draw_one_xq2_from_qa(iy: int,
                         h_meas: ROOT.TH2,
                         h_truth: ROOT.TH2,
                         outdir: str):
    """
    Draw phi vs counts for all (z, pt2) bins for a single Y-bin index 'iy'
    (one xQ2 bin). Uses h_meas (measured) and h_truth (truth) from qa_response.root.

    Truth is kept with its original normalization.
    Measured is rescaled per (z,pt2,xQ2) cell to match the truth integral.
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)

    yaxis = h_meas.GetYaxis()
    xq2_center = yaxis.GetBinCenter(iy)
    xq2_bin_int = int(round(xq2_center))

    width  = max(1200, int(N_PT2_BINS_WITH_OVERFLOW * 140))
    height = max(900,  int(N_ZBINS * 140))
    c = ROOT.TCanvas(f"c_xq2_{xq2_bin_int}",
                     f"QA: xQ^2 bin {xq2_bin_int}",
                     width, height)
    c.Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS)

    latex_objs = []
    graphs_meas = []
    graphs_truth = []

    for iz in range(1, N_ZBINS+1):
        for ipt in range(1, N_PT2_BINS_WITH_OVERFLOW+1):
            padnum = (iz-1)*N_PT2_BINS_WITH_OVERFLOW + ipt
            c.cd(padnum)
            pad = ROOT.gPad

            pad.SetFillColor(0)
            pad.SetFrameFillColor(0)
            pad.SetTicks(1,1)
            pad.SetLeftMargin(0.18 if ipt==1 else 0.06)
            pad.SetBottomMargin(0.22 if iz==N_ZBINS else 0.08)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBit(ROOT.TPad.kClipFrame, True)

            # Extract phi distributions for this (z, pt2, xq2)
            pts_meas  = get_phi_points_for_cell(h_meas,  iy, iz, ipt)
            pts_truth = get_phi_points_for_cell(h_truth, iy, iz, ipt)

            # --- rescale measured to have same area as truth (per cell) ---
            pts_meas_scaled = scale_meas_to_truth(pts_meas, pts_truth)

            # Compute y-max from truth + scaled measured
            ymax = 1.0
            if pts_truth:
                ymax = max(ymax, max(y for (_, y, _) in pts_truth))
            if pts_meas_scaled:
                ymax = max(ymax, max(y for (_, y, _) in pts_meas_scaled))
            if ymax <= 0:
                ymax = 1.0

            frame = pad.DrawFrame(0.0, 0.0, 360.0, 1.2*ymax)
            frame.GetXaxis().SetTitle("#phi_{Trento} [deg]")
            frame.GetYaxis().SetTitle("counts (meas scaled to truth)")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)

            if iz != N_ZBINS:
                frame.GetXaxis().SetLabelSize(0)
                frame.GetXaxis().SetTitleSize(0)
            if ipt != 1:
                frame.GetYaxis().SetLabelSize(0)
                frame.GetYaxis().SetTitleSize(0)

            # Measured (blue, scaled)
            if pts_meas_scaled:
                gM = ROOT.TGraphErrors(len(pts_meas_scaled))
                for i,(x,y,ye) in enumerate(pts_meas_scaled):
                    gM.SetPoint(i, x, y)
                    gM.SetPointError(i, 0.0, ye)
                gM.SetMarkerStyle(24)
                gM.SetMarkerSize(0.7)
                gM.SetLineWidth(1)
                gM.SetLineColor(ROOT.kBlue+1)
                gM.SetMarkerColor(ROOT.kBlue+1)
                gM.Draw("PC SAME")
                graphs_meas.append(gM)

            # Truth (black) – unchanged, no error bars ("PCX")
            if pts_truth:
                gT = ROOT.TGraphErrors(len(pts_truth))
                for i,(x,y,ye) in enumerate(pts_truth):
                    gT.SetPoint(i, x, y)
                    gT.SetPointError(i, 0.0, ye)
                gT.SetMarkerStyle(20)
                gT.SetMarkerSize(0.7)
                gT.SetLineWidth(1)
                gT.SetLineColor(ROOT.kBlack)
                gT.SetMarkerColor(ROOT.kBlack)
                gT.Draw("PCX SAME")
                graphs_truth.append(gT)

            # top row pt label
            if iz == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.18, 0.92,
                            f"p_{{T}}^{{2}} {bin_label(PT2_EDGES, ipt)}")
                latex_objs.append(t)

            # first column z label
            if ipt == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.20, 0.85,
                            f"z {bin_label(Z_EDGES, iz)}")
                latex_objs.append(t)

    # global title
    c.cd(0)
    st = ROOT.TLatex()
    st.SetNDC(); st.SetTextSize(0.04)
    st.DrawLatex(0.02, 0.98, f"xQ^2 bin (center): {xq2_bin_int}")

    # legend in pad 1
    c.cd(1)
    leg = ROOT.TLegend(0.10, 0.55, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    dummy_meas  = ROOT.TH1F("dummy_meas",  "", 1, 0, 1)
    dummy_truth = ROOT.TH1F("dummy_truth", "", 1, 0, 1)

    dummy_meas.SetLineColor(ROOT.kBlue+1)
    dummy_meas.SetMarkerColor(ROOT.kBlue+1)
    dummy_meas.SetMarkerStyle(24)

    dummy_truth.SetLineColor(ROOT.kBlack)
    dummy_truth.SetMarkerColor(ROOT.kBlack)
    dummy_truth.SetMarkerStyle(20)

    leg.AddEntry(dummy_truth, "truth MC (trained)",            "lep")
    leg.AddEntry(dummy_meas,  "measured MC (scaled to truth)", "lep")
    leg.Draw()

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"xq2_{xq2_bin_int:02d}.png")
    pdf = os.path.join(outdir, f"xq2_{xq2_bin_int:02d}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"  saved {png}")


# --------------------------------------------------------------------
# NEW: generic pair-drawing for extra 2D histograms
# --------------------------------------------------------------------
def draw_one_xq2_pair(iy: int,
                      h_first: ROOT.TH2,
                      h_second: ROOT.TH2,
                      outdir: str,
                      prefix: str,
                      label_first: str,
                      label_second: str,
                      scale_first_to_second: bool = True):
    """
    Generic 2-hist version: draw phi vs counts for all (z,pt2) bins for one xQ^2 bin.
    h_first is drawn in blue, optionally scaled so that integral(first)=integral(second)
    per (z,pt2,xQ^2) cell. h_second is drawn in black.
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)

    yaxis = h_first.GetYaxis()
    xq2_center = yaxis.GetBinCenter(iy)
    xq2_bin_int = int(round(xq2_center))

    width  = max(1200, int(N_PT2_BINS_WITH_OVERFLOW * 140))
    height = max(900,  int(N_ZBINS * 140))
    c = ROOT.TCanvas(f"c_{prefix}_xq2_{xq2_bin_int}",
                     f"{prefix}: xQ^2 bin {xq2_bin_int}",
                     width, height)
    c.Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS)

    latex_objs = []
    graphs_first = []
    graphs_second = []

    for iz in range(1, N_ZBINS+1):
        for ipt in range(1, N_PT2_BINS_WITH_OVERFLOW+1):
            padnum = (iz-1)*N_PT2_BINS_WITH_OVERFLOW + ipt
            c.cd(padnum)
            pad = ROOT.gPad

            pad.SetFillColor(0)
            pad.SetFrameFillColor(0)
            pad.SetTicks(1,1)
            pad.SetLeftMargin(0.18 if ipt==1 else 0.06)
            pad.SetBottomMargin(0.22 if iz==N_ZBINS else 0.08)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBit(ROOT.TPad.kClipFrame, True)

            pts_first  = get_phi_points_for_cell(h_first,  iy, iz, ipt)
            pts_second = get_phi_points_for_cell(h_second, iy, iz, ipt)

            if scale_first_to_second:
                pts_first_scaled = scale_meas_to_truth(pts_first, pts_second)
            else:
                pts_first_scaled = pts_first

            ymax = 1.0
            if pts_second:
                ymax = max(ymax, max(y for (_, y, _) in pts_second))
            if pts_first_scaled:
                ymax = max(ymax, max(y for (_, y, _) in pts_first_scaled))
            if ymax <= 0:
                ymax = 1.0

            ytitle = (f"counts ({label_first} scaled to {label_second})"
                      if scale_first_to_second else "counts")

            frame = pad.DrawFrame(0.0, 0.0, 360.0, 1.2*ymax)
            frame.GetXaxis().SetTitle("#phi_{Trento} [deg]")
            frame.GetYaxis().SetTitle(ytitle)
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)

            if iz != N_ZBINS:
                frame.GetXaxis().SetLabelSize(0)
                frame.GetXaxis().SetTitleSize(0)
            if ipt != 1:
                frame.GetYaxis().SetLabelSize(0)
                frame.GetYaxis().SetTitleSize(0)

            # First (blue)
            if pts_first_scaled:
                g1 = ROOT.TGraphErrors(len(pts_first_scaled))
                for i,(x,y,ye) in enumerate(pts_first_scaled):
                    g1.SetPoint(i, x, y)
                    g1.SetPointError(i, 0.0, ye)
                g1.SetMarkerStyle(24)
                g1.SetMarkerSize(0.7)
                g1.SetLineWidth(1)
                g1.SetLineColor(ROOT.kBlue+1)
                g1.SetMarkerColor(ROOT.kBlue+1)
                g1.Draw("PC SAME")
                graphs_first.append(g1)

            # Second (black)
            if pts_second:
                g2 = ROOT.TGraphErrors(len(pts_second))
                for i,(x,y,ye) in enumerate(pts_second):
                    g2.SetPoint(i, x, y)
                    g2.SetPointError(i, 0.0, ye)
                g2.SetMarkerStyle(20)
                g2.SetMarkerSize(0.7)
                g2.SetLineWidth(1)
                g2.SetLineColor(ROOT.kBlack)
                g2.SetMarkerColor(ROOT.kBlack)
                g2.Draw("PCX SAME")
                graphs_second.append(g2)

            # top row pt label
            if iz == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.18, 0.92,
                            f"p_{{T}}^{{2}} {bin_label(PT2_EDGES, ipt)}")
                latex_objs.append(t)

            # first column z label
            if ipt == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.20, 0.85,
                            f"z {bin_label(Z_EDGES, iz)}")
                latex_objs.append(t)

    # global title
    c.cd(0)
    st = ROOT.TLatex()
    st.SetNDC(); st.SetTextSize(0.04)
    st.DrawLatex(0.02, 0.98,
                 f"{prefix}, xQ^2 bin (center): {xq2_bin_int}")

    # legend in pad 1
    c.cd(1)
    leg = ROOT.TLegend(0.10, 0.55, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    dummy_first  = ROOT.TH1F(f"dummy_first_{prefix}_{xq2_bin_int}",  "", 1, 0, 1)
    dummy_second = ROOT.TH1F(f"dummy_second_{prefix}_{xq2_bin_int}", "", 1, 0, 1)

    dummy_first.SetLineColor(ROOT.kBlue+1)
    dummy_first.SetMarkerColor(ROOT.kBlue+1)
    dummy_first.SetMarkerStyle(24)

    dummy_second.SetLineColor(ROOT.kBlack)
    dummy_second.SetMarkerColor(ROOT.kBlack)
    dummy_second.SetMarkerStyle(20)

    text_first = label_first + (" (scaled)" if scale_first_to_second else "")
    leg.AddEntry(dummy_second, label_second, "lep")
    leg.AddEntry(dummy_first,  text_first,   "lep")
    leg.Draw()

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{prefix}_xq2_{xq2_bin_int:02d}.png")
    pdf = os.path.join(outdir, f"{prefix}_xq2_{xq2_bin_int:02d}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"  saved {png}")


# --------------------------------------------------------------------
# NEW: single-hist version (e.g. for h_check, residual map)
# --------------------------------------------------------------------
def draw_one_xq2_single(iy: int,
                        h2: ROOT.TH2,
                        outdir: str,
                        prefix: str,
                        label: str,
                        signed: bool = False):
    """
    Draw phi vs counts for all (z,pt2) bins for a single Y-bin index 'iy'
    using only one TH2 histogram.

    If signed=True, the y-axis is symmetric around 0 to show positive and
    negative contents (useful for residual maps).
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)

    yaxis = h2.GetYaxis()
    xq2_center = yaxis.GetBinCenter(iy)
    xq2_bin_int = int(round(xq2_center))

    width  = max(1200, int(N_PT2_BINS_WITH_OVERFLOW * 140))
    height = max(900,  int(N_ZBINS * 140))
    c = ROOT.TCanvas(f"c_{prefix}_xq2_{xq2_bin_int}",
                     f"{prefix}: xQ^2 bin {xq2_bin_int}",
                     width, height)
    c.Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS)

    latex_objs = []
    graphs = []

    for iz in range(1, N_ZBINS+1):
        for ipt in range(1, N_PT2_BINS_WITH_OVERFLOW+1):
            padnum = (iz-1)*N_PT2_BINS_WITH_OVERFLOW + ipt
            c.cd(padnum)
            pad = ROOT.gPad

            pad.SetFillColor(0)
            pad.SetFrameFillColor(0)
            pad.SetTicks(1,1)
            pad.SetLeftMargin(0.18 if ipt==1 else 0.06)
            pad.SetBottomMargin(0.22 if iz==N_ZBINS else 0.08)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBit(ROOT.TPad.kClipFrame, True)

            pts = get_phi_points_for_cell(h2, iy, iz, ipt)

            if signed:
                if pts:
                    ys = [y for (_, y, _) in pts]
                    ymax = max(abs(min(ys)), abs(max(ys)))
                    if ymax <= 0:
                        ymax = 1.0
                    ymin = -1.2 * ymax
                    ymax =  1.2 * ymax
                else:
                    ymin, ymax = -1.0, 1.0
            else:
                ymax = 1.0
                if pts:
                    ymax = max(ymax, max(y for (_, y, _) in pts))
                if ymax <= 0:
                    ymax = 1.0
                ymin = 0.0
                ymax = 1.2 * ymax

            frame = pad.DrawFrame(0.0, ymin, 360.0, ymax)
            frame.GetXaxis().SetTitle("#phi_{Trento} [deg]")
            frame.GetYaxis().SetTitle("counts" + (" (signed)" if signed else ""))
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)

            if iz != N_ZBINS:
                frame.GetXaxis().SetLabelSize(0)
                frame.GetXaxis().SetTitleSize(0)
            if ipt != 1:
                frame.GetYaxis().SetLabelSize(0)
                frame.GetYaxis().SetTitleSize(0)

            if pts:
                g = ROOT.TGraphErrors(len(pts))
                for i,(x,y,ye) in enumerate(pts):
                    g.SetPoint(i, x, y)
                    g.SetPointError(i, 0.0, ye)
                g.SetMarkerStyle(20)
                g.SetMarkerSize(0.7)
                g.SetLineWidth(1)
                g.SetLineColor(ROOT.kBlack)
                g.SetMarkerColor(ROOT.kBlack)
                g.Draw("PCX SAME")
                graphs.append(g)

            # top row pt label
            if iz == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.18, 0.92,
                            f"p_{{T}}^{{2}} {bin_label(PT2_EDGES, ipt)}")
                latex_objs.append(t)

            # first column z label
            if ipt == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.20, 0.85,
                            f"z {bin_label(Z_EDGES, iz)}")
                latex_objs.append(t)

    # global title
    c.cd(0)
    st = ROOT.TLatex()
    st.SetNDC(); st.SetTextSize(0.04)
    st.DrawLatex(0.02, 0.98,
                 f"{prefix} ({label}), xQ^2 bin (center): {xq2_bin_int}")

    # legend in pad 1
    c.cd(1)
    leg = ROOT.TLegend(0.10, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    dummy = ROOT.TH1F(f"dummy_single_{prefix}_{xq2_bin_int}", "", 1, 0, 1)
    dummy.SetLineColor(ROOT.kBlack)
    dummy.SetMarkerColor(ROOT.kBlack)
    dummy.SetMarkerStyle(20)

    leg.AddEntry(dummy, label, "lep")
    leg.Draw()

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{prefix}_xq2_{xq2_bin_int:02d}.png")
    pdf = os.path.join(outdir, f"{prefix}_xq2_{xq2_bin_int:02d}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"  saved {png}")


def main():
    ap = argparse.ArgumentParser(
        description="Draw φ-dependence (truth vs measured, and extra QA histos) "
                    "for all z–pt2 bins per xQ^2 bin, from qa_response.root."
    )
    ap.add_argument("--qa_root", default="../rooUnfold_bbb/qa_response.root",
                    help="qa_response.root produced by ResponseQA(pack, ...)")
    ap.add_argument("--outdir", default="plots_phi_qa",
                    help="Output directory for PNG/PDF grids")
    ap.add_argument("--max_ix", type=int, default=None,
                    help="Optional maximum xQ2 bin (center value) to draw")
    args = ap.parse_args()

    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    f = ROOT.TFile.Open(args.qa_root, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open QA file {args.qa_root}")

    h_meas  = f.Get("h_measured_trained")
    h_truth = f.Get("h_truth_trained")
    if not h_meas or not h_truth:
        raise RuntimeError("h_measured_trained and/or h_truth_trained not found in qa_response.root")

    # NEW: extra histograms (if present)
    h_check  = f.Get("h_check_debug")
    h_resid  = f.Get("h_truth_minus_rec_minus_hcheck_debug")
    h_rec    = f.Get("h_rec_vs_zrec_debug")
    h_gen    = f.Get("h_gen_vs_zgen_debug")

    if not h_check:
        print("NOTE: h_check_debug not found; skipping those plots.")
    if not h_resid:
        print("NOTE: h_truth_minus_rec_minus_hcheck_debug not found; skipping those plots.")
    if not h_rec or not h_gen:
        print("NOTE: h_rec_vs_zrec_debug and/or h_gen_vs_zgen_debug not found; "
              "skipping resp_fills plots.")

    # Sanity: shared axes as expected (base hist)
    if (h_meas.GetNbinsX() != h_truth.GetNbinsX() or
        h_meas.GetNbinsY() != h_truth.GetNbinsY()):
        print("WARNING: h_measured_trained and h_truth_trained have different binning!")

    ny = h_meas.GetNbinsY()
    yaxis = h_meas.GetYaxis()

    # Loop over all xQ2 Y bins
    for iy in range(1, ny+1):
        xq2_center = yaxis.GetBinCenter(iy)
        xq2_int = int(round(xq2_center))
        if args.max_ix is not None and xq2_int > args.max_ix:
            continue
        # Optionally skip xQ2 bins with no content at all in the base histos
        if h_meas.ProjectionX("_px_tmp", iy, iy).Integral() == 0 and \
           h_truth.ProjectionX("_px_tmp2", iy, iy).Integral() == 0:
            continue

        # Original QA: truth vs measured
        draw_one_xq2_from_qa(iy, h_meas, h_truth, args.outdir)

        # NEW: response fills (zrec,xrec) vs (zgen,xgen) used in resp->Fill
        # NEW: response reco fills (zrec,xrec) drawn alone
        if h_rec:
            out_rec = os.path.join(args.outdir, "resp_rec")
            draw_one_xq2_single(
                iy, h_rec,
                out_rec, "respReco",
                label="response reco (zrec,xrec)",
                signed=False
            )

        # NEW: response truth fills (zgen,xgen) drawn alone
        if h_gen:
            out_gen = os.path.join(args.outdir, "resp_gen")
            draw_one_xq2_single(
                iy, h_gen,
                out_gen, "respTruth",
                label="response truth (zgen,xgen)",
                signed=False
            )

        # NEW: h_check (zgen,xgen) alone
        if h_check:
            out_hchk = os.path.join(args.outdir, "h_check")
            draw_one_xq2_single(
                iy, h_check,
                out_hchk, "h_check",
                label="h_check (zgen,xgen)",
                signed=False
            )

        # NEW: residual map (truth - reco - h_check) alone, with signed y
        if h_resid:
            out_resid = os.path.join(args.outdir, "truth_minus_rec_minus_hcheck")
            draw_one_xq2_single(
                iy, h_resid,
                out_resid, "truthMinusRecMinusHcheck",
                label="truth minus reco minus h_check (residual)",
                signed=True
            )

    f.Close()


if __name__ == "__main__":
    main()
