
# default="/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/data_cover.root"

# File: plot_phi_by_xq2_root_clean.py
# Usage:
#   python plot_phi_by_xq2_root_clean.py /path/to/augmented.root \
#          --tree h22_fit --outdir plots_phi_by_xq2 --max_xq2 16

import argparse, os
import ROOT

# ---------------- binning (match your analysis) ----------------
N_ZBINS = 8
N_PT2_BINS_WITH_OVERFLOW = 11  # 10 + overflow
N_PHI = 9
Z_EDGES  = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
PT2_EDGES= [0, 0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1.0,1.5]

def phi_center_from_bin(phi_bin_1b: int) -> float:
    w = 360.0 / N_PHI
    return (phi_bin_1b - 0.5) * w

def bin_label(edges, i1b: int) -> str:
    if i1b < 1 or i1b >= len(edges): return "out"
    a, b = edges[i1b-1], edges[i1b]
    return f"{a:.2f}–{b:.2f}"

def gather_points_per_xq2(tree: ROOT.TTree):
    has_phi     = bool(tree.GetBranch("phi"))
    has_phi_bin = bool(tree.GetBranch("phi_bin"))
    if not has_phi and not has_phi_bin:
        raise RuntimeError("Need 'phi' or 'phi_bin' in the tree.")

    data = {}  # xq2 -> {(z,pt2): [(phi,y,yerr), ...]}
    nentries = tree.GetEntries()
    for i in range(nentries):
        tree.GetEntry(i)
        xq2  = int(getattr(tree, "xq2bin"))
        zb   = int(getattr(tree, "z_bin"))
        pb   = int(getattr(tree, "pt2_bin"))
        if zb < 1 or zb > N_ZBINS or pb < 1 or pb > N_PT2_BINS_WITH_OVERFLOW:
            continue
        y    = float(getattr(tree, "nPions"))
        yerr = float(getattr(tree, "errPions"))

        if has_phi:
            phi_deg = float(getattr(tree, "phi"))
        else:
            phib = int(getattr(tree, "phi_bin"))
            if phib < 1 or phib > N_PHI: continue
            phi_deg = phi_center_from_bin(phib)

        if not (0.0 <= phi_deg <= 360.0): continue

        data.setdefault(xq2, {}).setdefault((zb, pb), []).append((phi_deg, y, yerr))

    # sort by phi
    for mp in data.values():
        for k in mp:
            mp[k].sort(key=lambda t: t[0])
    return data

def draw_one_xq2(xq2v, binmap, outdir):
    # Global style: kill borders/stats, keep tick marks
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)

    width  = max(1200, int(N_PT2_BINS_WITH_OVERFLOW * 140))
    height = max(900,  int(N_ZBINS * 140))
    c = ROOT.TCanvas(f"c_xq2_{xq2v}", f"xQ^2 bin {xq2v}", width, height)
    c.Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS)  # nx, ny (cols, rows)

    # keep refs
    graphs = []
    labels = []

    for iz in range(1, N_ZBINS+1):
        for ip in range(1, N_PT2_BINS_WITH_OVERFLOW+1):
            padnum = (iz-1)*N_PT2_BINS_WITH_OVERFLOW + ip
            c.cd(padnum)
            pad = ROOT.gPad
            # clean pad look; no grid = no “fence” artifacts
            pad.SetFillColor(0)
            pad.SetFrameFillColor(0)
            pad.SetTicks(1,1)
            pad.SetLeftMargin(0.18 if ip==1 else 0.06)
            pad.SetBottomMargin(0.22 if iz==N_ZBINS else 0.08)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBit(ROOT.TPad.kClipFrame, True)  # enforce clipping

            pts = binmap.get((iz, ip), [])
            ymax = max((p[1] for p in pts), default=0.0)
            if ymax <= 0: ymax = 1.0

            # axis via DrawFrame (lighter and safer than TH2 per pad)
            frame = pad.DrawFrame(0.0, 0.0, 360.0, 1.2*ymax)
            frame.GetXaxis().SetTitle("#phi [deg]")
            frame.GetYaxis().SetTitle("nPions")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)
            if iz != N_ZBINS:
                frame.GetXaxis().SetLabelSize(0)
                frame.GetXaxis().SetTitleSize(0)
            if ip != 1:
                frame.GetYaxis().SetLabelSize(0)
                frame.GetYaxis().SetTitleSize(0)

            if pts:
                g = ROOT.TGraphErrors(len(pts))
                for i,(x,y,yerr) in enumerate(pts):
                    g.SetPoint(i, x, y)
                    g.SetPointError(i, 0.0, yerr)
                g.SetMarkerStyle(20)
                g.SetMarkerSize(0.7)
                g.SetLineWidth(1)
                g.SetLineColor(ROOT.kBlack)
                g.SetMarkerColor(ROOT.kBlack)
                g.Draw("P SAME")
                graphs.append(g)

            # top row: pt2 label; first col: z label
            if iz == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.18, 0.92, f"p_{{T}}^{{2}} {bin_label(PT2_EDGES, ip)}")
                labels.append(t)
            if ip == 1:
                t = ROOT.TLatex()
                t.SetNDC(); t.SetTextSize(0.04)
                t.DrawLatex(0.20, 0.85, f"z {bin_label(Z_EDGES, iz)}")
                labels.append(t)

    # super title
    c.cd(0)
    st = ROOT.TLatex()
    st.SetNDC(); st.SetTextSize(0.04)
    st.DrawLatex(0.02, 0.98, f"xQ^2 bin: {xq2v}")

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"xq2_{xq2v}.png")
    pdf = os.path.join(outdir, f"xq2_{xq2v}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"  saved {png}")

def main():
    ap = argparse.ArgumentParser(description="Clean ROOT-only φ-plots (no grid artifacts).")
    ap.add_argument("--rootfile", default="/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/data_cover.root")
    ap.add_argument("--tree", default="h22_fit")
    ap.add_argument("--outdir", default="plots_phi_by_xq2_clean")
    ap.add_argument("--max_xq2", type=int, default=None)
    args = ap.parse_args()

    f = ROOT.TFile.Open(args.rootfile, "READ")
    if not f or f.IsZombie(): raise RuntimeError(f"Cannot open {args.rootfile}")
    t = f.Get(args.tree)
    if not t: raise RuntimeError(f"Tree '{args.tree}' not found")

    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    data = gather_points_per_xq2(t)
    xq2_vals = sorted(data.keys())
    if args.max_xq2 is not None:
        xq2_vals = [v for v in xq2_vals if v <= args.max_xq2]

    for xv in xq2_vals:
        draw_one_xq2(xv, data[xv], args.outdir)

    f.Close()

if __name__ == "__main__":
    main()
