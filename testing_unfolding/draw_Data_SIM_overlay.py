#!/usr/bin/env python3
import os
import argparse
import ROOT

# ---------------- binning (must match your analysis) ----------------
N_ZBINS = 8
N_PT2_BINS_WITH_OVERFLOW = 11  # 10 + overflow
N_PHI = 8

Z_EDGES  = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
PT2_EDGES= [0, 0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1.0,1.5]


def phi_center_from_bin(phi_bin_1b: int) -> float:
    w = 360.0 / N_PHI
    return (phi_bin_1b - 0.5) * w


def bin_label(edges, i1b: int) -> str:
    if i1b < 1 or i1b >= len(edges):
        return "out"
    a, b = edges[i1b-1], edges[i1b]
    return f"{a:.2f}–{b:.2f}"


# ---------------------------------------------------------
# Load φ histograms from phi_slices-like file:
#   dirs: ix##/z##/phi_ix##_z##_pt##
# Returns: hmap[ix][(iz, ipt)] = TH1D (cloned, in memory)
# ---------------------------------------------------------
def load_phi_hists(rootfile: str):
    f = ROOT.TFile.Open(rootfile, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open %rootfile")

    hmap = {}
    keep_refs = [f]  # keep file alive if needed

    for key_ix in f.GetListOfKeys():
        obj_ix = key_ix.ReadObj()
        if not isinstance(obj_ix, ROOT.TDirectory):
            continue
        name_ix = obj_ix.GetName()  # "ix01"
        if not name_ix.startswith("ix"):
            continue
        try:
            ix = int(name_ix[2:])
        except ValueError:
            continue

        hmap.setdefault(ix, {})

        for key_z in obj_ix.GetListOfKeys():
            obj_z = key_z.ReadObj()
            if not isinstance(obj_z, ROOT.TDirectory):
                continue
            name_z = obj_z.GetName()  # "z03"
            if not name_z.startswith("z"):
                continue
            try:
                iz = int(name_z[1:])
            except ValueError:
                continue

            for key_h in obj_z.GetListOfKeys():
                h = key_h.ReadObj()
                if not isinstance(h, ROOT.TH1):
                    continue
                hname = h.GetName()      # "phi_ix01_z03_pt04"
                parts = hname.split("_")
                ipt = None
                for part in parts:
                    if part.startswith("pt"):
                        try:
                            ipt = int(part[2:])
                        except ValueError:
                            pass
                if ipt is None:
                    continue

                hc = h.Clone(f"{hname}_{os.path.basename(rootfile)}")
                hc.SetDirectory(0)
                hmap[ix][(iz, ipt)] = hc

    return hmap, keep_refs


# ---------------------------------------------------------
# Load fit parameters (p0,p1,p2) from phi_fit_results.root
# Tree must have branches: ix, iz, iPt, p0
# p1,p2 are optional; if missing, fits are skipped.
# Returns: params[(ix,iz,iPt)] = (p0,p1,p2), has_p1p2 (bool)
# ---------------------------------------------------------
def load_fit_params(rootfile: str):
    f = ROOT.TFile.Open(rootfile, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open %s" % rootfile)
    t = f.Get("phi_fit_p0")
    if not t:
        raise RuntimeError("Tree 'phi_fit_p0' not found in fit file")

    has_p1 = bool(t.GetBranch("p1"))
    has_p2 = bool(t.GetBranch("p2"))
    has_full = has_p1 and has_p2

    if not has_full:
        print("WARNING: p1 and/or p2 branches not found in fit tree. "
              "Will NOT draw fit curves, only use histograms + data.")

    params = {}
    n = t.GetEntries()
    for i in range(n):
        t.GetEntry(i)
        ix  = int(getattr(t, "ix"))
        iz  = int(getattr(t, "iz"))
        ipt = int(getattr(t, "iPt"))
        p0  = float(getattr(t, "p0"))
        p1  = float(getattr(t, "p1")) if has_p1 else 0.0
        p2  = float(getattr(t, "p2")) if has_p2 else 0.0
        params[(ix, iz, ipt)] = (p0, p1, p2)

    return params, has_full, [f]


# ---------------------------------------------------------
# Gather DATA points from tree (your original script logic)
# Returns: data[xq2][(z_bin, pt2_bin)] = [(phi,y,yerr), ...]
# ---------------------------------------------------------
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
            if phib < 1 or phib > N_PHI:
                continue
            phi_deg = phi_center_from_bin(phib)

        if not (0.0 <= phi_deg <= 360.0):
            continue

        data.setdefault(xq2, {}).setdefault((zb, pb), []).append((phi_deg, y, yerr))

    # sort by phi
    for mp in data.values():
        for k in mp:
            mp[k].sort(key=lambda t: t[0])
    return data


def gather_sim_points_per_xq2(files, tree_name="h22_fit"):
    """
    Aggregate simulation like the C++ draw_phi_grids_sim_impl:
      - chain all files
      - for each (xq2bin, z_pt2_phi_bin) sum nPions
        and add errPions in quadrature
      - return: sim[xq2][(z,pt2)] = [(phi_deg, y_sum, yerr_combined), ...]
    """
    if not files:
        return {}

    ch = ROOT.TChain(tree_name)
    for fn in files:
        ch.Add(fn)

    # tmp[xq2][(z_bin, pt2_bin)][phi_bin] = (sumY, sumErr2)
    tmp = {}

    nentries = ch.GetEntries()
    for i in range(nentries):
        ch.GetEntry(i)

        xq2  = int(getattr(ch, "xq2bin"))
        comp = int(getattr(ch, "z_pt2_phi_bin"))
        if comp < 1:
            continue

        # Decode composite index (same as C++ bins::decode_composite_1based)
        phi_bin  = ((comp - 1) % N_PHI) + 1
        zpt2_bin = ((comp - 1) // N_PHI) + 1
        z_bin    = ((zpt2_bin - 1) // N_PT2_BINS_WITH_OVERFLOW) + 1
        pt2_bin  = ((zpt2_bin - 1) %  N_PT2_BINS_WITH_OVERFLOW) + 1

        if (z_bin < 1 or z_bin > N_ZBINS or
            pt2_bin < 1 or pt2_bin > N_PT2_BINS_WITH_OVERFLOW or
            phi_bin < 1 or phi_bin > N_PHI):
            continue

        nPions   = float(getattr(ch, "nPions"))
        errPions = float(getattr(ch, "errPions"))

        cellmap = tmp.setdefault(xq2, {}).setdefault((z_bin, pt2_bin), {})
        sumY, sumE2 = cellmap.get(phi_bin, (0.0, 0.0))
        sumY  += nPions
        sumE2 += errPions * errPions
        cellmap[phi_bin] = (sumY, sumE2)

    # Convert aggregated dict into the same structure as data_all
    sim_all = {}
    for xq2, cellmap in tmp.items():
        sim_all[xq2] = {}
        for (z_bin, pt2_bin), phimap in cellmap.items():
            pts = []
            for phi_bin in sorted(phimap.keys()):
                y, e2 = phimap[phi_bin]
                phi_deg = phi_center_from_bin(phi_bin)
                pts.append((phi_deg, y, e2**0.5))
            sim_all[xq2][(z_bin, pt2_bin)] = pts

    return sim_all


# ---------------------------------------------------------
# Draw one xQ² bin as N_Z × N_PT2 grid with:
#  - unfolded φ hist (blue)
#  - data points (black)
#  - simulation points (green, globally normalized per xQ²)
#  - fit curve p0+p1cos+p2cos2 (red), if available
# ---------------------------------------------------------
def draw_one_xq2(ix, hmap, params, has_full_fit,
                 data_map_for_ix, sim_map_for_ix, outdir):
    """
    Draw one xQ² bin as N_Z × N_PT² grid with:
      - data points (black, with error bars)
      - simulation points (red, scaled globally per xQ²)
    Also dumps all data+sim points for this xQ² bin to a CSV file:
      outdir/csv/xq2_<ix>_data_sim.csv

    Note: hmap, params, has_full_fit are kept in the signature for
          backwards compatibility but are not used (no hist/fit drawn).
    """
    # unused inputs (kept for compatibility with existing main())
    _ = hmap
    _ = params
    _ = has_full_fit

    # ----- global style -----
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.055, "XYZ")
    ROOT.gStyle.SetLabelSize(0.045, "XYZ")

    width  = max(1200, int(N_PT2_BINS_WITH_OVERFLOW * 150))
    height = max(900,  int(N_ZBINS * 150))
    c = ROOT.TCanvas(f"c_xq2_{ix}", f"xQ^2 bin {ix}", width, height)
    c.Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS)

    latex_objs = []
    graphs_data = []
    graphs_sim  = []

    # ---- global normalization factor for simulation in this xQ² ----
    data_sum = 0.0
    sim_sum  = 0.0
    for pts in data_map_for_ix.values():
        for (_, y, _) in pts:
            data_sum += y
    for pts in sim_map_for_ix.values():
        for (_, y, _) in pts:
            sim_sum += y

    if sim_sum > 0 and data_sum > 0:
        sim_scale = data_sum / sim_sum
    else:
        sim_scale = 1.0  # nothing to scale, or one set is empty

    # ---- global y-range for this xQ² bin (makes grid nicer) ----
    ymax_global = 0.0
    for pts in data_map_for_ix.values():
        for (_, y, _) in pts:
            if y > ymax_global:
                ymax_global = y
    for pts in sim_map_for_ix.values():
        for (_, y, _) in pts:
            val = y * sim_scale
            if val > ymax_global:
                ymax_global = val
    if ymax_global <= 0:
        ymax_global = 1.0
    ymax_global *= 1.25  # small headroom

    # ---- container for CSV rows ----
    # each row: (xq2, z_bin, pt2_bin, phi_deg, y, yerr, source)
    csv_rows = []

    # ---------------------------------------------------------
    # loop over z, pT² pads
    # ---------------------------------------------------------
    for iz in range(1, N_ZBINS + 1):
        for ipt in range(1, N_PT2_BINS_WITH_OVERFLOW + 1):
            padnum = (iz - 1) * N_PT2_BINS_WITH_OVERFLOW + ipt
            c.cd(padnum)
            pad = ROOT.gPad

            pad.SetFillColor(0)
            pad.SetFrameFillColor(0)
            pad.SetTicks(1, 1)
            pad.SetGrid(1, 1)
            pad.SetLeftMargin(0.20 if ipt == 1 else 0.08)
            pad.SetBottomMargin(0.24 if iz == N_ZBINS else 0.08)
            pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.06)
            pad.SetBit(ROOT.TPad.kClipFrame, True)

            pts_data = data_map_for_ix.get((iz, ipt), [])
            pts_sim  = sim_map_for_ix.get((iz, ipt), [])

            # frame with shared y-range for this xQ²
            frame = pad.DrawFrame(0.0, 0.0, 360.0, ymax_global)
            frame.GetXaxis().SetTitle("#phi_{Trento} [deg]")
            frame.GetYaxis().SetTitle("Counts")
            frame.GetXaxis().CenterTitle(True)
            frame.GetYaxis().CenterTitle(True)
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)
            frame.GetXaxis().SetTitleOffset(1.1)
            frame.GetYaxis().SetTitleOffset(1.2)

            # hide inner axis labels
            if iz != N_ZBINS:
                frame.GetXaxis().SetLabelSize(0)
                frame.GetXaxis().SetTitleSize(0)
            if ipt != 1:
                frame.GetYaxis().SetLabelSize(0)
                frame.GetYaxis().SetTitleSize(0)

            # ---- simulation (red, scaled) ----
            if pts_sim:
                gs = ROOT.TGraphErrors(len(pts_sim))
                for i, (x, y, yerr) in enumerate(pts_sim):
                    y_scaled = y * sim_scale
                    yerr_scaled = yerr * sim_scale
                    gs.SetPoint(i, x, y_scaled)
                    gs.SetPointError(i, 0.0, yerr_scaled)

                    # add to CSV rows
                    csv_rows.append((ix, iz, ipt, x, y_scaled, yerr_scaled, "sim_scaled"))

                gs.SetMarkerStyle(21)
                gs.SetMarkerSize(0.9)
                gs.SetLineWidth(2)
                gs.SetLineStyle(2)
                gs.SetLineColor(ROOT.kRed + 1)
                gs.SetMarkerColor(ROOT.kRed + 1)
                gs.Draw("PC SAME")  # line + markers
                graphs_sim.append(gs)

            # ---- data (black, with error bars) ----
            if pts_data:
                g = ROOT.TGraphErrors(len(pts_data))
                for i, (x, y, yerr) in enumerate(pts_data):
                    g.SetPoint(i, x, y)
                    g.SetPointError(i, 0.0, yerr)

                    # add to CSV rows
                    csv_rows.append((ix, iz, ipt, x, y, yerr, "data"))

                g.SetMarkerStyle(20)
                g.SetMarkerSize(1.0)
                g.SetLineWidth(2)
                g.SetLineColor(ROOT.kBlack)
                g.SetMarkerColor(ROOT.kBlack)
                g.Draw("PE1 SAME")  # markers + error bars
                graphs_data.append(g)

            # top row pT² label
            if iz == 1:
                t = ROOT.TLatex()
                t.SetNDC()
                t.SetTextFont(42)
                t.SetTextSize(0.045)
                t.DrawLatex(0.18, 0.92,
                            f"p_{{T}}^{{2}} {bin_label(PT2_EDGES, ipt)}")
                latex_objs.append(t)

            # first column z label
            if ipt == 1:
                t = ROOT.TLatex()
                t.SetNDC()
                t.SetTextFont(42)
                t.SetTextSize(0.045)
                t.DrawLatex(0.20, 0.85,
                            f"z {bin_label(Z_EDGES, iz)}")
                latex_objs.append(t)

    # global title
    c.cd(0)
    st = ROOT.TLatex()
    st.SetNDC()
    st.SetTextFont(42)
    st.SetTextSize(0.05)
    st.DrawLatex(0.02, 0.98, f"xQ^2 bin: {ix}")

    # legend in pad 1
    c.cd(1)
    leg = ROOT.TLegend(0.10, 0.60, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)

    dummy_data = ROOT.TGraph()
    dummy_data.SetMarkerStyle(20)
    dummy_data.SetMarkerSize(1.0)
    dummy_data.SetLineWidth(2)
    dummy_data.SetLineColor(ROOT.kBlack)
    dummy_data.SetMarkerColor(ROOT.kBlack)

    dummy_sim = ROOT.TGraph()
    dummy_sim.SetMarkerStyle(21)
    dummy_sim.SetMarkerSize(0.9)
    dummy_sim.SetLineWidth(2)
    dummy_sim.SetLineStyle(2)
    dummy_sim.SetLineColor(ROOT.kRed + 1)
    dummy_sim.SetMarkerColor(ROOT.kRed + 1)

    leg.AddEntry(dummy_data, "data", "lep")
    leg.AddEntry(dummy_sim,  "simulation (scaled)", "lep")
    leg.Draw()

    # ----- save plots -----
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"xq2_{ix}.png")
    pdf = os.path.join(outdir, f"xq2_{ix}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)

    # ----- save CSV with points for this xQ² bin -----
    csv_dir = os.path.join(outdir, "csv")
    os.makedirs(csv_dir, exist_ok=True)
    csv_path = os.path.join(csv_dir, f"xq2_{ix}_data_sim.csv")

    with open(csv_path, "w") as f:
        f.write("xq2,z_bin,pt2_bin,phi_deg,y,yerr,source\n")
        for row in csv_rows:
            f.write("{:d},{:d},{:d},{:.6f},{:.6f},{:.6f},{}\n".format(*row))

    print(f"  saved {png}, {pdf}, and {csv_path}")




def main():
    ap = argparse.ArgumentParser(
        description="Overlay unfolded φ-slices + fit + data + simulation per xQ^2 bin."
    )
    ap.add_argument("--phi_slices", required=False,
                    default="../phi_slices.root",
                    help="phi_slices ROOT file (with ix##/z## dirs)")
    ap.add_argument("--fit_results", required=False,
                    default="../phi_fit_results.root",
                    help="phi_fit_results ROOT file (with tree phi_fit_p0)")
    ap.add_argument("--data_root",
                    default="/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/testing_unfolding/data_cover.root",
                    help="ROOT file with data tree (e.g. data_cover.root)")
    ap.add_argument("--data_tree", default="h22_fit",
                    help="Tree name in data_root")

    # NEW: simulation inputs
    ap.add_argument("--sim_root", default="../unfolding_rec_true/h3_bin_xBQ2_Valerii_20_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_7_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_6_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_12_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_19_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_11_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_5_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_15_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_16_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_14_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_4_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_3_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_10_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_13_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_9_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_8_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_18_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_1_fitted.root,../unfolding_rec_true/h3_bin_xBQ2_Valerii_2_fitted.root",
                    help="Comma-separated list of sim ROOT files (h3*_fitted.root, etc.)")
    ap.add_argument("--sim_tree", default="h22_fit",
                    help="Tree name in simulation files")
    ap.add_argument("--sim_use_gen", action="store_true",
                    help="Use *_gen branches in sim (xq2bin_gen, z_pt2_phi_bin_gen)")

    ap.add_argument("--outdir", default="plots_phi_unfold_fit_data_sim")
    ap.add_argument("--max_ix", type=int, default=None)
    args = ap.parse_args()

    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    # unfolded histograms
    hmap, refs1 = load_phi_hists(args.phi_slices)

    # fit parameters
    params, has_full_fit, refs2 = load_fit_params(args.fit_results)

    # data points
    fdata = ROOT.TFile.Open(args.data_root, "READ")
    if not fdata or fdata.IsZombie():
        raise RuntimeError(f"Cannot open data file {args.data_root}")
    tdata = fdata.Get(args.data_tree)
    if not tdata:
        raise RuntimeError(f"Tree '{args.data_tree}' not found in data file")
    data_all = gather_points_per_xq2(tdata)

    if args.sim_root:
        sim_files = [s for s in args.sim_root.split(",") if s.strip()]
        sim_all = gather_sim_points_per_xq2(sim_files, args.sim_tree)
    else:
        sim_all = {}

    # xQ² bins to draw = union of keys
    ix_vals = set(hmap.keys()) | set(data_all.keys()) | set(sim_all.keys())
    ix_vals = sorted(ix_vals)
    if args.max_ix is not None:
        ix_vals = [ix for ix in ix_vals if ix <= args.max_ix]

    for ix in ix_vals:
        data_map_for_ix = data_all.get(ix, {})
        sim_map_for_ix  = sim_all.get(ix, {})
        draw_one_xq2(ix, hmap, params, has_full_fit,
                     data_map_for_ix, sim_map_for_ix, args.outdir)

    fdata.Close()


if __name__ == "__main__":
    main()
