#!/usr/bin/env python3
# make_plots_rdf_mt_debug.py
import os
import glob
import argparse
import ROOT

REQUIRED_BRANCHES = ["xB", "Q2", "z", "pi0_sidis_PT2", "isEventINbins"]

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def warn(msg):
    print(f"[WARN] {msg}", flush=True)

def err(msg):
    print(f"[ERROR] {msg}", flush=True)

def find_root_files(input_dir: str):
    files = glob.glob(os.path.join(input_dir, "**", "*.root"), recursive=True)
    if not files:
        raise FileNotFoundError(f"No .root files found under: {input_dir}")
    return files

def verify_tree_and_branches(files, tree, sample_files=8):
    """Open a few files to ensure the tree and branches exist."""
    ok_files = 0
    first_good_file = None
    missing_tree = []
    unreadable = []
    for f in files[:sample_files]:
        tf = ROOT.TFile.Open(f)
        if not tf or tf.IsZombie():
            unreadable.append(f)
            continue
        obj = tf.Get(tree)
        if not obj:
            missing_tree.append(f)
            tf.Close()
            continue
        if obj.InheritsFrom("TTree"):
            if first_good_file is None:
                first_good_file = f
            ok_files += 1
        tf.Close()

    if unreadable:
        warn(f"{len(unreadable)} file(s) unreadable (first: {unreadable[0]})")
    if missing_tree:
        warn(f"Tree '{tree}' missing in {len(missing_tree)} sample file(s) (first: {missing_tree[0]})")

    if ok_files == 0:
        raise RuntimeError(f"Could not find TTree '{tree}' in the first {sample_files} files.")

    # Check required branches in one good file
    tf = ROOT.TFile.Open(first_good_file)
    t = tf.Get(tree)
    existing = {b.GetName() for b in t.GetListOfBranches()}
    tf.Close()

    missing_branches = [b for b in REQUIRED_BRANCHES if b not in existing]
    if missing_branches:
        warn(f"Missing branches in '{tree}': {missing_branches}")
        # Suggest close matches
        suggestions = {}
        for need in missing_branches:
            close = [b for b in existing if need.lower() in b.lower()]
            if close:
                suggestions[need] = close[:5]
        if suggestions:
            warn(f"Similar branch name suggestions: {suggestions}")

    log(f"Verified tree '{tree}' and inspected branches in: {first_good_file}")
    return missing_branches

def pad_range(xmin, xmax, frac=0.02):
    if xmin is None or xmax is None:
        return (0.0, 1.0)
    if xmin == xmax:
        eps = 1e-9 if xmin == 0 else abs(xmin) * 1e-6
        return (xmin - eps, xmax + eps)
    width = xmax - xmin
    pad = width * frac
    return (xmin - pad, xmax + pad)

def save_hist2d_png(hist, title, outpath):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    # Safety: zero-entry hist still produces a plot, but call out explicitly
    entries = int(hist.GetEntries())
    integral = hist.Integral()
    log(f"Saving '{outpath}' | entries={entries} integral={integral:.0f}")
    c = ROOT.TCanvas("c", "c", 900, 700)
    hist.SetTitle(title)
    hist.GetZaxis().SetTitle("Counts")
    hist.Draw("COLZ")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    c.SaveAs(outpath)
    c.Close()

def main():
    parser = argparse.ArgumentParser(
        description="RDataFrame (multithreaded) 2D plots with diagnostics: Q2 vs xB and z vs pi0_sidis_PT2, split by isEventINbins."
    )
    parser.add_argument("--input_dir", required=True, help="Folder containing .root files (searched recursively)")
    parser.add_argument("--tree", default="h22", help="TTree name (default: h22)")
    parser.add_argument("--out_dir", default="plots_rdf", help="Output directory for PNGs")
    parser.add_argument("--bins_q2_xb", type=int, nargs=2, default=(100, 100), help="Bins (nx ny) for xB vs Q2")
    parser.add_argument("--bins_z_pt2", type=int, nargs=2, default=(100, 100), help="Bins (nx ny) for z vs pi0_sidis_PT2")
    parser.add_argument("-j", "--threads", type=int, default=0, help="Threads for ROOT implicit MT (0=auto)")
    parser.add_argument("--cutvar", default="isEventINbins", help="Cut variable name (default: isEventINbins)")
    parser.add_argument("--debug_list_branches", action="store_true", help="List all branches (from first file with the tree) and exit")
    args = parser.parse_args()

    # --- Multithreading ---
    if args.threads is not None and args.threads >= 0:
        if args.threads == 0:
            ROOT.ROOT.EnableImplicitMT()  # auto
        else:
            ROOT.ROOT.EnableImplicitMT(args.threads)
    # Report MT status
    mt_on = bool(getattr(ROOT.ROOT, "IsImplicitMTEnabled", lambda: False)())
    try:
        pool = getattr(ROOT.ROOT, "GetThreadPoolSize", lambda: None)()
    except Exception:
        pool = None
    log(f"ROOT {ROOT.gROOT.GetVersion()} | ImplicitMT={'ON' if mt_on else 'OFF'}"
        + (f" (threads={pool})" if pool is not None else ""))

    # --- Files ---
    try:
        files = find_root_files(args.input_dir)
    except Exception as e:
        err(str(e))
        return
    log(f"Found {len(files)} .root file(s) under: {args.input_dir}")
    log(f"First 3 files: {files[:3]}")

    # --- Quick tree/branch sanity checks (before building RDF) ---
    missing_branches = verify_tree_and_branches(files, args.tree)
    if args.debug_list_branches:
        # Already printed warnings; list all branches from first good file and exit
        tf = ROOT.TFile.Open(files[0])
        t = tf.Get(args.tree)
        names = [b.GetName() for b in t.GetListOfBranches()]
        tf.Close()
        log(f"Branches in '{args.tree}' (first file): {names}")
        return
    if missing_branches:
        warn("Some required branches are missing. The script will proceed, but plots that rely on them may be empty.")

    # --- Build RDF ---
    df = ROOT.RDataFrame(args.tree, files)
    # Ensure cut column exists
    if not df.HasColumn(args.cutvar):
        warn(f"Cut column '{args.cutvar}' not found in dataframe columns. Available columns (sample):")
        log(list(df.GetColumnNames())[:50])
        # Try common alternates
        for alt in ["isEventINbins_m", "isEventINbins", "isEventINbins_true", "isEventINbins_reco"]:
            if df.HasColumn(alt):
                warn(f"Found alternate cut column: '{alt}'. Use: --cutvar {alt}")
        return

    # Count entries
    total_count = df.Count()
    n_total = int(total_count.GetValue())
    log(f"Total entries in '{args.tree}': {n_total}")
    if n_total == 0:
        warn("No entries found. Check that the tree name is correct and not empty in your files.")
        return

    # Ranges (these actions run in one event loop together)
    xB_min = df.Min("xB"); xB_max = df.Max("xB")
    Q2_min = df.Min("Q2"); Q2_max = df.Max("Q2")
    z_min = df.Min("z");   z_max = df.Max("z")
    pt2_min = df.Min("pi0_sidis_PT2"); pt2_max = df.Max("pi0_sidis_PT2")

    xB_min, xB_max = float(xB_min.GetValue()), float(xB_max.GetValue())
    Q2_min, Q2_max = float(Q2_min.GetValue()), float(Q2_max.GetValue())
    z_min, z_max = float(z_min.GetValue()), float(z_max.GetValue())
    pt2_min, pt2_max = float(pt2_min.GetValue()), float(pt2_max.GetValue())

    log(f"xB range: [{xB_min:.6g}, {xB_max:.6g}]")
    log(f"Q2 range: [{Q2_min:.6g}, {Q2_max:.6g}]")
    log(f"z  range: [{z_min:.6g}, {z_max:.6g}]")
    log(f"pi0_sidis_PT2 range: [{pt2_min:.6g}, {pt2_max:.6g}]")

    xB_min, xB_max = pad_range(xB_min, xB_max)
    Q2_min, Q2_max = pad_range(Q2_min, Q2_max)
    z_min, z_max = pad_range(z_min, z_max)
    pt2_min, pt2_max = pad_range(pt2_min, pt2_max)

    # --- Filters with names so df.Report() is useful ---
    cutvar = args.cutvar
    df0 = df.Filter(f"{cutvar} == 0", "isEventINbins_eq_0")
    df1 = df.Filter(f"{cutvar} == 1", "isEventINbins_eq_1")

    c0 = df0.Count()
    c1 = df1.Count()
    n0 = int(c0.GetValue())
    n1 = int(c1.GetValue())
    log(f"Entries passing {cutvar}=0: {n0}")
    log(f"Entries passing {cutvar}=1: {n1}")

    if n0 == 0 and n1 == 0:
        warn(f"Both filters are empty. Verify '{cutvar}' values in your files (maybe try --cutvar isEventINbins_m).")
        report = df.Report()
        report.Print()  # prints nothing if no named filters yet
        return

    # --- Histograms ---
    nx_xb, ny_q2 = args.bins_q2_xb
    nx_z, ny_pt2 = args.bins_z_pt2

    h_q2xb_0 = df0.Histo2D(("h_q2xb_0", "Q2 vs xB; xB; Q2", nx_xb, xB_min, xB_max, ny_q2, Q2_min, Q2_max), "xB", "Q2")
    h_q2xb_1 = df1.Histo2D(("h_q2xb_1", "Q2 vs xB; xB; Q2", nx_xb, xB_min, xB_max, ny_q2, Q2_min, Q2_max), "xB", "Q2")

    h_pt2z_0 = df0.Histo2D(("h_pt2z_0", "pi0_sidis_PT2 vs z; z; pi0_sidis_PT2",
                             nx_z, z_min, z_max, ny_pt2, pt2_min, pt2_max),
                            "z", "pi0_sidis_PT2")
    h_pt2z_1 = df1.Histo2D(("h_pt2z_1", "pi0_sidis_PT2 vs z; z; pi0_sidis_PT2",
                             nx_z, z_min, z_max, ny_pt2, pt2_min, pt2_max),
                            "z", "pi0_sidis_PT2")

    # Materialize and report
    out = args.out_dir
    os.makedirs(out, exist_ok=True)

    # Print filter report (now we have named filters)
    report = df.Report()
    log("Filter report:")
    report.Print()

    # Save plots; call GetValue() to trigger computation of each histogram
    save_hist2d_png(h_q2xb_0.GetValue(), f"Q2 vs xB ({cutvar}=0)", os.path.join(out, "Q2_vs_xB_isEventINbins0.png"))
    save_hist2d_png(h_q2xb_1.GetValue(), f"Q2 vs xB ({cutvar}=1)", os.path.join(out, "Q2_vs_xB_isEventINbins1.png"))
    save_hist2d_png(h_pt2z_0.GetValue(),  f"pi0_sidis_PT2 vs z ({cutvar}=0)", os.path.join(out, "pi0_sidis_PT2_vs_z_isEventINbins0.png"))
    save_hist2d_png(h_pt2z_1.GetValue(),  f"pi0_sidis_PT2 vs z ({cutvar}=1)", os.path.join(out, "pi0_sidis_PT2_vs_z_isEventINbins1.png"))

    log(f"Saved plots in: {out}")

if __name__ == "__main__":
    main()
