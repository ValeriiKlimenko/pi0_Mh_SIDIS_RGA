// ratio_mc_data.cxx
// Requires: fit_pi0_mass.cxx
//
// Changes:
// - Ratio is now MC/Data.
// - Drop (do not fill) bins where the relative uncertainty on the ratio exceeds 40%.
//   rel_unc = sqrt( (σ_MC/MC)^2 + (σ_Data/Data)^2 )
// - Additionally require BOTH Data and MC to have > 20 nPions in the bin before plotting the ratio.

#include <filesystem>
#include <memory>
#include <string>
#include <iostream>
#include <cmath>
#include <system_error>
#include <map>
#include <utility>
#include <algorithm>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLine.h"
#include "TLatex.h"

#include "fit_pi0_mass.cxx"

namespace fs = std::filesystem;

// ---------- PNG saver for mass-slice PNGs ----------
void SaveHistPNG(TH1* h, const std::string& outdir) {
    if (!h) return;
    fs::create_directories(outdir);

    gStyle->SetOptFit(111); // show fit box

    TCanvas c((std::string("c_") + h->GetName()).c_str(), "", 900, 700);
    c.SetGrid();
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->Draw("E1");

    TLegend leg(0.60, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h, h->GetName(), "lep");
    if (auto* fcry = h->GetFunction((std::string("fcry_") + h->GetName()).c_str()))
        leg.AddEntry(fcry, "Gaussian prefit", "l");
    if (auto* bkg  = h->GetFunction((std::string("bkg_") + h->GetName()).c_str()))
        leg.AddEntry(bkg,  "Background (pol3)", "l");
    if (auto* fc   = h->GetFunction((std::string("fc_") + h->GetName()).c_str()))
        leg.AddEntry(fc,   "Gaus + pol3 (final)", "l");
    leg.Draw();

    const std::string png = outdir + "/" + std::string(h->GetName()) + ".png";
    c.SaveAs(png.c_str());
}

// ---------- 3D → Z-slices ----------
struct ZSlice {
    int ix, iy;               // 1-based raw bin indices
    int cx, cy;               // integerized bin centers (rounded)
    std::unique_ptr<TH1D> h;  // owned slice
};

std::vector<ZSlice> makeZSlices(TH3D* h3) {
    std::vector<ZSlice> out;
    if (!h3) return out;

    const int nx = h3->GetNbinsX();
    const int ny = h3->GetNbinsY();

    if (h3->GetSumw2N() == 0) h3->Sumw2(kTRUE);

    // exclude under/overflow
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int cx = static_cast<int>(std::lround(h3->GetXaxis()->GetBinCenter(ix)));
            const int cy = static_cast<int>(std::lround(h3->GetYaxis()->GetBinCenter(iy)));

            std::string name  = std::string(h3->GetName()) +
                                "_Z_x" + std::to_string(ix) +
                                "_y" + std::to_string(iy) +
                                "_xy_centers_" + std::to_string(cx) + "_" + std::to_string(cy);
            std::string title = std::string(h3->GetTitle()) + "; " + h3->GetZaxis()->GetTitle();

            TH1D* raw = h3->ProjectionZ(name.c_str(), ix, ix, iy, iy);
            raw->SetDirectory(nullptr);
            raw->SetTitle(title.c_str());

            out.push_back(ZSlice{ix, iy, cx, cy, std::unique_ptr<TH1D>(raw)});
        }
    }
    return out;
}

// ---------- helper: get desired TH3 ----------
TH3D* GetSingleTH3D(TFile* f,
                    const std::string& preferred = "h3_binX_zpt_pi0m") {
    if (!f) return nullptr;

    // 1) try preferred name
    if (!preferred.empty()) {
        if (auto* h = dynamic_cast<TH3D*>(f->Get(preferred.c_str())))
            return h;
    }

    // 2) fallback: first TH3D in the file (top level)
    TIter nextKey(f->GetListOfKeys());
    while (auto* key = static_cast<TKey*>(nextKey())) {
        if (auto* cls = TClass::GetClass(key->GetClassName())) {
            if (cls->InheritsFrom(TH3D::Class())) {
                return dynamic_cast<TH3D*>(key->ReadObj());
            }
        }
    }
    return nullptr;
}

// --- helper: directory of this script, fallback to CWD ---
static fs::path script_dir() {
    fs::path here = fs::path(__FILE__).parent_path();
    if (!here.empty() && fs::exists(here)) return here;
    return fs::current_path();
}

// ---------------- fitter ----------------
// save_pngs=true will save PNGs under "<script_dir>/plots/"
void split_and_fit_data(const std::string& path_to_h3_file, bool save_pngs=true) {
    fs::path in_path{path_to_h3_file};
    const std::string stem   = in_path.stem().string();       // input file w/o .root
    const std::string folder = in_path.parent_path().string();

    std::unique_ptr<TFile> f_in(TFile::Open(in_path.string().c_str(), "READ"));
    if (!f_in || f_in->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + in_path.string());
    }

    TH3D* h3 = GetSingleTH3D(f_in.get()); // try "h3_binX_zpt_pi0m", then first TH3D
    if (!h3) {
        throw std::runtime_error("No TH3D found in file: " + in_path.string());
    }

    // Prepare output ROOT (tree only)
    const std::string out_root = folder + '/' + stem + "_fitted.root";
    TFile fout(out_root.c_str(), "RECREATE");

    // PNGs go to <script_dir>/plots/
    const std::string out_pngs_folder = (script_dir() / "plots").string();
    if (save_pngs) {
        std::error_code ec;
        fs::create_directories(out_pngs_folder, ec);
        if (ec) {
            std::cerr << "Failed to create directory " << out_pngs_folder
                      << ": " << ec.message() << '\n';
        }
    }

    // Output tree (NO *_gen branches)
    fout.cd();
    TTree tout("h22_fit", "Gaussian fit per (xq2bin, zpt2phi)");
    int    out_xq2{};                // from X-axis center (bin_xBQ2_Valerii)
    int    out_zpt2phi{};            // from Y-axis center (zpt2phit_8x8x9)
    int    out_zpt2phi_bin{};        // raw Y bin index (iy)
    double out_nPions{}, out_errPions{};

    tout.Branch("xq2bin",             &out_xq2);
    tout.Branch("z_pt2_phi_bin",      &out_zpt2phi);
    tout.Branch("z_pt2_phi_hist_bin", &out_zpt2phi_bin);
    tout.Branch("nPions",             &out_nPions);
    tout.Branch("errPions",           &out_errPions);

    // Slice along Z and fit each (x,y) slice
    auto slices = makeZSlices(h3);

    int png_counter = 0;
    for (auto& s : slices) {
        if (!s.h) continue;

        // quick stat filter
        if (SkipHist_lowStat(s.h.get())) continue;

        auto [nPions, errPions] = GetN_pions(s.h.get());
        if (nPions < 0 || errPions < 0) continue;

        out_nPions        = nPions;
        out_errPions      = errPions;
        out_xq2           = s.cx;     // X-axis center (bin_xBQ2_Valerii)
        out_zpt2phi       = s.cy;     // Y-axis center (zpt2phit_8x8x9)
        out_zpt2phi_bin   = s.iy;     // raw Y bin index

        tout.Fill();

        // Save PNGs under ./plots/
        if (save_pngs && (png_counter % 100 == 0)) SaveHistPNG(s.h.get(), out_pngs_folder);
        ++png_counter;
    }

    fout.cd();
    tout.Write();
    fout.Close();

    std::cout << "Wrote: " << out_root << "\n";
}

// --- read xq2bin -> (ybin -> (nPions, errPions)) from a fitted file ---
static std::map<int, std::map<int, std::pair<double,double>>>
read_fit_tree(const fs::path& fitted_root) {
    std::map<int, std::map<int, std::pair<double,double>>> out;
    std::unique_ptr<TFile> f(TFile::Open(fitted_root.string().c_str(), "READ"));
    if (!f || f->IsZombie()) {
        throw std::runtime_error("Cannot open fitted file: " + fitted_root.string());
    }
    TTree* t = nullptr;
    f->GetObject("h22_fit", t);
    if (!t) throw std::runtime_error("Tree 'h22_fit' not found in " + fitted_root.string());

    int    xq2bin{}, zpt2phi_hist_bin{};
    double nPions{}, errPions{};
    t->SetBranchAddress("xq2bin",             &xq2bin);
    t->SetBranchAddress("z_pt2_phi_hist_bin", &zpt2phi_hist_bin);
    t->SetBranchAddress("nPions",             &nPions);
    t->SetBranchAddress("errPions",           &errPions);

    const Long64_t nent = t->GetEntries();
    for (Long64_t i = 0; i < nent; ++i) {
        t->GetEntry(i);
        out[xq2bin][zpt2phi_hist_bin] = {nPions, errPions};
    }
    return out;
}

// ---------- draw helpers ----------
static void SaveRatioPlot(const TH1D& hRatio, int xq2bin, const std::string& outdir, double yMin=0.0, double yMax=2.0) {
    fs::create_directories(outdir);
    TCanvas c(("c_ratio_x"+std::to_string(xq2bin)).c_str(), "", 900, 700);
    c.SetGrid();
    TH1D h = hRatio; // local copy to tweak style
    h.SetLineWidth(2);
    h.SetMarkerStyle(20);
    h.SetMarkerSize(0.9);
    h.SetMinimum(yMin);
    h.SetMaximum(yMax);
    h.GetYaxis()->SetTitle("MC / Data");
    h.Draw("E1");

    TLine l( h.GetXaxis()->GetXmin(), 1.0, h.GetXaxis()->GetXmax(), 1.0);
    l.SetLineStyle(2);
    l.Draw();

    TLatex tx; tx.SetNDC(true); tx.SetTextSize(0.035);
    tx.DrawLatex(0.14,0.92,Form("xq2bin = %d  (bins with rel. unc. \\le 40%% and nPions>20 in Data & MC)", xq2bin));

    const std::string png = outdir + "/ratio_mc_over_data_xq2bin_" + std::to_string(xq2bin) + ".png";
    c.SaveAs(png.c_str());
}

static void SaveDataMCPlot(const TH1D& hDataIn, const TH1D& hMCIn, int xq2bin, const std::string& outdir) {
    fs::create_directories(outdir);
    TCanvas c(("c_datamc_x"+std::to_string(xq2bin)).c_str(), "", 900, 700);
    c.SetGrid();

    TH1D hD = hDataIn, hM = hMCIn; // copies to style safely
    hD.SetLineWidth(2);
    hD.SetMarkerStyle(20);
    hD.SetMarkerSize(0.9);

    hM.SetLineWidth(2);
    hM.SetLineColor(kRed);
    hM.SetMarkerColor(kRed);

    // autoscale y
    double ymax = std::max(hD.GetMaximum(), hM.GetMaximum());
    hD.SetMinimum(0);
    hD.SetMaximum(ymax * 1.35);

    hD.GetYaxis()->SetTitle("nPions");
    hD.GetXaxis()->SetTitle("z_{pt^{2}}#phi histogram bin");

    hD.Draw("E1");
    hM.Draw("HIST SAME");

    TLegend leg(0.60, 0.72, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&hD, "Data", "lep");
    leg.AddEntry(&hM, "MC",   "l");
    leg.Draw();

    TLatex tx; tx.SetNDC(true); tx.SetTextSize(0.035);
    tx.DrawLatex(0.14,0.92,Form("xq2bin = %d", xq2bin));

    const std::string png = outdir + "/datamc_xq2bin_" + std::to_string(xq2bin) + ".png";
    c.SaveAs(png.c_str());
}

// --- main driver: run fits on both h3 files, then save ratio histograms + PNGs ---
void ratio_mc_data() {
    const fs::path basedir = script_dir();
    const fs::path data_h3 = basedir / "h3_data.root";
    const fs::path mc_h3   = basedir / "h3_mc.root";

    // 1) Fit both and save PNGs under ./plots/
    split_and_fit_data(data_h3.string(), /*save_pngs=*/true);
    split_and_fit_data(mc_h3.string(),   /*save_pngs=*/true);

    // 2) Read back the fitted trees
    const fs::path data_fit = basedir / "h3_data_fitted.root";
    const fs::path mc_fit   = basedir / "h3_mc_fitted.root";
    auto D = read_fit_tree(data_fit);
    auto M = read_fit_tree(mc_fit);

    // 3) Create ratio histograms per xq2bin (MC/Data), dropping bins if:
    //    - relative uncertainty on the ratio exceeds 40%, OR
    //    - either Data or MC has ≤ 20 nPions in that bin.
    TFile fout((basedir / "pi0_nPions_ratio.root").string().c_str(), "RECREATE");

    const std::string plots_dir_ratio  = (basedir / "plots" / "ratios").string();
    const std::string plots_dir_datamc = (basedir / "plots" / "datamc").string();
    fs::create_directories(plots_dir_ratio);
    fs::create_directories(plots_dir_datamc);

    // union of xq2bin keys
    std::map<int, int> all_x;
    for (auto& kv : D) all_x[kv.first] = 1;
    for (auto& kv : M) all_x[kv.first] = 1;

    const double relUncMax = 0.4;   // 40%
    const double minNPions = 20.0;  // require > 20 in both Data and MC

    for (auto& xkv : all_x) {
        const int xb = xkv.first;
        const auto& Dy = D[xb];
        const auto& My = M[xb];

        // determine max y bin across data & MC for this x
        int maxY = 0;
        for (auto& kv : Dy) maxY = std::max(maxY, kv.first);
        for (auto& kv : My) maxY = std::max(maxY, kv.first);
        if (maxY <= 0) continue;

        // Histograms for ratio and for Data/MC counts overlay
        const std::string hnameRatio  = "ratio_mc_over_data_xq2bin_" + std::to_string(xb);
        const std::string htitleRatio = "MC/Data; z_{pt^{2}}#phi histogram bin; ratio";
        TH1D hRatio(hnameRatio.c_str(), htitleRatio.c_str(), maxY, 0.5, maxY + 0.5);
        hRatio.Sumw2(true);

        const std::string hnameD = "data_nPions_xq2bin_" + std::to_string(xb);
        const std::string hnameM = "mc_nPions_xq2bin_"   + std::to_string(xb);
        TH1D hData(hnameD.c_str(), "; z_{pt^{2}}#phi histogram bin; nPions", maxY, 0.5, maxY + 0.5);
        TH1D hMC  (hnameM.c_str(), "; z_{pt^{2}}#phi histogram bin; nPions", maxY, 0.5, maxY + 0.5);
        hData.Sumw2(true);
        hMC.Sumw2(true);

        for (int yb = 1; yb <= maxY; ++yb) {
            const auto Dit = Dy.find(yb);
            const auto Mit = My.find(yb);
            const double Dval = (Dit != Dy.end()) ? Dit->second.first  : 0.0;
            const double Derr = (Dit != Dy.end()) ? Dit->second.second : 0.0;
            const double Mval = (Mit != My.end()) ? Mit->second.first  : 0.0;
            const double Merr = (Mit != My.end()) ? Mit->second.second : 0.0;

            // fill Data/MC yield hists
            hData.SetBinContent(yb, Dval);
            hData.SetBinError  (yb, Derr);
            hMC.SetBinContent  (yb, Mval);
            hMC.SetBinError    (yb, Merr);

            // Ratio MC/Data: fill only if both have > 20 nPions and relative uncertainty <= 40%
            if (Mval > minNPions && Dval > minNPions) {
                const double ratio = Mval / Dval;
                const double rel2  = (Merr > 0.0 ? (Merr*Merr)/(Mval*Mval) : 0.0)
                                   + (Derr > 0.0 ? (Derr*Derr)/(Dval*Dval) : 0.0);
                const double rel   = std::sqrt(rel2);      // relative uncertainty
                const double er    = ratio * rel;          // absolute uncertainty
                if (rel <= relUncMax) {
                    hRatio.SetBinContent(yb, ratio);
                    hRatio.SetBinError  (yb, er);
                } else {
                    hRatio.SetBinContent(yb, 0.0);
                    hRatio.SetBinError  (yb, 0.0);
                }
            } else {
                hRatio.SetBinContent(yb, 0.0);
                hRatio.SetBinError  (yb, 0.0);
            }
        }

        // Save ratio PNG and Data-vs-MC PNG
        SaveRatioPlot(hRatio, xb, plots_dir_ratio, 0.0, 2.0);
        SaveDataMCPlot(hData, hMC, xb, plots_dir_datamc);

        // Store ratio histogram to ROOT file
        hRatio.Write();
    }

    fout.Close();
    std::cout << "[OK] Wrote ratio histograms to: "
              << (basedir / "pi0_nPions_ratio.root") << "\n";
    std::cout << "[OK] Saved PNGs to: " << (basedir / "plots") << "\n";
}

// To run in ROOT:
//   root -l -q ratio_mc_data.cxx
// then in the ROOT prompt:
//   ratio_mc_data();
