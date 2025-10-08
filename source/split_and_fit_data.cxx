// It takes the single TH3D from define_bin_migr_hists_data output,
// splits it into TH1D Z-slices, fits pi0 mass, and stores results in a TTree.
// Branches *do not* include any *_gen variables.

#include <filesystem>
#include <memory>
#include <string>
#include <iostream>
#include <cmath>
#include <system_error>

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

#include "fit_pi0_mass.cxx"

namespace fs = std::filesystem;

// ---------- PNG saver ----------
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

// ---------------- main driver ----------------
void split_and_fit_data(const std::string& path_to_h3_file) {
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

    // Prepare output ROOT + PNG directory
    const std::string out_root = folder + '/' + stem + "_fitted.root";
    TFile fout(out_root.c_str(), "RECREATE");

    std::string out_pngs_folder =
        "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/pi0_mass_fits/" +
        in_path.parent_path().filename().string() + '/';
    std::error_code ec;
    fs::create_directories(out_pngs_folder, ec);
    if (ec) {
        std::cerr << "Failed to create directory " << out_pngs_folder
                  << ": " << ec.message() << '\n';
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

        // save every 3rd PNG to limit I/O (adjust as you like)
        if (png_counter % 5 == 0) SaveHistPNG(s.h.get(), out_pngs_folder);
        ++png_counter;
    }

    fout.cd();
    tout.Write();
    fout.Close();

    std::cout << "Wrote: " << out_root << "\n";
}
