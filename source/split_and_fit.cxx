// It takes TH3D from define_... splits them into TH1D
// Fits pi0 mass with bkg + gaus and extracts sigma per (x,y) bin
// Designed for response files; can be reused for Data/Rec
// All bins are fitted; later migration uses bins 17–20

#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <regex>
#include <cmath>
#include <limits>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TCollection.h"   // TIter
#include "TH3D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"


#include "fit_pi0_mass.cxx"

namespace fs = std::filesystem;


// ---------------- utilities ----------------
std::vector<std::string> list_TH3D_in_current_dir() {
    std::vector<std::string> names;
    TDirectory* dir = gDirectory;
    if (!dir) return names;

    TIter next(dir->GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        if (auto* cls = TClass::GetClass(key->GetClassName())) {
            if (cls->InheritsFrom(TH3D::Class())) {
                names.emplace_back(key->GetName());
            }
        }
    }
    return names;
}

// extract 2nd & 3rd integers from a string (handles negatives)
std::pair<int,int> second_third_int_ignoreQ(const std::string& s) {
    // Match an integer that is either at start (^) or preceded by a non-Q/q char.
    // Capture group 1 holds the integer text we want.
    static const std::regex re(R"((?:^|[^Qq])(-?\d+))");

    std::sregex_iterator it(s.begin(), s.end(), re), end;
    int idx = 0;                 // 0 = first, 1 = second, 2 = third (after ignoring Q* numbers)
    int second = 0, third = 0;

    for (; it != end; ++it, ++idx) {
        const std::string num = (*it)[1].str();   // <-- use capture group 1
        if (idx == 1) second = std::stoi(num);
        if (idx == 2) { third = std::stoi(num); return {second, third}; }
    }
    throw std::runtime_error("Not enough integers (ignoring Q* numbers) in: " + s);
}

// --------- NEW: PNG saver ----------
void SaveHistPNG(TH1* h, const std::string& outdir) {
    if (!h) return;
    fs::create_directories(outdir);

    // (Optional) show fit params box on the canvas
    gStyle->SetOptFit(111);

    // Canvas per hist
    TCanvas c((std::string("c_") + h->GetName()).c_str(), "", 900, 700);
    c.SetGrid();

    // Style for points/bars
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);

    h->Draw("E1");       // draws the hist; ROOT will also paint attached functions

    // If you want a legend:
    TLegend leg(0.60, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h, h->GetName(), "lep");
    // Try to find your functions by name pattern:
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


// ---------------- 3D → Z-slices ----------------
struct ZSlice {
    int ix, iy;                 // 1-based bin indices (0/N+1 if you include UF/OF)
    int cx, cy;              // bin centers (can be outside axis if UF/OF)
    std::unique_ptr<TH1D> h;    // owned slice (nullptr if saved to a directory)
};

std::vector<ZSlice> makeZSlices(TH3D* h3, TDirectory* outDir = nullptr) {
    std::vector<ZSlice> out;
    if (!h3) return out;

    const int nx = h3->GetNbinsX();
    const int ny = h3->GetNbinsY();

    if (h3->GetSumw2N() == 0) h3->Sumw2(kTRUE);

    // choose whether to include under/overflow (0 and N+1). Here: EXclude them.
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
            raw->SetTitle(title.c_str());

            if (outDir) {
                raw->SetDirectory(outDir);     // file/dir owns it
                out.push_back(ZSlice{ix, iy, cx, cy, nullptr});
            } else {
                raw->SetDirectory(nullptr);    // we own it
                out.push_back(ZSlice{ix, iy, cx, cy, std::unique_ptr<TH1D>(raw)});
            }
        }
    }
    return out;
}

// ---------------- main driver ----------------
void split_and_fit(const std::string& path_to_xQ2bin_th3d) {
    fs::path file_in_path{path_to_xQ2bin_th3d};
    std::string stem = file_in_path.stem().string();
    std::string folder_in = file_in_path.parent_path().string();
  
    const std::string dir_with_hists = (stem.size() > 3) ? stem.substr(3) : std::string{};

    std::unique_ptr<TFile> f_in(TFile::Open(file_in_path.string().c_str(), "READ"));
    if (!f_in || f_in->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + file_in_path.string());
    }

    // locate and cd to the directory that holds the TH3Ds
    TDirectory* dir = f_in->GetDirectory(dir_with_hists.c_str(), /*printError=*/false);
    if (!dir) throw std::runtime_error("no dir: " + dir_with_hists);
    dir->cd();
    TDirectory* inDir = dir; 

    const auto th3_names = list_TH3D_in_current_dir();
    // output file + png directory
    const std::string outPath = folder_in + '/' + stem + "_fitted.root";
    TFile fout(outPath.c_str(), "RECREATE");
    fout.cd();
    TTree tout("h22_fit", "Gaussian fit per (xq2bin,xq2bin_gen,zpt2phi,zpt2phi_gen)");
    int    out_xq2{}, out_xq2g{}, out_zpt2phi{}, out_zpt2phi_gen{};
    int    out_zpt2phi_bin{}, out_zpt2phi_bin_gen{};
    double out_nPions{}, out_errPions{};
  
    tout.Branch("xq2bin",                 &out_xq2);
    tout.Branch("xq2bin_gen",             &out_xq2g);
    tout.Branch("z_pt2_phi_bin",          &out_zpt2phi);
    tout.Branch("z_pt2_phi_bin_gen",      &out_zpt2phi_gen);
    tout.Branch("z_pt2_phi_hist_bin",     &out_zpt2phi_bin);
    tout.Branch("z_pt2_phi_hist_bin_gen", &out_zpt2phi_bin_gen);
    tout.Branch("nPions",             &out_nPions);
    tout.Branch("errPions",             &out_errPions);

    // Optional: group PNGs per input TH3 name to keep things tidy
    // (e.g., png_only_out/<th3name>/...)

    // path to PNGs:
    std::string out_pngs_folder = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/pi0_mass_fits/" + file_in_path.parent_path().filename().string() + '/';


    std::error_code ec;
    fs::create_directories(out_pngs_folder, ec); // creates parents too
    if (ec) {
        std::cerr << "Failed to create directory " << out_pngs_folder
                  << ": " << ec.message() << '\n';
    }
  
    for (const auto& name : th3_names) {
        TH3D* h3 = nullptr;
        inDir->GetObject(name.c_str(), h3);
        if (!h3) { std::cerr << "Could not read TH3D: " << name << '\n'; continue; }

        auto [rec_xq2, gen_xq2] = second_third_int_ignoreQ(name);
        if (rec_xq2 == -5 || gen_xq2 == -5) continue;

        auto slices = makeZSlices(h3);  // we own the TH1Ds here



      
        const std::string pngSubdir = out_pngs_folder + name;  // per-TH3 subfolder
        fs::create_directories(pngSubdir);

        int nCounter = 0;
        for (auto& s : slices) {
            if (!s.h) continue;

            // check integrals is there is enough stat to fit:
            if (SkipHist_lowStat(s.h.get())) continue;
  

            auto [nPions, errPions] = GetN_pions(s.h.get());  // attaches functions to s.h

            if (nPions< 0 || errPions < 0) continue;

            out_nPions = nPions;
            out_errPions = errPions;

            // centers (may be outside axis for UF/OF)
            out_zpt2phi     = static_cast<int>(std::lround(s.cx));
            out_zpt2phi_gen = static_cast<int>(std::lround(s.cy));

            // raw bin indices
            out_zpt2phi_bin     = s.ix;
            out_zpt2phi_bin_gen = s.iy;

            out_xq2  = rec_xq2;
            out_xq2g = gen_xq2;

            tout.Fill();

            // ---- save PNG for this slice (includes attached functions) ----
            cout << rec_xq2 << " "  << gen_xq2 << " "<< name<<endl;
            if (nCounter%2 == 0) SaveHistPNG(s.h.get(), pngSubdir);
            nCounter++;
        }
    }

    fout.cd();
    tout.Write();
    fout.Close();
}
