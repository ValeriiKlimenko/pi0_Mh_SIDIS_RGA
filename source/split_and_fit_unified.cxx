// split_and_fit_unified.cxx
//
// Unified driver that merges the "data" workflow from split_and_fit_data.cxx
// and the "simulation" workflow from split_and_fit.cxx, while using the
// latest fit utilities from fit_pi0_mass.cxx.
//
// - Logic::Data reproduces the I/O and behavior of split_and_fit_data.cxx
// - Logic::Sim  reproduces the I/O and behavior of split_and_fit.cxx,
//               but updated to use PrecheckHistogram + FitPi0Mass and now
//               saves failed fits into a per-TH3 "failed/" subfolder.
//
// Notes for plots:
//   • Before saving, we re-synchronize the drawables so pictures show the
//     final combined-fit gaussian (first 3 params of "fc_*") and the
//     final background (params 3..N of "fc_*").
//   • A small text box is drawn with the synchronized gaussian width σ.
//
// Backward-compat wrappers are provided:
//   void split_and_fit_data(const std::string& path);
//   void split_and_fit     (const std::string& path);   // simulation

#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <regex>
#include <cmath>
#include <limits>
#include <system_error>
#include <iomanip>
#include <sstream>
#include <map>              // NEW: per-xQ2 counters

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TCollection.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TPaveText.h"

#include "fit_pi0_mass.cxx"  // SkipDecision, PrecheckHistogram, FitResult, FitPi0Mass



// --- Z⊗pT2⊗phi binning constants from  scheme ---
static constexpr int kNphi = 8;                  // N_phiTrbins
static constexpr int kNpt2_with_overflow = 11;   // N_pTbins_with_overflow (= 10 + overflow)

// Decode Y-axis combined bin (cy = 1..(N_zpt2*kNphi)) into (zbin, pt2bin, phibin)
static inline void DecodeZpt2Phi(int cy, int& zbin, int& pt2bin, int& phibin) {
    if (cy <= 0) { zbin = pt2bin = phibin = 0; return; }
    const int zpt2 = ((cy - 1) / kNphi) + 1;   // 1..(N_Zbins * kNpt2_with_overflow)
    phibin         = ((cy - 1) % kNphi) + 1;   // 1..kNphi
    zbin           = ((zpt2 - 1) / kNpt2_with_overflow) + 1; // 1..N_Zbins (=8)
    pt2bin         = ((zpt2 - 1) % kNpt2_with_overflow) + 1; // 1..kNpt2_with_overflow (=11)
}




constexpr int kMaxFails = std::numeric_limits<int>::max();

namespace fs = std::filesystem;

// ---------------------- PNG saver ----------------------
// Ensures the drawn gaussian/background reflect the *final synchronized* fit:
//   - Copies fc_*[0..2] → fcry_*[0..2]
//   - Copies fc_*[3..]  → bkg_*[0..]
// Also overlays a text box with σ from the synchronized gaussian.
// ---------------------- PNG saver ----------------------
static void SaveHistPNG(TH1* h, const std::string& outdir) {
    if (!h) return;
    std::error_code ec;
    fs::create_directories(outdir, ec);

    const std::string hname = h->GetName();

    // Prefer gaussian extracted from fc; fall back to legacy "fcry_*"
    TF1* fgaus = h->GetFunction((std::string("gaus_from_fc_") + hname).c_str());
    if (!fgaus) fgaus = h->GetFunction((std::string("fcry_") + hname).c_str());
    TF1* bkg  = h->GetFunction((std::string("bkg_")  + hname).c_str());
    TF1* fc   = h->GetFunction((std::string("fc_")   + hname).c_str());

    // --- Synchronize drawables with the final combined fit (no refits) ---
    if (fc) {
        if (fgaus && fgaus->GetNpar() >= 3 && fc->GetNpar() >= 3) {
            for (int i = 0; i < 3; ++i) fgaus->SetParameter(i, fc->GetParameter(i));
            fgaus->SetParName(2, "Sigma"); // label for clarity
            fgaus->SetRange(fc->GetXmin(), fc->GetXmax());
        }
        if (bkg) {
            const int npar_bkg = bkg->GetNpar();
            for (int i = 0; i < npar_bkg; ++i)
                if (3 + i < fc->GetNpar()) bkg->SetParameter(i, fc->GetParameter(3 + i));
            bkg->SetRange(fc->GetXmin(), fc->GetXmax());
        }
    }

    // Draw
    TCanvas c((std::string("c_") + hname).c_str(), "", 900, 700);
    c.SetGrid();

    gStyle->SetOptStat(0);
    h->SetStats(kFALSE);

    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->Draw("E1");

    if (bkg)  bkg->Draw("SAME");
    if (fc)   fc->Draw("SAME");
    if (fgaus) fgaus->Draw("SAME");

    // Legend
    TLegend leg(0.60, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h, hname.c_str(), "lep");
    if (fgaus) leg.AddEntry(fgaus, "Gaussian (from gaus(0)+)", "l");
    if (bkg) {
        const int npar = bkg->GetNpar();
        const int order = std::max(0, npar - 1);
        const std::string lbl = (order > 0)
          ? ("Background (pol" + std::to_string(order) + ", synced)")
          : "Background (synced)";
        leg.AddEntry(bkg, lbl.c_str(), "l");
    }
    if (fc) {
        const int order = std::max(0, fc->GetNpar() - 3);
        leg.AddEntry(fc, ("Gaus + pol" + std::to_string(order) + " (final)").c_str(), "l");
    }
    leg.Draw();

    // ---------- Draw ONLY σ from the final gaus(0) ----------
    double sigma = std::numeric_limits<double>::quiet_NaN();
    if (fc && fc->GetNpar() >= 3) {
        sigma = std::abs(fc->GetParameter(2));     // preferred
    } else if (fgaus && fgaus->GetNpar() >= 3) {
        sigma = std::abs(fgaus->GetParameter(2));  // fallback
    }

    if (std::isfinite(sigma)) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4) << "#sigma = " << sigma;

        // Place left to avoid the legend
        TPaveText* pave = new TPaveText(0.12, 0.78, 0.50, 0.92, "NDC");
        pave->SetName((std::string("sigma_box_") + hname).c_str());
        pave->SetFillColorAlpha(kWhite, 0.90);
        pave->SetFillStyle(1001);
        pave->SetBorderSize(1);
        pave->SetLineColor(kBlack);
        pave->SetTextAlign(12);
        pave->SetTextFont(42);
        pave->SetTextColor(kBlack);
        pave->SetTextSize(0.055);
        pave->AddText(ss.str().c_str());
        pave->Draw();
    }

    c.Modified();
    c.Update();
    const std::string png = outdir + "/" + hname + ".png";
    c.SaveAs(png.c_str());
}

// ---------------------- slice helper ----------------------
struct ZSlice {
    int ix, iy;               // 1-based raw bin indices
    int cx, cy;               // integerized bin centers (rounded)
    std::unique_ptr<TH1D> h;  // owned slice
};

static std::vector<ZSlice> makeZSlices(TH3D* h3) {
    std::vector<ZSlice> out;
    if (!h3) return out;

    const int nx = h3->GetNbinsX();
    const int ny = h3->GetNbinsY();

    if (h3->GetSumw2N() == 0) h3->Sumw2(kTRUE);

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

// ---------------------- TH3D discovery helpers ----------------------
static TH3D* GetSingleTH3D(TFile* f,
                           const std::string& preferred = "h3_binX_zpt_pi0m") {
    if (!f) return nullptr;

    if (!preferred.empty()) {
        if (auto* h = dynamic_cast<TH3D*>(f->Get(preferred.c_str())))
            return h;
    }

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

static std::vector<std::string> list_TH3D_in_current_dir() {
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

// extract 2nd & 3rd integers from a string (ignores numbers in Q* tokens)
static std::pair<int,int> second_third_int_ignoreQ(const std::string& s) {
    static const std::regex re(R"((?:^|[^Qq])(-?\d+))");
    std::sregex_iterator it(s.begin(), s.end(), re), end;
    int idx = 0;
    int second = 0, third = 0;
    for (; it != end; ++it, ++idx) {
        const std::string num = (*it)[1].str();
        if (idx == 1) second = std::stoi(num);
        if (idx == 2) { third = std::stoi(num); return {second, third}; }
    }
    throw std::runtime_error("Not enough integers (ignoring Q* numbers) in: " + s);
}

// ---------------------- counters ----------------------
struct FitCounters {
    long long success{0};
    long long failed{0};
    long long skipped{0};
};

// ---------------------- unified driver ----------------------
enum class Logic { Data, Sim };

static void split_and_fit_unified(const std::string& path_to_root, Logic logic,
                                  int png_every = -1) {
    const bool isData = (logic == Logic::Data);
    fs::path in_path{path_to_root};
    const std::string stem   = in_path.stem().string();
    const std::string folder = in_path.parent_path().string();

    std::unique_ptr<TFile> f_in(TFile::Open(in_path.string().c_str(), "READ"));
    if (!f_in || f_in->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + in_path.string());
    }

    const std::string out_root = folder + '/' + stem + "_fitted.root";
    TFile fout(out_root.c_str(), "RECREATE");

    std::string out_pngs_base =
        "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/pi0_mass_fits/" +
        in_path.parent_path().filename().string() + '/';

    std::error_code ec;
    fs::create_directories(out_pngs_base, ec);
    if (ec) {
        std::cerr << "Failed to create directory " << out_pngs_base
                  << ": " << ec.message() << '\n';
    }

    fout.cd();
    TTree tout("h22_fit", isData
        ? "Gaussian fit per (xq2bin, zpt2phi) — DATA"
        : "Gaussian fit per (xq2bin,xq2bin_gen,zpt2phi,zpt2phi_gen) — SIM");

    // Common outputs
    double out_nPions{}, out_errPions{};

    // ====================== DATA path ======================
    if (isData) {
        int out_xq2{}, out_zpt2phi{}, out_zpt2phi_hist_bin{};
      
        // NEW: decoded bins from cy
        int out_zbin{}, out_pt2bin{};

        // Per-xQ2 counters
        std::map<int, FitCounters> per_xq2;
        long long total_success = 0, total_failed = 0, total_skipped = 0;

        tout.Branch("xq2bin",             &out_xq2);
        tout.Branch("z_pt2_phi_bin",      &out_zpt2phi);
        tout.Branch("z_pt2_phi_hist_bin", &out_zpt2phi_hist_bin);

        // NEW: write decoded bins to tree
        tout.Branch("zbin",               &out_zbin);
        tout.Branch("pt2bin",             &out_pt2bin);

      
        tout.Branch("nPions",             &out_nPions);
        tout.Branch("errPions",           &out_errPions);

        std::string out_pngs_failed_folder  = out_pngs_base + "failed/";
        std::string out_pngs_nopeak_folder  = out_pngs_base + "no_peak/";
        fs::create_directories(out_pngs_failed_folder, ec);
        fs::create_directories(out_pngs_nopeak_folder, ec);

        if (png_every <= 0) png_every = 30;

        TH3D* h3 = GetSingleTH3D(f_in.get());
        if (!h3) {
            throw std::runtime_error("No TH3D found in file: " + in_path.string());
        }

        auto slices = makeZSlices(h3);
        int png_counter   = 0;
        int failed_count  = 0; // global cap helper

        for (auto& s : slices) {
            if (!s.h) continue;

            const int xq2_key = s.cx; // use X-axis integerized center as xQ2 key

            SkipDecision dec = PrecheckHistogram(s.h.get());
            if (dec.skip) {
                // Count as SKIPPED (empty/low-stat/no-peak, etc.)
                ++per_xq2[xq2_key].skipped;
                ++total_skipped;

                std::cerr << "[skip] " << s.h->GetName() << " — " << dec.reason << '\n';
                if (dec.is_no_peak) {
                    SaveHistPNG(s.h.get(), out_pngs_nopeak_folder);
                } else if (!dec.is_empty) {
                    SaveHistPNG(s.h.get(), out_pngs_failed_folder);
                    ++failed_count;
                    if (failed_count >= kMaxFails) {
                        std::cerr << "[unified:data] Reached " << kMaxFails
                                  << " failed/low-stat fits. Aborting loop.\n";
                        break;
                    }
                }
                continue;
            }

            FitResult fr = FitPi0Mass(s.h.get());
            if (!fr.ok) {
                ++per_xq2[xq2_key].failed;
                ++total_failed;
                std::cerr << "[failed fit] " << s.h->GetName() << " — " << fr.reason << '\n';
                SaveHistPNG(s.h.get(), out_pngs_failed_folder);
                ++failed_count;
                if (failed_count >= kMaxFails) {
                    std::cerr << "[unified:data] Reached " << kMaxFails
                              << " failed/low-stat fits. Aborting loop.\n";
                    break;
                }
                continue;
            }

            // Success
            ++per_xq2[xq2_key].success;
            ++total_success;
            
            out_nPions           = fr.nPions;
            out_errPions         = fr.errPions;
            out_xq2              = s.cx;
            out_zpt2phi          = s.cy;
            out_zpt2phi_hist_bin = s.iy;
            
            // NEW: decode (zbin, pt2bin, phibin) from cy
            int phibin = 0;
            DecodeZpt2Phi(out_zpt2phi, out_zbin, out_pt2bin, phibin);
            
            // Fill tree with zbin/pt2bin
            tout.Fill();
            
            // NEW: save ALL successful φ fits, grouped by xq2 → (z,pt2)
            {
                const std::string xq2_dir   = out_pngs_base + "xq2_" + std::to_string(out_xq2) + "/";
                const std::string zpt2_dir  = xq2_dir + "z" + std::to_string(out_zbin)
                                              + "_pt2" + std::to_string(out_pt2bin) + "/";
                SaveHistPNG(s.h.get(), zpt2_dir);
            }
        }

        // --- Summary (DATA) ---
        std::cout << "\n=== Fit summary by xQ2 (DATA) ===\n";
        std::cout << "xQ2    success   failed   skipped   total\n";
        std::cout << "-----------------------------------------\n";
        for (const auto& kv : per_xq2) {
            const int x = kv.first;
            const auto& c = kv.second;
            std::cout << std::setw(3) << x << std::setw(11) << c.success
                      << std::setw(9)  << c.failed  << std::setw(10) << c.skipped
                      << std::setw(8)  << (c.success + c.failed + c.skipped) << "\n";
        }
        std::cout << "-----------------------------------------\n";
        std::cout << "TOT" << std::setw(10) << total_success
                  << std::setw(9)  << total_failed
                  << std::setw(10) << total_skipped
                  << std::setw(8)  << (total_success + total_failed + total_skipped) << "\n";

        fout.cd(); tout.Write(); fout.Close();
        std::cout << "Wrote: " << out_root << "\n";
        return;
    }

    // ====================== SIMULATION path (with global fail cap) ======================
    if (png_every <= 0) png_every = 10;

    const std::string dir_with_hists = (stem.size() > 3) ? stem.substr(3) : std::string{};
    TDirectory* dir = f_in->GetDirectory(dir_with_hists.c_str(), /*printError=*/false);
    if (!dir) throw std::runtime_error("no dir: " + dir_with_hists);
    dir->cd();
    TDirectory* inDir = dir;

    int out_xq2{}, out_xq2g{};
    int out_zpt2phi{}, out_zpt2phi_gen{};
    int out_zpt2phi_hist_bin{}, out_zpt2phi_hist_bin_gen{};

    tout.Branch("xq2bin",                 &out_xq2);
    tout.Branch("xq2bin_gen",             &out_xq2g);
    tout.Branch("z_pt2_phi_bin",          &out_zpt2phi);
    tout.Branch("z_pt2_phi_bin_gen",      &out_zpt2phi_gen);
    tout.Branch("z_pt2_phi_hist_bin",     &out_zpt2phi_hist_bin);
    tout.Branch("z_pt2_phi_hist_bin_gen", &out_zpt2phi_hist_bin_gen);
    tout.Branch("nPions",                 &out_nPions);
    tout.Branch("errPions",               &out_errPions);

    const auto th3_names = list_TH3D_in_current_dir();

    // Per-xQ2 counters (reconstructed xQ2)
    std::map<int, FitCounters> per_xq2;
    long long total_success = 0, total_failed = 0, total_skipped = 0;

    int failed_count = 0;

    for (const auto& name : th3_names) {
        if (failed_count >= kMaxFails) {
            std::cerr << "[unified:sim] Stop processing further TH3s due to fail cap.\n";
            break;
        }

        TH3D* h3 = nullptr;
        inDir->GetObject(name.c_str(), h3);
        if (!h3) { std::cerr << "Could not read TH3D: " << name << '\n'; continue; }

        int rec_xq2 = 0, gen_xq2 = 0;
        try {
            auto p = second_third_int_ignoreQ(name);
            rec_xq2 = p.first;
            gen_xq2 = p.second;
        } catch (const std::exception& e) {
            std::cerr << "[warn] name parse failed for " << name << ": " << e.what() << "\n";
            continue;
        }
        if (rec_xq2 == -5 || gen_xq2 == -5) continue;

        // Per-TH3 OK/FAILED subfolders
        const std::string pngSubdirOK     = out_pngs_base + name + "/";
        const std::string pngSubdirFailed = pngSubdirOK + "failed/";
        fs::create_directories(pngSubdirOK, ec);
        fs::create_directories(pngSubdirFailed, ec);

        auto slices = makeZSlices(h3);
        int png_counter = 0;

        for (auto& s : slices) {
            if (failed_count >= kMaxFails) {
                std::cerr << "[unified:sim] Reached " << kMaxFails
                          << " failed/low-stat fits. Aborting loop.\n";
                break;
            }
            if (!s.h) continue;

            const int xq2_key = rec_xq2;   // count by reconstructed xQ2

            // Precheck: save non-empty skips to failed/
            SkipDecision dec = PrecheckHistogram(s.h.get());
            if (dec.skip) {
                ++per_xq2[xq2_key].skipped;
                ++total_skipped;

                if (!dec.is_empty) {
                    SaveHistPNG(s.h.get(), pngSubdirFailed);
                    ++failed_count;
                    if (failed_count >= kMaxFails) {
                        std::cerr << "[unified:sim] Reached " << kMaxFails
                                  << " failed/low-stat fits. Aborting loop.\n";
                        break;
                    }
                }
                continue;
            }

            // Fit: save failures to failed/
            FitResult fr = FitPi0Mass(s.h.get());
            if (!fr.ok) {
                ++per_xq2[xq2_key].failed;
                ++total_failed;

                std::cerr << "[failed fit] " << s.h->GetName() << " — " << fr.reason << '\n';
                SaveHistPNG(s.h.get(), pngSubdirFailed);
                ++failed_count;
                if (failed_count >= kMaxFails) {
                    std::cerr << "[unified:sim] Reached " << kMaxFails
                              << " failed/low-stat fits. Aborting loop.\n";
                    break;
                }
                continue;
            }

            // Success → fill tree, save "ok" PNGs sparsely
            ++per_xq2[xq2_key].success;
            ++total_success;

            out_nPions               = fr.nPions;
            out_errPions             = fr.errPions;
            out_zpt2phi              = static_cast<int>(std::lround(s.cx));
            out_zpt2phi_gen          = static_cast<int>(std::lround(s.cy));
            out_zpt2phi_hist_bin     = s.ix;
            out_zpt2phi_hist_bin_gen = s.iy;
            out_xq2                  = rec_xq2;
            out_xq2g                 = gen_xq2;

            tout.Fill();

            if (png_counter % png_every == 0) SaveHistPNG(s.h.get(), pngSubdirOK);
            ++png_counter;
        }

        if (failed_count >= kMaxFails) {
            std::cerr << "[unified:sim] Stop processing further TH3s due to fail cap.\n";
            break;
        }

        // optional debug (kept from original)
        std::cout << rec_xq2 << " "  << gen_xq2 << " " << name << std::endl;
    }

    // --- Summary (SIM, reconstructed xQ2) ---
    std::cout << "\n=== Fit summary by xQ2 (SIM, reconstructed) ===\n";
    std::cout << "xQ2    success   failed   skipped   total\n";
    std::cout << "-----------------------------------------\n";
    for (const auto& kv : per_xq2) {
        const int x = kv.first;
        const auto& c = kv.second;
        std::cout << std::setw(3) << x << std::setw(11) << c.success
                  << std::setw(9)  << c.failed  << std::setw(10) << c.skipped
                  << std::setw(8)  << (c.success + c.failed + c.skipped) << "\n";
    }
    std::cout << "-----------------------------------------\n";
    std::cout << "TOT" << std::setw(10) << total_success
              << std::setw(9)  << total_failed
              << std::setw(10) << total_skipped
              << std::setw(8)  << (total_success + total_failed + total_skipped) << "\n";

    fout.cd();
    tout.Write();
    fout.Close();
    std::cout << "Wrote: " << out_root << "\n";
}

// ---------------------- public wrappers ----------------------
void split_and_fit_data(const std::string& path_to_h3_file) {
    split_and_fit_unified(path_to_h3_file, Logic::Data);
}

void split_and_fit(const std::string& path_to_xQ2bin_th3d) {
    split_and_fit_unified(path_to_xQ2bin_th3d, Logic::Sim);
}

// Optional switch with explicit flag; png_every<=0 uses defaults (30 data, 10 sim)
void split_and_fit_switch(const std::string& path, const std::string& logic, int png_every = -1) {
    if (logic == "data" || logic == "DATA")
        split_and_fit_unified(path, Logic::Data, png_every);
    else if (logic == "sim" || logic == "SIM" || logic == "mc" || logic == "MC")
        split_and_fit_unified(path, Logic::Sim, png_every);
    else
        throw std::runtime_error("Unknown logic flag: " + logic + " (use \"data\" or \"sim\")");
}
