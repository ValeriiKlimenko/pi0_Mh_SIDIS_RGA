// slice_and_fit_selected_hists.cxx
//
// Default: runs on ALL three files in:
//   /lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/
//
// Files:
//   rec_data_hists.root       (TH3D -> slice Z -> FitPi0Mass)
//   rec_true_hists.root       (TH3D -> slice Z -> FitPi0Mass)
//   gen_binning_hists.root    (TH2D -> slice projections; no pi0 fit)
//
// Build:
//   g++ -O2 -std=c++17 `root-config --cflags --libs` slice_and_fit_selected_hists.cxx -o slice_and_fit_selected
//
// ROOT:
//   root -l
//   root [0] .L slice_and_fit_selected_hists.cxx+
//   root [1] run_slice_and_fit_all_selected();   // default: all 3
//
// Output (per input):
//   <same folder>/<stem>_sliced_fitted.root
//
// Requires:
//   #include "fit_pi0_mass.cxx"

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
#include <map>
#include <cstring>
#include <cstdlib>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TCollection.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TSystem.h"

#include "fit_pi0_mass.cxx"  // SkipDecision, PrecheckHistogram, FitResult, FitPi0Mass

namespace fs = std::filesystem;

// ---------------- HARD-CODED DEFAULT INPUTS ----------------
static const char* kDefaultHistsDir =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/";

static const char* kRecDataHists =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/rec_data_hists.root";

static const char* kRecTrueHists =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/rec_true_hists.root";

static const char* kGenBinningHists =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/gen_binning_hists.root";

// ---- scheme constants for decoding zpt2phi bin ----
static constexpr int kNphi = 8;
static constexpr int kNpt2_with_overflow = 11;

static inline void DecodeZpt2Phi(int cy, int& zbin, int& pt2bin, int& phibin) {
  if (cy <= 0) { zbin = pt2bin = phibin = 0; return; }
  const int zpt2 = ((cy - 1) / kNphi) + 1;
  phibin         = ((cy - 1) % kNphi) + 1;
  zbin           = ((zpt2 - 1) / kNpt2_with_overflow) + 1;
  pt2bin         = ((zpt2 - 1) % kNpt2_with_overflow) + 1;
}

static std::string BaseNameFromPath(const std::string& p) {
  auto pos = p.find_last_of('/');
  return (pos == std::string::npos) ? p : p.substr(pos + 1);
}

static TDirectory* MkdirP(TDirectory* top, const std::string& path) {
  if (!top) return nullptr;
  if (path.empty()) return top;
  TDirectory* cur = top;
  std::stringstream ss(path);
  std::string item;
  while (std::getline(ss, item, '/')) {
    if (item.empty()) continue;
    TDirectory* next = cur->GetDirectory(item.c_str());
    if (!next) next = cur->mkdir(item.c_str());
    cur = next;
  }
  return cur;
}

// ---------------------- PNG saver (optional) ----------------------
static void SaveHistPNG(TH1* h, const std::string& outdir) {
  if (!h) return;
  std::error_code ec;
  fs::create_directories(outdir, ec);

  const std::string hname = h->GetName();

  TF1* fgaus = h->GetFunction((std::string("gaus_from_fc_") + hname).c_str());
  if (!fgaus) fgaus = h->GetFunction((std::string("fcry_") + hname).c_str());
  TF1* bkg  = h->GetFunction((std::string("bkg_")  + hname).c_str());
  TF1* fc   = h->GetFunction((std::string("fc_")   + hname).c_str());

  if (fc) {
    if (fgaus && fgaus->GetNpar() >= 3 && fc->GetNpar() >= 3) {
      for (int i = 0; i < 3; ++i) fgaus->SetParameter(i, fc->GetParameter(i));
      fgaus->SetRange(fc->GetXmin(), fc->GetXmax());
    }
    if (bkg) {
      const int npar_bkg = bkg->GetNpar();
      for (int i = 0; i < npar_bkg; ++i)
        if (3 + i < fc->GetNpar()) bkg->SetParameter(i, fc->GetParameter(3 + i));
      bkg->SetRange(fc->GetXmin(), fc->GetXmax());
    }
  }

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

  double sigma = std::numeric_limits<double>::quiet_NaN();
  if (fc && fc->GetNpar() >= 3) sigma = std::abs(fc->GetParameter(2));
  else if (fgaus && fgaus->GetNpar() >= 3) sigma = std::abs(fgaus->GetParameter(2));

  if (std::isfinite(sigma)) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << "#sigma = " << sigma;
    TPaveText* pave = new TPaveText(0.12, 0.78, 0.50, 0.92, "NDC");
    pave->SetFillColorAlpha(kWhite, 0.90);
    pave->SetBorderSize(1);
    pave->SetTextFont(42);
    pave->SetTextSize(0.055);
    pave->AddText(ss.str().c_str());
    pave->Draw();
  }

  c.Modified(); c.Update();
  const std::string png = outdir + "/" + hname + ".png";
  c.SaveAs(png.c_str());
}

// ---------------------- slice structs ----------------------
struct ZSlice {
  int ix, iy;
  int cx, cy;
  std::unique_ptr<TH1D> h;
};

static std::vector<ZSlice> makeZSlices(TH3D* h3) {
  std::vector<ZSlice> out;
  if (!h3) return out;

  const int nx = h3->GetNbinsX();
  const int ny = h3->GetNbinsY();
  if (h3->GetSumw2N() == 0) h3->Sumw2(kTRUE);

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      const int cx = (int)std::lround(h3->GetXaxis()->GetBinCenter(ix));
      const int cy = (int)std::lround(h3->GetYaxis()->GetBinCenter(iy));

      std::string name =
        std::string(h3->GetName()) +
        "_Z_ix" + std::to_string(ix) +
        "_iy" + std::to_string(iy) +
        "_centers_" + std::to_string(cx) + "_" + std::to_string(cy);

      TH1D* raw = h3->ProjectionZ(name.c_str(), ix, ix, iy, iy);
      raw->SetDirectory(nullptr);
      out.push_back(ZSlice{ix, iy, cx, cy, std::unique_ptr<TH1D>(raw)});
    }
  }
  return out;
}

struct H1Slice {
  int ifix;
  int cfix;
  bool projX;
  std::unique_ptr<TH1D> h;
};

static std::vector<H1Slice> make2DSlices(TH2D* h2, bool doProjX_byY=true, bool doProjY_byX=true) {
  std::vector<H1Slice> out;
  if (!h2) return out;

  if (doProjX_byY) {
    const int ny = h2->GetNbinsY();
    for (int iy = 1; iy <= ny; ++iy) {
      const int cy = (int)std::lround(h2->GetYaxis()->GetBinCenter(iy));
      std::string name = std::string(h2->GetName()) + "_PX_iy" + std::to_string(iy) + "_cy" + std::to_string(cy);
      TH1D* raw = h2->ProjectionX(name.c_str(), iy, iy);
      raw->SetDirectory(nullptr);
      out.push_back(H1Slice{iy, cy, true, std::unique_ptr<TH1D>(raw)});
    }
  }

  if (doProjY_byX) {
    const int nx = h2->GetNbinsX();
    for (int ix = 1; ix <= nx; ++ix) {
      const int cx = (int)std::lround(h2->GetXaxis()->GetBinCenter(ix));
      std::string name = std::string(h2->GetName()) + "_PY_ix" + std::to_string(ix) + "_cx" + std::to_string(cx);
      TH1D* raw = h2->ProjectionY(name.c_str(), ix, ix);
      raw->SetDirectory(nullptr);
      out.push_back(H1Slice{ix, cx, false, std::unique_ptr<TH1D>(raw)});
    }
  }

  return out;
}

// ---------------------- recursive collectors ----------------------
struct ObjPath {
  std::string dir;
  std::string name;
  std::string full() const { return dir.empty() ? name : (dir + "/" + name); }
};

static void collect_TH3_TH2(TDirectory* dir,
                            const std::string& dirpath,
                            std::vector<ObjPath>& out_th3,
                            std::vector<ObjPath>& out_th2) {
  if (!dir) return;

  TIter nextKey(dir->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(nextKey())) {
    auto* cls = TClass::GetClass(key->GetClassName());
    if (!cls) continue;

    if (cls->InheritsFrom(TDirectory::Class())) {
      auto* sub = dynamic_cast<TDirectory*>(key->ReadObj());
      if (!sub) continue;
      const std::string subpath = dirpath.empty() ? sub->GetName() : (dirpath + "/" + sub->GetName());
      collect_TH3_TH2(sub, subpath, out_th3, out_th2);
      continue;
    }

    if (cls->InheritsFrom(TH3D::Class())) out_th3.push_back(ObjPath{dirpath, key->GetName()});
    else if (cls->InheritsFrom(TH2D::Class())) out_th2.push_back(ObjPath{dirpath, key->GetName()});
  }
}

static TH3D* getTH3(TFile* f, const ObjPath& p) {
  if (!f) return nullptr;
  return dynamic_cast<TH3D*>(f->Get(p.full().c_str()));
}
static TH2D* getTH2(TFile* f, const ObjPath& p) {
  if (!f) return nullptr;
  return dynamic_cast<TH2D*>(f->Get(p.full().c_str()));
}

// ---------------------- per-file worker ----------------------
int slice_and_fit_selected(const std::string& in_root,
                           int png_every = -1,
                           bool save_pngs = false) {
  fs::path in_path{in_root};
  const std::string stem   = in_path.stem().string();
  const std::string folder = in_path.parent_path().string();

  std::unique_ptr<TFile> f_in(TFile::Open(in_path.string().c_str(), "READ"));
  if (!f_in || f_in->IsZombie()) {
    std::cerr << "[ERROR] Failed to open: " << in_path << "\n";
    return 2;
  }

  const std::string out_root = folder + "/" + stem + "_sliced_fitted.root";
  TFile fout(out_root.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Failed to create output: " << out_root << "\n";
    return 3;
  }

  std::string out_pngs_base = folder + "/pi0_mass_fits/" + stem + "/";
  std::error_code ec;
  if (save_pngs) fs::create_directories(out_pngs_base, ec);

  std::vector<ObjPath> th3_paths, th2_paths;
  collect_TH3_TH2(f_in.get(), "", th3_paths, th2_paths);

  std::cout << "\n[FILE] " << in_root << "\n";
  std::cout << "[INFO] Found TH3D: " << th3_paths.size() << "  TH2D: " << th2_paths.size() << "\n";

  // --- Fit tree (TH3 only) ---
  TTree tfit("h22_fit", "FitPi0Mass per (xq2bin, zpt2phi) slice from TH3D");

  int    out_ix=0, out_iy=0;
  int    out_xq2=0, out_zpt2phi=0;
  int    out_zbin=0, out_pt2bin=0, out_phibin=0;
  double out_nPions=0, out_errPions=0;
  int    out_ok=0;
  int    out_bkg_order=-1;
  double out_chi2ndf=1e9;
  char   out_reason[256] = {0};
  char   out_src[256]    = {0};

  tfit.Branch("src",        out_src,     "src/C");
  tfit.Branch("ix",         &out_ix,     "ix/I");
  tfit.Branch("iy",         &out_iy,     "iy/I");
  tfit.Branch("xq2bin",     &out_xq2,    "xq2bin/I");
  tfit.Branch("z_pt2_phi",  &out_zpt2phi,"z_pt2_phi/I");
  tfit.Branch("zbin",       &out_zbin,   "zbin/I");
  tfit.Branch("pt2bin",     &out_pt2bin, "pt2bin/I");
  tfit.Branch("phibin",     &out_phibin, "phibin/I");
  tfit.Branch("ok",         &out_ok,     "ok/I");
  tfit.Branch("nPions",     &out_nPions, "nPions/D");
  tfit.Branch("errPions",   &out_errPions,"errPions/D");
  tfit.Branch("bkg_order",  &out_bkg_order,"bkg_order/I");
  tfit.Branch("chi2ndf",    &out_chi2ndf,"chi2ndf/D");
  tfit.Branch("reason",     out_reason,  "reason/C");

  // --- Gen tree (TH2 only; no fit) ---
  TTree tgen("gen_slices", "1D projections from TH2D (no pi0 mass fit)");
  char   gen_src[256] = {0};
  int    gen_projX=0;
  int    gen_ifix=0;
  int    gen_cfix=0;
  double gen_integral=0;
  double gen_mean=0;
  double gen_rms=0;

  tgen.Branch("src",      gen_src,     "src/C");
  tgen.Branch("projX",    &gen_projX,  "projX/I");
  tgen.Branch("ifix",     &gen_ifix,   "ifix/I");
  tgen.Branch("cfix",     &gen_cfix,   "cfix/I");
  tgen.Branch("integral", &gen_integral,"integral/D");
  tgen.Branch("mean",     &gen_mean,   "mean/D");
  tgen.Branch("rms",      &gen_rms,    "rms/D");

  if (png_every <= 0) png_every = 30;

  // ---- TH3: slice Z and fit ----
  for (const auto& p : th3_paths) {
    TH3D* h3 = getTH3(f_in.get(), p);
    if (!h3) continue;

    const std::string src = p.full();
    std::snprintf(out_src, sizeof(out_src), "%s", src.c_str());

    auto slices = makeZSlices(h3);
    long long counter = 0;

    for (auto& s : slices) {
      if (!s.h) continue;

      out_ix = s.ix; out_iy = s.iy;
      out_xq2 = s.cx;
      out_zpt2phi = s.cy;
      DecodeZpt2Phi(out_zpt2phi, out_zbin, out_pt2bin, out_phibin);

      SkipDecision dec = PrecheckHistogram(s.h.get());
      if (dec.skip) {
        out_ok = 0;
        out_nPions = -1; out_errPions = -1;
        out_bkg_order = -1;
        out_chi2ndf = 1e9;
        std::snprintf(out_reason, sizeof(out_reason), "%s", dec.reason.c_str());

        fout.cd();
        std::string base = "TH3_slices/" + BaseNameFromPath(src) + "/skipped/";
        TDirectory* d = MkdirP(&fout, base);
        d->cd();
        s.h->Write();

        if (save_pngs && !dec.is_empty) SaveHistPNG(s.h.get(), out_pngs_base + "/skipped/");
        tfit.Fill();
        continue;
      }

      FitResult fr = FitPi0Mass(s.h.get());
      out_ok = fr.ok ? 1 : 0;
      out_nPions = fr.nPions;
      out_errPions = fr.errPions;
      out_bkg_order = fr.bkg_order;
      out_chi2ndf = fr.chi2ndf;
      std::snprintf(out_reason, sizeof(out_reason), "%s", fr.ok ? "OK" : fr.reason.c_str());

      fout.cd();
      std::string base = "TH3_slices/" + BaseNameFromPath(src)
        + "/xq2_" + std::to_string(out_xq2)
        + "/z" + std::to_string(out_zbin) + "_pt2" + std::to_string(out_pt2bin)
        + "/phi" + std::to_string(out_phibin);

      if (!fr.ok) base = "TH3_slices/" + BaseNameFromPath(src) + "/failed/";
      TDirectory* d = MkdirP(&fout, base);
      d->cd();
      s.h->Write();

      if (save_pngs && fr.ok && (counter % png_every) == 0)
        SaveHistPNG(s.h.get(), out_pngs_base + "/" + BaseNameFromPath(src) + "/");

      ++counter;
      tfit.Fill();
    }
  }

  // ---- TH2: slice to 1D projections (no fit) ----
  for (const auto& p : th2_paths) {
    TH2D* h2 = getTH2(f_in.get(), p);
    if (!h2) continue;

    const std::string src = p.full();
    std::snprintf(gen_src, sizeof(gen_src), "%s", src.c_str());

    auto slices = make2DSlices(h2, true, true);
    for (auto& s : slices) {
      if (!s.h) continue;

      gen_projX    = s.projX ? 1 : 0;
      gen_ifix     = s.ifix;
      gen_cfix     = s.cfix;
      gen_integral = s.h->Integral();
      gen_mean     = s.h->GetMean();
      gen_rms      = s.h->GetRMS();

      fout.cd();
      std::string base = "TH2_slices/" + BaseNameFromPath(src)
        + (s.projX ? "/projX_byY/" : "/projY_byX/")
        + (s.projX ? ("y_" + std::to_string(gen_cfix)) : ("x_" + std::to_string(gen_cfix)));

      TDirectory* d = MkdirP(&fout, base);
      d->cd();
      s.h->Write();

      tgen.Fill();
    }
  }

  fout.cd();
  tfit.Write();
  tgen.Write();
  fout.Write();
  fout.Close();

  std::cout << "[OK] Wrote: " << out_root << "\n";
  return 0;
}

// ---------------- RUN ONE FILE ----------------
void run_slice_and_fit_selected(const char* in_root,
                               int png_every = -1,
                               bool save_pngs = false) {
  if (!in_root || std::strlen(in_root) == 0) {
    std::cerr << "[ERROR] run_slice_and_fit_selected: in_root is empty\n";
    return;
  }
  slice_and_fit_selected(std::string(in_root), png_every, save_pngs);
}

// ---------------- RUN ALL THREE (DEFAULT) ----------------
void run_slice_and_fit_all_selected(int png_every = -1, bool save_pngs = false) {
  std::cout << "[RUN ALL] Using default folder:\n  " << kDefaultHistsDir << "\n";
  std::cout << "[RUN ALL] Files:\n"
            << "  " << kRecDataHists << "\n"
            << "  " << kRecTrueHists << "\n"
            << "  " << kGenBinningHists << "\n";

  slice_and_fit_selected(std::string(kRecDataHists),  png_every, save_pngs);
  slice_and_fit_selected(std::string(kRecTrueHists),  png_every, save_pngs);
  slice_and_fit_selected(std::string(kGenBinningHists), png_every, save_pngs);
}

#ifndef __CLING__
int main(int argc, char** argv) {
  // Default behavior: run all three files.
  // Optional: if you pass one argument, treat it as a single input file to process.
  // Optional: you can also pass png_every and save_pngs.
  //
  // Examples:
  //   ./slice_and_fit_selected
  //   ./slice_and_fit_selected rec_data_hists.root
  //   ./slice_and_fit_selected rec_data_hists.root 30 1
  //   ./slice_and_fit_selected ALL 30 0

  int png_every = (argc >= 3) ? std::atoi(argv[2]) : -1;
  bool save_pngs = (argc >= 4) ? (std::atoi(argv[3]) != 0) : false;

  if (argc < 2) {
    run_slice_and_fit_all_selected(png_every, save_pngs);
    return 0;
  }

  const std::string arg1 = argv[1];
  if (arg1 == "ALL" || arg1 == "all") {
    run_slice_and_fit_all_selected(png_every, save_pngs);
    return 0;
  }

  // single file
  return slice_and_fit_selected(arg1, png_every, save_pngs);
}
#endif
