// unfold_onepass_1D.cxx
// 1D unfolding version.
// We flatten the original 2D bin pair (z_pt2_phi_bin, xq2bin) into a single
// global 1D bin: global_bin = xq2bin * RESP::nZ + z_pt2_phi_bin.

#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TCanvas.h"

// RooUnfold
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldResponse.h"

#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TString.h>
#include <TROOT.h>
#include <TKey.h>

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <limits>
#include <climits>
#include <fstream>   // <-- added for chi2 text output

#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins, etc.

// --------------------------- Bundle we pass around ---------------------------
struct ResponseBundle {
  std::unique_ptr<RooUnfoldResponse> resp;  // trained response
  std::unique_ptr<TH1D> h_meas;             // MC reco (measured space), 1D global bin
  std::unique_ptr<TH1D> h_true;             // MC truth (truth space), 1D global bin
  std::unique_ptr<TH1D> h_meas_data;        // measured data (preferred input), 1D global bin
};

// ---------------- Fixed response geometry (as before) ----------------
namespace RESP {
  static const int    nX   = 21;   // xq2bin: 0..20
  static const double x_lo = -0.5;
  static const double x_hi = x_lo + nX;

  // z_pt2_phi_bin: N_Zbins * N_pTbins_with_overflow * N_phiTrbins bins, plus +1
  static const int    nZ   = N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1;
  static const double z_lo = -0.5;
  static const double z_hi = z_lo + nZ;

  // Flattened global bin for (z, x)
  static const int    nTot = nZ * nX;
  static const double g_lo = -0.5;
  static const double g_hi = g_lo + nTot;
}

// Map (z_bin, x_bin) -> global 1D bin index
static inline int GlobalBinIndex(int z_bin, int x_bin) {
  return x_bin * RESP::nZ + z_bin;   // z is the fast index
}

// Check that two 1D histograms have identical binning
static inline bool axes_identical(const TH1& a, const TH1& b) {
  const TAxis* A = a.GetXaxis();
  const TAxis* B = b.GetXaxis();
  return (A->GetNbins() == B->GetNbins()) &&
         (std::fabs(A->GetXmin() - B->GetXmin()) < 1e-9) &&
         (std::fabs(A->GetXmax() - B->GetXmax()) < 1e-9);
}

// ---------------- Load X-vs-Z histogram for "Misses" helper ----------------
//
// File expected to contain a TH2: X = bin_xBQ2_Valerii (xq2bin), Y = zpt2phit_8x8x9 (Z).
// We transpose into (Z on X-axis, X on Y-axis) with RESP binning.
static std::unique_ptr<TH2D>
LoadTruthLikeFromGen2D(const char* filename,
                       const char* hname, // try this first; fall back to first TH2 in file
                       int nZ, double z_lo, double z_hi,
                       int nX, double x_lo, double x_hi)
{
  if (gSystem->AccessPathName(filename)) {
    Warning("LoadTruthLikeFromGen2D","File '%s' not found.", filename);
    return nullptr;
  }
  TFile f(filename, "READ");
  if (f.IsZombie()) {
    Warning("LoadTruthLikeFromGen2D","Cannot open '%s'.", filename);
    return nullptr;
  }

  TH2* src = nullptr;
  if (hname && *hname) f.GetObject(hname, src);
  if (!src) {
    // fallback: first TH2 in the file
    TIter next(f.GetListOfKeys());
    while (TObject* ok = next()) {
      auto* key = dynamic_cast<TKey*>(ok);
      if (!key) continue;
      TObject* obj = key->ReadObj();
      if (obj && obj->InheritsFrom(TH2::Class())) {
        src = dynamic_cast<TH2*>(obj);
        break;
      }
    }
  }
  if (!src) {
    Warning("LoadTruthLikeFromGen2D","No TH2 found in '%s'.", filename);
    return nullptr;
  }

  // Build output with RESP geometry: X-axis = Z bins, Y-axis = X bins
  auto out = std::make_unique<TH2D>("h_truth_like",
    "Truth-like;zpt2phit_8x8x9;bin_xBQ2_Valerii",
    nZ, z_lo, z_hi, nX, x_lo, x_hi);
  out->Sumw2(); out->SetDirectory(nullptr);

  // src is (X=bin_xBQ2_Valerii, Y=zpt2phit_8x8x9). Transpose into out (Z,X).
  for (int ix=1; ix<=src->GetNbinsX(); ++ix) {
    const double xcen = src->GetXaxis()->GetBinCenter(ix); // X (bin_xBQ2_Valerii)
    for (int iy=1; iy<=src->GetNbinsY(); ++iy) {
      const double zcen = src->GetYaxis()->GetBinCenter(iy); // Z (zpt2phit_8x8x9)
      const double w    = src->GetBinContent(ix, iy);
      if (w==0.0) continue;
      if (zcen < z_lo || zcen >= z_hi || xcen < x_lo || xcen >= x_hi) continue;

      const int bx = out->GetXaxis()->FindBin(zcen); // -> Z axis
      const int by = out->GetYaxis()->FindBin(xcen); // -> X axis
      out->AddBinContent(out->GetBin(bx,by), w);
    }
  }
  return out;
}

// -------------------------- Build response in memory (1D) -------------------------
static bool build_response_in_memory(ResponseBundle& out,
                                     bool include_fakes = false,
                                     bool also_write_response_file = false,
                                     const char* optional_resp_file = "response_out.root")
{
  const char* treename     = "h22_fit";
  const char* treename_gen = "h22";

  // Measured data file (for h_meas_data)
  const char* data_file =
    "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_data/"
    "h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root";

  ROOT::DisableImplicitMT();

  // Gather rec_true and (optionally) rec_fake files
  std::vector<std::string> files_rt, files_fk;
  for (int i=1;i<=20;++i) {
    std::string f_rt = Form("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_true/"
                            "h3_bin_xBQ2_Valerii_%d_fitted.root", i);
    if (!gSystem->AccessPathName(f_rt.c_str())) {
      TFile tf(f_rt.c_str(),"READ");
      TTree* t=nullptr; tf.GetObject(treename,t);
      bool ok = t && t->GetBranch("nPions") && t->GetBranch("errPions")
                 && t->GetBranch("z_pt2_phi_bin") && t->GetBranch("z_pt2_phi_bin_gen")
                 && t->GetBranch("xq2bin") && t->GetBranch("xq2bin_gen");
      if (ok) files_rt.push_back(f_rt);
    }

    if (include_fakes) {
      std::string f_fk = Form("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_fake/"
                              "h3_bin_xBQ2_Valerii_%d_fitted.root", i);
      if (!gSystem->AccessPathName(f_fk.c_str())) {
        TFile tf(f_fk.c_str(),"READ");
        TTree* t=nullptr; tf.GetObject(treename,t);
        bool ok = t && t->GetBranch("nPions") && t->GetBranch("errPions")
                   && t->GetBranch("z_pt2_phi_bin") && t->GetBranch("xq2bin");
        if (ok) files_fk.push_back(f_fk);
      }
    }
  }

  if (files_rt.empty()) {
    Error("build_response_in_memory","No usable rec_true files.");
    return false;
  }
  if (files_fk.empty()) {
    if (include_fakes)
      Warning("build_response_in_memory","No usable rec_fake files; proceeding without fakes.");
    include_fakes = false;
  }

  ROOT::RDataFrame df_rt(treename, files_rt);
  std::unique_ptr<ROOT::RDataFrame> df_fk_ptr;
  if (include_fakes && !files_fk.empty())
    df_fk_ptr.reset(new ROOT::RDataFrame(treename, files_fk));

  // Optional: measured DATA frame
  std::unique_ptr<ROOT::RDataFrame> df_data_ptr;
  bool have_data = (!gSystem->AccessPathName(data_file));
  if (have_data) {
    TFile tfd(data_file,"READ");
    TTree* td=nullptr; tfd.GetObject(treename, td);
    have_data = td && td->GetBranch("xq2bin") && td->GetBranch("z_pt2_phi_bin")
                   && td->GetBranch("nPions") && td->GetBranch("errPions");
  }
  if (have_data) df_data_ptr.reset(new ROOT::RDataFrame(treename, data_file));
  else Warning("build_response_in_memory","Measured data file missing or branches not found: %s", data_file);

  // Filter finite positive weights
  auto flt_rt  = df_rt.Filter(
    [](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; },
    {"nPions","errPions"}
  );

  if (include_fakes && df_fk_ptr) {
    // just to instantiate; we'll filter again later
  }

  const int    nZ   = RESP::nZ;
  const int    nX   = RESP::nX;
  const int    nTot = RESP::nTot;
  const double z_lo = RESP::z_lo, z_hi = RESP::z_hi;
  const double x_lo = RESP::x_lo, x_hi = RESP::x_hi;
  const double g_lo = RESP::g_lo, g_hi = RESP::g_hi;

  // --- create 1D histograms (global bin index) ---
  out.h_meas.reset(new TH1D("h_meas","Measured;(z,x) global index;counts",
                            nTot, g_lo, g_hi));
  out.h_true.reset(new TH1D("h_true","Truth;(z,x) global index;counts",
                            nTot, g_lo, g_hi));
  out.h_meas->Sumw2(true);
  out.h_true->Sumw2(true);

  out.h_meas_data.reset(new TH1D("h_meas_data","Measured data (to unfold);(z,x) global index;counts",
                                 nTot, g_lo, g_hi));
  out.h_meas_data->Sumw2(true);

  // Accumulate Σerr^2
  TH1D h_meas_err2_acc     ("h_meas_err2_acc",     "", nTot, g_lo, g_hi);
  TH1D h_true_err2_acc     ("h_true_err2_acc",     "", nTot, g_lo, g_hi);
  TH1D h_meas_data_err2_acc("h_meas_data_err2_acc","", nTot, g_lo, g_hi);
  h_meas_err2_acc.SetDirectory(nullptr);
  h_true_err2_acc.SetDirectory(nullptr);
  h_meas_data_err2_acc.SetDirectory(nullptr);

  // 2D helper for Miss (truth occupancy reconstructed in MC)
  TH2D h_check2D("h_check2D","Reco-matched truth occupancy;z_pt2_phi_bin_gen;xq2bin_gen",
                 nZ, z_lo, z_hi, nX, x_lo, x_hi);
  h_check2D.Sumw2();
  h_check2D.SetDirectory(nullptr);

  // In-memory response object (1D)
  out.resp.reset(new RooUnfoldResponse(out.h_meas.get(), out.h_true.get(),
                                       "response", "response"));
  out.resp->UseOverflow(false);

  auto in_range = [](double v, double lo, double hi){
    return std::isfinite(v) && v >= lo && v < hi;
  };

  // ---- Fill measured DATA (sum contents and Σerr^2) ----
  if (df_data_ptr) {
    auto flt_data_fill = df_data_ptr->Filter(
      [](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; },
      {"nPions","errPions"}
    );
    Long64_t skipped = 0;
    flt_data_fill.Foreach(
      [&](int zbin, int xbin, double wgt, double egt){
        const double z = (double)zbin;
        const double x = (double)xbin;
        if (in_range(z, z_lo, z_hi) && in_range(x, x_lo, x_hi)) {
          const int g    = GlobalBinIndex(zbin, xbin);
          const double gcoord = (double)g;
          const int ibin = out.h_meas_data->GetXaxis()->FindBin(gcoord);
          out.h_meas_data->AddBinContent(ibin, wgt);
          h_meas_data_err2_acc.AddBinContent(ibin, egt*egt);
        } else {
          ++skipped;
        }
      },
      {"z_pt2_phi_bin","xq2bin","nPions","errPions"}
    );
    if (skipped>0)
      Warning("build_response_in_memory","Measured data: skipped %lld out-of-range entries.",
              (long long)skipped);

    // finalize data bin errors
    for (int ibin=1; ibin<=out.h_meas_data->GetNbinsX(); ++ibin) {
      const double e2 = h_meas_data_err2_acc.GetBinContent(ibin);
      out.h_meas_data->SetBinError(ibin, (e2>0.0 ? std::sqrt(e2) : 0.0));
    }
  }

  // ---- Matched rec<->truth events ----
  flt_rt.Foreach(
    [&](int z_gen, int z_rec, int x_gen, int x_rec, double ww, double ew){
      const double zrec = (double)z_rec;
      const double xrec = (double)x_rec;
      const double zgen = (double)z_gen;
      const double xgen = (double)x_gen;

      if (!(in_range(zrec,z_lo,z_hi) && in_range(xrec,x_lo,x_hi) &&
            in_range(zgen,z_lo,z_hi) && in_range(xgen,x_lo,x_hi)))
        return;

      const int gReco  = GlobalBinIndex(z_rec, x_rec);
      const int gTruth = GlobalBinIndex(z_gen, x_gen);
      const double coordReco  = (double)gReco;
      const double coordTruth = (double)gTruth;

      out.resp->Fill(coordReco, coordTruth, ww);

      // Accumulate errors for measured and truth MC
      const int ibinReco  = h_meas_err2_acc.GetXaxis()->FindBin(coordReco);
      const int ibinTruth = h_true_err2_acc.GetXaxis()->FindBin(coordTruth);
      h_meas_err2_acc.AddBinContent(ibinReco,  ew*ew);
      h_true_err2_acc.AddBinContent(ibinTruth, ew*ew);

      // 2D truth occupancy helper for Miss
      h_check2D.Fill(zgen, xgen, ww);
    },
    {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions","errPions"}
  );

  // ---- Reco-only fakes (optional) ----
  if (include_fakes && df_fk_ptr) {
    auto flt_fk = df_fk_ptr->Filter(
      [](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; },
      {"nPions","errPions"}
    );
    flt_fk.Foreach(
      [&](int /*z_gen*/, int z_rec, int /*x_gen*/, int x_rec, double ww, double ew){
        const double zrec = (double)z_rec;
        const double xrec = (double)x_rec;
        if (!(in_range(zrec,z_lo,z_hi) && in_range(xrec,x_lo,x_hi))) return;

        const int gReco  = GlobalBinIndex(z_rec, x_rec);
        const double coordReco = (double)gReco;
        out.resp->Fake(coordReco, ww);

        const int ibinReco = h_meas_err2_acc.GetXaxis()->FindBin(coordReco);
        h_meas_err2_acc.AddBinContent(ibinReco, ew*ew);
      },
      {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions","errPions"}
    );
  }

  // ---- Truth-only Misses (from gen_binning_2D.root) ----
  {
    const char* miss2d_file  = "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/gen_binning_2D.root";
    const char* miss2d_hname = "h2_binX_vs_z"; // change if your object name differs

    auto h_truth_like = LoadTruthLikeFromGen2D(
      miss2d_file, miss2d_hname,
      nZ, z_lo, z_hi,
      nX, x_lo, x_hi);

    if (!h_truth_like) {
      Warning("build_response_in_memory",
              "Misses: could not load truth-like 2D from '%s' (skipping Miss step).",
              miss2d_file);
    } else {
      // subtract reconstructed occupancy to get truth-only misses
      TH2D h_truth_minus_rec(*h_truth_like);
      h_truth_minus_rec.SetDirectory(nullptr);

      if (h_truth_minus_rec.GetNbinsX()==h_check2D.GetNbinsX() &&
          h_truth_minus_rec.GetNbinsY()==h_check2D.GetNbinsY()) {
        h_truth_minus_rec.Add(&h_check2D, -1.0);
      }

      for (int iz=1; iz<=h_truth_minus_rec.GetNbinsX(); ++iz) {
        const double zc = h_truth_minus_rec.GetXaxis()->GetBinCenter(iz);
        const int z_int = static_cast<int>(std::round(zc));
        for (int ix=1; ix<=h_truth_minus_rec.GetNbinsY(); ++ix) {
          const double xc    = h_truth_minus_rec.GetYaxis()->GetBinCenter(ix);
          const int x_int    = static_cast<int>(std::round(xc));
          const double wmiss = h_truth_minus_rec.GetBinContent(iz, ix);
          if (wmiss <= 0.0) continue;

          if (!in_range(z_int, z_lo, z_hi) || !in_range(x_int, x_lo, x_hi)) continue;

          const int gTruth      = GlobalBinIndex(z_int, x_int);
          const double coordTruth = (double)gTruth;
          out.resp->Miss(coordTruth, wmiss);
        }
      }
    }
  }

  // >>> COPY INTERNAL ROOUNFOLD HISTOS INTO YOUR MC MAPS <<<
  if (out.resp) {
    if (auto* hm_int = dynamic_cast<TH1*>(out.resp->Hmeasured()))
      out.h_meas->Add(hm_int);
    if (auto* ht_int = dynamic_cast<TH1*>(out.resp->Htruth()))
      out.h_true->Add(ht_int);
  }

  // ---- Overwrite bin errors with sqrt(Σ err^2) ----
  for (int ibin=1; ibin<=out.h_meas->GetNbinsX(); ++ibin) {
    const double e2m = h_meas_err2_acc.GetBinContent(ibin);
    out.h_meas->SetBinError(ibin, (e2m>0.0 ? std::sqrt(e2m) : 0.0));
  }
  for (int ibin=1; ibin<=out.h_true->GetNbinsX(); ++ibin) {
    const double e2t = h_true_err2_acc.GetBinContent(ibin);
    out.h_true->SetBinError(ibin, (e2t>0.0 ? std::sqrt(e2t) : 0.0));
  }

  // Optional: also save the response histograms to a file
  if (also_write_response_file) {
    TFile fout(optional_resp_file, "RECREATE");
    if (!fout.IsZombie()) {
      if (out.h_meas) {
        TH1* h = (TH1*)out.h_meas->Clone("h_meas_mc");
        h->SetDirectory(&fout);
        h->Write();
      }
      if (out.h_true) {
        TH1* h = (TH1*)out.h_true->Clone("h_true_mc");
        h->SetDirectory(&fout);
        h->Write();
      }
      if (out.h_meas_data) {
        TH1* h = (TH1*)out.h_meas_data->Clone("h_meas_data");
        h->SetDirectory(&fout);
        h->Write();
      }
      fout.Write();
      fout.Close();
    } else {
      Error("build_response_in_memory","Cannot create '%s'.", optional_resp_file);
    }
  }

  return true;
}

//------- QA: simple 1D QA + migration writeout -------------------
static void ResponseQA(const ResponseBundle& pack, const char* out="qa_response.root") {
  // pull what we can
  TH1* hm = (pack.resp ? dynamic_cast<TH1*>(pack.resp->Hmeasured()) : nullptr);
  TH1* ht = (pack.resp ? dynamic_cast<TH1*>(pack.resp->Htruth())    : nullptr);
  TH2* M  = (pack.resp ? dynamic_cast<TH2*>(pack.resp->Hresponse()) : nullptr); // migration matrix (reco X, truth Y)

  TString out_full = gSystem->ConcatFileName(gSystem->WorkingDirectory(), out);
  TFile f(out_full, "RECREATE");
  if (f.IsZombie()) {
    Error("ResponseQA","Cannot create '%s' (cwd=%s).", out, gSystem->WorkingDirectory());
    return;
  }

  auto write_if = [&](TObject* o, const char* name=nullptr){
    if (!o) return;
    if (name && *name) o->Write(name);
    else               o->Write();
  };

  // ---- Totals ----
  std::cout << "QA totals:\n";
  if (hm) {
    const double Imeas = hm->Integral(0,hm->GetNbinsX()+1);
    std::cout << "  Integral(Hmeasured) = " << Imeas << "\n";
  } else {
    std::cout << "  Hmeasured: (missing)\n";
  }
  if (ht) {
    const double Itruth = ht->Integral(0,ht->GetNbinsX()+1);
    std::cout << "  Integral(Htruth)    = " << Itruth << "\n";
  } else {
    std::cout << "  Htruth: (missing)\n";
  }
  if (M) {
    const double sumM_rows = M->Integral(); // total matched weight
    std::cout << "  Sum matched (Mij)   = " << sumM_rows << "\n";
  } else {
    std::cout << "  Migration M: (missing)\n";
  }

  // ---- Always write the base histos we have ----
  write_if(hm, "h_measured_trained");
  write_if(ht, "h_truth_trained");

  // ---- Migration (may be gigantic) ----
  if (M) {
    const long long nBins = 1LL * M->GetNbinsX() * M->GetNbinsY();
    if (nBins <= 50000000LL) {
      write_if(M, "h_migration");
    } else {
      std::cerr << "QA: skip writing migration matrix: "
                << M->GetNbinsX() << "x" << M->GetNbinsY()
                << " bins (" << nBins << " total).\n";
    }
  }

  f.Write();
  f.Close();
  std::cout << "QA wrote: " << out_full << "\n";
}

// Small helper for writing 1D hists
static void WriteHist1D(TH1* h, const char* name, TFile& fout) {
  if (!h) return;
  h->SetName(name);
  h->SetDirectory(&fout);
  h->Write();
}

// ---------------------- Unfold (Bayes or BinByBin, 1D) --------------------------
static int unfold_from_memory(const ResponseBundle& pack,
                              bool use_bayes = true,
                              int  nIter    = 5,
                              const char* out_file ="unfold_out.root")
{
  if (!pack.resp || !pack.h_true || !pack.h_meas) {
    std::cerr << "ERROR: response / histos not built.\n"; return 1;
  }
  TH1* h_input = pack.h_meas_data ? (TH1*)pack.h_meas_data.get()
                                  : (TH1*)pack.h_meas.get();
  if (!h_input) {
    std::cerr << "ERROR: no measured input.\n"; return 2;
  }

  if (!axes_identical(*h_input, *pack.h_meas)) {
    std::cerr << "WARNING: input measured binning != response measured binning.\n";
  }

  std::unique_ptr<RooUnfold> unfold;
  TString tag;
  if (use_bayes) {
    unfold.reset(new RooUnfoldBayes(pack.resp.get(), h_input, nIter));
    tag = TString::Format("Bayes_iter%d", nIter);
  } else {
    unfold.reset(new RooUnfoldBinByBin(pack.resp.get(), h_input));
    tag = "BinByBin";
  }

  unfold->SetVerbose(0);
  auto errModeHist = RooUnfold::kErrors;
  TH1* h_unfold_raw = dynamic_cast<TH1*>(unfold->Hunfold(errModeHist));
  if (!h_unfold_raw) {
    std::cerr << "ERROR: Hunfold returned null.\n"; return 3;
  }

  TString hname = TString::Format("unfold_%s", tag.Data());
  TH1* h_unfold = (TH1*)h_unfold_raw->Clone(hname);
  h_unfold->SetTitle(TString::Format("Unfolded spectrum (%s);global_bin;entries",
                     tag.Data()));

  // ---- Write outputs ----
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) { std::cerr << "ERROR: cannot create " << out_file << "\n"; return 4; }

  if (pack.h_meas_data) {
    TH1* h_data_clone = (TH1*)pack.h_meas_data->Clone("h_meas_data_input");
    WriteHist1D(h_data_clone, "h_meas_data_input", fout);
  }
  {
    TH1* h_meas_clone = (TH1*)pack.h_meas->Clone("h_meas_mc");
    WriteHist1D(h_meas_clone, "h_meas_mc", fout);
  }
  {
    TH1* h_true_clone = (TH1*)pack.h_true->Clone("h_true_mc");
    WriteHist1D(h_true_clone, "h_true_mc", fout);
  }

  WriteHist1D(h_unfold, hname, fout);

  // --- Chi2 vs truth for this (single) unfolding ---
  if (use_bayes && pack.h_true) {
    double chi2_truth = unfold->Chi2(pack.h_true.get(), RooUnfold::kErrors);
    std::cout << "Bayes unfolding (nIter=" << nIter
              << ")  chi2(unfolded vs MC truth) = " << chi2_truth << std::endl;

    // store in the ROOT file as a parameter
    TParameter<double> pChi2("chi2_vs_truth", chi2_truth);
    pChi2.Write();
  }

  fout.Write(); fout.Close();

  std::cout << "Unfolding done (" << (use_bayes ? "Bayes" : "BinByBin")
            << "). Wrote: " << out_file << "\n";
  return 0;
}

// Save Bayes unfolded output after *every* iteration (1..maxIter).
// Still 1D, flattened global bin.
static int unfold_bayes_save_each_iter(const ResponseBundle& pack,
                                       int  maxIter          = 6,
                                       const char* out_base  = "unfold_bayes",
                                       bool put_all_in_one   = true,
                                       bool /*make_pngs*/    = false)
{
  if (!pack.resp || !pack.h_true || !pack.h_meas) {
    std::cerr << "ERROR: need trained response + MC histos.\n";
    return 1;
  }
  TH1* h_input = pack.h_meas_data ? (TH1*)pack.h_meas_data.get()
                                  : (TH1*)pack.h_meas.get();

  if (!h_input) { std::cerr << "ERROR: no measured input.\n"; return 2; }
  if (!axes_identical(*h_input, *pack.h_meas)) {
    std::cerr << "WARNING: input measured binning != response measured binning.\n";
  }

  auto write_header = [&](TFile& f){
    if (pack.h_meas_data) {
      TH1* h = (TH1*)pack.h_meas_data->Clone("h_meas_data_input");
      WriteHist1D(h, "h_meas_data_input", f);
    }
    {
      TH1* h = (TH1*)pack.h_meas->Clone("h_meas_mc");
      WriteHist1D(h, "h_meas_mc", f);
    }
    {
      TH1* h = (TH1*)pack.h_true->Clone("h_true_mc");
      WriteHist1D(h, "h_true_mc", f);
    }
  };

  // one iteration, with chi2 bookkeeping
  auto one_iter = [&](int i,
                      TFile* f_external,
                      TH1D* hChi2Truth,
                      TH1D* hChi2Change,
                      TH1*& hPrevUnfold,
                      std::ofstream& chi2txt)
  {
    RooUnfoldBayes u(pack.resp.get(), h_input, i);
    u.SetVerbose(0);  // set to 1 if you also want RooUnfold's own console "Chi^2 of change"

    TH1* h_unfold_raw = dynamic_cast<TH1*>(u.Hunfold(RooUnfold::kErrors));
    if (!h_unfold_raw) {
      std::cerr << "ERROR: Hunfold returned null at iteration " << i << ".\n";
      return 3;
    }

    TString hname = Form("unfold_Bayes_iter%02d", i);
    TH1* hC = (TH1*)h_unfold_raw->Clone(hname);
    hC->SetTitle(Form("Unfolded spectrum (Bayes, %d iter);global_bin;entries", i));

    // --- Chi2 vs MC truth ---
    double chi2_truth = -1.0;
    if (pack.h_true) {
      chi2_truth = u.Chi2(pack.h_true.get(), RooUnfold::kErrors);
      if (hChi2Truth) hChi2Truth->SetBinContent(i, chi2_truth);
    }

    // --- Chi2 "change" between iterations: current vs previous unfolded ---
    double chi2_change = -1.0;
    if (hPrevUnfold) {
      chi2_change = u.Chi2(hPrevUnfold, RooUnfold::kErrors);
      if (hChi2Change) hChi2Change->SetBinContent(i, chi2_change);
    }

    std::cout << "Iter " << i
              << "  chi2_truth="  << chi2_truth
              << "  chi2_change=" << chi2_change << std::endl;

    if (chi2txt.is_open()) {
      chi2txt << i << " " << chi2_change << " " << chi2_truth << "\n";
    }

    // update previous-unfolded spectrum
    if (hPrevUnfold) {
      delete hPrevUnfold;
      hPrevUnfold = nullptr;
    }
    hPrevUnfold = (TH1*)hC->Clone("prev_unfold");
    hPrevUnfold->SetDirectory(nullptr);

    // write unfolded histogram
    if (f_external) {
      WriteHist1D(hC, hname, *f_external);
    } else {
      TString ofn = Form("%s_iter%02d.root", out_base, i);
      TFile fout(ofn, "RECREATE");
      if (fout.IsZombie()) { std::cerr << "ERROR: cannot create " << ofn << "\n"; return 4; }
      write_header(fout);
      WriteHist1D(hC, hname, fout);

      // also store chi2 values as TParameters in per-iter files
      TParameter<double> pChi2Truth("chi2_vs_truth", chi2_truth);
      TParameter<double> pChi2Change("chi2_change_vs_prev", chi2_change);
      pChi2Truth.Write();
      pChi2Change.Write();

      fout.Write(); fout.Close();
      std::cout << "Wrote: " << ofn << "\n";
    }
    return 0;
  };

  if (put_all_in_one) {
    TString ofn = Form("%s_allIters.root", out_base);
    TFile fout(ofn, "RECREATE");
    if (fout.IsZombie()) { std::cerr << "ERROR: cannot create " << ofn << "\n"; return 5; }
    write_header(fout);

    // histograms of chi2 vs iteration
    TH1D hChi2Truth("hChi2Truth",
                    "Bayes #chi^{2} vs MC truth;iteration;#chi^{2}",
                    maxIter, 0.5, maxIter+0.5);

    TH1D hChi2Change("hChi2Change",
                     "Bayes #chi^{2} change between iterations;iteration;#chi^{2}",
                     maxIter, 0.5, maxIter+0.5);

    TH1* hPrevUnfold = nullptr;

    std::ofstream chi2txt("bayes_chi2_per_iter.txt");
    if (chi2txt.is_open()) {
      chi2txt << "# iter  chi2_change  chi2_vs_truth\n";
    }

    int rc = 0;
    for (int i=1; i<=maxIter; ++i) {
      rc = one_iter(i, &fout, &hChi2Truth, &hChi2Change, hPrevUnfold, chi2txt);
      if (rc) {
        std::cerr << "Stopping at iteration " << i << " due to error.\n";
        break;
      }
    }

    if (hPrevUnfold) {
      delete hPrevUnfold;
      hPrevUnfold = nullptr;
    }

    // write chi2 histograms into same ROOT file
    hChi2Truth.Write("hChi2Truth");
    hChi2Change.Write("hChi2Change");

    fout.Write(); fout.Close();
    if (!rc) std::cout << "Bayes (1.." << maxIter << ") wrote: " << ofn << "\n";

    if (chi2txt.is_open()) chi2txt.close();
    return rc;
  } else {
    // one file per iteration; still collect chi2 vs iter and save them
    TH1D hChi2Truth("hChi2Truth",
                    "Bayes #chi^{2} vs MC truth;iteration;#chi^{2}",
                    maxIter, 0.5, maxIter+0.5);

    TH1D hChi2Change("hChi2Change",
                     "Bayes #chi^{2} change between iterations;iteration;#chi^{2}",
                     maxIter, 0.5, maxIter+0.5);

    TH1* hPrevUnfold = nullptr;

    std::ofstream chi2txt("bayes_chi2_per_iter.txt");
    if (chi2txt.is_open()) {
      chi2txt << "# iter  chi2_change  chi2_vs_truth\n";
    }

    int rc = 0;
    for (int i=1; i<=maxIter; ++i) {
      rc = one_iter(i, /*f_external=*/nullptr,
                    &hChi2Truth, &hChi2Change,
                    hPrevUnfold, chi2txt);
      if (rc) break;
    }

    if (hPrevUnfold) {
      delete hPrevUnfold;
      hPrevUnfold = nullptr;
    }

    // store chi2-vs-iter histograms in a separate file
    TString ofn_chi2 = Form("%s_chi2.root", out_base);
    TFile fchi2(ofn_chi2, "RECREATE");
    if (!fchi2.IsZombie()) {
      hChi2Truth.Write("hChi2Truth");
      hChi2Change.Write("hChi2Change");
      fchi2.Write();
      fchi2.Close();
      std::cout << "Wrote chi2 histograms: " << ofn_chi2 << "\n";
    } else {
      std::cerr << "ERROR: cannot create " << ofn_chi2 << " for chi2 histograms.\n";
    }

    if (chi2txt.is_open()) chi2txt.close();
    return rc;
  }
}

// ------------------------------- Driver -------------------------------------
// Call from ROOT:
//   .x unfold_onepass_1D.cxx+("bayes",5)
//   .x unfold_onepass_1D.cxx+("bbb")
// or similar.

int onepass_unfold(const char* method = "",
                   int nIter = 9,
                   const char* out_bayes = "unfold_out_bayes.root",
                   const char* out_bbb   = "unfold_out_bbb.root",
                   bool write_response_snapshot = false)
{
  cout<<"starting unf: "<<method<<endl;
  ResponseBundle pack;
  if (!build_response_in_memory(pack, /*include_fakes=*/false,
                                write_response_snapshot, /*optional file*/"response_out.root"))
    return 10;

  cout <<"Performing QA"<<endl;
  ResponseQA(pack, "qa_response.root");
  cout <<"QA is done"<<endl;

  TString m(method); m.ToLower();
  if (m=="bayes" || m=="bayes_iter") {
    return unfold_from_memory(pack, true, nIter, out_bayes);

  } else if (m=="bbb" || m=="binbybin" || m=="bbb_roo" || m=="roo_bbb" || m=="roobbb") {
    // RooUnfold Bin-by-bin on the flattened 1D global bin
    return unfold_from_memory(pack, false, 0, out_bbb);

  } else if (m=="both") {
    int rc1 = unfold_from_memory(pack, true, nIter, out_bayes);
    int rc2 = unfold_from_memory(pack, false, 0,    out_bbb);
    return rc1 ? rc1 : rc2;

  } else if (m=="bayes_sweep" || m=="bayes_scan" || m=="bayes_all") {
    // writes unfold_Bayes_iter01..iterN into one ROOT file
    return unfold_bayes_save_each_iter(pack, nIter, /*out_base=*/"unfold_bayes",
                                       /*put_all_in_one=*/true,  /*make_pngs=*/false);

  } else {
    std::cerr << "Unknown method: " << method
              << " (use 'bayes', 'bbb', 'both', or 'bayes_sweep')\n";
    return 11;
  }
}
