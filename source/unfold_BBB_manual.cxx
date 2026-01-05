// make_scaled_hist.cxx
// Build h_meas, h_true and h_meas_data the same way as in unfold_onepass.cxx,
// then produce h_out = h_meas_data * h_true / h_meas and save everything.

// Standard
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TH2.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TKey.h"
#include "TError.h"

#include <ROOT/RDataFrame.hxx>

// Your binning
#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins, etc.

// ---------------- Fixed response geometry (same as your code) ----------------
namespace RESP {
  static const int    nX   = 21;   // xq2bin: 0..20 (20 bins, edges -0.5..20.5)
  static const double x_lo = -0.5;
  static const double x_hi = x_lo + nX;

  // z_pt2_phi_bin: N_Zbins * N_pTbins_with_overflow * N_phiTrbins bins, plus +1
  static const int    nZ   = N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1;
  static const double z_lo = -0.5;
  static const double z_hi = z_lo + nZ;
}

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

// Compare axes (same binning ranges, same number of bins)
static bool axes_identical(const TH2& a, const TH2& b) {
  auto same_axis = [](const TAxis* A, const TAxis* B){
    return A->GetNbins()==B->GetNbins()
        && std::fabs(A->GetXmin()-B->GetXmin())<1e-9
        && std::fabs(A->GetXmax()-B->GetXmax())<1e-9;
  };
  return same_axis(a.GetXaxis(), b.GetXaxis()) && same_axis(a.GetYaxis(), b.GetYaxis());
}

// Subtract common logical bins even if ranges differ
static void subtract_common_bins(TH2& dest, const TH2& sub) {
  auto same_axis = [](const TAxis* A, const TAxis* B){
    return A->GetNbins()==B->GetNbins()
        && std::fabs(A->GetXmin()-B->GetXmin())<1e-9
        && std::fabs(A->GetXmax()-B->GetXmax())<1e-9;
  };
  if (same_axis(dest.GetXaxis(), sub.GetXaxis()) && same_axis(dest.GetYaxis(), sub.GetYaxis())) {
    dest.Add(&sub, -1.0);
    return;
  }
  const TAxis* dx = dest.GetXaxis();
  const TAxis* dy = dest.GetYaxis();
  const TAxis* sx = sub.GetXaxis();
  const TAxis* sy = sub.GetYaxis();

  for (int iz=1; iz<=dx->GetNbins(); ++iz) {
    const double zc = dx->GetBinCenter(iz);
    const int iz_s = sx->FindBin(zc);
    if (iz_s<1 || iz_s>sx->GetNbins()) continue;
    if (std::fabs(sx->GetBinCenter(iz_s) - zc) > 0.25) continue;

    for (int ix=1; ix<=dy->GetNbins(); ++ix) {
      const double xc = dy->GetBinCenter(ix);
      const int ix_s = sy->FindBin(xc);
      if (ix_s<1 || ix_s>sy->GetNbins()) continue;
      if (std::fabs(sy->GetBinCenter(ix_s) - xc) > 0.25) continue;
      const double v = dest.GetBinContent(iz, ix) - sub.GetBinContent(iz_s, ix_s);
      dest.SetBinContent(iz, ix, v);
    }
  }
}

// Load X-vs-Z histogram from file and convert to a "truth-like" TH2D with
// axes (Z on X-axis, X on Y-axis) and RESP binning.
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
  out->Sumw2();
  out->SetDirectory(nullptr);

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

// -----------------------------------------------------------------------------
// Bundle to hold our three histograms
// -----------------------------------------------------------------------------
struct HistBundle {
  std::unique_ptr<TH2D> h_meas;      // MC reco (measured space), with errPions errors
  std::unique_ptr<TH2D> h_true;      // MC truth (truth space), with errPions errors + Misses
  std::unique_ptr<TH2D> h_meas_data; // measured data, with errPions errors
};

// -----------------------------------------------------------------------------
// Build h_meas, h_true, h_meas_data in the same way as your code
// -----------------------------------------------------------------------------
static bool build_histograms(HistBundle& out)
{
  const char* treename     = "h22_fit";

  // Measured data file (for h_meas_data)
  const char* data_file =
    "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_data_Mxcut_1/"
    "h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root";

  // rec_true MC files (same pattern as your code)
  std::vector<std::string> files_rt;
  for (int i=1;i<=20;++i) {
    std::string f_rt = Form(
      "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/"
      "unfolding_rec_Mxcut_1/h3_bin_xBQ2_Valerii_%d_fitted.root", i);
    if (!gSystem->AccessPathName(f_rt.c_str())) {
      TFile tf(f_rt.c_str(),"READ");
      TTree* t=nullptr;
      tf.GetObject(treename,t);
      bool ok = t && t->GetBranch("nPions") && t->GetBranch("errPions")
                 && t->GetBranch("z_pt2_phi_bin") && t->GetBranch("z_pt2_phi_bin_gen")
                 && t->GetBranch("xq2bin") && t->GetBranch("xq2bin_gen");
      if (ok) files_rt.push_back(f_rt);
    }
  }
  if (files_rt.empty()) {
    Warning("build_histograms","No usable rec_true files found.");
    return false;
  }

  ROOT::DisableImplicitMT();
  ROOT::RDataFrame df_rt(treename, files_rt);

  // Optional: measured DATA frame
  std::unique_ptr<ROOT::RDataFrame> df_data_ptr;
  bool have_data = (!gSystem->AccessPathName(data_file));
  if (have_data) {
    TFile tfd(data_file,"READ");
    TTree* td=nullptr;
    tfd.GetObject(treename, td);
    have_data = td && td->GetBranch("xq2bin") && td->GetBranch("z_pt2_phi_bin")
                   && td->GetBranch("nPions") && td->GetBranch("errPions");
  }
  if (have_data) {
    df_data_ptr.reset(new ROOT::RDataFrame(treename, data_file));
  } else {
    Warning("build_histograms","Measured data file missing or invalid: %s", data_file);
  }

  const Int_t nZ = RESP::nZ, nX = RESP::nX;
  const double z_lo=RESP::z_lo, z_hi=RESP::z_hi;
  const double x_lo=RESP::x_lo, x_hi=RESP::x_hi;

  // --- Create histograms ---
  out.h_meas.reset(new TH2D("h_meas","Measured MC;z_pt2_phi_bin;xq2bin",
                            nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_true.reset(new TH2D("h_true","Truth MC;z_pt2_phi_bin_gen;xq2bin_gen",
                            nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_meas_data.reset(new TH2D("h_meas_data","Measured data;z_pt2_phi_bin;xq2bin",
                                 nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_meas->Sumw2(true);
  out.h_true->Sumw2(true);
  out.h_meas_data->Sumw2(true);
  out.h_meas->SetDirectory(nullptr);
  out.h_true->SetDirectory(nullptr);
  out.h_meas_data->SetDirectory(nullptr);

  // h_check: reco-like map in truth-space to subtract from generator map
  std::unique_ptr<TH2D> h_check(new TH2D("h_check",
                                         "Reco-like occupancy;z_pt2_phi_bin;xq2bin",
                                         nZ, z_lo, z_hi, nX, x_lo, x_hi));
  h_check->Sumw2(true);
  h_check->SetDirectory(nullptr);

  // Accumulate Σerr^2
  TH2D h_meas_err2_acc("h_meas_err2_acc","",
                       nZ,z_lo,z_hi, nX,x_lo,x_hi);
  TH2D h_true_err2_acc("h_true_err2_acc","",
                       nZ,z_lo,z_hi, nX,x_lo,x_hi);
  TH2D h_meas_data_err2_acc("h_meas_data_err2_acc","",
                            nZ,z_lo,z_hi, nX,x_lo,x_hi);
  h_meas_err2_acc.SetDirectory(nullptr);
  h_true_err2_acc.SetDirectory(nullptr);
  h_meas_data_err2_acc.SetDirectory(nullptr);

  auto in_range = [&](double v, double lo, double hi){
    return std::isfinite(v) && v >= lo && v < hi;
  };

  // ---- Fill measured DATA (sum contents and Σerr^2) ----
  if (df_data_ptr) {
    auto flt_data_fill = df_data_ptr->Filter(
      [](double w,double ew){
        return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0;
      },
      {"nPions","errPions"}
    );
    long long skipped = 0;
    flt_data_fill.Foreach(
      [&](int zbin, int xbin, double wgt, double egt){
        const double z = (double)zbin;
        const double x = (double)xbin;
        if (in_range(z, z_lo, z_hi) && in_range(x, x_lo, x_hi)) {
          out.h_meas_data->Fill(z, x, wgt);
          h_meas_data_err2_acc.Fill(z, x, egt*egt);
        } else {
          ++skipped;
        }
      },
      {"z_pt2_phi_bin","xq2bin","nPions","errPions"}
    );
    if (skipped>0) {
      Warning("build_histograms",
              "Measured data: skipped %lld out-of-range entries.",
              (long long)skipped);
    }
  }

  // ---- Matched rec<->truth events (MC) ----
  auto flt_rt = df_rt.Filter(
    [](double w,double ew){
      return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0;
    },
    {"nPions","errPions"}
  );

  flt_rt.Foreach(
    [&](int z_gen, int z_rec, int x_gen, int x_rec, double ww, double ew){
      const double zrec = (double)z_rec;
      const double xrec = (double)x_rec;
      const double zgen = (double)z_gen;
      const double xgen = (double)x_gen;

      if (!(in_range(zrec,z_lo,z_hi) && in_range(xrec,x_lo,x_hi)
            && in_range(zgen,z_lo,z_hi) && in_range(xgen,x_lo,x_hi))) {
        return;
      }

      // Measured MC
      out.h_meas->Fill(zrec, xrec, ww);
      h_meas_err2_acc.Fill(zrec, xrec, ew*ew);

      // Truth MC (matched part)
      out.h_true->Fill(zgen, xgen, ww);
      h_true_err2_acc.Fill(zgen, xgen, ew*ew);

      // For Miss subtraction: this is the part of truth that *was* reconstructed
      h_check->Fill(zgen, xgen, ww);
    },
    {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions","errPions"}
  );

  // ---- Truth-only Misses ---- (read from gen_binning_2D.root)
  {
    const char* miss2d_file  = "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/gen_binning_2D_MxCut_1.root";
    const char* miss2d_hname = "h2_binX_vs_z"; // adjust if your object name differs

    auto h_truth_like = LoadTruthLikeFromGen2D(
        miss2d_file, miss2d_hname,
        nZ, z_lo, z_hi,
        nX, x_lo, x_hi);

    if (!h_truth_like) {
      Warning("build_histograms",
              "Misses: could not load truth-like 2D from '%s' (skipping Miss step).",
              miss2d_file);
    } else {
      // subtract what you actually reconstructed (reco-like occupancy)
      TH2D h_truth_minus_rec(*h_truth_like);
      h_truth_minus_rec.SetDirectory(nullptr);
      subtract_common_bins(h_truth_minus_rec, *h_check);

      // Add positive residuals as truth-only Miss contributions
      for (int ix=1; ix<=h_truth_minus_rec.GetNbinsX(); ++ix) {
        const double zc = h_truth_minus_rec.GetXaxis()->GetBinCenter(ix);
        for (int iy=1; iy<=h_truth_minus_rec.GetNbinsY(); ++iy) {
          const double xc    = h_truth_minus_rec.GetYaxis()->GetBinCenter(iy);
          const double wmiss = h_truth_minus_rec.GetBinContent(ix, iy);
          if (wmiss > 0.0 &&
              in_range(zc,z_lo,z_hi) && in_range(xc,x_lo,x_hi)) {
            // Miss entries are added to truth but have no errPions → no extra error
            out.h_true->Fill(zc, xc, wmiss);
          }
        }
      }
    }
  }

  // ---- Finalize bin errors from Σerr^2 accumulators ----
  for (int ix=1; ix<=out.h_meas->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=out.h_meas->GetNbinsY(); ++iy) {
      const double e2m = h_meas_err2_acc.GetBinContent(ix,iy);
      out.h_meas->SetBinError(ix,iy, (e2m>0.0 ? std::sqrt(e2m) : 0.0));
    }
  }
  for (int ix=1; ix<=out.h_true->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=out.h_true->GetNbinsY(); ++iy) {
      const double e2t = h_true_err2_acc.GetBinContent(ix,iy);
      out.h_true->SetBinError(ix,iy, (e2t>0.0 ? std::sqrt(e2t) : 0.0));
    }
  }
  for (int ix=1; ix<=out.h_meas_data->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=out.h_meas_data->GetNbinsY(); ++iy) {
      const double e2d = h_meas_data_err2_acc.GetBinContent(ix,iy);
      out.h_meas_data->SetBinError(ix,iy, (e2d>0.0 ? std::sqrt(e2d) : 0.0));
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Build h_out = h_meas_data * h_true / h_meas and save everything
// -----------------------------------------------------------------------------
static int make_scaled_hist(const char* out_file = "measData_times_truth_over_meas.root")
{
  HistBundle pack;
  if (!build_histograms(pack)) {
    std::cerr << "ERROR: could not build histograms.\n";
    return 1;
  }

  if (!pack.h_meas || !pack.h_true || !pack.h_meas_data) {
    std::cerr << "ERROR: missing one or more histograms.\n";
    return 2;
  }

  if (!axes_identical(*pack.h_meas, *pack.h_true) ||
      !axes_identical(*pack.h_meas, *pack.h_meas_data)) {
    Warning("make_scaled_hist",
            "Binning mismatch between h_meas, h_true and/or h_meas_data.");
  }

  // Create output histogram with same binning as truth
  std::unique_ptr<TH2D> h_out(
    (TH2D*)pack.h_true->Clone("h_measData_times_truth_over_meas"));
  h_out->Reset();
  h_out->Sumw2(true);
  h_out->SetTitle("h_meas_data * h_true / h_meas;z_pt2_phi_bin;xq2bin");
  h_out->SetDirectory(nullptr);

  const int nx = h_out->GetNbinsX();
  const int ny = h_out->GetNbinsY();

  const double eps = 0.0;             // threshold for M=0
  const bool include_mc_stat = false; // match your manual BBB default

  for (int ix=1; ix<=nx; ++ix) {
    for (int iy=1; iy<=ny; ++iy) {
      const double D  = pack.h_meas_data->GetBinContent(ix, iy);
      const double eD = pack.h_meas_data->GetBinError(ix, iy);

      const double T  = pack.h_true->GetBinContent(ix, iy);
      const double eT = pack.h_true->GetBinError(ix, iy);

      const double M  = pack.h_meas->GetBinContent(ix, iy);
      const double eM = pack.h_meas->GetBinError(ix, iy);

      if (std::fabs(M) <= eps) {
        h_out->SetBinContent(ix,iy, 0.0);
        h_out->SetBinError(ix,iy,   0.0);
        continue;
      }

      const double SF  = T / M;          // scale factor T_MC / M_MC
      const double val = D * SF;         // = D * T / M

      double err = 0.0;
      if (!include_mc_stat) {
        // Only propagate data stat (like include_mc_stat = false in your BBB)
        err = std::fabs(SF) * eD;
      } else {
        // Full error propagation (if you ever want it)
        const double termD = SF * eD;                 // d f / dD
        const double termT = (D / M) * eT;            // d f / dT
        const double termM = (T * D / (M*M)) * eM;    // d f / dM
        err = std::sqrt(termD*termD + termT*termT + termM*termM);
      }

      h_out->SetBinContent(ix,iy, val);
      h_out->SetBinError(ix,iy,   err);
    }
  }

  // ---- Write everything to file ----
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: cannot create output file " << out_file << "\n";
    return 3;
  }

  if (pack.h_meas)      pack.h_meas->Write("h_meas_mc");
  if (pack.h_true)      pack.h_true->Write("h_true_mc");
  if (pack.h_meas_data) pack.h_meas_data->Write("h_meas_data");
  if (h_out)            h_out->Write("h_measData_times_truth_over_meas");

  fout.Write();
  fout.Close();

  std::cout << "Wrote histograms to " << out_file << "\n";
  return 0;
}


