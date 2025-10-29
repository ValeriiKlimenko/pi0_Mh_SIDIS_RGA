// create_response_obj_sparse.cxx
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <limits>


// there is no true_matching to implement it properly
// do not turn on, unless new MC is avaliable and previous steps are updated
const bool include_fakes = false;

// ---------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------
bool has_all_branches(const std::string& file,
                      const char* treename,
                      const std::vector<std::string>& cols) {
  TFile f(file.c_str(), "READ");
  if (f.IsZombie()) { std::cerr << "Bad file: " << file << "\n"; return false; }
  TTree* t=nullptr; f.GetObject(treename, t);
  if (!t) { std::cerr << "Missing tree '" << treename << "' in: " << file << "\n"; return false; }
  for (auto& c: cols) {
    if (!t->GetBranch(c.c_str())) {
      std::cerr << "Missing branch '" << c << "' in: " << file << "\n";
      return false;
    }
  }
  return true;
}

static inline UInt_t flat_index_2d(UInt_t iz, UInt_t ix, UInt_t nX) {
  // (z,x) -> i = z*nX + x ; 0-based, matches how bins were filled
  return iz * nX + ix;
}

static inline bool axes_identical(const TH2& a, const TH2& b) {
  auto same_axis = [](const TAxis* A, const TAxis* B){
    return A->GetNbins()==B->GetNbins()
        && std::fabs(A->GetXmin()-B->GetXmin())<1e-9
        && std::fabs(A->GetXmax()-B->GetXmax())<1e-9;
  };
  return same_axis(a.GetXaxis(), b.GetXaxis()) && same_axis(a.GetYaxis(), b.GetYaxis());
}

// Subtract 'sub' from 'dest' matching the *same* logical bins even if axis ranges differ.
// Assumes both histos use 1-wide integer-coded bins (centers at integers).
static void subtract_common_bins(TH2& dest, const TH2& sub) {
  if (axes_identical(dest, sub)) {
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

// List all ROOT files in a directory and keep only those that contain required branches.
static std::vector<std::string>
list_root_files_with_branches(const char* dirpath,
                              const char* treename,
                              const std::vector<std::string>& needed_cols)
{
  std::vector<std::string> out;
  TSystemDirectory dir("miss_dir", dirpath);
  TList* flist = dir.GetListOfFiles();
  if (!flist) return out;

  TIter next(flist);
  while (TSystemFile* f = (TSystemFile*)next()) {
    if (f->IsDirectory()) continue;
    TString nm = f->GetName();
    if (!nm.EndsWith(".root", TString::kIgnoreCase)) continue;
    std::string full = std::string(dirpath) + "/" + nm.Data();
    if (gSystem->AccessPathName(full.c_str())) continue;
    if (has_all_branches(full, treename, needed_cols))
      out.push_back(full);
  }
  return out;
}

// ---------------------------------------------------------------
// PART A: build and SAVE sparse response
// ---------------------------------------------------------------
void create_response_obj() {
  const char* treename = "h22_fit";
  const char* treename_gen = "h22";

  // --- measured data input (single file) ---
  const char* data_file =
    "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_data/"
    "h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root";

  const std::vector<std::string> needed = {
    "z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions"
  };
  const std::vector<std::string> miss_cols = {
    "bin_xBQ2_Valerii","zpt2phit_8x8x9"
  };
  const char* miss_dir = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning";

  // --- collect files with needed branches ---
  std::vector<std::string> files_rt; // rec_true
  std::vector<std::string> files_fk; // rec_fake

  for (int i=1;i<=20;++i) {
    std::string f_rt = Form("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_true/h3_bin_xBQ2_Valerii_%d_fitted.root", i);
    if (!gSystem->AccessPathName(f_rt.c_str()) && has_all_branches(f_rt, treename, needed))
      files_rt.push_back(f_rt);

    std::string f_fk = Form("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_fake/h3_bin_xBQ2_Valerii_%d_fitted.root", i);
    if (!gSystem->AccessPathName(f_fk.c_str()) && has_all_branches(f_fk, treename, needed))
      files_fk.push_back(f_fk);
  }

  if (files_rt.empty()) { Error("create_response_obj","No usable rec_true files."); return; }
  if (files_fk.empty()) { Warning("create_response_obj","No usable rec_fake files found; proceeding without fakes."); }

  std::vector<std::string> files_all = files_rt;
  files_all.insert(files_all.end(), files_fk.begin(), files_fk.end());

  // --- collect files for MISSES from the directory (all ROOTs with required branches)
  std::vector<std::string> files_miss = list_root_files_with_branches(miss_dir, treename_gen, miss_cols);
  if (files_miss.empty()) {
    Warning("create_response_obj","No miss files found in %s with required branches; Miss() will be skipped.", miss_dir);
  } else {
    std::cout << "Miss: found " << files_miss.size() << " ROOT files in " << miss_dir << "\n";
  }

  ROOT::DisableImplicitMT(); // response.Fill/Fake + manual TTree filling not MT-safe
  ROOT::RDataFrame df_all(treename, files_all);
  ROOT::RDataFrame df_rt (treename, files_rt);
  std::unique_ptr<ROOT::RDataFrame> df_fk_ptr;
  if (!files_fk.empty()) df_fk_ptr.reset(new ROOT::RDataFrame(treename, files_fk));
  std::unique_ptr<ROOT::RDataFrame> df_miss_ptr;
  if (!files_miss.empty()) df_miss_ptr.reset(new ROOT::RDataFrame(treename_gen, files_miss));

  // Optional: measured DATA frame (for maxima + filling)
  std::unique_ptr<ROOT::RDataFrame> df_data_ptr;
  bool have_data = (!gSystem->AccessPathName(data_file) &&
                    has_all_branches(data_file, treename, {"xq2bin","z_pt2_phi_bin","nPions"}));
  if (have_data) df_data_ptr.reset(new ROOT::RDataFrame(treename, data_file));
  else Warning("create_response_obj","Measured data file missing or branches not found: %s", data_file);

  // Filter finite positive weights only
  auto flt_all = df_all.Filter([](double w){ return std::isfinite(w) && w>0.0; }, {"nPions"});
  auto flt_rt  = df_rt .Filter([](double w){ return std::isfinite(w) && w>0.0; }, {"nPions"});

  // --- derive binning from rec_true + rec_fake + DATA (take the largest) ---
  auto max_z_rec_all = flt_all.Max<int>("z_pt2_phi_bin");
  auto max_z_gen_all = flt_all.Max<int>("z_pt2_phi_bin_gen");
  auto max_x_rec_all = flt_all.Max<int>("xq2bin");
  auto max_x_gen_all = flt_all.Max<int>("xq2bin_gen");

  int max_z_data_val = -1;
  int max_x_data_val = -1;
  if (df_data_ptr) {
    auto flt_data_for_max = df_data_ptr->Filter([](double w){ return std::isfinite(w) && w>0.0; }, {"nPions"});
    auto max_z_data = flt_data_for_max.Max<int>("z_pt2_phi_bin");
    auto max_x_data = flt_data_for_max.Max<int>("xq2bin");
    max_z_data_val = *max_z_data;
    max_x_data_val = *max_x_data;
  }

  const Long64_t maxZ = std::max<Long64_t>(
                          std::max<Long64_t>(*max_z_rec_all, *max_z_gen_all),
                          static_cast<Long64_t>(std::max(0, max_z_data_val)));
  const Long64_t maxX = std::max<Long64_t>(
                          std::max<Long64_t>(*max_x_rec_all, *max_x_gen_all),
                          static_cast<Long64_t>(std::max(0, max_x_data_val)));

  const Long64_t nZ_ll = maxZ + 1;
  const Long64_t nX_ll = maxX + 1;

  auto saneN = [](const char* what, Long64_t n){
    if (!(n > 0 && n <= 2'000'000)) {
      Error("create_response_obj","Bad %s: n=%lld", what, (long long)n);
      return false;
    }
    return true;
  };
  if (!saneN("Z*phi*pt bin count", nZ_ll) || !saneN("xQ2 bin count", nX_ll)) return;

  const Int_t nZ = static_cast<Int_t>(nZ_ll);
  const Int_t nX = static_cast<Int_t>(nX_ll);

  const double z_lo = -0.5, z_hi = z_lo + nZ;
  const double x_lo = -0.5, x_hi = x_lo + nX;

  // --- create histograms (not attached to any file yet) ---
  TH2D h_meas("h_meas","Measured;z_pt2_phi_bin;xq2bin",
              nZ, z_lo, z_hi, nX, x_lo, x_hi);
  TH2D h_true("h_true","Truth;z_pt2_phi_bin_gen;xq2bin_gen",
              nZ, z_lo, z_hi, nX, x_lo, x_hi);
  TH2D h_check("filled_pi0_true_fake",";z_pt2_phi_bin;xq2bin",
               nZ, z_lo, z_hi, nX, x_lo, x_hi);

  // --- measured data (to be unfolded): X=z, Y=xQ2 ---
  TH2F h_meas_data("h_meas_data",
                   "Measured data (to unfold);z_pt2_phi_bin;xq2bin",
                   nZ, z_lo, z_hi,
                   nX, x_lo, x_hi);
  h_meas_data.Sumw2(true);

  RooUnfoldResponse response(&h_meas, &h_true, "response", "response");
  response.UseOverflow(false);

  // Prepare output file and ensure it's the current directory
  TFile fout("response_out.root","RECREATE");
  if (fout.IsZombie()) { Error("create_response_obj","Cannot create output file"); return; }
  fout.cd(); // <<< IMPORTANT: make 'fout' the current gDirectory

  // meta (binning sizes)
  TParameter<int> p_nZ("nZ", nZ);
  TParameter<int> p_nX("nX", nX);
  p_nZ.Write(); p_nX.Write();

  UInt_t irec=0, itruth=0; Float_t w=0.f;
  TTree spr("resp_sparse","sparse response (irec, itruth, w) -- only matched rec<->truth");
  spr.SetDirectory(&fout); // <<< ensure the TTree belongs to the file
  spr.Branch("irec",   &irec,   "irec/i");
  spr.Branch("itruth", &itruth, "itruth/i");
  spr.Branch("w",      &w,      "w/F");

  auto in_range = [&](double v, double lo, double hi){ return std::isfinite(v) && v >= lo && v < hi; };

  // Fill measured data histogram from the single experimental file
  if (df_data_ptr) {
    auto flt_data_fill = df_data_ptr->Filter([](double w){ return std::isfinite(w) && w>0.0; }, {"nPions"});
    Long64_t skipped = 0;
    flt_data_fill.Foreach(
      [&](int zbin, int xbin, double wgt){
        const double z = static_cast<double>(zbin);
        const double x = static_cast<double>(xbin);
        if (in_range(z, z_lo, z_hi) && in_range(x, x_lo, x_hi)) {
          h_meas_data.Fill(z, x, wgt);
        } else {
          ++skipped;
        }
      },
      {"z_pt2_phi_bin","xq2bin","nPions"}
    );
    if (skipped > 0) Warning("create_response_obj","Measured data: skipped %lld out-of-range entries.", (long long)skipped);
  }

  // 1) Matched rec<->truth events
  flt_rt.Foreach(
    [&](int z_gen, int z_rec, int x_gen, int x_rec, double ww){
      const double zrec = static_cast<double>(z_rec);
      const double xrec = static_cast<double>(x_rec);
      const double zgen = static_cast<double>(z_gen);
      const double xgen = static_cast<double>(x_gen);
      if (!(in_range(zrec, z_lo, z_hi) && in_range(xrec, x_lo, x_hi) &&
            in_range(zgen, z_lo, z_hi) && in_range(xgen, x_lo, x_hi))) return;

      response.Fill(zrec, xrec, zgen, xgen, ww);
      h_check.Fill(zrec, xrec, ww);

      const UInt_t uzrec  = static_cast<UInt_t>(z_rec);
      const UInt_t uxrec  = static_cast<UInt_t>(x_rec);
      const UInt_t uzgen  = static_cast<UInt_t>(z_gen);
      const UInt_t uxgen  = static_cast<UInt_t>(x_gen);

      irec   = flat_index_2d(uzrec, uxrec, nX);
      itruth = flat_index_2d(uzgen, uxgen, nX);
      w      = static_cast<Float_t>(ww);
      spr.Fill();
    },
    {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions"}
  );

  // 2) Reco-only fakes
  if (df_fk_ptr && include_fakes) {
    auto flt_fk = df_fk_ptr->Filter([](double w){ return std::isfinite(w) && w>0.0; }, {"nPions"});
    flt_fk.Foreach(
      [&](int /*z_gen*/, int z_rec, int /*x_gen*/, int x_rec, double ww){
        const double zrec = static_cast<double>(z_rec);
        const double xrec = static_cast<double>(x_rec);
        if (!(in_range(zrec, z_lo, z_hi) && in_range(xrec, x_lo, x_hi))) return;

        response.Fake(zrec, xrec, ww);
        h_check.Fill(zrec, xrec, ww);
      },
      {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions"}
    );
  }

  // 3) Truth-only Misses
  if (df_miss_ptr) {
    TH2D h_truth_like("h_truth_like","Truth-like;zpt2phit_8x8x9;bin_xBQ2_Valerii",
                      nZ, z_lo, z_hi, nX, x_lo, x_hi);
    h_truth_like.SetDirectory(nullptr);

    df_miss_ptr->Foreach(
      [&](int z_truth_like, int x_truth_like){
        const double zt = static_cast<double>(z_truth_like);
        const double xt = static_cast<double>(x_truth_like);
        if (!(in_range(zt, z_lo, z_hi) && in_range(xt, x_lo, x_hi))) return;
        h_truth_like.Fill(zt, xt, 1.0);
      },
      {"zpt2phit_8x8x9","bin_xBQ2_Valerii"}
    );

    TH2D h_truth_minus_rec(h_truth_like);
    h_truth_minus_rec.SetName("h_truth_minus_rec_temp");
    h_truth_minus_rec.SetDirectory(nullptr);
    subtract_common_bins(h_truth_minus_rec, h_check);

    for (int iz=1; iz<=h_truth_minus_rec.GetXaxis()->GetNbins(); ++iz) {
      const double zc = h_truth_minus_rec.GetXaxis()->GetBinCenter(iz);
      for (int ix=1; ix<=h_truth_minus_rec.GetYaxis()->GetNbins(); ++ix) {
        const double xc = h_truth_minus_rec.GetYaxis()->GetBinCenter(ix);
        const double wmiss = h_truth_minus_rec.GetBinContent(iz, ix);
        if (wmiss > 0.0) response.Miss(zc, xc, wmiss);
      }
    }
  } else {
    Warning("create_response_obj","Miss() step skipped: no files in %s with (bin_xBQ2_Valerii, zpt2phit_8x8x9).", miss_dir);
  }

  // --- Writeout (ensure 'fout' is the current directory) ---
  fout.cd(); // <<< IMPORTANT
  h_meas.SetDirectory(&fout);
  h_true.SetDirectory(&fout);
  h_check.SetDirectory(&fout);
  h_meas_data.SetDirectory(&fout);

  h_meas.Write();
  h_true.Write();
  h_check.Write();
  h_meas_data.Write();

  spr.Write();

  // Optional: write dense migration if truly tiny. Otherwise rely on resp_sparse.
  const long long nm = 1LL * nZ * nX;
  std::cout << "nZ=" << nZ << " nX=" << nX << " nm=" << nm << "\n";
  
  // If nm is big, skip the dense matrix entirely.
  const long long max_nm_for_dense = 4000; // 4000x4000 ~ 16M bins → comfortably < 1 GiB even with some overhead
  if (nm > max_nm_for_dense) {
    Warning("create_response_obj",
            "Skipping dense migration (nm=%lld → %lldx%lld bins). Using resp_sparse only.",
            nm, nm, nm);
  } else {
    // If you insist on writing it when small:
    TH2F h_migF("h_response_migration_f",";irec;itruth",
                (Int_t)nm, -0.5, (double)nm-0.5,
                (Int_t)nm, -0.5, (double)nm-0.5);
    // Do NOT call Sumw2 here—keep it tiny.
    UInt_t irec2=0, itruth2=0; Float_t w2=0.f;
    spr.SetBranchAddress("irec",   &irec2);
    spr.SetBranchAddress("itruth", &itruth2);
    spr.SetBranchAddress("w",      &w2);
    for (Long64_t i=0, n=spr.GetEntries(); i<n; ++i) {
      spr.GetEntry(i);
      h_migF.Fill((double)irec2 + 0.5, (double)itruth2 + 0.5, (double)w2);
    }
    h_migF.Write();
  }


  fout.Write();
  fout.Close();
  std::cout << "Saved response_out.root (h_meas, h_true, h_meas_data, resp_sparse, h_check, [h_response_migration_f if small])\n";
}

// ---------------------------------------------------------------
// PART B: RECONSTRUCT from the sparse TTree (also loads h_meas_data)
// ---------------------------------------------------------------
bool load_response_from_ttree(const char* infile,
                              RooUnfoldResponse*& outResp,
                              TH2*& outHmeas,
                              TH2*& outHtrue,
                              TH2*& outHmig,
                              TH2*& outHmeasData)
{
  outResp = nullptr; outHmeas = nullptr; outHtrue = nullptr; outHmig = nullptr; outHmeasData = nullptr;

  TFile f(infile, "READ");
  if (f.IsZombie()) { Error("load_response_from_ttree","Bad file: %s", infile); return false; }

  TH2* h_meas = nullptr;  f.GetObject("h_meas", h_meas);
  TH2* h_true = nullptr;  f.GetObject("h_true", h_true);
  if (!h_meas || !h_true) {
    Error("load_response_from_ttree","Missing h_meas and/or h_true in %s", infile);
    return false;
  }

  TH2* h_meas_data = nullptr; f.GetObject("h_meas_data", h_meas_data);

  outHmeas = dynamic_cast<TH2*>(h_meas->Clone("h_meas_recon"));
  outHtrue = dynamic_cast<TH2*>(h_true->Clone("h_true_recon"));
  outHmeas->SetDirectory(nullptr);
  outHtrue->SetDirectory(nullptr);
  if (h_meas_data) {
    outHmeasData = dynamic_cast<TH2*>(h_meas_data->Clone("h_meas_data_recon"));
    outHmeasData->SetDirectory(nullptr);
  } else {
    Warning("load_response_from_ttree","h_meas_data not found (continuing without it).");
  }

  const Long64_t nm = 1LL * outHmeas->GetNbinsX() * outHmeas->GetNbinsY();
  const Long64_t nt = 1LL * outHtrue->GetNbinsX() * outHtrue->GetNbinsY();
  if (nm <= 0 || nt <= 0) {
    Error("load_response_from_ttree","Bad dimensions: nm=%lld nt=%lld", (long long)nm, (long long)nt);
    return false;
  }
  if (nm > INT_MAX || nt > INT_MAX) {
    Error("load_response_from_ttree","Too many bins for a TH2 axis (nm=%lld, nt=%lld).", (long long)nm, (long long)nt);
    return false;
  }

  TH2* h_mig_file = nullptr; f.GetObject("h_response_migration_f", h_mig_file);
  if (h_mig_file) {
    outHmig = dynamic_cast<TH2*>(h_mig_file->Clone("h_response_migration_f_recon"));
    outHmig->SetDirectory(nullptr);
  } else {
    TTree* spr = nullptr; f.GetObject("resp_sparse", spr);
    if (!spr) { Error("load_response_from_ttree","Missing TTree 'resp_sparse' and no 'h_response_migration_f'."); return false; }
    outHmig = new TH2F("h_response_migration_f_recon",";irec;itruth",
                       (Int_t)nm, -0.5, (double)nm-0.5,
                       (Int_t)nt, -0.5, (double)nt-0.5);
    outHmig->Sumw2(true);

    UInt_t irec=0, itruth=0; Float_t w=0.f;
    spr->SetBranchAddress("irec",   &irec);
    spr->SetBranchAddress("itruth", &itruth);
    spr->SetBranchAddress("w",      &w);

    const Long64_t nentries = spr->GetEntries();
    for (Long64_t i=0; i<nentries; ++i) {
      spr->GetEntry(i);
      if (irec < (UInt_t)nm && itruth < (UInt_t)nt) {
        outHmig->Fill((double)irec + 0.5, (double)itruth + 0.5, (double)w);
      }
    }
  }

  outResp = new RooUnfoldResponse(outHmeas, outHtrue, outHmig, "response", "response");
  outResp->UseOverflow(false);
  f.Close();
  return true;
}
