// unfold_onepass.cxx
// Build response in-memory and immediately unfold: no intermediate save/read.

#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH2F.h"
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
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TKey.h>

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <limits>
#include <climits>
#include <TLegend.h>

#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins, etc.




// forward declarations so helpers compile before the full defs appear
// --------------------------- Bundle we pass around ---------------------------
struct ResponseBundle {
  std::unique_ptr<RooUnfoldResponse> resp;  // trained response
  std::unique_ptr<TH2D> h_meas;             // MC reco (measured space), with errPions errors
  std::unique_ptr<TH2D> h_true;             // MC truth (truth space), with errPions errors
  std::unique_ptr<TH2D> h_meas_data;        // measured data (preferred input), with errPions errors
  std::unique_ptr<TH2D> h_check;            // diagnostics (reco-like)
  // migration histogram is intentionally omitted to avoid huge allocations
};


// ------- Make TH2 like another (same variable/regular binning) -------
static TH2D* MakeEmptyTH2DLike(const TH2* tmpl, const char* name, const char* title) {
  const TAxis *ax = tmpl->GetXaxis(), *ay = tmpl->GetYaxis();
  const Int_t nx = ax->GetNbins(), ny = ay->GetNbins();
  const TArrayD* xb = ax->GetXbins();
  const TArrayD* yb = ay->GetXbins();
  TH2D* h = nullptr;
  if (xb && xb->GetSize()==nx+1 && yb && yb->GetSize()==ny+1) {
    h = new TH2D(name, title, nx, xb->GetArray(), ny, yb->GetArray());
  } else if (xb && xb->GetSize()==nx+1) {
    h = new TH2D(name, title, nx, xb->GetArray(), ny, ay->GetXmin(), ay->GetXmax());
  } else if (yb && yb->GetSize()==ny+1) {
    h = new TH2D(name, title, nx, ax->GetXmin(), ax->GetXmax(), ny, yb->GetArray());
  } else {
    h = new TH2D(name, title, nx, ax->GetXmin(), ax->GetXmax(), ny, ay->GetXmin(), ay->GetXmax());
  }
  h->Sumw2();
  h->SetDirectory(nullptr);
  return h;
}


// --- tiny style + drawing helpers ---
static void QA_SetStyle() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
}

static void Save2D(TH2* h, const char* outdir, const char* base, bool logz=false) {
  if (!h) return;
  gSystem->mkdir(outdir, /*recursive=*/kTRUE);
  TCanvas c("c", "c", 1200, 1000);
  c.SetRightMargin(0.15);
  if (logz) c.SetLogz();
  TString ttl = h->GetTitle(); // keep user titles
  h->SetTitle(ttl);
  h->Draw("COLZ");
  // draw integral
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
  tx.DrawLatex(0.12, 0.93,
               Form("Integral = %.6g", h->Integral(0,h->GetNbinsX()+1,0,h->GetNbinsY()+1)));
  c.SaveAs(Form("%s/%s%s.png", outdir, base, logz ? "_logz" : ""));
}

static TH2D* MakeRatio2D(const TH2* num, const TH2* den,
                         const char* name, const char* title, double eps=0.0) {
  if (!num || !den) return nullptr;
  TH2D* r = MakeEmptyTH2DLike(den, name, title);
  const int nx = den->GetNbinsX(), ny = den->GetNbinsY();
  for (int ix=1; ix<=nx; ++ix) for (int iy=1; iy<=ny; ++iy) {
    const double d = den->GetBinContent(ix,iy);
    const double n = num->GetBinContent(ix,iy);
    r->SetBinContent(ix,iy, (std::fabs(d) > eps) ? (n/d) : 0.0);
  }
  return r;
}

static TH1D* ProjYSum(const TH2* h, const char* name, const char* title) {
  if (!h) return nullptr;
  // Sum all X-bins to get distribution vs Y (xq2bin)
  TH1D* p = h->ProjectionY(name, 1, h->GetNbinsX(), "e");
  p->SetTitle(title);
  return p;
}

static void SaveOverlayY(const TH2* h_true, const TH2* h_meas, const TH2* h_unf,
                         const char* outdir, const char* base) {
  if (!h_meas || !h_true || !h_unf) return;
  gSystem->mkdir(outdir, /*recursive=*/kTRUE);

  std::unique_ptr<TH1D> t(ProjYSum(h_true, "py_true",  "Sum over Z; xq2bin; counts"));
  std::unique_ptr<TH1D> m(ProjYSum(h_meas, "py_meas",  "Sum over Z; xq2bin; counts"));
  std::unique_ptr<TH1D> u(ProjYSum(h_unf,  "py_unfold","Sum over Z; xq2bin; counts"));

  TCanvas c("c_py","c_py", 1100, 850);
  gStyle->SetOptStat(0);
  t->SetLineColor(kBlack); t->SetLineWidth(2);
  m->SetLineColor(kRed);   m->SetLineWidth(2);
  u->SetLineColor(kBlue);  u->SetLineWidth(2);
  double ymax = std::max({t->GetMaximum(), m->GetMaximum(), u->GetMaximum()}) * 1.15;
  t->SetMaximum(ymax);
  t->Draw("HIST");
  m->Draw("HIST SAME");
  u->Draw("HIST SAME");

  TLegend leg(0.62,0.75,0.88,0.90);
  leg.AddEntry(t.get(),"Truth MC","l");
  leg.AddEntry(m.get(),"Measured MC (or data)","l");
  leg.AddEntry(u.get(),"Unfolded (BBB)","l");
  leg.Draw();

  c.SaveAs(Form("%s/%s.png", outdir, base));
}

// --- one-stop QA for BBB ---
static void QA_BBB_SavePNGs(const ResponseBundle& pack,
                            TH2* h_unfold,
                            TH2* h_sf,
                            TH2* h_zero,
                            const char* outdir = "png_bbb") {
  QA_SetStyle();
  // originals / inputs
  Save2D(pack.h_true.get(),      outdir, "true_mc");
  Save2D(pack.h_true.get(),      outdir, "true_mc", true);

  Save2D(pack.h_meas.get(),      outdir, "meas_mc");
  Save2D(pack.h_meas.get(),      outdir, "meas_mc", true);

  if (pack.h_meas_data) {
    Save2D(pack.h_meas_data.get(), outdir, "meas_data");
    Save2D(pack.h_meas_data.get(), outdir, "meas_data", true);
  }

  // results from BBB
  Save2D(h_unfold, outdir, "unfold_bbb");
  Save2D(h_unfold, outdir, "unfold_bbb", true);

  Save2D(h_sf,     outdir, "bbb_scale_factor");
  Save2D(h_zero,   outdir, "bbb_zero_mask");

  // diagnostics: ratio (unfold / truth) and difference
  std::unique_ptr<TH2D> h_ratio(
    MakeRatio2D(h_unfold, pack.h_true.get(),
                "h_ratio_unfold_true", "Unfold / Truth;Z index; xq2bin", /*eps=*/0.0));
  if (h_ratio) {
    Save2D(h_ratio.get(), outdir, "ratio_unfold_over_truth");
  }

  std::unique_ptr<TH2D> h_diff(
    MakeEmptyTH2DLike(h_unfold, "h_diff_unfold_minus_truth", "Unfold - Truth"));
  if (h_diff) {
    h_diff->Add(h_unfold, 1.0);
    h_diff->Add(pack.h_true.get(), -1.0);
    Save2D(h_diff.get(), outdir, "diff_unfold_minus_truth");
  }

  // migration (if available in response)
  if (pack.resp) {
    if (auto* M = dynamic_cast<TH2*>(pack.resp->Hresponse())) {
      Save2D(M, outdir, "migration_matrix");
      Save2D(M, outdir, "migration_matrix", true);
    }
  }

  // simple 1D overlay vs xq2bin
  const TH2* h_meas_like = pack.h_meas_data ? (TH2*)pack.h_meas_data.get()
                                            : (TH2*)pack.h_meas.get();
  SaveOverlayY(pack.h_true.get(), h_meas_like, h_unfold, outdir, "projY_xq2bin_overlay");
}


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

static inline bool axes_identical(const TH2& a, const TH2& b) {
  auto same_axis = [](const TAxis* A, const TAxis* B){
    return A->GetNbins()==B->GetNbins()
        && std::fabs(A->GetXmin()-B->GetXmin())<1e-9
        && std::fabs(A->GetXmax()-B->GetXmax())<1e-9;
  };
  return same_axis(a.GetXaxis(), b.GetXaxis()) && same_axis(a.GetYaxis(), b.GetYaxis());
}

// Subtract common logical bins even if ranges differ (your helper)
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

// List ROOTs in a directory that contain required branches.
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

    TFile tf(full.c_str(), "READ");
    if (tf.IsZombie()) continue;
    TTree* t=nullptr; tf.GetObject(treename, t);
    if (!t) continue;

    bool ok=true;
    for (auto& c: needed_cols) {
      if (!t->GetBranch(c.c_str())) { ok=false; break; }
    }
    if (ok) out.push_back(full);
  }
  return out;
}

// ------------ sparse/dense write helpers (unchanged) ------------
static bool IsHuge(const TH2* h, long long maxBins = 20000000LL) {
  const long long nb = 1LL * h->GetNbinsX() * h->GetNbinsY();
  return nb > maxBins;
}
static void WriteAxisEdges(const TAxis& ax, const char* arrName, TDirectory* dir) {
  dir->cd();
  const Int_t nb = ax.GetNbins();
  const TArrayD* xb = ax.GetXbins();
  TVectorD edges(nb+1);
  if (xb && xb->GetSize() == nb+1) {
    for (int i=0; i<=nb; ++i) edges[i] = (*xb)[i];
  } else {
    const double x0 = ax.GetXmin(), x1 = ax.GetXmax();
    const double w  = (x1-x0)/nb;
    for (int i=0;i<=nb;++i) edges[i] = x0 + i*w;
  }
  edges.Write(arrName);
}
static void WriteTH2AsSparseTree(const TH2* h, const char* treeName, TDirectory* dir,
                                 double eps=0.0, bool onlyNonZero=true) {
  dir->cd();
  auto* tr = new TTree(treeName, "Sparse representation of a huge 2D histogram");
  Int_t   ix=0, iy=0;
  Float_t val=0.f, err=0.f;
  tr->Branch("ix",  &ix);
  tr->Branch("iy",  &iy);
  tr->Branch("val", &val);
  tr->Branch("err", &err);
  tr->SetAutoSave(300*1024*1024);
  tr->SetAutoFlush(1000000);

  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  WriteAxisEdges(*h->GetXaxis(), TString::Format("%s_xedges", treeName).Data(), dir);
  WriteAxisEdges(*h->GetYaxis(), TString::Format("%s_yedges", treeName).Data(), dir);
  { TParameter<int>("nx", nx).Write(TString::Format("%s_nx", treeName));
    TParameter<int>("ny", ny).Write(TString::Format("%s_ny", treeName)); }

  for (int jy=1; jy<=ny; ++jy) {
    for (int ix1=1; ix1<=nx; ++ix1) {
      const double v = h->GetBinContent(ix1, jy);
      const double e = h->GetBinError(ix1,   jy);
      if (onlyNonZero && std::fabs(v) <= eps && std::fabs(e) <= eps) continue;
      ix  = ix1; iy = jy;
      val = static_cast<Float_t>(v);
      err = static_cast<Float_t>(e);
      tr->Fill();
    }
  }
  tr->Write();
}
static void WriteSmart(TH2* h, const char* name, TFile& fout,
                       double eps=0.0, bool forceSparse=false) {
  if (!h) return;
  h->SetName(name);
  if (forceSparse || IsHuge(h)) {
    auto* sub = fout.mkdir(name);
    WriteTH2AsSparseTree(h, "bins", sub, eps, /*onlyNonZero=*/true);
  } else {
    h->SetDirectory(&fout);
    h->Write();
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
      if (obj && obj->InheritsFrom(TH2::Class())) { src = dynamic_cast<TH2*>(obj); break; }
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

// -------------------------- Build response in memory -------------------------
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

  // Required branches
  const std::vector<std::string> needed = {
    "z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions","errPions"
  };
  const std::vector<std::string> miss_cols = {
    "bin_xBQ2_Valerii","zpt2phit_8x8x9"
  };
  const char* miss_dir = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning";

  // Gather inputs (same logic as your code)
  std::vector<std::string> files_rt, files_fk;
  for (int i=1;i<=3;++i) {
    std::string f_rt = Form("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_true/h3_bin_xBQ2_Valerii_%d_fitted.root", i);
    if (!gSystem->AccessPathName(f_rt.c_str())) {
      TFile tf(f_rt.c_str(),"READ"); TTree* t=nullptr; tf.GetObject(treename,t);
      bool ok = t && t->GetBranch("nPions") && t->GetBranch("errPions")
                 && t->GetBranch("z_pt2_phi_bin") && t->GetBranch("z_pt2_phi_bin_gen")
                 && t->GetBranch("xq2bin") && t->GetBranch("xq2bin_gen");
      if (ok) files_rt.push_back(f_rt);
    }
    /*
    std::string f_fk = Form("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_fake/h3_bin_xBQ2_Valerii_%d_fitted.root", i);
    if (!gSystem->AccessPathName(f_fk.c_str())) {
      TFile tf(f_fk.c_str(),"READ"); TTree* t=nullptr; tf.GetObject(treename,t);
      bool ok = t && t->GetBranch("nPions") && t->GetBranch("errPions")
                 && t->GetBranch("z_pt2_phi_bin") && t->GetBranch("xq2bin");
      if (ok) files_fk.push_back(f_fk);
    }
    */
  }
  if (files_rt.empty()) { Error("build_response_in_memory","No usable rec_true files."); return false; }
  if (files_fk.empty()) { Warning("build_response_in_memory","No usable rec_fake files; proceeding without fakes."); include_fakes=false; }

  std::vector<std::string> files_all = files_rt;
  files_all.insert(files_all.end(), files_fk.begin(), files_fk.end());

  std::vector<std::string> files_miss = list_root_files_with_branches(miss_dir, treename_gen, miss_cols);
  if (files_miss.empty()) {
    Warning("build_response_in_memory","No Miss files found in %s (Miss step skipped).", miss_dir);
  }

  ROOT::DisableImplicitMT();
  ROOT::RDataFrame df_all(treename, files_all);
  ROOT::RDataFrame df_rt (treename, files_rt);
  std::unique_ptr<ROOT::RDataFrame> df_fk_ptr;
  if (!files_fk.empty()) df_fk_ptr.reset(new ROOT::RDataFrame(treename, files_fk));
  std::unique_ptr<ROOT::RDataFrame> df_miss_ptr;
  if (!files_miss.empty()) df_miss_ptr.reset(new ROOT::RDataFrame(treename_gen, files_miss));

  // Optional: measured DATA frame
  std::unique_ptr<ROOT::RDataFrame> df_data_ptr;
  bool have_data = (!gSystem->AccessPathName(data_file));
  if (have_data) {
    TFile tfd(data_file,"READ"); TTree* td=nullptr; tfd.GetObject(treename, td);
    have_data = td && td->GetBranch("xq2bin") && td->GetBranch("z_pt2_phi_bin")
                   && td->GetBranch("nPions") && td->GetBranch("errPions");
  }
  if (have_data) df_data_ptr.reset(new ROOT::RDataFrame(treename, data_file));
  else Warning("build_response_in_memory","Measured data file missing or branches not found: %s", data_file);

  // Filter finite positive weights
  auto flt_all = df_all.Filter([](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; }, {"nPions","errPions"});
  auto flt_rt  = df_rt .Filter([](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; }, {"nPions","errPions"});

  const Int_t nZ = RESP::nZ, nX = RESP::nX;
  const double z_lo=RESP::z_lo, z_hi=RESP::z_hi, x_lo=RESP::x_lo, x_hi=RESP::x_hi;

  // --- create histograms (on heap; response will reference them) ---
  out.h_meas.reset(new TH2D("h_meas","Measured;z_pt2_phi_bin;xq2bin", nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_true.reset(new TH2D("h_true","Truth;z_pt2_phi_bin_gen;xq2bin_gen", nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_check.reset(new TH2D("filled_pi0_true_fake",";z_pt2_phi_bin;xq2bin", nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_meas->Sumw2(true); out.h_true->Sumw2(true); out.h_check->Sumw2(true);

  out.h_meas_data.reset(new TH2D("h_meas_data","Measured data (to unfold);z_pt2_phi_bin;xq2bin", nZ, z_lo, z_hi, nX, x_lo, x_hi));
  out.h_meas_data->Sumw2(true);

  // Accumulate Σerr^2
  TH2D h_meas_err2_acc("h_meas_err2_acc","", nZ,z_lo,z_hi, nX,x_lo,x_hi);
  TH2D h_true_err2_acc("h_true_err2_acc","", nZ,z_lo,z_hi, nX,x_lo,x_hi);
  TH2D h_check_err2_acc("h_check_err2_acc","", nZ,z_lo,z_hi, nX,x_lo,x_hi);
  TH2D h_meas_data_err2_acc("h_meas_data_err2_acc","", nZ,z_lo,z_hi, nX,x_lo,x_hi);
  h_meas_err2_acc.SetDirectory(nullptr);
  h_true_err2_acc.SetDirectory(nullptr);
  h_check_err2_acc.SetDirectory(nullptr);
  h_meas_data_err2_acc.SetDirectory(nullptr);

  // In-memory response object
  out.resp.reset(new RooUnfoldResponse(out.h_meas.get(), out.h_true.get(), "response", "response"));
  out.resp->UseOverflow(false);

  auto in_range = [&](double v, double lo, double hi){ return std::isfinite(v) && v >= lo && v < hi; };


  // QA:
  // --- AUDIT: everything we feed to RooUnfoldResponse (per-xQ2 slices) ---
  TFile* fAudit = new TFile("response_audit.root","RECREATE");
  fAudit->cd();

  // Row-per-action trace of what got filled
  TTree* tfill = new TTree("resp_fills","Everything fed into RooUnfoldResponse");
  Int_t   a_zgen=-1, a_xgen=-1, a_zrec=-1, a_xrec=-1, a_kind=-1;  // kind: 0=match,1=fake,2=miss
  Double_t a_w=0.0;
  tfill->Branch("z_gen",&a_zgen); tfill->Branch("x_gen",&a_xgen);
  tfill->Branch("z_rec",&a_zrec); tfill->Branch("x_rec",&a_xrec);
  tfill->Branch("w",&a_w);
  tfill->Branch("kind",&a_kind);

  // Per-xQ2 “slices” of migration in z: z_rec (X) vs z_gen (Y)
  std::vector<TH2D*> M_given_truthX, M_given_recoX;
  M_given_truthX.reserve(nX);
  M_given_recoX .reserve(nX);
  for (int ix=0; ix<nX; ++ix) {
    auto* ht = new TH2D(Form("M_zrec_vs_zgen_given_xgen_%02d", ix),
                        Form("Matched;z_{rec};z_{gen}  [x_{gen}=%d]", ix),
                        nZ, z_lo, z_hi, nZ, z_lo, z_hi);
    ht->Sumw2(); ht->SetDirectory(fAudit);
    M_given_truthX.push_back(ht);

    auto* hr = new TH2D(Form("M_zrec_vs_zgen_given_xrec_%02d", ix),
                        Form("Matched;z_{rec};z_{gen}  [x_{rec}=%d]", ix),
                        nZ, z_lo, z_hi, nZ, z_lo, z_hi);
    hr->Sumw2(); hr->SetDirectory(fAudit);
    M_given_recoX.push_back(hr);
  }

  // Compact fake/miss maps (handy to project later per xQ2)
  TH2D* h_fakes_zrec_xrec = new TH2D("fakes_zrec_vs_xrec","Fakes;z_{rec};x_{rec}",
                                     nZ,z_lo,z_hi, nX,x_lo,x_hi);
  TH2D* h_miss_zgen_xgen  = new TH2D("miss_zgen_vs_xgen","Misses;z_{gen};x_{gen}",
                                     nZ,z_lo,z_hi, nX,x_lo,x_hi);
  h_fakes_zrec_xrec->Sumw2(); h_fakes_zrec_xrec->SetDirectory(fAudit);
  h_miss_zgen_xgen ->Sumw2(); h_miss_zgen_xgen ->SetDirectory(fAudit);

  ////
  

  // ---- Fill measured DATA (sum contents and Σerr^2) ----
  if (df_data_ptr) {
    auto flt_data_fill = df_data_ptr->Filter(
      [](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; },
      {"nPions","errPions"}
    );
    Long64_t skipped = 0;
    flt_data_fill.Foreach(
      [&](int zbin, int xbin, double wgt, double egt){
        const double z = (double)zbin, x=(double)xbin;
        if (in_range(z, z_lo, z_hi) && in_range(x, x_lo, x_hi)) {
          const int bx = out.h_meas_data->GetXaxis()->FindBin(z);
          const int by = out.h_meas_data->GetYaxis()->FindBin(x);
          out.h_meas_data->AddBinContent(out.h_meas_data->GetBin(bx,by), wgt);
          h_meas_data_err2_acc.AddBinContent(h_meas_data_err2_acc.GetBin(bx,by), egt*egt);
        } else ++skipped;
      },
      {"z_pt2_phi_bin","xq2bin","nPions","errPions"}
    );
    if (skipped>0) Warning("build_response_in_memory","Measured data: skipped %lld out-of-range entries.", (long long)skipped);
    // finalize data bin errors
    for (int ix=1; ix<=out.h_meas_data->GetNbinsX(); ++ix) {
      for (int iy=1; iy<=out.h_meas_data->GetNbinsY(); ++iy) {
        const double e2 = h_meas_data_err2_acc.GetBinContent(ix,iy);
        out.h_meas_data->SetBinError(ix,iy, (e2>0.0 ? std::sqrt(e2) : 0.0));
      }
    }
  }

  // ---- Matched rec<->truth events ----
  flt_rt.Foreach(
    [&](int z_gen, int z_rec, int x_gen, int x_rec, double ww, double ew){
      const double zrec = (double)z_rec, xrec = (double)x_rec;
      const double zgen = (double)z_gen, xgen = (double)x_gen;
      if (!(in_range(zrec,z_lo,z_hi) && in_range(xrec,x_lo,x_hi)
            && in_range(zgen,z_lo,z_hi) && in_range(xgen,x_lo,x_hi))) return;

      out.resp->Fill(zrec, xrec, zgen, xgen, ww);

      a_zgen = z_gen; a_xgen = x_gen; a_zrec = z_rec; a_xrec = x_rec;
      a_w = ww; a_kind = 0; tfill->Fill();
      if (x_gen>=0 && x_gen<nX) M_given_truthX[x_gen]->Fill(zrec, zgen, ww);
      if (x_rec>=0 && x_rec<nX) M_given_recoX [x_rec]->Fill(zrec, zgen, ww);


      
      out.h_check->Fill(zrec, xrec, ww);

      const int bx_m = out.h_meas->GetXaxis()->FindBin(zrec);
      const int by_m = out.h_meas->GetYaxis()->FindBin(xrec);
      h_meas_err2_acc.AddBinContent(h_meas_err2_acc.GetBin(bx_m,by_m), ew*ew);
      h_check_err2_acc.AddBinContent(h_check_err2_acc.GetBin(bx_m,by_m), ew*ew);

      const int bx_t = out.h_true->GetXaxis()->FindBin(zgen);
      const int by_t = out.h_true->GetYaxis()->FindBin(xgen);
      h_true_err2_acc.AddBinContent(h_true_err2_acc.GetBin(bx_t,by_t), ew*ew);
    },
    {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions","errPions"}
  );

  // ---- Reco-only fakes (off by default; leave logic here) ----
  if (include_fakes && df_fk_ptr) {
    auto flt_fk = df_fk_ptr->Filter([](double w,double ew){ return std::isfinite(w) && w>0.0 && std::isfinite(ew) && ew>=0.0; }, {"nPions","errPions"});
    flt_fk.Foreach(
      [&](int /*z_gen*/, int z_rec, int /*x_gen*/, int x_rec, double ww, double ew){
        const double zrec = (double)z_rec, xrec = (double)x_rec;
        if (!(in_range(zrec,z_lo,z_hi) && in_range(xrec,x_lo,x_hi))) return;
        out.resp->Fake(zrec, xrec, ww);
        out.h_check->Fill(zrec, xrec, ww);
        const int bx_m = out.h_meas->GetXaxis()->FindBin(zrec);
        const int by_m = out.h_meas->GetYaxis()->FindBin(xrec);
        h_meas_err2_acc.AddBinContent(h_meas_err2_acc.GetBin(bx_m,by_m), ew*ew);
        h_check_err2_acc.AddBinContent(h_check_err2_acc.GetBin(bx_m,by_m), ew*ew);
      },
      {"z_pt2_phi_bin_gen","z_pt2_phi_bin","xq2bin_gen","xq2bin","nPions","errPions"}
    );
  }

// ---- Truth-only Misses ---- (read from gen_binning_2D.root)
// File should contain a TH2: X=bin_xBQ2_Valerii, Y=zpt2phit_8x8x9
{
  const char* miss2d_file  = "gen_binning_2D.root";
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
    // subtract what you actually reconstructed (reco-like occupancy) to get truth-only misses
    TH2D h_truth_minus_rec(*h_truth_like);
    h_truth_minus_rec.SetDirectory(nullptr);
    subtract_common_bins(h_truth_minus_rec, *out.h_check);

    // --- Save per-xQ2 zpt2phi (Z) distributions BEFORE and AFTER subtraction
    // Put them neatly into response_audit.root/miss_slices
    TDirectory* dMiss = nullptr;
    if (fAudit) {
      fAudit->cd();
      dMiss = dynamic_cast<TDirectory*>(fAudit->Get("miss_slices"));
      if (!dMiss) dMiss = fAudit->mkdir("miss_slices");
    }
    if (dMiss) dMiss->cd();

    // Y-axis of h_truth_like/h_truth_minus_rec is X (xQ2 index), X-axis is Z
    // We project Z for each fixed X bin (1..nX). We name by 0-based x-index for clarity.
    for (int ybin=1; ybin<=nX; ++ybin) {
      const int xidx = ybin - 1; // your convention: centers are 0..nX-1

      // BEFORE: truth-like Z given x
      TH1D* hZ_before = h_truth_like->ProjectionX(
        Form("z_truthLike_given_x%02d", xidx), ybin, ybin, "e");
      if (hZ_before) {
        hZ_before->SetTitle(Form("Truth-like Z | x=%d; zpt2phit_8x8x9; weight", xidx));
        if (dMiss) hZ_before->SetDirectory(dMiss);
        else       hZ_before->SetDirectory(fAudit);
      }

      // AFTER: residual (truth - reco) Z given x (signed residual)
      TH1D* hZ_after = h_truth_minus_rec.ProjectionX(
        Form("z_truthResidual_given_x%02d", xidx), ybin, ybin, "e");
      if (hZ_after) {
        hZ_after->SetTitle(Form("Residual (truth-reco) Z | x=%d; zpt2phit_8x8x9; weight", xidx));
        if (dMiss) hZ_after->SetDirectory(dMiss);
        else       hZ_after->SetDirectory(fAudit);
      }
    }

    // --- Feed Miss(...) for positive residuals; also fill your audit products
    for (int iz=1; iz<=h_truth_minus_rec.GetNbinsX(); ++iz) {
      const double zc = h_truth_minus_rec.GetXaxis()->GetBinCenter(iz);
      for (int ix=1; ix<=h_truth_minus_rec.GetNbinsY(); ++ix) {
        const double xc    = h_truth_minus_rec.GetYaxis()->GetBinCenter(ix);
        const double wmiss = h_truth_minus_rec.GetBinContent(iz, ix);
        if (wmiss > 0.0) {
          out.resp->Miss(zc, xc, wmiss);

          // --- AUDIT: record
          a_zgen = static_cast<int>(zc);
          a_xgen = static_cast<int>(xc);
          a_zrec = -1; a_xrec = -1;
          a_w = wmiss; a_kind = 2;
          tfill->Fill();
          h_miss_zgen_xgen->Fill(a_zgen, a_xgen, wmiss);
        }
      }
    }
  }
}


  // ---- Overwrite bin errors with sqrt(Σ err^2) ----
  for (int ix=1; ix<=out.h_meas->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=out.h_meas->GetNbinsY(); ++iy) {
      const double e2m = h_meas_err2_acc.GetBinContent(ix,iy);
      const double e2c = h_check_err2_acc.GetBinContent(ix,iy);
      out.h_meas->SetBinError(ix,iy, (e2m>0.0 ? std::sqrt(e2m) : 0.0));
      out.h_check->SetBinError(ix,iy, (e2c>0.0 ? std::sqrt(e2c) : 0.0));
    }
  }
  for (int ix=1; ix<=out.h_true->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=out.h_true->GetNbinsY(); ++iy) {
      const double e2t = h_true_err2_acc.GetBinContent(ix,iy);
      out.h_true->SetBinError(ix,iy, (e2t>0.0 ? std::sqrt(e2t) : 0.0));
    }
  }

  // Optional: also save the response to a file (off by default)
  if (also_write_response_file) {
    TFile fout(optional_resp_file, "RECREATE");
    out.h_meas->SetDirectory(&fout);
    out.h_true->SetDirectory(&fout);
    out.h_check->SetDirectory(&fout);
    if (out.h_meas_data) out.h_meas_data->SetDirectory(&fout);
    out.h_meas->Write(); out.h_true->Write(); out.h_check->Write();
    if (out.h_meas_data) out.h_meas_data->Write();
    fout.Write(); fout.Close();
  }

  // --- finalize audit outputs
  fAudit->cd();
  tfill->Write();  // tree + all histos (already attached) into response_audit.root
  fAudit->Write();
  fAudit->Close();

  return true;
}


//------- QA:
static void ResponseQA(const ResponseBundle& pack, const char* out="qa_response.root") {
  // pull what we can
  TH2* hm = (pack.resp ? dynamic_cast<TH2*>(pack.resp->Hmeasured())  : nullptr);
  TH2* ht = (pack.resp ? dynamic_cast<TH2*>(pack.resp->Htruth())     : nullptr);
  TH2* M  = (pack.resp ? dynamic_cast<TH2*>(pack.resp->Hresponse())  : nullptr); // reco=X, truth=Y

  // where are we writing?
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

  // ---- Totals (print what we can) ----
  std::cout << "QA totals:\n";
  if (hm) {
    const double Imeas = hm->Integral(0,hm->GetNbinsX()+1,0,hm->GetNbinsY()+1);
    std::cout << "  Integral(Hmeasured) = " << Imeas << "\n";
  } else {
    std::cout << "  Hmeasured: (missing)\n";
  }
  if (ht) {
    const double Itruth = ht->Integral(0,ht->GetNbinsX()+1,0,ht->GetNbinsY()+1);
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
bool wrote_migration = false;
if (M) {
  // treat anything above ~5e7 bins as "too big" for safe TH2 writing
  const bool M_huge = IsHuge(M, /*maxBins*/ 50000000LL);
  if (M_huge) {
    std::cerr << "QA: skip writing migration matrix: "
              << M->GetNbinsX() << "x" << M->GetNbinsY()
              << " bins (" << 1LL*M->GetNbinsX()*M->GetNbinsY() << " total).\n";
  } else {
    write_if(M, "h_migration");
    wrote_migration = true;
  }

  // Only create normalized views if the base matrix is not huge
  if (!M_huge) {
    TH2* M_col = (TH2*) M->Clone("h_mig_colNorm");
    TH2* M_row = (TH2*) M->Clone("h_mig_rowNorm");
    const int nx = M->GetNbinsX(), ny = M->GetNbinsY();

    for (int jy=1; jy<=ny; ++jy) {
      double s=0; for (int ix=1; ix<=nx; ++ix) s += M->GetBinContent(ix,jy);
      for (int ix=1; ix<=nx; ++ix)
        M_col->SetBinContent(ix,jy, s>0 ? M->GetBinContent(ix,jy)/s : 0.0);
    }
    for (int ix=1; ix<=nx; ++ix) {
      double s=0; for (int jy=1; jy<=ny; ++jy) s += M->GetBinContent(ix,jy);
      for (int jy=1; jy<=ny; ++jy)
        M_row->SetBinContent(ix,jy, s>0 ? M->GetBinContent(ix,jy)/s : 0.0);
    }

    write_if(M_col);
    write_if(M_row);
  } else {
    std::cerr << "QA: skip M_col/M_row (matrix is too large).\n";
  }
} else {
  std::cout << "  Migration M: (missing)\n";
}

  // ---- Efficiency / fake / purity / stability (need M + hm/ht accordingly) ----
  TH2 *hEff=nullptr, *hMiss=nullptr, *hMF=nullptr, *hFake=nullptr, *hPur=nullptr, *hSta=nullptr;
  if (M && hm && ht) {
    hEff  = (TH2*) ht->Clone("h_efficiency");   hEff->Reset("ICES");
    hMiss = (TH2*) ht->Clone("h_miss");         hMiss->Reset("ICES");
    hMF   = (TH2*) hm->Clone("h_matchedFrac");  hMF->Reset("ICES");
    hFake = (TH2*) hm->Clone("h_fakeFrac");     hFake->Reset("ICES");
    hPur  = (TH2*) hm->Clone("h_purity");       hPur->Reset("ICES");
    hSta  = (TH2*) ht->Clone("h_stability");    hSta->Reset("ICES");

    const int nx = M->GetNbinsX(), ny = M->GetNbinsY();
    for (int jy=1; jy<=ny; ++jy) {
      for (int ix=1; ix<=nx; ++ix) {
        double colSum=0, rowSum=0;
        for (int ix2=1; ix2<=nx; ++ix2) colSum += M->GetBinContent(ix2, jy); // sum over reco (fixed truth)
        for (int jy2=1; jy2<=ny; ++jy2) rowSum += M->GetBinContent(ix,  jy2); // sum over truth (fixed reco)

        const double T   = ht->GetBinContent(ix,jy);
        const double R   = hm->GetBinContent(ix,jy);
        const double diag= M->GetBinContent(ix,jy);

        const double eff = (T>0)     ? colSum/T    : 0.0;
        const double mf  = (R>0)     ? rowSum/R    : 0.0;
        const double pur = (rowSum>0)? diag/rowSum : 0.0;
        const double sta = (colSum>0)? diag/colSum : 0.0;

        hEff ->SetBinContent(ix,jy, eff);
        hMiss->SetBinContent(ix,jy, 1.0-eff);
        hMF  ->SetBinContent(ix,jy, mf);
        hFake->SetBinContent(ix,jy, 1.0-mf);
        hPur ->SetBinContent(ix,jy, pur);
        hSta ->SetBinContent(ix,jy, sta);
      }
    }
    write_if(hEff); write_if(hMiss);
    write_if(hMF);  write_if(hFake);
    write_if(hPur); write_if(hSta);
  } else {
    if (!M)  std::cerr << "QA: skip efficiency/purity/etc (missing M)\n";
    if (!hm) std::cerr << "QA: skip matchedFrac/fake/purity (missing Hmeasured)\n";
    if (!ht) std::cerr << "QA: skip efficiency/stability (missing Htruth)\n";
  }

  // ---- MC closure (need response + hm; ratio needs ht) ----
  TH2* hU = nullptr;
  TH2* hRatio = nullptr;
  if (pack.resp && hm && false) {
    RooUnfoldBayes cl(pack.resp.get(), hm, /*iter=*/5);
    // cl.SetVerbose(0); // optional if available
    hU = dynamic_cast<TH2*>(cl.Hunfold(RooUnfold::kErrors));
    write_if(hU, "h_unfold_closure");

    if (hU && ht) {
      hRatio = (TH2*) ht->Clone("h_closure_ratio"); hRatio->Reset("ICES");
      const int nx = hU->GetNbinsX(), ny = hU->GetNbinsY();
      for (int ix=1; ix<=nx; ++ix) for (int jy=1; jy<=ny; ++jy) {
        const double t = ht->GetBinContent(ix,jy);
        const double u = hU->GetBinContent(ix,jy);
        hRatio->SetBinContent(ix,jy, (t>0) ? (u/t) : 0.0);
      }
      write_if(hRatio);
    } else if (!ht) {
      std::cerr << "QA: skip closure ratio (missing Htruth)\n";
    }
  } else {
    if (!pack.resp) std::cerr << "QA: skip closure (missing response)\n";
    if (!hm)        std::cerr << "QA: skip closure (missing Hmeasured)\n";
  }

  f.Write();
  f.Close();
  std::cout << "QA wrote: " << out_full << "\n";
}




// ---------------------- Unfold (Bayes or BinByBin) --------------------------
static int unfold_from_memory(const ResponseBundle& pack,
                              bool use_bayes=true,
                              int  nIter=5,
                              const char* out_file ="unfold_out.root")
{
  if (!pack.resp || !pack.h_true || !pack.h_meas) {
    std::cerr << "ERROR: response / histos not built.\n"; return 1;
  }
  TH2* h_input = pack.h_meas_data ? (TH2*)pack.h_meas_data.get() : (TH2*)pack.h_meas.get();
  if (!h_input) { std::cerr << "ERROR: no measured input.\n"; return 2; }

  if (!axes_identical(*h_input, *pack.h_meas)) {
    std::cerr << "WARNING: input measured binning != response measured binning.\n";
  }

  std::unique_ptr<RooUnfold> unfold;
  TString tag;
  if (use_bayes) { unfold.reset(new RooUnfoldBayes(pack.resp.get(), h_input, nIter)); tag = TString::Format("Bayes_iter%d", nIter); }
  else           { unfold.reset(new RooUnfoldBinByBin(pack.resp.get(), h_input));     tag = "BinByBin"; }

  unfold->SetVerbose(0); 
  auto errModeHist = RooUnfold::kErrors;
  TH2* h_unfold = dynamic_cast<TH2*>(unfold->Hunfold(errModeHist));
  if (!h_unfold) { std::cerr << "ERROR: Hunfold returned null.\n"; return 3; }

  h_unfold->SetName(TString::Format("unfold_%s", tag.Data()));
  h_unfold->SetTitle(TString::Format("Unfolded spectrum (%s);%s;%s",
                     tag.Data(),
                     pack.h_true->GetXaxis()->GetTitle(),
                     pack.h_true->GetYaxis()->GetTitle()));

  // ---- Write outputs (only final things; no response saving here) ----
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) { std::cerr << "ERROR: cannot create " << out_file << "\n"; return 4; }
  if (pack.h_meas_data) WriteSmart((TH2*)pack.h_meas_data->Clone(), "h_meas_data_input", fout);
  WriteSmart((TH2*)pack.h_meas->Clone(), "h_meas_mc", fout);
  WriteSmart((TH2*)pack.h_true->Clone(), "h_true_mc", fout);
  WriteSmart(h_unfold, TString::Format("unfold_%s", tag.Data()), fout, /*eps=*/0.0, /*forceSparse=*/true);
  fout.Write(); fout.Close();

  std::cout << "Unfolding done (" << (use_bayes ? "Bayes" : "BinByBin") << "). Wrote: " << out_file << "\n";
  return 0;
}

// ---------------------- Manual 2D Bin-by-Bin from memory --------------------
static int unfold_manual_bbb_from_memory(const ResponseBundle& pack,
                                         const char* out_file="unfold_out_manual_bbb.root",
                                         bool include_mc_stat=false,
                                         double eps=0.0)
{
  if (!pack.h_meas || !pack.h_true) {
    std::cerr << "ERROR: need MC measured and truth histos.\n"; return 1;
  }
  TH2* h_data = pack.h_meas_data ? (TH2*)pack.h_meas_data.get() : (TH2*)pack.h_meas.get();

  // Preconditions for simple bin-by-bin (your checks)
  if (!axes_identical(*h_data, *pack.h_meas) || !axes_identical(*pack.h_true, *pack.h_meas)) {
    std::cerr << "ERROR: bin-by-bin requires identical binning across data/measMC/trueMC.\n";
    return 2;
  }

  std::unique_ptr<TH2D> h_unfold ( MakeEmptyTH2DLike(pack.h_true.get(), "unfold_ManualBinByBin", "Unfolded spectrum (Manual Bin-by-Bin)") );
  std::unique_ptr<TH2D> h_sf     ( MakeEmptyTH2DLike(pack.h_true.get(), "h_bbb_scale_factor",   "Scale factor T_MC/M_MC per bin") );
  std::unique_ptr<TH2D> h_zero   ( MakeEmptyTH2DLike(pack.h_true.get(), "h_bbb_zeroeff",        "Mask: 1 where M_MC<=eps, else 0") );

  const int nx = pack.h_true->GetNbinsX();
  const int ny = pack.h_true->GetNbinsY();
  for (int jy=1; jy<=ny; ++jy) for (int ix=1; ix<=nx; ++ix) {
    const double T  = pack.h_true->GetBinContent(ix, jy);
    const double eT = pack.h_true->GetBinError  (ix, jy);
    const double M  = pack.h_meas->GetBinContent(ix, jy);
    const double eM = pack.h_meas->GetBinError  (ix, jy);
    const double D  = h_data    ->GetBinContent(ix, jy);
    const double eD = h_data    ->GetBinError  (ix, jy);

    if (std::fabs(M) <= eps) {
      h_sf->SetBinContent(ix,jy,0); h_sf->SetBinError(ix,jy,0);
      h_unfold->SetBinContent(ix,jy,0); h_unfold->SetBinError(ix,jy,0);
      h_zero->SetBinContent(ix,jy,1); continue;
    }
    const double SF = T/M;
    const double TS = SF*D;
    double eTS = 0.0;
    if (!include_mc_stat) eTS = std::fabs(SF)*eD;
    else {
      const double termD = (SF*eD);
      const double termT = (D / M) * eT;
      const double termM = (T * D / (M*M)) * eM;
      eTS = std::sqrt(termD*termD + termT*termT + termM*termM);
    }
    h_sf->SetBinContent(ix,jy, SF);
    if (include_mc_stat) {
      const double relT = (std::fabs(T)>eps) ? eT/std::fabs(T) : 0.0;
      const double relM = (std::fabs(M)>eps) ? eM/std::fabs(M) : 0.0;
      h_sf->SetBinError(ix,jy, std::fabs(SF)*std::sqrt(relT*relT + relM*relM));
    }
    h_unfold->SetBinContent(ix,jy, TS);
    h_unfold->SetBinError  (ix,jy, eTS);
    h_zero->SetBinContent(ix,jy, 0.0);
  }

  QA_BBB_SavePNGs(pack, h_unfold.get(), h_sf.get(), h_zero.get(), "png_bbb");

  
  // write
  TFile fout(out_file,"RECREATE");
  if (fout.IsZombie()) { std::cerr << "ERROR: cannot create " << out_file << "\n"; return 3; }
  if (pack.h_meas_data) WriteSmart((TH2*)pack.h_meas_data->Clone(), "h_meas_data_input", fout);
  WriteSmart((TH2*)pack.h_meas->Clone(), "h_meas_mc", fout);
  WriteSmart((TH2*)pack.h_true->Clone(), "h_true_mc", fout);
  WriteSmart(h_unfold.release(), "unfold_ManualBinByBin", fout, /*eps=*/0.0, /*forceSparse=*/true);
  WriteSmart(h_sf.release(),     "h_bbb_scale_factor",   fout);
  WriteSmart(h_zero.release(),   "h_bbb_zeroeff",        fout);
  fout.Write(); fout.Close();

  std::cout << "Manual bin-by-bin done. Wrote: " << out_file << "\n";
  return 0;
}

// ------------------------------- Driver -------------------------------------
// Call this from ROOT:  .x unfold_onepass.cxx+("bayes",5)
// or compile w/ ACLiC or your build system.
int onepass_unfold(const char* method = "",
                   int nIter = 1,
                   const char* out_bayes = "unfold_out.root",
                   const char* out_bbb   = "unfold_out_manual_bbb.root",
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

  return unfold_manual_bbb_from_memory(pack, out_bbb, false, 0.0);

  /*
  TString m(method); m.ToLower();
  if (m=="bayes" || m=="bayes_iter") {
    cout<<"Performing Bayes"<<endl;
    return unfold_from_memory(pack, true, nIter, out_bayes);
  } else if (m=="bbb" || m=="binbybin") {
    cout<<"Performing BBB"<<endl;
    return unfold_manual_bbb_from_memory(pack, out_bbb, false, 0.0);
  } else if (m=="both") {
    int rc1 = unfold_from_memory(pack, true, nIter, out_bayes);
    int rc2 = unfold_manual_bbb_from_memory(pack, out_bbb, false, 0.0);
    return rc1 ? rc1 : rc2;
  } else {
    std::cerr << "Unknown method: " << method << " (use 'bayes', 'bbb', or 'both')\n";
    return 11;
  }
  */
}


