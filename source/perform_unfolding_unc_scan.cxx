// perform_unfolding_sparse_with_toys.cxx
// Safe for huge 2D grids; optional stat uncertainties from Data + MC toys.
// Compile/run in ROOT (Cling) or build into your analysis as needed.

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom3.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"

// RooUnfold headers
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldResponse.h"
#endif

// Ensure iostream is present in interactive too
#include <iostream>

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <cassert>
#include <iomanip>

#include "TError.h"   // Info/Warning/Error
#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TParameter.h"
#include "TString.h"
#include "TSystem.h" 

// ---------------- Debug helpers & toggles ----------------
static bool gVerboseDefault = true;
#define DBG_IF(v) if (v)
#define DBG(v, ...) do { if (v) Info("unfold_debug", __VA_ARGS__); } while(0)

static void DumpHistSummary(const char* tag, const TH2* h, bool verbose=true) {
  if (!h) { DBG(verbose, "[%s] nullptr", tag); return; }
  const int nx=h->GetNbinsX(), ny=h->GetNbinsY();
  double sum=0, sumw2=0, minv=+std::numeric_limits<double>::infinity(), maxv=-minv;
  long long nNaN=0, nInf=0, nNeg=0;
  for (int j=1;j<=ny;++j){
    for (int i=1;i<=nx;++i){
      const double v = h->GetBinContent(i,j);
      const double e = h->GetBinError(i,j);
      if (std::isnan(v) || std::isnan(e)) ++nNaN;
      if (!std::isfinite(v) || !std::isfinite(e)) ++nInf;
      if (v < 0) ++nNeg;
      sum += v; sumw2 += e*e;
      if (v<minv) minv=v; if (v>maxv) maxv=v;
    }
  }
  DBG(verbose, "[%s] %s: %dx%d; sum=%.6g; sumw2=%.6g; min=%.6g; max=%.6g; NaN=%lld; non-finite=%lld; negatives=%lld",
      tag, h->GetName(), nx, ny, sum, sumw2, minv, maxv, nNaN, nInf, nNeg);
}

static bool CheckFiniteHistogram(const char* tag, const TH2* h, bool verbose=true) {
  if (!h) { Error(tag, "Histogram is null"); return false; }
  const int nx=h->GetNbinsX(), ny=h->GetNbinsY();
  for (int j=1;j<=ny;++j)
    for (int i=1;i<=nx;++i){
      const double v=h->GetBinContent(i,j), e=h->GetBinError(i,j);
      if (!std::isfinite(v) || !std::isfinite(e)){
        Error(tag, "Found non-finite at bin (%d,%d): v=%.6g e=%.6g", i,j,v,e);
        return false;
      }
    }
  return true;
}

// If you have the loader compiled elsewhere, you can include it instead.
// Here we forward-declare and re-implement it so this file is self-contained.
static inline bool axes_identical(const TH2& a, const TH2& b) {
  auto same_axis = [](const TAxis* A, const TAxis* B){
    return A->GetNbins()==B->GetNbins()
        && std::fabs(A->GetXmin()-B->GetXmin())<1e-9
        && std::fabs(A->GetXmax()-B->GetXmax())<1e-9;
  };
  return same_axis(a.GetXaxis(), b.GetXaxis()) && same_axis(a.GetYaxis(), b.GetYaxis());
}

static inline UInt_t flat_index_2d(UInt_t iz, UInt_t ix, UInt_t nX) {
  return iz * nX + ix; // (z,x) -> i
}
static inline void unflat_index_2d(UInt_t i, UInt_t nX, UInt_t& iz, UInt_t& ix) {
  iz = i / nX; ix = i % nX;
}

// ---------- Utilities for sparse writing ----------
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

static bool IsHuge(const TH2* h, long long maxBins = 20000000LL) {
  const long long nb = 1LL * h->GetNbinsX() * h->GetNbinsY();
  return nb > maxBins;
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

// ---------- Online stats (diagonal & banded) ----------
struct OnlineDiag {
  std::vector<double> mean, M2;
  Long64_t n=0;
  OnlineDiag() {}
  explicit OnlineDiag(size_t m): mean(m,0.0), M2(m,0.0) {}
  void reset(size_t m){ mean.assign(m,0.0); M2.assign(m,0.0); n=0; }
  void add(const std::vector<double>& x){
    ++n; for (size_t i=0;i<x.size();++i){ double d=x[i]-mean[i]; mean[i]+=d/n; M2[i]+=d*(x[i]-mean[i]); }
  }
  void finish(std::vector<double>& var) const {
    var.resize(mean.size());
    for (size_t i=0;i<mean.size();++i) var[i] = (n>1? M2[i]/(n-1): 0.0);
  }
};

struct OnlineBand {
  int B=0; Long64_t n=0;
  std::vector<double> sumx, sumx2;
  std::vector<std::vector<double>> sumxx; // [Δ-1][i]
  OnlineBand() {}
  OnlineBand(size_t m, int B_): B(B_), sumx(m,0.0), sumx2(m,0.0), sumxx(B_, std::vector<double>(m,0.0)) {}
  void reset(size_t m, int B_){ B=B_; n=0; sumx.assign(m,0.0); sumx2.assign(m,0.0); sumxx.assign(B, std::vector<double>(m,0.0)); }
  void add(const std::vector<double>& x){
    ++n;
    const size_t m=x.size();
    for (size_t i=0;i<m;++i){ sumx[i]+=x[i]; sumx2[i]+=x[i]*x[i]; }
    for (int d=1; d<=B; ++d){
      const size_t stop = m - d;
      for (size_t i=0;i<stop; ++i) sumxx[d-1][i] += x[i]*x[i+d];
    }
  }
  void finish(std::vector<double>& mean, std::vector<double>& var, std::vector<std::vector<double>>& covBands) const {
    const size_t m=sumx.size();
    mean.resize(m); var.resize(m);
    for (size_t i=0;i<m;++i){
      mean[i]=sumx[i]/n;
      double ex2 = sumx2[i]/n;
      var[i] = (n>1? (ex2 - mean[i]*mean[i])*n/(n-1) : 0.0);
    }
    covBands.assign(B, std::vector<double>(m,0.0));
    for (int d=1; d<=B; ++d){
      const size_t stop = m - d;
      for (size_t i=0;i<stop; ++i){
        double exy = sumxx[d-1][i]/n;
        covBands[d-1][i] = (n>1? (exy - mean[i]*mean[i+d])*n/(n-1) : 0.0);
      }
    }
  }
};

// ---------- Flatten helpers ----------
static void flatten(const TH2* h, std::vector<double>& v) {
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  v.resize( size_t(nx)*size_t(ny) );
  size_t k=0;
  for (int jy=1; jy<=ny; ++jy)
    for (int ix=1; ix<=nx; ++ix)
      v[k++] = h->GetBinContent(ix,jy);
}

// ---------- Load response from response_out.root (sparse->dense-safe) ----------
bool load_response_from_ttree(const char* infile,
                              RooUnfoldResponse*& outResp,
                              TH2*& outHmeas,
                              TH2*& outHtrue,
                              TH2*& outHmig,
                              TH2*& outHmeasData)
{
  outResp = nullptr; outHmeas = nullptr; outHtrue = nullptr; outHmig = nullptr; outHmeasData = nullptr;

  DBG(gVerboseDefault, "[LOAD] Opening file: %s", infile);
  TFile f(infile, "READ");
  if (f.IsZombie()) { Error("load_response_from_ttree","Bad file: %s", infile); return false; }

  TH2* h_meas = nullptr;  f.GetObject("h_meas", h_meas);
  TH2* h_true = nullptr;  f.GetObject("h_true", h_true);
  if (!h_meas || !h_true) {
    Error("load_response_from_ttree","Missing h_meas and/or h_true in %s", infile);
    return false;
  }

  TH2* h_meas_data = nullptr; f.GetObject("h_meas_data", h_meas_data);

  // Debug summaries
  DumpHistSummary("LOAD/h_meas", h_meas);
  DumpHistSummary("LOAD/h_true", h_true);
  if (h_meas_data) DumpHistSummary("LOAD/h_meas_data", h_meas_data);

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

  // Finite checks
  if (!CheckFiniteHistogram("LOAD/h_meas", outHmeas) || !CheckFiniteHistogram("LOAD/h_true", outHtrue)) {
    Error("load_response_from_ttree","Non-finite bins in inputs.");
    return false;
  }

  // Rebuild migration TH2F from sparse, but ONLY if small enough (avoid huge allocations).
  TH2* h_mig_file = nullptr; f.GetObject("h_response_migration_f", h_mig_file);
  if (h_mig_file) {
    outHmig = dynamic_cast<TH2*>(h_mig_file->Clone("h_response_migration_f_recon"));
    outHmig->SetDirectory(nullptr);
  } else {
    TTree* spr = nullptr; f.GetObject("resp_sparse", spr);
    if (spr) {
      if (!spr->GetBranch("irec") || !spr->GetBranch("itruth") || !spr->GetBranch("w")) {
        Error("load_response_from_ttree","resp_sparse missing branches (irec, itruth, w).");
        spr = nullptr;
      } else {
        const Long64_t nm = 1LL * outHmeas->GetNbinsX() * outHmeas->GetNbinsY();
        const Long64_t nt = 1LL * outHtrue->GetNbinsX() * outHtrue->GetNbinsY();
        const long double estBins = (long double)nm * (long double)nt;
        DBG(gVerboseDefault, "[LOAD] resp_sparse present. nm=%lld nt=%lld nm*nt≈%.3Le bins", nm, nt, estBins);
        if (estBins <= 16e6L) { // ~16M bins is safe
          outHmig = new TH2F("h_response_migration_f_recon",";irec;itruth",
                             (Int_t)nm, -0.5, (double)nm-0.5,
                             (Int_t)nt, -0.5, (double)nt-0.5);
          outHmig->Sumw2(true);
          UInt_t irec=0, itruth=0; Float_t w=0.f;
          spr->SetBranchAddress("irec",   &irec);
          spr->SetBranchAddress("itruth", &itruth);
          spr->SetBranchAddress("w",      &w);
          const Long64_t nentries = spr->GetEntries();
          DBG(gVerboseDefault, "[LOAD] Filling dense migration from %lld entries", nentries);
          for (Long64_t i=0; i<nentries; ++i) {
            spr->GetEntry(i);
            if (irec < (UInt_t)nm && itruth < (UInt_t)nt)
              outHmig->Fill((double)irec + 0.5, (double)itruth + 0.5, (double)w);
          }
        } else {
          DBG(gVerboseDefault, "[LOAD] Dense migration skipped (too large).");
          outHmig = nullptr; // too big; not needed for unfolding object
        }
      }
    }
  }

  outResp = new RooUnfoldResponse(outHmeas, outHtrue, outHmig, "response", "response");
  outResp->UseOverflow(false);
  f.Close();
  DBG(gVerboseDefault, "[LOAD] Response built. Overflow off.");
  return true;
}

// ---------- Build compact MC-toy "means" from response_out.root ----------
struct MCMeans {
  // Sums over matched pairs (unique cells) and per-axis totals
  std::vector<UInt_t> irec, itruth;        // unique cell indices (flat)
  std::vector<double> mu_pair;             // expected weight per (irec, itruth)
  std::vector<double> mu_rec_total;        // meas axis totals per irec (from h_meas)
  std::vector<double> mu_truth_total;      // truth axis totals per itruth (from h_true)
  std::vector<double> mu_rec_matched;      // matched sum per irec
  std::vector<double> mu_truth_matched;    // matched sum per itruth
  UInt_t nX_meas=0, nX_truth=0;            // for unflatten
  int nZ_meas=0, nZ_truth=0;
};

static MCMeans MakeMCMeans(const char* in_file, const TH2* h_meas, const TH2* h_true)
{
  MCMeans M;
  // Geometry
  M.nZ_meas   = h_meas->GetNbinsX();
  M.nX_meas   = h_meas->GetYaxis()->GetNbins();
  M.nZ_truth  = h_true->GetNbinsX();
  M.nX_truth  = h_true->GetYaxis()->GetNbins();
  const size_t nm_meas  = (size_t)M.nZ_meas  * (size_t)M.nX_meas;
  const size_t nm_truth = (size_t)M.nZ_truth * (size_t)M.nX_truth;
  M.mu_rec_total.resize(nm_meas,  0.0);
  M.mu_truth_total.resize(nm_truth,0.0);
  M.mu_rec_matched.resize(nm_meas, 0.0);
  M.mu_truth_matched.resize(nm_truth, 0.0);

  // Fill totals from histograms
  {
    size_t k=0;
    for (int jy=1; jy<=h_meas->GetNbinsY(); ++jy)
      for (int ix=1; ix<=h_meas->GetNbinsX(); ++ix, ++k)
        M.mu_rec_total[k] = h_meas->GetBinContent(ix,jy);
  }
  {
    size_t k=0;
    for (int jy=1; jy<=h_true->GetNbinsY(); ++jy)
      for (int ix=1; ix<=h_true->GetNbinsX(); ++ix, ++k)
        M.mu_truth_total[k] = h_true->GetBinContent(ix,jy);
  }

  // Read resp_sparse and compress to unique (irec, itruth)
  TFile f(in_file,"READ");
  TTree* spr=nullptr; f.GetObject("resp_sparse", spr);
  if (!spr) {
    Warning("MakeMCMeans","No TTree 'resp_sparse' in %s. MC toys will be disabled.", in_file);
    return M;
  }
  UInt_t irec=0, itruth=0; Float_t w=0.f;
  spr->SetBranchAddress("irec", &irec);
  spr->SetBranchAddress("itruth",&itruth);
  spr->SetBranchAddress("w",    &w);

  // accumulate in hash map
  std::unordered_map<unsigned long long, double> acc;
  acc.reserve( 2u * (unsigned) spr->GetEntries() / 3u );

  const Long64_t nentries = spr->GetEntries();
  for (Long64_t i=0; i<nentries; ++i) {
    spr->GetEntry(i);
    const unsigned long long key = ( (unsigned long long)itruth << 32ull ) | (unsigned long long)irec;
    acc[key] += (double)w;
    // per-axis matched sums
    if (irec < nm_meas)  M.mu_rec_matched  [irec]  += (double)w;
    if (itruth< nm_truth)M.mu_truth_matched[itruth]+= (double)w;
  }
  f.Close();

  // Move to vectors for fast iteration
  M.irec.reserve(acc.size());
  M.itruth.reserve(acc.size());
  M.mu_pair.reserve(acc.size());
  for (auto& kv: acc) {
    const unsigned long long key = kv.first;
    const UInt_t irec0   = (UInt_t)( key         & 0xffffffffull);
    const UInt_t itruth0 = (UInt_t)((key >> 32) & 0xffffffffull);
    M.irec.push_back(irec0);
    M.itruth.push_back(itruth0);
    M.mu_pair.push_back(kv.second);
  }
  return M;
}

// Build one MC toy response from the compressed means
static std::unique_ptr<RooUnfoldResponse>
MakeMCToyResponse(const MCMeans& M, TRandom3& rng,
                  const TH2* h_meas_proto, const TH2* h_true_proto)
{
  // IMPORTANT: keep prototypes alive (avoid use-after-free with some RooUnfold builds)
  TH2* hM = static_cast<TH2*>(h_meas_proto->Clone("hM_proto"));
  hM->Reset(); hM->SetDirectory(nullptr);
  
  TH2* hT = static_cast<TH2*>(h_true_proto->Clone("hT_proto"));
  hT->Reset(); hT->SetDirectory(nullptr);

  auto resp = std::make_unique<RooUnfoldResponse>(hM, hT, "toyresp","toyresp");
  resp->UseOverflow(false);

  // 1) Matched pairs: Poisson fluctuate each (irec, itruth) cell
  const size_t NP = M.mu_pair.size();
  for (size_t j=0; j<NP; ++j) {
    const double mu = M.mu_pair[j];
    if (mu<=0) continue;
    const double draw = rng.PoissonD(mu);
    if (draw<=0) continue;

    // measured pair
    UInt_t zr, xr; unflat_index_2d(M.irec[j],   M.nX_meas,  zr, xr);
    double xrec = hM->GetXaxis()->GetBinCenter(zr+1);
    double yrec = hM->GetYaxis()->GetBinCenter(xr+1);
    
    // truth pair
    UInt_t zt, xt; unflat_index_2d(M.itruth[j], M.nX_truth, zt, xt);
    double xtru = hT->GetXaxis()->GetBinCenter(zt+1);
    double ytru = hT->GetYaxis()->GetBinCenter(xt+1);
    
    resp->Fill(xrec, yrec, xtru, ytru, draw);
  }

  // 2) Fakes (reco-only)
  const size_t nm_meas = M.mu_rec_total.size();
  for (size_t i=0; i<nm_meas; ++i) {
    double mu_fake = M.mu_rec_total[i] - M.mu_rec_matched[i];
    if (mu_fake <= 0) continue;
    const double draw = rng.PoissonD(mu_fake);
    if (draw <= 0) continue;
    UInt_t zr, xr; unflat_index_2d((UInt_t)i, M.nX_meas, zr, xr);
    double xrec = hM->GetXaxis()->GetBinCenter(zr+1);
    double yrec = hM->GetYaxis()->GetBinCenter(xr+1);
    resp->Fake(xrec, yrec, draw);
  }

  // 3) Misses (truth-only)
  const size_t nm_truth = M.mu_truth_total.size();
  for (size_t i=0; i<nm_truth; ++i) {
    double mu_miss = M.mu_truth_total[i] - M.mu_truth_matched[i];
    if (mu_miss <= 0) continue;
    const double draw = rng.PoissonD(mu_miss);
    if (draw <= 0) continue;
    UInt_t zt, xt; unflat_index_2d((UInt_t)i, M.nX_truth, zt, xt);
    double xtru = hT->GetXaxis()->GetBinCenter(zt+1);
    double ytru = hT->GetYaxis()->GetBinCenter(xt+1);
    resp->Miss(xtru, ytru, draw);
  }

  return resp;
}

// -----------------------------------------------------------------------------
// Main entry: perform_unfolding
//   - use_bayes: true=Bayes (D'Agostini), false=BinByBin
//   - nIter:     Bayes iterations (ignored for BinByBin)
//   - include_unc: if true, run toys for Data+MC and attach per-bin errors
//   - ntoys_data, ntoys_mc: number of toys for each source
//   - bandWidth: store banded covariance up to Δ bins (0 disables)
//   - rng_seed:  RNG seed
//   - in_file:   response file (from create_response_obj)
//   - out_file:  base output name; if include_unc, "_wStatDataMC" is appended
int perform_unfolding_unc(bool use_bayes=true,
                      int  nIter=3,
                      bool include_unc=true,
                      int  ntoys_data=600,
                      int  ntoys_mc=400,
                      int  bandWidth=2,
                      int  rng_seed=12345,
                      const char* in_file  ="response_out.root",
                      const char* out_file ="unfold_out.root",
                      bool verbose = gVerboseDefault)
{
  DBG(verbose, "[PHASE] LOAD");
  TTree* covTree = nullptr;
  RooUnfoldResponse* resp = nullptr;
  TH2* h_meas = nullptr;     // MC reco (diagnostic)
  TH2* h_true = nullptr;     // MC truth (diagnostic)
  TH2* h_mig  = nullptr;     // migration (if rebuilt or stored; not required)
  TH2* h_meas_data = nullptr;// measured data (preferred)

  if (!load_response_from_ttree(in_file, resp, h_meas, h_true, h_mig, h_meas_data)) {
    std::cerr << "ERROR: Failed to load response from: " << in_file << "\n";
    return 1;
  }

  TH2* h_input = h_meas_data ? h_meas_data : h_meas;
  if (!h_input) {
    std::cerr << "ERROR: No measured histogram available (neither h_meas_data nor h_meas).\n";
    return 2;
  }
  DumpHistSummary("INPUT/h_input", h_input, verbose);
  if (!CheckFiniteHistogram("INPUT/h_input", h_input, verbose)) return 2;

  if (h_meas && !axes_identical(*h_input, *h_meas)) {
    std::cerr << "WARNING: Input measured histogram binning != response measured binning.\n";
  }

  // Build chosen unfold object
  DBG(verbose, "[PHASE] CENTRAL: building unfold object (use_bayes=%d, nIter=%d)", (int)use_bayes, nIter);
  std::unique_ptr<RooUnfold> unfold;
  TString tag;
  if (use_bayes) {
    unfold.reset(new RooUnfoldBayes(resp, h_input, nIter));
    tag = TString::Format("Bayes_iter%d", nIter);
  } else {
    unfold.reset(new RooUnfoldBinByBin(resp, h_input));
    tag = "BinByBin";
  }

  // Unfold once (central), per-bin (diagonal) errors only for now
  auto errModeHist = RooUnfolding::kErrors;
  TH2* h_unfold = dynamic_cast<TH2*>(unfold->Hunfold(errModeHist));
  if (!h_unfold) {
    std::cerr << "ERROR: Unfolding failed (null output from Hunfold).\n";
    return 3;
  }

  DumpHistSummary("CENTRAL/h_unfold", h_unfold, verbose);
  if (!CheckFiniteHistogram("CENTRAL/h_unfold", h_unfold, verbose)) { return 3; }

  h_unfold->SetName(TString::Format("unfold_%s", tag.Data()));
  h_unfold->SetTitle(TString::Format("Unfolded spectrum (%s);%s;%s",
                          tag.Data(),
                          h_true ? h_true->GetXaxis()->GetTitle() : "true X",
                          h_true ? h_true->GetYaxis()->GetTitle() : "true Y"));

  // --- Optional: statistical uncertainties via toys (Data + MC) ---
  std::vector<double> var_data, var_mc; // flattened variances
  std::vector<std::vector<double>> covBands_data, covBands_mc;
  std::vector<double> dummy_mean, dummy_varB;
  TRandom3 rng(rng_seed);

  const int nx = h_unfold->GetNbinsX();
  const int ny = h_unfold->GetNbinsY();
  const size_t m = size_t(nx) * size_t(ny);

  if (include_unc) {
    // Prepare accumulators
    OnlineDiag diag_data(m), diag_mc(m);
    OnlineBand band_data, band_mc;
    if (bandWidth > 0) {
      band_data.reset(m, bandWidth);
      band_mc.reset(m, bandWidth);
    }

    // ---------- DATA toys ----------
    DBG(verbose, "[PHASE] DATA_TOYS: N=%d bandWidth=%d", ntoys_data, bandWidth);
    if (ntoys_data > 0) {
      std::vector<double> xvec(m);
      for (int t=0; t<ntoys_data; ++t){
        // 1) Poisson toy of measured histogram
        std::unique_ptr<TH2> h_toy((TH2*)h_input->Clone("h_meas_toy"));
        h_toy->SetDirectory(nullptr);
        for (int jy=1; jy<=h_toy->GetNbinsY(); ++jy){
          for (int ix=1; ix<=h_toy->GetNbinsX(); ++ix){
            double mu = std::max(0.0, h_input->GetBinContent(ix,jy));
            double val = rng.PoissonD(mu);
            h_toy->SetBinContent(ix,jy,val);
            h_toy->SetBinError(ix,jy, std::sqrt(std::max(0.0,val)));
          }
        }
        // 2) Unfold quickly (no per-toy error propagation)
        std::unique_ptr<RooUnfold> ut( use_bayes
                ? (RooUnfold*) new RooUnfoldBayes(resp, h_toy.get(), nIter)
                : (RooUnfold*) new RooUnfoldBinByBin(resp, h_toy.get()) );
        TH2* h_u = dynamic_cast<TH2*>( ut->Hunfold(RooUnfolding::kNoError) );
        if (!h_u) continue;
        flatten(h_u, xvec);
        diag_data.add(xvec);
        if (bandWidth>0) band_data.add(xvec);
        if ((t+1)%50==0) DBG(verbose, "  data toy %d/%d", t+1, ntoys_data);
      }
      diag_data.finish(var_data);
      if (bandWidth>0) {
        band_data.finish(dummy_mean, dummy_varB, covBands_data); // we only need covBands_data
      }
    } else {
      var_data.assign(m, 0.0);
      covBands_data.assign(bandWidth, std::vector<double>(m,0.0));
    }

    // ---------- MC toys (response variation) ----------
    DBG(verbose, "[PHASE] MC_TOYS: N=%d", ntoys_mc);
    if (ntoys_mc > 0) {
      // Precompute response means from 'resp_sparse', 'h_meas', 'h_true'
      MCMeans means = MakeMCMeans(in_file, h_meas, h_true);
      if (means.mu_pair.empty() && (means.mu_rec_total.empty() || means.mu_truth_total.empty())) {
        Warning("perform_unfolding","MC toys disabled: couldn't build response means.");
        var_mc.assign(m, 0.0);
        covBands_mc.assign(bandWidth, std::vector<double>(m,0.0));
      } else {
        std::vector<double> xvec(m);
        for (int t=0; t<ntoys_mc; ++t){
          auto toyResp = MakeMCToyResponse(means, rng, h_meas, h_true);
          std::unique_ptr<RooUnfold> ut( use_bayes
                  ? (RooUnfold*) new RooUnfoldBayes(toyResp.get(), h_input, nIter)
                  : (RooUnfold*) new RooUnfoldBinByBin(toyResp.get(), h_input) );
          TH2* h_u = dynamic_cast<TH2*>( ut->Hunfold(RooUnfolding::kNoError) );
          if (!h_u) continue;
          flatten(h_u, xvec);
          diag_mc.add(xvec);
          if (bandWidth>0) band_mc.add(xvec);
          if ((t+1)%5==0) DBG(verbose, "  MC toy %d/%d", t+1, ntoys_mc);
        }
        diag_mc.finish(var_mc);
        if (bandWidth>0) {
          band_mc.finish(dummy_mean, dummy_varB, covBands_mc);
        }
      }
    } else {
      var_mc.assign(m, 0.0);
      covBands_mc.assign(bandWidth, std::vector<double>(m,0.0));
    }

    // ---------- Combine (Data ⊕ MC) & attach per-bin σ ----------
    std::vector<double> var_tot(m, 0.0);
    for (size_t i=0;i<m;++i) var_tot[i] = var_data[i] + var_mc[i];

    // Put σ into h_unfold bin errors
    size_t k=0;
    for (int jy=1; jy<=ny; ++jy)
      for (int ix=1; ix<=nx; ++ix, ++k)
        h_unfold->SetBinError(ix, jy, std::sqrt(std::max(0.0, var_tot[k])));

    // We'll also write banded covariances (Δ=1..bandWidth) if requested
    // by summing the data and MC bands element-wise.
    if (bandWidth>0) {
      TTree* ct = new TTree("cov_bands","Banded covariance (delta,i)->cov (Data+MC)");
      covTree = ct;
      Int_t delta=0, ix=0, iy=0; Float_t cov=0.f;
      ct->Branch("delta",&delta); ct->Branch("ix",&ix); ct->Branch("iy",&iy); ct->Branch("cov",&cov);

      for (int d=1; d<=bandWidth; ++d){
        const auto& Bd = covBands_data[d-1];
        const auto& Md = covBands_mc  [d-1];
        size_t idx=0;
        for (int jy=1; jy<=ny; ++jy){
          for (int ix1=1; ix1<=nx; ++ix1, ++idx){
            delta = d; ix = ix1; iy = jy;
            cov = (Float_t)( (idx < Bd.size()? Bd[idx]:0.0) + (idx < Md.size()? Md[idx]:0.0) );
            ct->Fill();
          }
        }
      }
      DBG(verbose, "[WRITE] Prepared cov_bands tree with Δ up to %d", bandWidth);
    }
  }

  // ---------------- Write outputs (sparse where appropriate) ----------------
  // Adjust output name if uncertainties requested
  TString outName(out_file);
  if (include_unc) {
    if (outName.EndsWith(".root")) outName.ReplaceAll(".root","_wStatDataMC.root");
    else outName += "_wStatDataMC.root";
  }

  DBG(verbose, "[PHASE] WRITE: out=%s (include_unc=%d)", outName.Data(), (int)include_unc);
  TFile fout(outName, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: Cannot create output file: " << outName << "\n";
    return 4;
  }
  fout.cd();
  if (covTree) { covTree->SetDirectory(&fout); covTree->Write(); }

  // Metadata
  TParameter<int>("use_bayes", (int)use_bayes).Write("use_bayes");
  TParameter<int>("nIter", nIter).Write("nIter");
  TParameter<int>("include_unc", (int)include_unc).Write("include_unc");
  TParameter<int>("ntoys_data", ntoys_data).Write("ntoys_data");
  TParameter<int>("ntoys_mc", ntoys_mc).Write("ntoys_mc");
  TParameter<int>("bandWidth", bandWidth).Write("bandWidth");
  TParameter<int>("rng_seed", rng_seed).Write("rng_seed");

  if (h_meas_data) WriteSmart((TH2*)h_meas_data->Clone(), "h_meas_data_input", fout);
  if (h_meas)      WriteSmart((TH2*)h_meas->Clone(),      "h_meas_mc",         fout);
  if (h_true)      WriteSmart((TH2*)h_true->Clone(),      "h_true_mc",         fout);
  if (h_mig)       WriteSmart((TH2*)h_mig->Clone(),       "h_mig_dense",       fout);

  // Always write unfolded sparsely
  WriteSmart(h_unfold, TString::Format("unfold_%s", tag.Data()), fout,
             /*eps=*/0.0, /*forceSparse=*/true);

  // If we created a banded covariance tree, write it now into file root
  TObject* obj = gDirectory->Get("cov_bands");
  if (obj) obj->Write();

  fout.Write();
  fout.Close();

  std::cout << "Unfolding done (" << (use_bayes ? "Bayes" : "BinByBin")
            << "). Wrote: " << outName << "\n"
            << "  - unfolded (sparse) under directory: unfold_" << tag << "/\n"
            << (include_unc ? "  - per-bin σ from Data⊕MC toys written in bin errors\n" : "")
            << (include_unc && bandWidth>0 ? "  - banded covariance tree: cov_bands\n" : "");
  return 0;
}


int perform_unfolding_unc_scan(bool use_bayes=true,
                               int  max_iter=8,           // go up to 8 iterations
                               bool include_unc=true,
                               int  ntoys_data=300,
                               int  ntoys_mc=300,
                               int  bandWidth=2,
                               int  seed_base=12345,      // will use seed_base+iter
                               const char* in_file="response_out.root",
                               const char* out_dir="unfold_iter_scan",
                               bool verbose = gVerboseDefault)
{
  if (max_iter < 1) max_iter = 1;
  // make the output folder (OS directory)
  gSystem->mkdir(out_dir, /*recursive=*/true);

  int status = 0;
  for (int iter=1; iter<=max_iter; ++iter) {
    DBG(verbose, "============================================================");
    DBG(verbose, "[SCAN] iter=%d/%d  seed=%d", iter, max_iter, seed_base+iter);
    TString outfile = TString::Format("%s/unfold_iter%02d.root", out_dir, iter);

    // Bump seed so toys are independent per iteration (helps comparisons)
    int rc = perform_unfolding_unc(use_bayes,
                                   /*nIter=*/iter,
                                   include_unc,
                                   ntoys_data,
                                   ntoys_mc,
                                   bandWidth,
                                   /*rng_seed=*/seed_base + iter,
                                   in_file,
                                   outfile.Data(),
                                   /*verbose=*/verbose);

    if (rc != 0) {
      Error("perform_unfolding_unc_scan",
            "Iteration %d failed with code %d (continuing).", iter, rc);
      status = rc; // remember last non-zero
    } else {
      DBG(verbose, "[SCAN] iter=%d OK → %s", iter, outfile.Data());
    }
  }
  return status;
}
