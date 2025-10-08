// perform_unfolding_sparse.cxx  (safe for huge 2D binning like 16000x16000)

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"

// RooUnfold headers
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"  // <-- added
#include "RooUnfoldResponse.h"
#endif

#include <memory>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"

#include "create_response_obj.cxx"

// Forward-declare the loader
bool load_response_from_ttree(const char* infile,
                              RooUnfoldResponse*& outResp,
                              TH2*& outHmeas,
                              TH2*& outHtrue,
                              TH2*& outHmig,
                              TH2*& outHmeasData);

// -------- Helpers (unchanged) --------------------------------
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

static bool IsHuge(const TH2* h, long long maxBins = 20000000LL) {
  const long long nb = 1LL * h->GetNbinsX() * h->GetNbinsY();
  return nb > maxBins;
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

// -----------------------------------------------------------------------------
// Now with a flag: use_bayes=true (Bayes), false -> BinByBin
int perform_unfolding(bool use_bayes=true,
                      int  nIter=3,
                      const char* in_file  ="response_out.root",
                      const char* out_file ="unfold_out.root")
{
  RooUnfoldResponse* resp = nullptr;
  TH2* h_meas = nullptr;     // MC reco (diagnostic)
  TH2* h_true = nullptr;     // MC truth (diagnostic)
  TH2* h_mig  = nullptr;     // migration (if rebuilt or stored)
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

  if (h_meas && !axes_identical(*h_input, *h_meas)) {
    std::cerr << "WARNING: Input measured histogram binning != response measured binning.\n";
  }

  // Build chosen unfold object
  std::unique_ptr<RooUnfold> unfold;
  TString tag;
  if (use_bayes) {
    unfold.reset(new RooUnfoldBayes(resp, h_input, nIter));
    tag = TString::Format("Bayes_iter%d", nIter);
  } else {
    unfold.reset(new RooUnfoldBinByBin(resp, h_input));
    tag = "BinByBin";
  }

  // Diagonal (per-bin) errors only
  auto errModeHist = RooUnfolding::kErrors;
  TH2* h_unfold = dynamic_cast<TH2*>(unfold->Hunfold(errModeHist));
  if (!h_unfold) {
    std::cerr << "ERROR: Unfolding failed (null output from Hunfold).\n";
    return 3;
  }

  h_unfold->SetName(TString::Format("unfold_%s", tag.Data()));
  h_unfold->SetTitle(TString::Format("Unfolded spectrum (%s);%s;%s",
                          tag.Data(),
                          h_true ? h_true->GetXaxis()->GetTitle() : "true X",
                          h_true ? h_true->GetYaxis()->GetTitle() : "true Y"));

  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: Cannot create output file: " << out_file << "\n";
    return 4;
  }
  fout.cd();

  if (h_meas_data) WriteSmart((TH2*)h_meas_data->Clone(), "h_meas_data_input", fout);
  if (h_meas)      WriteSmart((TH2*)h_meas->Clone(),      "h_meas_mc",         fout);
  if (h_true)      WriteSmart((TH2*)h_true->Clone(),      "h_true_mc",         fout);
  if (h_mig)       WriteSmart((TH2*)h_mig->Clone(),       "h_mig_dense",       fout);

  // Always write unfolded sparsely (safe for huge grids)
  WriteSmart(h_unfold, TString::Format("unfold_%s", tag.Data()), fout,
             /*eps=*/0.0, /*forceSparse=*/true);

  fout.Write();
  fout.Close();

  std::cout << "Unfolding done (" << (use_bayes ? "Bayes" : "BinByBin")
            << "). Wrote: " << out_file << "\n"
            << "  - unfolded (sparse) under directory: unfold_" << tag << "/\n";
  return 0;
}
