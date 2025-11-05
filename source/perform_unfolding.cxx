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
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldResponse.h"
#endif

#include <memory>
#include <string>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"

#include "TSystem.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"

#include "create_response_obj.cxx"  // provides load_response_from_ttree, axes_identical, constants, etc.

// Forward-declare the loader (already provided in create_response_obj.cxx; redeclare to be safe)
bool load_response_from_ttree(const char* infile,
                              RooUnfoldResponse*& outResp,
                              TH2*& outHmeas,
                              TH2*& outHtrue,
                              TH2*& outHmig,
                              TH2*& outHmeasData);

// -------- Helpers --------------------------------
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

// Make an empty TH2D with same (possibly variable) binning
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

// Unflatten a TH1 (x-fastest packing) into a TH2 with the same binning as `shape`
static TH2D* UnflattenLike(const TH2* shape, const TH1* flat, const char* name){
  if (!shape || !flat) return nullptr;
  const int nx = shape->GetNbinsX();
  const int ny = shape->GetNbinsY();
  auto* out = MakeEmptyTH2DLike(shape, name, name);
  int k = 0;
  for (int iy=1; iy<=ny; ++iy) {
    for (int ix=1; ix<=nx; ++ix) {
      ++k; // x-fastest packing
      const double v = (k<=flat->GetNbinsX()) ? flat->GetBinContent(k) : 0.0;
      const double e = (k<=flat->GetNbinsX()) ? flat->GetBinError  (k) : 0.0;
      out->SetBinContent(ix, iy, v);
      out->SetBinError  (ix, iy, e);
    }
  }
  out->SetDirectory(nullptr);
  return out;
}

// -----------------------------------------------------------------------------
// Now with a flag: use_bayes=true (Bayes), false -> BinByBin
int perform_unfolding(bool use_bayes=true,
                      int  nIter=5,
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

  std::unique_ptr<RooUnfold> unfold;
  TString tag;
  if (use_bayes) {
    unfold.reset(new RooUnfoldBayes(resp, h_input, nIter));
    tag = TString::Format("Bayes_iter%d", nIter);
  } else {
    unfold.reset(new RooUnfoldBinByBin(resp, h_input));
    tag = "BinByBin";
  }

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

  WriteSmart(h_unfold, TString::Format("unfold_%s", tag.Data()), fout,
             /*eps=*/0.0, /*forceSparse=*/true);

  fout.Write();
  fout.Close();

  std::cout << "Unfolding done (" << (use_bayes ? "Bayes" : "BinByBin")
            << "). Wrote: " << out_file << "\n"
            << "  - unfolded (sparse) under directory: unfold_" << tag << "/\n";
  return 0;
}


// -----------------------------------------------------------------------------
// Save φ grids: one PNG per xQ2 bin; z×pT subpads, φ on x-axis.
// This version fixes the packed-index mapping (off-by-one) and supports both
// layouts: packed on X or packed on Y.
//
// Assumptions:
//  - Packed index is ordered with φ fastest, then pT^2, then z.
//  - Packed index values stored in TH2 axes are represented as 1-wide,
//    integer-coded bins (axis bins are 1..N, centers at 0..N-1 if axis is [-0.5, N-0.5]).
// -----------------------------------------------------------------------------
static void SavePhiGridsPerIx(const TH2* h2,
                              const char* outDir,
                              const char* tag,
                              int   max_ix     = 16,
                              int   nZ         = N_Zbins,
                              int   nPt        = N_pTbins_with_overflow,
                              int   nPhi       = N_phiTrbins,
                              double emptyEps  = 0.0)
{
  if (!h2) return;
  gSystem->mkdir(outDir, /*recursive*/true);

  const int packedNeed = nZ * nPt * nPhi;
  const int nx = h2->GetNbinsX();
  const int ny = h2->GetNbinsY();

  enum class Layout { PackOnX, PackOnY, Unknown } layout = Layout::Unknown;

  // Prefer exact matches first, then fall back to heuristics
  if (nx == packedNeed) layout = Layout::PackOnX;
  else if (ny == packedNeed) layout = Layout::PackOnY;
  else if (nx > ny && nx >= packedNeed) layout = Layout::PackOnX;
  else if (ny > nx && ny >= packedNeed) layout = Layout::PackOnY;

  if (layout == Layout::Unknown) {
    std::cerr << "SavePhiGridsPerIx(" << tag << "): cannot detect layout. "
              << "nbinsX=" << nx << " nbinsY=" << ny
              << " expected packed=" << packedNeed << "\n";
    return;
  }

  // Number of xQ2 bins lives on the non-packed axis
  const int nIxq2 = (layout == Layout::PackOnX ? ny : nx);
  const int ixMax = std::min(max_ix, nIxq2);

  // Helper: convert (z,pt,phi) -> packed 0-based index, then to TH2 bin number (1-based)
  auto packed_bin_1based = [&](int iz, int ipt, int iph) -> int {
    // Inputs are 1-based (iz in 1..nZ, ipt in 1..nPt, iph in 1..nPhi)
    const int comp0 = ((iz-1) * nPt + (ipt-1)) * nPhi + (iph-1); // 0..packedNeed-1
    return comp0 + 2; // TH2 axis bin numbering is 1..N
  };

  for (int ixq2 = 1; ixq2 <= ixMax; ++ixq2) {
    const int w = 260 * nPt;
    const int h = 240 * nZ;
    TCanvas c(Form("c_%s_ix%02d", tag, ixq2), Form("%s ix=%d", tag, ixq2), w, h);
    c.Divide(nPt, nZ, 0.001, 0.001);

    std::vector<std::unique_ptr<TH1D>> keep;
    keep.reserve(nZ * nPt);

    int pad = 0;
    for (int iz = 1; iz <= nZ; ++iz) {
      for (int ipt = 1; ipt <= nPt; ++ipt) {
        ++pad; c.cd(pad);
        gPad->SetTicks(1,1);
        gPad->SetLeftMargin(0.18);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.22);
        gPad->SetTopMargin(0.10);

        auto hphi = std::make_unique<TH1D>(
          Form("hphi_%s_ix%02d_z%02d_pt%02d", tag, ixq2, iz, ipt),
          ";#phi_{Trento} [deg];Counts", nPhi, 0.0, 360.0
        );
        hphi->Sumw2(); hphi->SetDirectory(nullptr);

        bool any = false;
        for (int iph = 1; iph <= nPhi; ++iph) {
          const int binPacked = packed_bin_1based(iz, ipt, iph); // 1..packedNeed

          double v = 0.0, e = 0.0;
          if (layout == Layout::PackOnX) {
            if (binPacked >= 1 && binPacked <= nx) {
              v = h2->GetBinContent(binPacked, ixq2);
              e = h2->GetBinError  (binPacked, ixq2);
            }
          } else { // PackOnY
            if (binPacked >= 1 && binPacked <= ny) {
              v = h2->GetBinContent(ixq2, binPacked);
              e = h2->GetBinError  (ixq2, binPacked);
            }
          }

          hphi->SetBinContent(iph, v);
          hphi->SetBinError  (iph, e);
          if (std::fabs(v) > emptyEps || std::fabs(e) > emptyEps) any = true;
        }

        double ymax = hphi->GetMaximum();
        if (!(ymax > 0)) ymax = 1.0;
        hphi->SetMaximum(1.15 * ymax);
        if (hphi->GetMinimum() < 0) hphi->SetMinimum(0);
        hphi->SetStats(0);
        hphi->GetXaxis()->SetLabelSize(0.08);
        hphi->GetYaxis()->SetLabelSize(0.08);
        hphi->GetXaxis()->SetTitleSize(0.09);
        hphi->GetYaxis()->SetTitleSize(0.09);
        hphi->GetYaxis()->SetTitleOffset(1.0);

        hphi->Draw("AXIS");
        hphi->Draw("E1 SAME");

        TLatex lab; lab.SetNDC(); lab.SetTextSize(0.10);
        lab.DrawLatex(0.15, 0.88, Form("z%02d pT%02d", iz, ipt));

        keep.emplace_back(std::move(hphi));
      }
    }

    c.cd(1);
    TLatex head; head.SetNDC(); head.SetTextSize(0.06);
    head.DrawLatex(0.02, 0.97, Form("%s  |  xQ^{2} bin = %d  |  layout = %s",
                                    tag, ixq2,
                                    (layout==Layout::PackOnX ? "PackOnX" : "PackOnY")));

    c.Update();
    c.SaveAs(Form("%s/%s_ix%02d.png", outDir, tag, ixq2));
  }
}





// -----------------------------------------------------------------------------
// Manual 2D Bin-by-Bin unfolding (uses MC truth/reco scale factors per bin).
int perform_unfolding_manual_bbb(const char* in_file          ="response_out.root",
                                 const char* out_file         ="unfold_out_manual_bbb.root",
                                 bool  include_mc_stat        = false,
                                 double eps                   = 0.0,
                                 // --- plotting options ---
                                 bool   make_phi_plots        = true,
                                 const char* plot_dir_base    = "phi_plots_manual_bbb",
                                 int    max_ix_plot           = 17)
{
  RooUnfoldResponse* resp = nullptr;
  TH2* h_meas_mc   = nullptr; // MC reco (measured space) prototype
  TH2* h_true_mc   = nullptr; // MC truth (truth space) prototype
  TH2* h_mig       = nullptr; // optional migration matrix (unused here)
  TH2* h_meas_data = nullptr; // data in measured space

  if (!load_response_from_ttree(in_file, resp, h_meas_mc, h_true_mc, h_mig, h_meas_data)) {
    std::cerr << "ERROR: Failed to load response from: " << in_file << "\n";
    return 1;
  }

  TH2* h_data = h_meas_data ? h_meas_data : h_meas_mc;
  if (!h_data || !h_meas_mc || !h_true_mc) {
    std::cerr << "ERROR: Need h_data, h_meas_mc, and h_true_mc.\n";
    return 2;
  }

  auto nonempty = [](TH2* h) -> bool {
    return h && (h->GetEntries() > 0 || h->Integral(0, h->GetNbinsX()+1, 0, h->GetNbinsY()+1) != 0.0);
  };

  std::cout << "DBG: prototypes  h_meas_mc int=" << (h_meas_mc? h_meas_mc->Integral(): -999)
            << "  h_true_mc int="              << (h_true_mc? h_true_mc->Integral(): -999) << "\n";

  // Pull trained histograms from the response; RooUnfold often stores them flattened (TH1)
  TH2* h_meas_trained = nullptr;
  TH2* h_true_trained = nullptr;
  if (resp) {
    TH1* hm1 = resp->Hmeasured();
    TH1* ht1 = resp->Htruth();
    std::cout << "DBG: resp Hmeasured class=" << (hm1?hm1->ClassName():"(null)")
              << "  nbins=" << (hm1?hm1->GetNbinsX():-1)
              << "  integral=" << (hm1?hm1->Integral():0) << "\n";
    std::cout << "DBG: resp Htruth    class=" << (ht1?ht1->ClassName():"(null)")
              << "  nbins=" << (ht1?ht1->GetNbinsX():-1)
              << "  integral=" << (ht1?ht1->Integral():0) << "\n";

    if (auto* h2m = dynamic_cast<TH2*>(hm1)) {
      h_meas_trained = (TH2*) h2m->Clone("h_meas_mc_trained");
      h_meas_trained->SetDirectory(nullptr);
    } else if (hm1) {
      h_meas_trained = UnflattenLike(h_meas_mc, hm1, "h_meas_mc_trained");
    }

    if (auto* h2t = dynamic_cast<TH2*>(ht1)) {
      h_true_trained = (TH2*) h2t->Clone("h_true_mc_trained");
      h_true_trained->SetDirectory(nullptr);
    } else if (ht1) {
      h_true_trained = UnflattenLike(h_true_mc, ht1, "h_true_mc_trained");
    }
  }

  // Choose the actual sources we will use
  TH2* h_meas_src = (nonempty(h_meas_mc) ? h_meas_mc : h_meas_trained);
  TH2* h_true_src = (nonempty(h_true_mc) ? h_true_mc : h_true_trained);

  if (!h_meas_src || !h_true_src) {
    std::cerr << "ERROR: No non-empty MC histograms available (meas/truth)." << std::endl;
    return 2;
  }

  // Enforce bin-by-bin preconditions on the *sources*
  if (!axes_identical(*h_data, *h_meas_src)) {
    std::cerr << "ERROR: Data measured binning must match MC measured binning for bin-by-bin.\n";
    return 3;
  }
  if (!axes_identical(*h_true_src, *h_meas_src)) {
    std::cerr << "ERROR: MC truth and MC measured binnings must be identical for bin-by-bin.\n";
    return 4;
  }

  // --- Create outputs ---
  std::unique_ptr<TH2D> h_unfold (
    MakeEmptyTH2DLike(h_true_src,
                      "unfold_ManualBinByBin",
                      "Unfolded spectrum (Manual Bin-by-Bin)")
  );

  std::unique_ptr<TH2D> h_sf (
    MakeEmptyTH2DLike(h_true_src,
                      "h_bbb_scale_factor",
                      "Scale factor T_MC/M_MC per bin")
  );

  std::unique_ptr<TH2D> h_zeroeff (
    MakeEmptyTH2DLike(h_true_src,
                      "h_bbb_zeroeff",
                      "Mask: 1 where M_MC<=eps, else 0")
  );

  // --- Compute per-bin
  const int nx = h_true_src->GetNbinsX();
  const int ny = h_true_src->GetNbinsY();

  for (int jy=1; jy<=ny; ++jy) {
    for (int ix=1; ix<=nx; ++ix) {
      const double T  = h_true_src->GetBinContent(ix, jy);
      const double eT = h_true_src->GetBinError  (ix, jy);
      const double M  = h_meas_src->GetBinContent(ix, jy);
      const double eM = h_meas_src->GetBinError  (ix, jy);
      const double D  = h_data    ->GetBinContent(ix, jy);
      const double eD = h_data    ->GetBinError  (ix, jy);

      if (std::fabs(M) <= eps) {
        h_sf->SetBinContent(ix,jy, 0.0);
        h_sf->SetBinError  (ix,jy, 0.0);
        h_unfold->SetBinContent(ix,jy, 0.0);
        h_unfold->SetBinError  (ix,jy, 0.0);
        h_zeroeff->SetBinContent(ix,jy, 1.0);
        continue;
      }

      const double SF  = T / M;
      const double TS  = SF * D;

      // Errors
      double eTS = 0.0;
      if (!include_mc_stat) {
        eTS = std::fabs(SF) * eD; // data-only
      } else {
        // Var(TS) = (SF^2)*eD^2 + (D/M)^2*eT^2 + (T*D/M^2)^2*eM^2 (no covariances)
        const double termD = (SF*eD);
        const double termT = (D / M) * eT;
        const double termM = (T * D / (M*M)) * eM;
        eTS = std::sqrt(termD*termD + termT*termT + termM*termM);
      }

      // Store
      h_sf->SetBinContent(ix,jy, SF);
      if (include_mc_stat) {
        const double relT = (std::fabs(T)>eps) ? eT/std::fabs(T) : 0.0;
        const double relM = (std::fabs(M)>eps) ? eM/std::fabs(M) : 0.0;
        const double eSF  = std::fabs(SF) * std::sqrt(relT*relT + relM*relM);
        h_sf->SetBinError(ix,jy, eSF);
      } else {
        h_sf->SetBinError(ix,jy, 0.0);
      }

      h_unfold->SetBinContent(ix,jy, TS);
      h_unfold->SetBinError  (ix,jy, eTS);
      h_zeroeff->SetBinContent(ix,jy, 0.0);
    }
  }

  // Titles/axes
  if (h_true_src->GetXaxis() && h_true_src->GetYaxis()) {
    h_unfold->GetXaxis()->SetTitle(h_true_src->GetXaxis()->GetTitle());
    h_unfold->GetYaxis()->SetTitle(h_true_src->GetYaxis()->GetTitle());
  }

  // --- Write results ---
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: Cannot create output file: " << out_file << "\n";
    return 5;
  }
  if (h_meas_data) WriteSmart((TH2*)h_meas_data->Clone(), "h_meas_data_input", fout);

  // IMPORTANT: write the sources we actually used (not the empty placeholders)
  // If you want to force dense write instead of sparse, replace WriteSmart(...) with the
  // small lambda below:
  // auto WriteDense = [&](TH2* h,const char* name){ if (!h) return; h->SetName(name); h->SetDirectory(&fout); h->Write(); };
  WriteSmart((TH2*)h_meas_src->Clone(), "h_meas_mc", fout);
  WriteSmart((TH2*)h_true_src->Clone(), "h_true_mc", fout);

  if (h_mig) WriteSmart((TH2*)h_mig->Clone(), "h_mig_dense", fout);

  WriteSmart(h_unfold.release(), "unfold_ManualBinByBin", fout, /*eps=*/0.0, /*forceSparse=*/true);
  WriteSmart(h_sf.release(),      "h_bbb_scale_factor",   fout);
  WriteSmart(h_zeroeff.release(), "h_bbb_zeroeff",        fout);

  std::cout << "DBG: write h_meas_mc int(all bins incl. u/o) = "
            << h_meas_src->Integral(0,h_meas_src->GetNbinsX()+1,0,h_meas_src->GetNbinsY()+1) << "\n";
  std::cout << "DBG: write h_true_mc int(all bins incl. u/o) = "
            << h_true_src->Integral(0,h_true_src->GetNbinsX()+1,0,h_true_src->GetNbinsY()+1) << "\n";

  fout.Write();
  fout.Close();

  std::cout << "Manual bin-by-bin unfolding done. Wrote: " << out_file << "\n"
            << "  - unfolded (sparse) under directory: unfold_ManualBinByBin/\n"
            << "  - also wrote h_bbb_scale_factor and h_bbb_zeroeff\n";

  // Optional φ-plots — use the sources we actually used so plots aren't empty
  if (make_phi_plots) {
    if (h_meas_data) {
      TString d = TString::Format("%s/meas_data", plot_dir_base);
      SavePhiGridsPerIx(h_meas_data, d.Data(), "meas_data", max_ix_plot,
                        N_Zbins, N_pTbins_with_overflow, N_phiTrbins, /*emptyEps=*/0.0);
    }
    if (h_meas_src) {
      TString d = TString::Format("%s/meas_mc", plot_dir_base);
      SavePhiGridsPerIx(h_meas_src, d.Data(), "meas_mc", max_ix_plot,
                        N_Zbins, N_pTbins_with_overflow, N_phiTrbins, /*emptyEps=*/0.0);
    }
    if (h_true_src) {
      TString d = TString::Format("%s/true_mc", plot_dir_base);
      SavePhiGridsPerIx(h_true_src, d.Data(), "true_mc", max_ix_plot,
                        N_Zbins, N_pTbins_with_overflow, N_phiTrbins, /*emptyEps=*/0.0);
    }
    std::cout << "Saved φ-plots under: " << plot_dir_base
              << "  (meas_data/, meas_mc/, true_mc/)\n";
  }

  return 0;
}
