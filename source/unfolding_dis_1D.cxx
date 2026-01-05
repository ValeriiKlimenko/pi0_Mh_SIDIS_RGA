// unfold_dis_xQ2_bayes.cxx
//
// Bayesian unfolding of the 1D xQ2 distribution using histograms produced by
// make_dis_histos_per_file.cxx.
//
// Typical workflow:
//   1) Run the histogram maker (standalone exe):
//        ./make_dis_histos_per_file 1 1 1
//      which produces (by default) in: unfolding_dis/
//        - dis_data_hists_all.root  (h_xQ2_data)
//        - dis_rec_hists_all.root   (h_xQ2_rec, h_xQ2_response)
//        - dis_gen_hists_all.root   (h_xQ2_gen)
//
//   2) From ROOT, load and run this macro:
//        .L unfold_dis_xQ2_bayes.cxx+
//        unfold_dis_xQ2_bayes_root(4, "unfolding_dis/dis_unfold_bayes.root");
//
// The core worker is:
//   int unfold_dis_xQ2_bayes(int nIter, const char* out_file);

#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TString.h"

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

using std::cout;
using std::cerr;
using std::endl;

// Simple axis check
static bool axes_identical_1D(const TH1* a, const TH1* b) {
  if (!a || !b) return false;
  const TAxis* A = a->GetXaxis();
  const TAxis* B = b->GetXaxis();
  return (A->GetNbins() == B->GetNbins()) &&
         (std::fabs(A->GetXmin() - B->GetXmin()) < 1e-9) &&
         (std::fabs(A->GetXmax() - B->GetXmax()) < 1e-9);
}

// Transpose the response matrix so that:
//   input:  h_in (X = true, Y = reco)
//   output: h_out (X = reco, Y = true)
// RooUnfold expects X = measured (reco), Y = truth.
static TH2D* transpose_trueReco_to_recoTrue(const TH2D* h_in,
                                            const char* new_name,
                                            const char* new_title)
{
  if (!h_in) return nullptr;

  int nTrue = h_in->GetNbinsX();
  int nReco = h_in->GetNbinsY();

  double true_lo = h_in->GetXaxis()->GetXmin();
  double true_hi = h_in->GetXaxis()->GetXmax();
  double reco_lo = h_in->GetYaxis()->GetXmin();
  double reco_hi = h_in->GetYaxis()->GetXmax();

  // New hist: X=reco, Y=true
  TH2D* h_out = new TH2D(new_name, new_title,
                         nReco, reco_lo, reco_hi,
                         nTrue, true_lo, true_hi);
  h_out->Sumw2();
  h_out->SetDirectory(nullptr);

  for (int ix_true = 1; ix_true <= nTrue; ++ix_true) {
    for (int iy_reco = 1; iy_reco <= nReco; ++iy_reco) {
      double c  = h_in->GetBinContent(ix_true, iy_reco);
      double ce = h_in->GetBinError  (ix_true, iy_reco);
      if (c == 0.0 && ce == 0.0) continue;

      int ix_reco_new = iy_reco; // X in output
      int iy_true_new = ix_true; // Y in output

      h_out->SetBinContent(ix_reco_new, iy_true_new, c);
      h_out->SetBinError  (ix_reco_new, iy_true_new, ce);
    }
  }

  return h_out;
}

// --- NEW: simple manual bin-by-bin unfolding: data * (gen / rec)
// h_unfolded = data * (gen/rec)
// h_correction = gen/rec
static void manual_bin_bin(const TH1D* h_data,
                           const TH1D* h_gen,
                           const TH1D* h_rec,
                           TH1D*& h_unfolded,
                           TH1D*& h_correction)
{
  h_unfolded  = nullptr;
  h_correction = nullptr;

  if (!h_data || !h_gen || !h_rec) {
    cerr << "ERROR [manual_bin_bin]: null input histogram" << endl;
    return;
  }

  if (!axes_identical_1D(h_data, h_gen) || !axes_identical_1D(h_data, h_rec)) {
    cerr << "WARNING [manual_bin_bin]: input histograms have different binning." << endl;
  }

  // Correction factor gen/rec
  TH1D* h_ratio = dynamic_cast<TH1D*>(h_gen->Clone("h_xQ2_gen_over_rec_mc"));
  if (!h_ratio) {
    cerr << "ERROR [manual_bin_bin]: failed to clone h_gen" << endl;
    return;
  }
  h_ratio->SetDirectory(nullptr);
  h_ratio->SetTitle("MC truth / MC reco;bin_xBQ2_Valerii;ratio");
  h_ratio->Sumw2();
  h_ratio->Divide(h_rec);  // bin-by-bin gen/rec

  // Manual unfolded: data * (gen/rec)
  TH1D* h_manual = dynamic_cast<TH1D*>(h_data->Clone("h_xQ2_manual_bin_bin"));
  if (!h_manual) {
    cerr << "ERROR [manual_bin_bin]: failed to clone h_data" << endl;
    delete h_ratio;
    return;
  }
  h_manual->SetDirectory(nullptr);
  h_manual->SetTitle("Manual bin-by-bin unfolded DIS xQ2;bin_xBQ2_Valerii;Events");
  h_manual->Sumw2();
  h_manual->Multiply(h_ratio);

  h_unfolded  = h_manual;
  h_correction = h_ratio;
}

// Core unfolding routine (can be called from ROOT or from another C++ file)
int unfold_dis_xQ2_bayes(int nIter = 4,
                         const char* out_file = "unfolding_dis/dis_unfold_bayes.root")
{
  cout << "=== DIS xQ2 Bayesian unfolding ===" << endl;
  cout << "  nIter    = " << nIter << endl;
  cout << "  out_file = " << out_file << endl;

  if (nIter < 1) {
    cerr << "ERROR: nIter must be >= 1" << endl;
    return 100;
  }

  // -------------------------------------------------------------------------
  // 1) Open input files & grab histograms
  // -------------------------------------------------------------------------
  const char* data_fname = "unfolding_dis/dis_data_hists_all.root";
  const char* rec_fname  = "unfolding_dis/dis_rec_hists_all.root";
  const char* gen_fname  = "unfolding_dis/dis_gen_hists_all.root";

  TFile fData(data_fname, "READ");
  TFile fRec (rec_fname,  "READ");
  TFile fGen (gen_fname,  "READ");

  if (fData.IsZombie() || fRec.IsZombie() || fGen.IsZombie()) {
    cerr << "ERROR: could not open one of the input files:" << endl;
    cerr << "  " << data_fname << endl;
    cerr << "  " << rec_fname  << endl;
    cerr << "  " << gen_fname  << endl;
    return 1;
  }

  TH1D* h_data = nullptr;    // measured data to unfold
  TH1D* h_rec  = nullptr;    // MC reco
  TH1D* h_gen  = nullptr;    // MC truth
  TH2D* h_resp_in = nullptr; // response as written by make_dis_histos_per_file

  fData.GetObject("h_xQ2_data", h_data);
  fRec .GetObject("h_xQ2_rec",  h_rec);
  fRec .GetObject("h_xQ2_response", h_resp_in);
  fGen .GetObject("h_xQ2_gen",  h_gen);

  if (!h_data || !h_rec || !h_gen || !h_resp_in) {
    cerr << "ERROR: missing histograms in input files." << endl;
    cerr << "  data: h_xQ2_data        = " << (void*)h_data    << endl;
    cerr << "  rec : h_xQ2_rec         = " << (void*)h_rec     << endl;
    cerr << "  rec : h_xQ2_response    = " << (void*)h_resp_in << endl;
    cerr << "  gen : h_xQ2_gen         = " << (void*)h_gen     << endl;
    return 2;
  }

  // Detach from files so they survive after closing
  h_data->SetDirectory(nullptr);
  h_rec ->SetDirectory(nullptr);
  h_gen ->SetDirectory(nullptr);
  h_resp_in->SetDirectory(nullptr);

  fData.Close();
  fRec.Close();
  fGen.Close();

  // -------------------------------------------------------------------------
  // 2) Basic consistency checks
  // -------------------------------------------------------------------------
  if (!axes_identical_1D(h_data, h_rec)) {
    cerr << "WARNING: h_xQ2_data and h_xQ2_rec have different X-axes!" << endl;
  }
  if (!axes_identical_1D(h_data, h_gen)) {
    cerr << "WARNING: h_xQ2_data and h_xQ2_gen have different X-axes!" << endl;
  }

  cout << "Input integrals:" << endl;
  cout << "  Data (measured) : " << h_data->Integral(0, h_data->GetNbinsX()+1) << endl;
  cout << "  MC reco         : " << h_rec ->Integral(0, h_rec ->GetNbinsX()+1) << endl;
  cout << "  MC truth        : " << h_gen ->Integral(0, h_gen ->GetNbinsX()+1) << endl;

  // -------------------------------------------------------------------------
  // 3) Transpose response so it becomes (reco X, truth Y) for RooUnfold
  //    Original: x-axis = true bin_xBQ2_Valerii
  //              y-axis = reco bin_xBQ2_Valeriim
  // -------------------------------------------------------------------------
  TH2D* h_resp_rt = transpose_trueReco_to_recoTrue(
      h_resp_in,
      "h_xQ2_response_reco_true",
      "DIS response;reco bin_xBQ2_Valeriim;true bin_xBQ2_Valerii");

  if (!h_resp_rt) {
    cerr << "ERROR: failed to transpose response matrix." << endl;
    return 3;
  }

  cout << "Response matrix (reco X, truth Y) integral = "
       << h_resp_rt->Integral() << endl;

  // -------------------------------------------------------------------------
  // 4) Build RooUnfoldResponse from MC histograms and migration matrix
  // -------------------------------------------------------------------------
  // Note: constructor expects X = measured (reco), Y = truth.
  RooUnfoldResponse response(h_rec, h_gen, h_resp_rt,
                             "response", "DIS xQ2 response");

  // -------------------------------------------------------------------------
  // 5) Bayes unfolding sweep over iterations 1..nIter
  // -------------------------------------------------------------------------

  // To store χ² vs truth and χ² change vs previous iteration
  TH1D hChi2Truth("hChi2Truth",
                  "Bayes #chi^{2} vs MC truth;iteration;#chi^{2}",
                  nIter, 0.5, nIter + 0.5);
  TH1D hChi2Change("hChi2Change",
                   "Bayes #chi^{2} change between iterations;iteration;#chi^{2}",
                   nIter, 0.5, nIter + 0.5);
  hChi2Truth.SetDirectory(nullptr);
  hChi2Change.SetDirectory(nullptr);

  TH1*  hPrevUnfold = nullptr;        // previous iteration's unfolded hist
  TH1D* h_final     = nullptr;        // unfolded spectrum at final iteration

  cout << "=== Bayes iteration scan ===" << endl;

  for (int i = 1; i <= nIter; ++i) {
    RooUnfoldBayes unfold(&response, h_data, i);
    unfold.SetVerbose(0);

    // Use Hunfold with explicit error treatment (your RooUnfold version uses RooUnfolding enum)
    TH1* h_unfold_raw = unfold.Hunfold(RooUnfolding::kErrors);
    if (!h_unfold_raw) {
      cerr << "ERROR: RooUnfoldBayes::Hunfold returned null at iteration " << i << endl;
      if (hPrevUnfold) { delete hPrevUnfold; }
      if (h_final)     { delete h_final;     }
      return 4;
    }

    TString hname = Form("h_xQ2_unfold_bayes_iter%d", i);
    TH1D* h_iter = dynamic_cast<TH1D*>(h_unfold_raw->Clone(hname));
    if (!h_iter) {
      cerr << "ERROR: could not cast unfolded spectrum to TH1D at iteration " << i << endl;
      if (hPrevUnfold) { delete hPrevUnfold; }
      if (h_final)     { delete h_final;     }
      return 5;
    }
    h_iter->SetDirectory(nullptr);
    h_iter->SetTitle(
      Form("Unfolded DIS xQ2 (Bayes, %d iterations);bin_xBQ2_Valerii;Events", i));

    // χ² vs MC truth for this iteration
    double chi2_truth = unfold.Chi2(h_gen, RooUnfolding::kCovariance);
    hChi2Truth.SetBinContent(i, chi2_truth);

    // χ² of change relative to previous iteration (if any)
    double chi2_change = -1.0;
    if (hPrevUnfold) {
      chi2_change = unfold.Chi2(hPrevUnfold, RooUnfolding::kCovariance);
      hChi2Change.SetBinContent(i, chi2_change);
    }

    cout << "  Iter " << i
         << " : chi2_vs_truth = " << chi2_truth;
    if (hPrevUnfold) {
      cout << " , chi2_change_vs_prev = " << chi2_change;
    }
    cout << endl;

    // Update previous-unfolded spectrum
    if (hPrevUnfold) {
      delete hPrevUnfold;
      hPrevUnfold = nullptr;
    }
    hPrevUnfold = (TH1*) h_iter->Clone("h_prev_unfold");
    hPrevUnfold->SetDirectory(nullptr);

    // Keep a copy of the final iteration spectrum for output
    if (h_final) {
      delete h_final;
      h_final = nullptr;
    }
    h_final = (TH1D*) h_iter->Clone("h_xQ2_unfold_bayes");
    h_final->SetDirectory(nullptr);

    // This local clone is no longer needed
    delete h_iter;
  }

  if (hPrevUnfold) {
    delete hPrevUnfold;
    hPrevUnfold = nullptr;
  }

  if (!h_final) {
    cerr << "ERROR: no final unfolded spectrum was created." << endl;
    return 6;
  }

  cout << "Final unfolded integral (iter " << nIter << ") : "
       << h_final->Integral(0, h_final->GetNbinsX()+1) << endl;

  // For convenience, also quote final iteration χ² vs truth:
  double chi2_final = hChi2Truth.GetBinContent(nIter);
  cout << "Final chi2 (iter " << nIter << " vs MC truth) = " << chi2_final << endl;

  // -------------------------------------------------------------------------
  // 5b) NEW: manual bin-by-bin unfolding and correction histograms
  // -------------------------------------------------------------------------
  TH1D* h_manual_binbin = nullptr;
  TH1D* h_gen_over_rec  = nullptr;
  manual_bin_bin(h_data, h_gen, h_rec, h_manual_binbin, h_gen_over_rec);

  if (h_manual_binbin) {
    cout << "Manual bin-by-bin unfolded integral : "
         << h_manual_binbin->Integral(0, h_manual_binbin->GetNbinsX()+1) << endl;
  }

  // -------------------------------------------------------------------------
  // 6) Write outputs
  // -------------------------------------------------------------------------
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    cerr << "ERROR: cannot create output file: " << out_file << endl;
    delete h_final;
    if (h_manual_binbin) delete h_manual_binbin;
    if (h_gen_over_rec)  delete h_gen_over_rec;
    return 7;
  }

  // Make clones with clean names for output
  TH1D* h_data_out = dynamic_cast<TH1D*>(h_data->Clone("h_xQ2_data_input"));
  TH1D* h_rec_out  = dynamic_cast<TH1D*>(h_rec ->Clone("h_xQ2_rec_mc"));
  TH1D* h_gen_out  = dynamic_cast<TH1D*>(h_gen ->Clone("h_xQ2_gen_mc"));
  TH2D* h_resp_out = dynamic_cast<TH2D*>(h_resp_rt->Clone("h_xQ2_response_reco_true"));

  if (h_data_out)      h_data_out     ->SetDirectory(&fout);
  if (h_rec_out)       h_rec_out      ->SetDirectory(&fout);
  if (h_gen_out)       h_gen_out      ->SetDirectory(&fout);
  if (h_resp_out)      h_resp_out     ->SetDirectory(&fout);
  h_final->SetDirectory(&fout);
  if (h_manual_binbin) h_manual_binbin->SetDirectory(&fout);
  if (h_gen_over_rec)  h_gen_over_rec ->SetDirectory(&fout);

  if (h_data_out)      h_data_out     ->Write();
  if (h_rec_out)       h_rec_out      ->Write();
  if (h_gen_out)       h_gen_out      ->Write();
  if (h_resp_out)      h_resp_out     ->Write();
  h_final->Write("h_xQ2_unfold_bayes");

  // --- NEW: write manual bin-by-bin and correction histograms
  if (h_manual_binbin) h_manual_binbin->Write("h_xQ2_manual_bin_bin");
  if (h_gen_over_rec)  h_gen_over_rec ->Write("h_xQ2_gen_over_rec_mc");

  // χ² vs iteration histograms
  hChi2Truth.Write("hChi2Truth");
  hChi2Change.Write("hChi2Change");

  // Store final-iter chi2 as a TParameter for quick access
  TParameter<double> pChi2("chi2_vs_truth_final", chi2_final);
  pChi2.Write();

  fout.Write();
  fout.Close();

  cout << "Wrote unfolding results to: " << out_file << endl;
  cout << "=== Done ===" << endl;

  // h_final and others are now owned by the file
  return 0;
}

// ROOT-friendly wrapper (no main)
void unfold_dis_xQ2_bayes_root(int nIter = 4,
                               const char* out_file = "unfolding_dis/dis_unfold_bayes.root")
{
  unfold_dis_xQ2_bayes(nIter, out_file);
}
