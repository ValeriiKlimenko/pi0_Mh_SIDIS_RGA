// File: draw_closure_1D.C
//
// Quick visualization + numeric QC for the 1D closure test.
//
// Expected input (from RunClosureTestBayes):
//   closure_test_bayes.root
//     h_meas_mc         (MC reco, used as pseudo-data)
//     h_true_mc         (MC truth)
//     h_unfold_closure  (unfolded MC reco)
//     h_closure_ratio   (unfolded / truth)
//     h_closure_diff    (unfolded - truth)
//
// or from ResponseQA:
//   qa_response.root
//     h_measured_trained
//     h_truth_trained
//     h_unfold_closure
//     h_closure_ratio
//
// What this macro does:
//   * Draw spectra: truth vs reco vs unfolded
//   * Draw ratio unfolded/truth with bands around 1
//   * Compute chi2/ndf (unfolded vs truth) and max |ratio - 1|
//   * Print a PASS/FAIL message based on user-defined tolerances
//
// Usage examples:
//   root -l -b -q 'draw_closure_1D.C("closure_test_bayes.root")'
//   root -l -b -q 'draw_closure_1D.C("qa_response.root","closure_plots_1D",0.05,2.0)'
//
// Arguments:
//   closureFile      : ROOT file from closure test
//   outDir           : directory for PNGs
//   ratioTolerance   : allowed |ratio-1| per bin (default 0.1 -> 10%)
//   chi2PerNdfMax    : allowed chi2/ndf (default 3.0)
//

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <cmath>
#include <algorithm>

using std::cout;
using std::cerr;
using std::endl;

// ----------------------------------------------------------------------
// Helper to fetch histograms by either "new" or "QA" names
// ----------------------------------------------------------------------
static TH1D* GetHist1D(TFile& f, const char* name1, const char* name2 = nullptr)
{
  TH1D* h = nullptr;
  if (name1 && *name1) {
    h = dynamic_cast<TH1D*>(f.Get(name1));
    if (h) return h;
  }
  if (name2 && *name2) {
    h = dynamic_cast<TH1D*>(f.Get(name2));
    if (h) return h;
  }
  return nullptr;
}

// ----------------------------------------------------------------------
// Main driver
// ----------------------------------------------------------------------
void draw_closure_1D(const char* closureFile   = "closure_test_bayes.root",
                     const char* outDir        = "closure_plots_1D",
                     double       ratioTolerance = 0.10,  // 10% by default
                     double       chi2PerNdfMax  = 3.0)
{
  cout << "draw_closure_1D: reading closure file: " << closureFile << endl;

  gSystem->mkdir(outDir, /*recursive=*/true);

  TFile f(closureFile, "READ");
  if (f.IsZombie()) {
    ::Error("draw_closure_1D", "Cannot open file '%s'", closureFile);
    return;
  }

  // Try both naming conventions (RunClosureTestBayes vs ResponseQA)
  TH1D* hTrue  = GetHist1D(f, "h_true_mc", "h_truth_trained");
  TH1D* hMeas  = GetHist1D(f, "h_meas_mc", "h_measured_trained");
  TH1D* hUnf   = dynamic_cast<TH1D*>(f.Get("h_unfold_closure"));
  TH1D* hRatio = dynamic_cast<TH1D*>(f.Get("h_closure_ratio"));
  TH1D* hDiff  = dynamic_cast<TH1D*>(f.Get("h_closure_diff"));

  if (!hTrue || !hUnf) {
    ::Error("draw_closure_1D",
            "Missing essential histograms (need at least truth + unfolded). "
            "hTrue=%p, hUnf=%p", (void*)hTrue, (void*)hUnf);
    f.Close();
    return;
  }

  // --------------------------------------------------------------------
  // 1) Draw spectra: truth, reco, unfolded
  // --------------------------------------------------------------------
  TCanvas cSpec("c_closure_spectra", "Closure spectra", 1000, 800);
  cSpec.Divide(1, 2);

  // Top pad: spectra
  cSpec.cd(1);
  gPad->SetGrid();

  // Style
  hTrue->SetLineColor(kBlack);
  hTrue->SetMarkerColor(kBlack);
  hTrue->SetMarkerStyle(21);

  if (hMeas) {
    hMeas->SetLineColor(kBlue+1);
    hMeas->SetMarkerColor(kBlue+1);
    hMeas->SetMarkerStyle(20);
  }

  hUnf->SetLineColor(kRed+1);
  hUnf->SetMarkerColor(kRed+1);
  hUnf->SetMarkerStyle(22);

  hTrue->SetTitle("MC closure: truth vs reco vs unfolded;global bin;entries");

  // Choose a reasonable maximum
  double ymax = hTrue->GetMaximum();
  if (hMeas) ymax = std::max(ymax, hMeas->GetMaximum());
  ymax = std::max(ymax, hUnf->GetMaximum());
  if (!(ymax > 0)) ymax = 1.0;
  hTrue->SetMaximum(1.2 * ymax);

  hTrue->Draw("E1");
  if (hMeas) hMeas->Draw("E1 SAME");
  hUnf->Draw("E1 SAME");

  TLegend leg(0.65, 0.70, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hTrue, "MC truth", "lep");
  if (hMeas) leg.AddEntry(hMeas, "MC reco (input)", "lep");
  leg.AddEntry(hUnf,  "Unfolded (closure)", "lep");
  leg.Draw();

  // Bottom pad: difference or ratio
  cSpec.cd(2);
  gPad->SetGrid();

  if (hDiff) {
    hDiff->SetLineColor(kMagenta+2);
    hDiff->SetMarkerColor(kMagenta+2);
    hDiff->SetMarkerStyle(20);
    hDiff->SetTitle("Closure difference: unfolded - truth;global bin;unfolded - truth");
    hDiff->Draw("E1");
  } else if (hRatio) {
    hRatio->SetLineColor(kBlue+1);
    hRatio->SetMarkerColor(kBlue+1);
    hRatio->SetMarkerStyle(20);
    hRatio->SetTitle("Closure ratio: unfolded / truth;global bin;ratio");
    hRatio->Draw("E1");

    double xmin = hRatio->GetXaxis()->GetXmin();
    double xmax = hRatio->GetXaxis()->GetXmax();

    // line at ratio = 1
    TLine line1(xmin, 1.0, xmax, 1.0);
    line1.SetLineStyle(2);
    line1.SetLineColor(kGray+2);
    line1.Draw();
  }

  cSpec.SaveAs(TString::Format("%s/closure_spectra.png", outDir));
  cSpec.SaveAs(TString::Format("%s/closure_spectra.pdf", outDir));

  // --------------------------------------------------------------------
  // 2) Dedicated ratio canvas + numeric metrics
  // --------------------------------------------------------------------
  int    nUsedRatioBins = 0;
  int    nBadBins       = 0;
  double maxAbsDev      = 0.0;

  if (hRatio) {
    TCanvas cRatio("c_closure_ratio", "Closure ratio", 1000, 600);
    gPad->SetGrid();

    hRatio->SetLineColor(kBlue+1);
    hRatio->SetMarkerColor(kBlue+1);
    hRatio->SetMarkerStyle(20);
    hRatio->SetTitle("Closure ratio: unfolded / truth;global bin;unfolded / truth");
    hRatio->Draw("E1");

    const int nb = hRatio->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
      double r = hRatio->GetBinContent(i);
      double e = hRatio->GetBinError(i);
      if (r == 0.0 && e == 0.0) continue; // skip completely empty
      ++nUsedRatioBins;

      double dev = std::fabs(r - 1.0);
      maxAbsDev = std::max(maxAbsDev, dev);
      if (dev > ratioTolerance) ++nBadBins;
    }

    double xmin = hRatio->GetXaxis()->GetXmin();
    double xmax = hRatio->GetXaxis()->GetXmax();

    // reference lines
    TLine line1(xmin, 1.0, xmax, 1.0);
    line1.SetLineStyle(2);
    line1.SetLineColor(kGray+2);
    line1.Draw();

    TLine lineUp(xmin, 1.0 + ratioTolerance, xmax, 1.0 + ratioTolerance);
    lineUp.SetLineStyle(3);
    lineUp.SetLineColor(kRed+1);
    lineUp.Draw();

    TLine lineDn(xmin, 1.0 - ratioTolerance, xmax, 1.0 - ratioTolerance);
    lineDn.SetLineStyle(3);
    lineDn.SetLineColor(kRed+1);
    lineDn.Draw();

    // Annotate
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.88, Form("Max |ratio-1| = %.3g", maxAbsDev));
    lat.DrawLatex(0.15, 0.84, Form("Bins outside #pm %.3g: %d / %d",
                                   ratioTolerance, nBadBins, nUsedRatioBins));

    cRatio.SaveAs(TString::Format("%s/closure_ratio.png", outDir));
    cRatio.SaveAs(TString::Format("%s/closure_ratio.pdf", outDir));
  }

  // --------------------------------------------------------------------
  // 3) Compute chi2/ndf between unfolded and truth
  // --------------------------------------------------------------------
  double chi2 = 0.0;
  int    ndf  = 0;

  if (hTrue && hUnf) {
    const int nb = hTrue->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
      double t  = hTrue->GetBinContent(i);
      double u  = hUnf->GetBinContent(i);
      double eT = hTrue->GetBinError(i);
      double eU = hUnf->GetBinError(i);
      double e2 = eT*eT + eU*eU;
      if (e2 <= 0.0) continue;

      double d = u - t;
      chi2 += (d*d) / e2;
      ++ndf;
    }
  }
  double chi2PerNdf = (ndf > 0 ? chi2 / ndf : 0.0);

  // --------------------------------------------------------------------
  // 4) Decide pass / fail and print summary
  // --------------------------------------------------------------------
  bool passChi2   = (chi2PerNdfMax <= 0.0) ? true : (chi2PerNdf <= chi2PerNdfMax);
  bool passRatio  = (hRatio ? (nBadBins == 0) : true);
  bool closureOK  = passChi2 && passRatio;

  cout << "\n============================================================\n";
  cout << "Closure QC summary for file: " << closureFile << "\n";
  cout << "  chi2/ndf = " << chi2 << " / " << ndf
       << " = " << chi2PerNdf << "  (limit: " << chi2PerNdfMax << ")\n";

  if (hRatio) {
    cout << "  ratio tolerance = +/- " << ratioTolerance << "\n";
    cout << "  bins outside tolerance = " << nBadBins
         << " / " << nUsedRatioBins
         << "  (max |ratio-1| = " << maxAbsDev << ")\n";
  } else {
    cout << "  ratio histogram 'h_closure_ratio' not found; "
            "skipping per-bin ratio check.\n";
  }

  cout << "  ==> Closure " << (closureOK ? "SUCCESS" : "FAILURE") << "\n";
  cout << "============================================================\n\n";

  f.Close();
}
