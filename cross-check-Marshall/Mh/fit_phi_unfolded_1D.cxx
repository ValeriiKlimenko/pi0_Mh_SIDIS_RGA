// ------------------------------------------------------------------------------------------
// File: fit_phi_unfolded_1D.C
//
// Usage examples (ROOT):
//   root -l -b -q 'fit_phi_unfolded_1D.C("unfold_out.root","unfold_Bayes_iter5")'
//   root -l -b -q 'fit_phi_unfolded_1D.C("unfold_out.root","unfold_Bayes_iter5",
//                                        "phi_slices_1D.root","phi_fit_results_1D.root",
//                                        "phi_fit_plots_1D",0.0,16)'
//
// This is the 1D analogue of make_phi_slices.C + fit_phi_unfolded.C:
//   * Input:  1D unfolded TH1D with global bin index g = xq2bin * nZ + z_pt2_phi_bin
//   * Output: phi_slices_1D.root with TH1D(phi) for each (xq2, z, pT)
//             phi_fit_results_1D.root with p0,p1,p2 fits
//             phi_fit_plots_1D/*.png with the fitted curves
//
// Geometry MUST match unfold_onepass_1D.cxx:
//   RESP::nX  = 21
//   RESP::nZ  = N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1
//   global_bin = xq2bin * nZ + z_pt2_phi_bin
//   TH1D axis: nTot = nX*nZ bins, x from -0.5 .. nTot-0.5, bin i center = (i-1)
// ------------------------------------------------------------------------------------------



#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TParameter.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <algorithm>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins

// ----------------------------------------------------------------------
//  Geometry helpers (must match unfold_onepass_1D.cxx)
// ----------------------------------------------------------------------

namespace RESP1D {
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

// Map (z_bin, x_bin) -> global 1D bin index (same as unfold_onepass_1D.cxx)

static inline int GlobalBinIndex1D(int z_bin, int x_bin) {
  return x_bin * RESP1D::nZ + z_bin;   // z is the fast index
}

// ----------------------------------------------------------------------
// Split 1D unfolded TH1D into phi histograms per (xq2bin,z,pT)
// Layout of output file matches the old 2D version:
//
//   phi_slices_1D.root
//     /ix01/z01/phi_ix01_z01_pt01
//     /ix01/z01/phi_ix01_z01_pt02
//     ...
//     /ix02/z01/...
// ----------------------------------------------------------------------

static void SplitPhiHists1D(const TH1D* h1,
                            const char* outFile       = "phi_slices.root",
                            double      emptyEps      = 0.0,
                            int         max_ix        = 16, // cap xQ2 to 16 by default
                            int         nZ            = N_Zbins,
                            int         nPt           = N_pTbins_with_overflow,
                            int         nPhi          = N_phiTrbins)
{
  if (!h1) {
    cerr << "SplitPhiHists1D: null TH1D\n";
    return;
  }

  // Check geometry vs expectation
  const int nComposite = nZ * nPt * nPhi;      // number of z⊗pT⊗phi bins (no +1)
  const int nZglobal   = nComposite + 1;       // includes extra "+1" bin (same as RESP1D::nZ)
  if (nZglobal != RESP1D::nZ) {
    cerr << "WARNING: nZglobal(" << nZglobal
         << ") != RESP1D::nZ(" << RESP1D::nZ
         << "). Using RESP1D::nZ for indexing.\n";
  }

  const int nX       = RESP1D::nX;
  const int nTot_exp = RESP1D::nTot;
  const int nTot_have = h1->GetNbinsX();

  if (nTot_have < nTot_exp) {
    cerr << "WARNING: unfolded TH1D has only " << nTot_have
         << " bins, but geometry expects at least " << nTot_exp << ".\n";
  }

  const int nX_limit = std::min(max_ix, nX);   // we will use xq2bin = 0..(nX_limit-1)

  long long nCandidates  = 0; // total (xQ2,z,pT) combos
  long long nWritten     = 0; // non-empty φ hists saved
  long long nEmpty       = 0; // empty (not saved)
  long long nOutOfRange  = 0; // global bins outside TH1D
  long long nClamped     = 0; // (unused, for symmetry with 2D version)

  TFile fout(outFile, "RECREATE");
  if (fout.IsZombie()) {
    cerr << "SplitPhiHists1D: ERROR: cannot create " << outFile << "\n";
    return;
  }

  // Loop over xQ2 bins
  //  - directory naming is 1-based ("ix01", "ix02", ...) to match old code
  for (int ix_bin=0; ix_bin<nX_limit; ++ix_bin) {
    const int ix_dir = ix_bin + 1;    // for directory and hist titles

    TDirectory* d_ix = fout.mkdir(Form("ix%02d", ix_dir));
    if (!d_ix) {
      cerr << "SplitPhiHists1D: failed to mkdir ix%02d\n" << ix_dir;
      continue;
    }

    for (int iz=1; iz<=nZ; ++iz) {
      TDirectory* d_z = d_ix->mkdir(Form("z%02d", iz));
      if (!d_z) {
        cerr << "SplitPhiHists1D: failed to mkdir z%02d in ix%02d\n" << iz << " " << ix_dir << "\n";
        continue;
      }

      for (int ipt=1; ipt<=nPt; ++ipt) {
        ++nCandidates;

        const int zpt2 = (iz-1)*nPt + ipt;   // 1..(nZ*nPt)
        const int j0   = zpt2 - 1;           // 0-based index over (z,pT)

        std::unique_ptr<TH1D> hphi(new TH1D(
          Form("phi_ix%02d_z%02d_pt%02d", ix_dir, iz, ipt),
          Form("#phi_{Trento};#phi_{Trento} [deg];Unfolded counts (xQ^{2}=%d, z=%d, pT=%d)",
               ix_dir, iz, ipt),
          nPhi, 0.0, 360.0));
        hphi->Sumw2();

        bool nonEmpty = false;

        for (int iph=1; iph<=nPhi; ++iph) {
          // Composite index 1..(nComposite); 0 reserved
          const int comp = j0*nPhi + iph;   // 1..nComposite
          if (comp < 1 || comp > nComposite) break;

          const int zbin = comp;           // this is exactly z_pt2_phi_bin used in 1D code

          // global_bin = xq2bin * nZ + z_bin
          const int g = GlobalBinIndex1D(zbin, ix_bin);   // 0..(nTot-1)
          const int ibin = g + 1;                         // ROOT bin number (1..nTot)

          if (ibin < 1 || ibin > nTot_have) {
            ++nOutOfRange;
            break;
          }

          const double v = h1->GetBinContent(ibin);
          const double e = h1->GetBinError  (ibin);

          hphi->SetBinContent(iph, v);
          hphi->SetBinError  (iph, e);

          if (std::fabs(v) > emptyEps || std::fabs(e) > emptyEps) nonEmpty = true;
        } // iph

        if (nonEmpty) {
          d_z->cd();
          hphi->Write();
          ++nWritten;
        } else {
          ++nEmpty;
        }
      } // ipt
    } // iz
  } // ix_bin

  fout.Write();
  fout.Close();

  cout << "SplitPhiHists1D: wrote phi slices to: " << outFile << "\n"
       << "  candidates (xQ2×z×pT): " << nCandidates << "\n"
       << "  written (non-empty):    " << nWritten << "\n"
       << "  empty (not saved):      " << nEmpty << "\n"
       << "  out-of-range globals:   " << nOutOfRange << "\n"
       << "  clamped packs (unused): " << nClamped << "\n"
       << "  nX_limit / nX = " << nX_limit << " / " << nX << "\n"
       << "  nZglobal / nComposite+1 = "
       << nZglobal << " / " << (nComposite+1) << "\n"
       << "  nTot_have / nTot_exp = "
       << nTot_have << " / " << nTot_exp << endl;
}

// ----------------------------------------------------------------------
// Fit every phi histogram with
//   p0 + p1 cos(phi) + p2 cos(2phi) (phi in degrees)
// and write:
//   * ROOT file with TTree(ix,iz,iPt,p0,p1,p2)
//   * PNG plots with data + fit
//
// This is almost identical to your 2D version, just renamed to *_1D.
// ----------------------------------------------------------------------
void fit_phi_slices_1D(const char* inPhiFile   = "phi_slices_1D.root",
                       const char* outRootFile = "phi_fit_results_1D.root",
                       const char* plotDir     = "phi_fit_plots_1D",
                       double      emptyEps    = 0.0,
                       int         max_ix      = 16,              // cap to 16 by default
                       int         nZ          = N_Zbins,
                       int         nPt         = N_pTbins_with_overflow,
                       int         minPoints   = 6)
{
  // Cuts:
  const double relUncMax = 0.90;  // currently NOT enforced (commented below)
  const double outlierK  = 1.9;   // currently NOT enforced (commented below)

  gSystem->mkdir(plotDir, /*recursive*/true);

  TFile fIn(inPhiFile, "READ");
  if (fIn.IsZombie()) {
    ::Error("fit_phi_slices_1D","Cannot open %s", inPhiFile);
    return;
  }

  TFile fOut(outRootFile, "RECREATE");
  if (fOut.IsZombie()) {
    ::Error("fit_phi_slices_1D","Cannot create %s", outRootFile);
    return;
  }

  Int_t   out_ix=0, out_iz=0, out_iPt=0;
  Double_t out_p0=0.0, out_p1=0.0, out_p2=0.0;

  TTree t("phi_fit_p0",
          "Fit results: p0,p1,p2 from p0 + p1*cos + p2*cos2 (1D unfolding)");
  t.Branch("ix",  &out_ix);
  t.Branch("iz",  &out_iz);
  t.Branch("iPt", &out_iPt);
  t.Branch("p0",  &out_p0);
  t.Branch("p1",  &out_p1);
  t.Branch("p2",  &out_p2);

  TF1 fphi("fphi",
           "[0] + [1]*cos(x*TMath::Pi()/180.) + [2]*cos(2.*x*TMath::Pi()/180.)",
           0.0, 360.0);

  long long nVisited=0, nFittedSaved=0, nTooFewPts=0, nMissing=0;

  for (int ix=1; ix<=max_ix; ++ix) {
    TDirectory* dix = fIn.GetDirectory(Form("ix%02d", ix));
    if (!dix) continue;

    for (int iz=1; iz<=nZ; ++iz) {
      TDirectory* dz = dix->GetDirectory(Form("z%02d", iz));
      if (!dz) continue;

      for (int ipt=1; ipt<=nPt; ++ipt) {
        ++nVisited;

        const TString hname = Form("phi_ix%02d_z%02d_pt%02d", ix, iz, ipt);
        TH1D* h = dynamic_cast<TH1D*>(dz->Get(hname));
        if (!h) { ++nMissing; continue; }

        const int nBins = h->GetNbinsX();

        // Histogram average over non-empty bins
        double sumNonZero = 0.0;
        int    nNonZero   = 0;
        for (int b = 1; b <= nBins; ++b) {
          const double v = h->GetBinContent(b);
          const double e = h->GetBinError(b);
          if (std::fabs(v) <= emptyEps && std::fabs(e) <= emptyEps) continue;
          sumNonZero += v;
          ++nNonZero;
        }
        const double histAvg = (nNonZero > 0) ? sumNonZero / nNonZero : 0.0;

        // Simple Fourier-like initial guesses for p1, p2
        double sumC1 = 0.0, sumC2 = 0.0;
        double sumC1C1 = 0.0, sumC2C2 = 0.0;
        double sumW = 0.0;

        for (int b = 1; b <= nBins; ++b) {
          const double v = h->GetBinContent(b);
          const double e = h->GetBinError(b);
          if (std::fabs(v) <= emptyEps && std::fabs(e) <= emptyEps) continue;

          const double phi_deg = h->GetXaxis()->GetBinCenter(b);
          const double phi     = phi_deg * TMath::Pi() / 180.0;

          const double w = (e > emptyEps) ? 1.0/(e*e) : 1.0;

          const double c1 = std::cos(phi);
          const double c2 = std::cos(2*phi);
          const double y  = v - histAvg;   // subtract mean

          sumW      += w;
          sumC1     += w * y * c1;
          sumC2     += w * y * c2;
          sumC1C1   += w * c1 * c1;
          sumC2C2   += w * c2 * c2;
        }

        double p0_init = std::max(1e-12, std::max(0.0, histAvg));
        double p1_init = 0.0;
        double p2_init = 0.0;
        if (sumC1C1 > 0) p1_init = sumC1 / sumC1C1;
        if (sumC2C2 > 0) p2_init = sumC2 / sumC2C2;

        fphi.SetParameters(p0_init, p1_init, p2_init);

        double scale = std::max(histAvg, 1.0);
        fphi.SetParLimits(0, 0.0,       20.0*scale);
        fphi.SetParLimits(1, -10.0*scale, 10.0*scale);
        fphi.SetParLimits(2, -10.0*scale, 10.0*scale);

        // Build TGraphErrors of points (with cuts)
        std::vector<double> x, y, ex, ey;
        x.reserve(nBins); y.reserve(nBins);
        ex.reserve(nBins); ey.reserve(nBins);
        std::vector<char> keep(nBins, 0);

        for (int b=1; b<=nBins; ++b) {
          const double v = h->GetBinContent(b);
          const double e = h->GetBinError(b);
          if (std::fabs(v) <= emptyEps && std::fabs(e) <= emptyEps) continue;

          const double denom  = std::fabs(v);
          const double relErr = (denom > emptyEps)
                                  ? std::fabs(e) / denom
                                  : std::numeric_limits<double>::infinity();
          // if (relErr > relUncMax) continue;      // keep disabled for now
          // if (histAvg > emptyEps && v > outlierK * histAvg) continue;

          x.push_back(h->GetXaxis()->GetBinCenter(b));
          y.push_back(v);
          ex.push_back(0.0);
          ey.push_back(e);
          keep[b-1] = 1;
        }

        const int npts = static_cast<int>(x.size());
        if (npts < minPoints) { ++nTooFewPts; continue; }

        TGraphErrors gr(npts, x.data(), y.data(), ex.data(), ey.data());
        gr.SetTitle(hname);
        gr.GetXaxis()->SetTitle("#phi_{Trento} [deg]");
        gr.GetYaxis()->SetTitle("Unfolded counts");

        TFitResultPtr r = gr.Fit(&fphi, "QS");   // quiet fit

        // Draw & annotate
        TCanvas c(Form("c_%s", hname.Data()), hname, 900, 700);

        TH1D* hplot = (TH1D*)h->Clone(Form("%s_clean", hname.Data()));
        hplot->SetStats(0);
        for (int b=1; b<=nBins; ++b) {
          if (!keep[b-1]) {
            hplot->SetBinContent(b, 0.0);
            hplot->SetBinError  (b, 0.0);
          }
        }

        double fmax = fphi.GetMaximum(0.0, 360.0);
        double ymax = std::max(hplot->GetMaximum(), fmax);
        if (!(ymax > 0)) ymax = 1.0;
        hplot->SetMaximum(1.2 * ymax);
        if (hplot->GetMinimum() < 0) hplot->SetMinimum(0);

        hplot->Draw("E1");
        gr.SetMarkerStyle(20);
        gr.SetMarkerSize(1.0);
        gr.Draw("P SAME");
        fphi.SetLineWidth(2);
        fphi.Draw("SAME");

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.15, 0.89, Form("histAvg = %.4g", histAvg));
        lat.DrawLatex(0.15, 0.85, Form("npts = %d", npts));
        lat.DrawLatex(0.15, 0.81, Form("p0 = %.4g", fphi.GetParameter(0)));
        lat.DrawLatex(0.15, 0.77, Form("p1 = %.4g", fphi.GetParameter(1)));
        lat.DrawLatex(0.15, 0.73, Form("p2 = %.4g", fphi.GetParameter(2)));

        out_ix  = ix;
        out_iz  = iz;
        out_iPt = ipt;
        out_p0  = fphi.GetParameter(0);
        out_p1  = fphi.GetParameter(1);
        out_p2  = fphi.GetParameter(2);
        t.Fill();

        c.SaveAs(TString::Format("%s/%s.png", plotDir, hname.Data()));

        delete hplot;
        ++nFittedSaved;
      } // ipt
    } // iz
  } // ix

  fOut.cd();
  t.Write();
  fOut.Close();
  fIn.Close();

  ::Info("fit_phi_slices_1D",
         "Visited %lld; fitted+saved %lld; too-few-points %lld; missing %lld. "
         "Plots in '%s', results in '%s'.",
         nVisited, nFittedSaved, nTooFewPts, nMissing, plotDir, outRootFile);
}

// ----------------------------------------------------------------------
// Convenience driver:
//   1) read 1D unfolded TH1D from unfoldFile (hist name 'hname')
//   2) split into phi slices with SplitPhiHists1D
//   3) fit all phi slices with fit_phi_slices_1D
// ----------------------------------------------------------------------
void fit_phi_unfolded_1D(const char* unfoldFile    = "unfold_bayes_iter01.root",
                         const char* hname         = "unfold_Bayes_iter01",
                         const char* phiSlicesFile = "phi_slices_1D.root",
                         const char* fitOutFile    = "phi_fit_results_1D.root",
                         const char* plotDir       = "phi_fit_plots_1D",
                         double       emptyEps     = 0.0,
                         int          max_ix       = 16)  // cap xQ2 to 16 by default
{
  TFile fin(unfoldFile, "READ");
  if (fin.IsZombie()) {
    cerr << "fit_phi_unfolded_1D: ERROR: cannot open " << unfoldFile << "\n";
    return;
  }

  TH1D* h1 = dynamic_cast<TH1D*>(fin.Get(hname));
  if (!h1) {
    cerr << "fit_phi_unfolded_1D: ERROR: cannot find TH1D '" << hname
         << "' in file " << unfoldFile << "\n";
    fin.Close();
    return;
  }

  // OPTIONAL: simple geometry sanity check
  if (h1->GetNbinsX() != RESP1D::nTot ||
      std::fabs(h1->GetXaxis()->GetXmin() - RESP1D::g_lo) > 1e-3 ||
      std::fabs(h1->GetXaxis()->GetXmax() - RESP1D::g_hi) > 1e-3) {
    cerr << "fit_phi_unfolded_1D: WARNING: TH1D geometry does not exactly match "
            "RESP1D (nTot=" << RESP1D::nTot
         << ", range=[" << RESP1D::g_lo << "," << RESP1D::g_hi << "]).\n";
  }

  // 1) build phi slices from the 1D unfolded spectrum
  SplitPhiHists1D(h1, phiSlicesFile, emptyEps,
                  max_ix, N_Zbins, N_pTbins_with_overflow, N_phiTrbins);

  fin.Close();

  // 2) fit each phi histogram and save plots + p0,p1,p2 table
  fit_phi_slices_1D(phiSlicesFile, fitOutFile, plotDir,
                    emptyEps, max_ix, N_Zbins, N_pTbins_with_overflow);
}
