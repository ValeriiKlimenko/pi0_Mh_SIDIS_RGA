// File: fit_phi_from_scaled.C
//
// Usage:
//   root -l -b -q 'fit_phi_from_scaled.C("phi_1d_scaled.root",
//                                        "phi_fit_results_scaled.root",
//                                        "phi_fit_plots_scaled",
//                                        0.0)'
//
// Expects histograms in inHistFile named:
//   h_phi_xq2_%02d_z%d_pt2%d
//
// These should be 1D TH1D with 8 φ bins between 0 and 360 degrees.
//
// For each xQ2 bin, creates a canvas with N_ZBINS rows and
// N_PT2_BINS_WITH_OVERFLOW columns, fits each φ spectrum with
//   f(φ) = p0 + p1*cos(φ) + p2*cos(2φ)
// and draws data+fit.
//
// Legend per pad:
//   line 1: p0
//   line 2: average of non-zero bins.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPad.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using std::cout;
using std::endl;

// -------- Binning: must match your analysis / Python script --------
static const int   N_ZBINS                  = 8;
static const int   N_PT2_BINS_WITH_OVERFLOW = 11; // 10 + overflow
static const int   N_PHI                    = 8;
static const int   MAX_XQ2_BIN              = 20; // xq2 in 0..20

// Optional nice bin labels (for z, pt2) if you want them on the plots
static const double Z_EDGES[9] = {
  0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0
};
static const double PT2_EDGES[12] = {
  0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50,
  0.65, 0.80, 1.00, 1.50
};

static double get_bin_width(int iz, int ipt)
{
  // iz: 1..N_ZBINS  (Z_EDGES has N_ZBINS+1 entries)
  // ipt: 1..N_PT2_BINS_WITH_OVERFLOW (PT2_EDGES has N_PT2_BINS_WITH_OVERFLOW+1 entries)
  if (iz < 1 || iz > N_ZBINS) return 0.0;
  if (ipt < 1 || ipt > N_PT2_BINS_WITH_OVERFLOW) return 0.0;

  double z_lo  = Z_EDGES[iz - 1];
  double z_hi  = Z_EDGES[iz];
  double pt_lo = PT2_EDGES[ipt - 1];
  double pt_hi = PT2_EDGES[ipt];

  return (z_hi - z_lo) * (pt_hi - pt_lo);
}

static TString bin_label(const double* edges, int nbins, int i1b) {
  if (i1b < 1 || i1b > nbins) return "out";
  double a = edges[i1b-1];
  double b = edges[i1b];
  return Form("%.2f–%.2f", a, b);
}

// --------------------------------------------------------------------
void fit_phi_from_scaled(const char* inHistFile  = "phi_1d_scaled.root",
                         const char* outRootFile = "phi_fit_results_scaled.root",
                         const char* plotDir     = "phi_fit_plots_scaled",
                         double       emptyEps   = 0.0)
{
  gSystem->mkdir(plotDir, kTRUE);

  TFile fIn(inHistFile, "READ");
  if (fIn.IsZombie()) {
    ::Error("fit_phi_from_scaled", "Cannot open input file %s", inHistFile);
    return;
  }

  TFile fOut(outRootFile, "RECREATE");
  if (fOut.IsZombie()) {
    ::Error("fit_phi_from_scaled", "Cannot create output file %s", outRootFile);
    return;
  }

  // Prepare CSV output: same basename as outRootFile but .csv
  TString csvName(outRootFile);
  csvName.ReplaceAll(".root", ".csv");
  std::ofstream csv(csvName.Data());
  if (!csv) {
    ::Error("fit_phi_from_scaled", "Cannot create CSV file %s", csvName.Data());
 } else {
    csv << "ixq2,iz,ipt,p0,mean_nonzero,no_limit_hit,bin_width\n";
  }

  // Tree with fit results
  Int_t   out_ixq2 = 0, out_iz = 0, out_ipt = 0;
  Double_t out_p0 = 0.0, out_p1 = 0.0, out_p2 = 0.0;
  Double_t out_mean_nonzero = 0.0;
  Int_t   out_no_limit_hit = 1; 

  TTree t("phi_fit_scaled", "Fits of phi slices from scaled histograms");
  t.Branch("ixq2", &out_ixq2,        "ixq2/I");
  t.Branch("iz",   &out_iz,          "iz/I");
  t.Branch("ipt",  &out_ipt,         "ipt/I");
  t.Branch("p0",   &out_p0,          "p0/D");
  t.Branch("p1",   &out_p1,          "p1/D");
  t.Branch("p2",   &out_p2,          "p2/D");
  t.Branch("mean_nonzero", &out_mean_nonzero, "mean_nonzero/D");
  t.Branch("no_limit_hit", &out_no_limit_hit, "no_limit_hit/I");


  long long nVisited = 0;
  long long nFitted  = 0;

  // Loop over xQ2 bins (using the center integer: 0..20)
  for (int ixq2 = 0; ixq2 <= MAX_XQ2_BIN; ++ixq2) {

    // Check if there is any histogram for this xQ2
    bool hasAny = false;
    for (int iz = 1; iz <= N_ZBINS && !hasAny; ++iz) {
      for (int ipt = 1; ipt <= N_PT2_BINS_WITH_OVERFLOW && !hasAny; ++ipt) {
        TString hname = Form("h_phi_xq2_%02d_z%d_pt2%d", ixq2, iz, ipt);
        if (fIn.Get(hname)) hasAny = true;
      }
    }
    if (!hasAny) continue;

    // Create canvas for this xQ2
    int width  = std::max(1400, N_PT2_BINS_WITH_OVERFLOW * 150);
    int height = std::max(900,  N_ZBINS * 150);
    TCanvas* c = new TCanvas(Form("c_scaled_xq2_%02d", ixq2),
                             Form("Scaled data + fit: xQ^{2} bin %d", ixq2),
                             width, height);
    c->Divide(N_PT2_BINS_WITH_OVERFLOW, N_ZBINS);

    // Loop over z and pT^2
    for (int iz = 1; iz <= N_ZBINS; ++iz) {
      for (int ipt = 1; ipt <= N_PT2_BINS_WITH_OVERFLOW; ++ipt) {
        ++nVisited;

        const int padnum = (iz - 1) * N_PT2_BINS_WITH_OVERFLOW + ipt;
        c->cd(padnum);
        TPad* pad = (TPad*)gPad;
        pad->SetFillColor(0);
        pad->SetFrameFillColor(0);
        pad->SetTicks(1,1);
        pad->SetLeftMargin(0.18);
        pad->SetBottomMargin((iz == N_ZBINS) ? 0.22 : 0.08);
        pad->SetRightMargin(0.05);
        pad->SetTopMargin(0.08);
        pad->SetBit(TPad::kClipFrame, kTRUE);

        TString hname = Form("h_phi_xq2_%02d_z%d_pt2%d", ixq2, iz, ipt);
        TH1D* h = dynamic_cast<TH1D*>(fIn.Get(hname));
        if (!h) {
          // Nothing for this cell – just draw axes so layout is consistent
          TH1D dummy("dummy","", N_PHI, 0.0, 360.0);
          dummy.SetStats(0);
          dummy.GetXaxis()->SetTitle("#phi_{Trento} [deg]");
          dummy.GetYaxis()->SetTitle("");
          dummy.GetXaxis()->SetNdivisions(505);
          dummy.GetYaxis()->SetNdivisions(505);

          if (iz == N_ZBINS) {
            dummy.GetXaxis()->SetLabelSize(0.08);
            dummy.GetXaxis()->SetTitleSize(0.08);
          } else {
            dummy.GetXaxis()->SetLabelSize(0.0);
            dummy.GetXaxis()->SetTitleSize(0.0);
          }
          dummy.GetYaxis()->SetLabelSize(0.08);
          dummy.GetYaxis()->SetTitleSize(0.01);
          dummy.Draw("AXIS");
          continue;
        }

        // Build vectors for TGraphErrors
        const int nBins = h->GetNbinsX();
        std::vector<double> vx, vy, vex, vey;
        vx.reserve(nBins);
        vy.reserve(nBins);
        vex.reserve(nBins);
        vey.reserve(nBins);

        // Track which phi bins are non-empty in the original histogram
        std::vector<bool> phiOccupied(nBins, false);

        // Compute mean over non-zero bins (same as before)
        double sumNonZero = 0.0;
        int    nNonZero   = 0;

        for (int b = 1; b <= nBins; ++b) {
          double y = h->GetBinContent(b);
          double e = h->GetBinError(b);

          if (std::fabs(y) <= emptyEps && std::fabs(e) <= emptyEps) continue;

          double x = h->GetXaxis()->GetBinCenter(b);
          vx.push_back(x);
          vy.push_back(y);
          vex.push_back(0.0);
          vey.push_back(e);

          sumNonZero += y;
          ++nNonZero;

          // Mark this phi-bin as occupied (b runs 1..nBins)
          phiOccupied[b - 1] = true;
        }

        if (nNonZero == 0) {
          // No non-zero bins – draw axes only
          TH1D dummy("dummy","", N_PHI, 0.0, 360.0);
          dummy.SetStats(0);
          dummy.GetXaxis()->SetTitle("#phi_{Trento} [deg]");
          dummy.GetYaxis()->SetTitle("");
          dummy.GetXaxis()->SetNdivisions(505);
          dummy.GetYaxis()->SetNdivisions(505);

          if (iz == N_ZBINS) {
            dummy.GetXaxis()->SetLabelSize(0.08);
            dummy.GetXaxis()->SetTitleSize(0.08);
          } else {
            dummy.GetXaxis()->SetLabelSize(0.0);
            dummy.GetXaxis()->SetTitleSize(0.0);
          }
          dummy.GetYaxis()->SetLabelSize(0.08);
          dummy.GetYaxis()->SetTitleSize(0.01);
          dummy.Draw("AXIS");
          continue;
        }

        double mean_nonzero = sumNonZero / nNonZero;


                // ------------------------------------------------------------
        // Check phi coverage:
        //   - allow first and last bin (b=1 and b=nBins) to be empty
        //   - require at least 3 interior bins (2..nBins-1) to be non-empty
        //   - and avoid large gaps: no run of >= 2 empty interior bins
        // ------------------------------------------------------------
        bool goodPhiCoverage = false;
        if (nBins >= 3) {
          int nInteriorOccupied = 0;
          int maxEmptyRun       = 0;
          int currentEmptyRun   = 0;

          // interior bins: ROOT indices 2..(nBins-1)
          for (int b = 2; b <= nBins - 1; ++b) {
            bool occ = phiOccupied[b - 1];
            if (occ) {
              ++nInteriorOccupied;
              currentEmptyRun = 0;
            } else {
              ++currentEmptyRun;
              if (currentEmptyRun > maxEmptyRun)
                maxEmptyRun = currentEmptyRun;
            }
          }

          // tweak these two numbers as you like:
          //   - require at least 3 interior bins with data
          //   - and no gap of 2 or more consecutive empty interior bins
          if (nInteriorOccupied >= 3 && maxEmptyRun <= 1) {
            goodPhiCoverage = true;
          }
        }


        // Make a TGraphErrors from these points (heap allocation so it survives)
        const int npts = (int)vx.size();
        TGraphErrors* gr = new TGraphErrors(npts,
                                            vx.data(),  vy.data(),
                                            vex.data(), vey.data());
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.7);
        gr->SetLineWidth(1);
        gr->SetLineColor(kBlue+1);
        gr->SetMarkerColor(kBlue+1);

        // Per-pad Y max for frame
        double ymax = 0.0;
        for (int i = 0; i < npts; ++i) {
          if (vy[i] > ymax) ymax = vy[i];
        }
        if (!(ymax > 0.0)) ymax = 1.0;
        ymax *= 1.2;

        // Draw initial frame (like in the Python script)
        TH1* frame = pad->DrawFrame(0.0, 0.0, 360.0, ymax);
        frame->GetXaxis()->SetTitle("#phi_{Trento} [deg]");
        frame->GetYaxis()->SetTitle("");  // no Y-label text
        frame->GetXaxis()->SetNdivisions(505);
        frame->GetYaxis()->SetNdivisions(505);

        if (iz == N_ZBINS) {
          frame->GetXaxis()->SetLabelSize(0.08);
          frame->GetXaxis()->SetTitleSize(0.08);
        } else {
          frame->GetXaxis()->SetLabelSize(0.0);
          frame->GetXaxis()->SetTitleSize(0.0);
        }
        frame->GetYaxis()->SetLabelSize(0.08);
        frame->GetYaxis()->SetTitleSize(0.01);

        // Fit if we have at least 3 points
        bool canFit = (npts >= 4 && goodPhiCoverage);
        double p0 = 0.0, p1 = 0.0, p2 = 0.0;
        TF1* fphi = nullptr;  // per-pad fit function

        if (canFit) {
          // Create a unique TF1 for this pad
          TString fname = Form("fphi_ixq2_%02d_z%d_pt2%d", ixq2, iz, ipt);
          fphi = new TF1(fname,
                         "[0] + [1]*cos(x*TMath::Pi()/180.) + [2]*cos(2.*x*TMath::Pi()/180.)",
                         0.0, 360.0);

          // Initial parameters
          fphi->SetParameters(mean_nonzero, 0.0, 0.0);
          fphi->SetRange(0.0, 360.0); // ensure fit function range is 0–360

          double scale = std::max(mean_nonzero, 1.0);
          fphi->SetParLimits(0, 0.0,          5.0 * scale);
          fphi->SetParLimits(1, -0.5 * scale, 0.5 * scale);
          fphi->SetParLimits(2, -0.3 * scale, 0.3 * scale);

          // "QSN": Quiet, Store result, No automatic draw
          gr->Fit(fphi, "QSN");

          p0 = fphi->GetParameter(0);
          p1 = fphi->GetParameter(1);
          p2 = fphi->GetParameter(2);
          

        // Check if any parameter hit its limits (0,1,2)
        bool hitLimit = false;
        const double relTol = 1e-3; // 0.1% of range
        
        for (int ipar = 1; ipar < 3; ++ipar) {
          double lo = 0.0, hi = 0.0;
          fphi->GetParLimits(ipar, lo, hi);
        
          // If no limits set, ROOT keeps lo=0, hi=0 (default); skip those
          if (lo == 0.0 && hi == 0.0) continue;
        
          double val   = fphi->GetParameter(ipar);
          double range = hi - lo;
          if (range <= 0.0) continue;
        
          double distLo = std::fabs(val - lo) / range;
          double distHi = std::fabs(val - hi) / range;
        
          if (distLo < relTol || distHi < relTol) {
            hitLimit = true;
            break;
          }
        }

// 1 if no limit was hit, 0 if any limit was hit
out_no_limit_hit = hitLimit ? 0 : 1;

          

          // Re-adjust y max to include fit over 0–360
          double fmax = fphi->GetMaximum(0.0, 360.0);
          double ymax2 = std::max(ymax, 1.2 * fmax);
          if (ymax2 <= 0.0) ymax2 = ymax;
          pad->Clear();
          frame = pad->DrawFrame(0.0, 0.0, 360.0, ymax2);
          frame->GetXaxis()->SetTitle("#phi_{Trento} [deg]");
          frame->GetYaxis()->SetTitle("");
          frame->GetXaxis()->SetNdivisions(505);
          frame->GetYaxis()->SetNdivisions(505);

          if (iz == N_ZBINS) {
            frame->GetXaxis()->SetLabelSize(0.08);
            frame->GetXaxis()->SetTitleSize(0.08);
          } else {
            frame->GetXaxis()->SetLabelSize(0.0);
            frame->GetXaxis()->SetTitleSize(0.0);
          }
          frame->GetYaxis()->SetLabelSize(0.08);
          frame->GetYaxis()->SetTitleSize(0.01);
        }

        // Draw the underlying histogram (so you clearly see the "histogram")
        h->SetStats(0);
        h->SetLineColor(kBlue+1);
        h->SetMarkerColor(kBlue+1);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.7);
        h->Draw("E1 SAME");

        // Draw points
        gr->Draw("P SAME");

        // Draw fit function explicitly over 0–360 if we have a fit
        if (fphi) {
          fphi->SetRange(0.0, 360.0);
          fphi->SetLineWidth(2);
          fphi->Draw("SAME");
        }

        // Top-row pt2 label
        if (iz == 1) {
          TLatex t;
          t.SetNDC();
          t.SetTextSize(0.05);
          t.DrawLatex(0.18, 0.92,
                      Form("p_{T}^{2} %s",
                           bin_label(PT2_EDGES, 11, ipt).Data()));
        }

        // First-column z label
        if (ipt == 1) {
          TLatex t;
          t.SetNDC();
          t.SetTextSize(0.05);
          t.DrawLatex(0.20, 0.85,
                      Form("z %s",
                           bin_label(Z_EDGES, 8, iz).Data()));
        }

        // Legend: only p0 and mean_nonzero
        TLegend* leg = new TLegend(0.15, 0.70, 0.95, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        if (fphi) {
          leg->AddEntry((TObject*)0,
                        Form("p_{0} = %.4g", p0),
                        "");
        } else {
          leg->AddEntry((TObject*)0,
                        "p_{0}: fit failed / few pts",
                        "");
        }
        leg->AddEntry((TObject*)0,
                      Form("<y>_{non-zero} = %.4g", mean_nonzero),
                      "");
        leg->Draw();
        
        // Store in TTree and CSV if fitted
        if (fphi) {
          out_ixq2         = ixq2;
          out_iz           = iz;
          out_ipt          = ipt;
          out_p0           = p0;
          out_p1           = p1;
          out_p2           = p2;
          out_mean_nonzero = mean_nonzero;
          // out_no_limit_hit was already set just after the fit
          t.Fill();
          ++nFitted;
        
          if (csv) {
            double bin_width = get_bin_width(iz, ipt);

            csv << ixq2 << ","
                << iz   << ","
                << ipt  << ","
                << p0   << ","
                << mean_nonzero << ","
                << out_no_limit_hit << ","
                << bin_width << "\n";
          }

        }

      } // ipt
    } // iz

    // Global title
    c->cd(0);
    TLatex st;
    st.SetNDC();
    st.SetTextSize(0.045);
    st.DrawLatex(0.02, 0.98, Form("xQ^{2} bin (center): %d", ixq2));

    // Save canvas
    TString pngName = Form("%s/xq2_%02d.png", plotDir, ixq2);
    TString pdfName = Form("%s/xq2_%02d.pdf", plotDir, ixq2);
    c->SaveAs(pngName);
    c->SaveAs(pdfName);

    delete c;
  } // ixq2

  fOut.cd();
  t.Write();
  fOut.Close();
  fIn.Close();

  if (csv) csv.close();

  ::Info("fit_phi_from_scaled",
         "Visited %lld (xQ2×z×pt2), fitted %lld. "
         "Results: '%s', plots: '%s', CSV: '%s'.",
         nVisited, nFitted, outRootFile, plotDir, csvName.Data());
}
