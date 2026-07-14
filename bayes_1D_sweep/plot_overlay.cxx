// overlay_phi_unfolded_iters_1D.C
//
// Overlay unfolded phi spectra for 4 Bayes iterations,
// using the same 1D geometry and mapping as unfold_onepass_1D.cxx
// and fit_phi_unfolded_1D.C.
//
// Default expects files:
//   unfold_bayes_iter01.root ... unfold_bayes_iter04.root
// with histograms:
//   unfold_Bayes_iter01 ... unfold_Bayes_iter04
//
// Usage example (from ROOT):
//   root -l -b -q 'overlay_phi_unfolded_iters_1D.C()'
// or with explicit args:
//   root -l -b -q \
//   'overlay_phi_unfolded_iters_1D.C("unfold_bayes_iter%02d.root","unfold_Bayes_iter%02d",4,"phi_overlay_iters_1D",0.0,16)'

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
#include "TStyle.h"

#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <algorithm>

// Your binning
#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins

// ----------------------------------------------------------------------
//  Geometry helpers (must match unfold_onepass_1D.cxx)
// ----------------------------------------------------------------------
namespace RESP1D_OL {
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
static inline int GlobalBinIndex1D_OL(int z_bin, int x_bin) {
  return x_bin * RESP1D_OL::nZ + z_bin;   // z is the fast index
}

// Simple binning check
static bool axes_identical_1D(const TH1* a, const TH1* b) {
  if (!a || !b) return false;
  const TAxis* A = a->GetXaxis();
  const TAxis* B = b->GetXaxis();
  return (A->GetNbins() == B->GetNbins()) &&
         (std::fabs(A->GetXmin() - B->GetXmin()) < 1e-9) &&
         (std::fabs(A->GetXmax() - B->GetXmax()) < 1e-9);
}

// ----------------------------------------------------------------------
// Build phi histogram (for given ix,z,pT) from a 1D unfolded spectrum
// ----------------------------------------------------------------------
static TH1D* MakePhiSliceFromGlobal(const TH1D* hGlobal,
                                    int ix_bin0, // 0-based xq2bin
                                    int iz,      // 1..N_Zbins
                                    int ipt,     // 1..N_pTbins_with_overflow
                                    const char* name,
                                    double emptyEps = 0.0,
                                    int nZ = N_Zbins,
                                    int nPt = N_pTbins_with_overflow,
                                    int nPhi = N_phiTrbins)
{
  if (!hGlobal) return nullptr;

  const int nTot_have = hGlobal->GetNbinsX();
  const int zpt2 = (iz-1)*nPt + ipt;   // 1..(nZ*nPt)
  const int j0   = zpt2 - 1;           // 0-based over (z,pT)

  TH1D* hphi = new TH1D(
    name,
    Form("#phi_{Trento};#phi_{Trento} [deg];Unfolded counts (xQ^{2}=%d, z=%d, pT=%d)",
         ix_bin0+1, iz, ipt),
    nPhi, 0.0, 360.0);
  hphi->Sumw2();

  bool nonEmpty = false;

  for (int iph=1; iph<=nPhi; ++iph) {
    // Composite index 1..(nZ*nPt*nPhi); 0 is reserved
    const int comp = j0*nPhi + iph;    // 1..(nZ*nPt*nPhi)
    const int zbin = comp;             // this is z_pt2_phi_bin

    // global_bin = xq2bin * nZ_global + zbin
    const int g    = GlobalBinIndex1D_OL(zbin, ix_bin0);
    const int ibin = g + 1;            // ROOT bin index (1..nTot)

    if (ibin < 1 || ibin > nTot_have) continue;

    const double v = hGlobal->GetBinContent(ibin);
    const double e = hGlobal->GetBinError  (ibin);

    hphi->SetBinContent(iph, v);
    hphi->SetBinError  (iph, e);

    if (std::fabs(v) > emptyEps || std::fabs(e) > emptyEps) nonEmpty = true;
  }

  if (!nonEmpty) {
    delete hphi;
    return nullptr;
  }

  return hphi;
}

// ----------------------------------------------------------------------
// Main overlay driver
// ----------------------------------------------------------------------
void overlay_phi_unfolded_iters_1D(const char* filePattern = "unfold_bayes_iter%02d.root",
                                   const char* histPattern = "unfold_Bayes_iter%02d",
                                   int         nIter       = 4,
                                   const char* outDir      = "phi_overlay_iters_1D",
                                   double      emptyEps    = 0.0,
                                   int         max_ix      = 16,
                                   int         nZ          = N_Zbins,
                                   int         nPt         = N_pTbins_with_overflow,
                                   int         nPhi        = N_phiTrbins)
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, /*recursive*/true);

  const int nX_limit = std::min(max_ix, RESP1D_OL::nX);

  // ------------------------------------------------------------------
  // Open all iteration files and load their unfolded 1D histograms
  // ------------------------------------------------------------------
  std::vector<TFile*> files(nIter, nullptr);
  std::vector<TH1D*>  globals(nIter, nullptr);

  for (int iIter = 0; iIter < nIter; ++iIter) {
    TString fname = Form(filePattern, iIter+1);
    files[iIter] = TFile::Open(fname, "READ");
    if (!files[iIter] || files[iIter]->IsZombie()) {
      std::cerr << "overlay_phi_unfolded_iters_1D: WARNING: cannot open " << fname << "\n";
      if (files[iIter]) { delete files[iIter]; files[iIter] = nullptr; }
      continue;
    }
    TString hname = Form(histPattern, iIter+1);
    globals[iIter] = dynamic_cast<TH1D*>(files[iIter]->Get(hname));
    if (!globals[iIter]) {
      std::cerr << "overlay_phi_unfolded_iters_1D: WARNING: cannot find " << hname
                << " in " << fname << "\n";
    }
  }

  bool anyGlobal = false;
  for (auto* h : globals) {
    if (h) { anyGlobal = true; break; }
  }
  if (!anyGlobal) {
    std::cerr << "overlay_phi_unfolded_iters_1D: ERROR: no unfolded histograms loaded.\n";
    return;
  }

  // Check binning consistency
  for (int i = 1; i < nIter; ++i) {
    if (globals[0] && globals[i] &&
        !axes_identical_1D(globals[0], globals[i])) {
      std::cerr << "overlay_phi_unfolded_iters_1D: WARNING: binning mismatch between iter 1 and iter "
                << (i+1) << "\n";
    }
  }

  std::vector<int> colors = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kCyan+1};

  long long nCand = 0;
  long long nPlotted = 0;
  long long nSkippedEmpty = 0;

  // ------------------------------------------------------------------
  // Loop over xQ², z, pT; build phi-hists for each iteration; overlay
  // ------------------------------------------------------------------
  for (int ix_bin0 = 0; ix_bin0 < nX_limit; ++ix_bin0) {
    int ix_dir = ix_bin0 + 1; // 1-based for labels

    for (int iz_idx = 1; iz_idx <= nZ; ++iz_idx) {
      for (int ipt = 1; ipt <= nPt; ++ipt) {
        ++nCand;

        std::vector<TH1D*> hs(nIter, nullptr);
        bool anyNonNull = false;

        for (int it = 0; it < nIter; ++it) {
          if (!globals[it]) continue;
          TString hnameSlice = Form("phi_iter%02d_ix%02d_z%02d_pt%02d",
                                    it+1, ix_dir, iz_idx, ipt);
          hs[it] = MakePhiSliceFromGlobal(globals[it],
                                          ix_bin0, iz_idx, ipt,
                                          hnameSlice.Data(),
                                          emptyEps, nZ, nPt, nPhi);
          if (hs[it]) anyNonNull = true;
        }

        if (!anyNonNull) {
          ++nSkippedEmpty;
          for (auto* h : hs) delete h;
          continue;
        }

        // Compute common y-max
        double ymax = 0.0;
        for (auto* h : hs) {
          if (!h) continue;
          ymax = std::max(ymax, h->GetMaximum());
        }
        if (!(ymax > 0.0)) ymax = 1.0;

        TString cname = Form("c_phi_ix%02d_z%02d_pt%02d", ix_dir, iz_idx, ipt);
        TString ctitle = Form("Unfolded #phi overlay: ix=%d, z=%d, pT=%d",
                              ix_dir, iz_idx, ipt);
        TCanvas c(cname, ctitle, 900, 700);

        TLegend leg(0.15, 0.65, 0.50, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.035);

        bool first = true;
        for (int it = 0; it < nIter; ++it) {
          TH1D* h = hs[it];
          if (!h) continue;

          int col = colors[it % (int)colors.size()];
          h->SetStats(0);
          h->SetLineColor(col);
          h->SetMarkerColor(col);
          h->SetMarkerStyle(20 + it);
          h->SetLineWidth(2);
          h->SetTitle(Form("Unfolded #phi;#phi_{Trento} [deg];Counts (ix=%d, z=%d, pT=%d)",
                           ix_dir, iz_idx, ipt));
          h->SetMaximum(1.2 * ymax);
          if (h->GetMinimum() < 0.0) h->SetMinimum(0.0);

          if (first) {
            h->Draw("E1");
            first = false;
          } else {
            h->Draw("E1 SAME");
          }

          leg.AddEntry(h, Form("Bayes iter %d", it+1), "lep");
        }

        leg.Draw();

        TString outName =
          Form("%s/phi_overlay_ix%02d_z%02d_pt%02d.png",
               outDir, ix_dir, iz_idx, ipt);
        c.SaveAs(outName);

        for (auto* h : hs) delete h;
        ++nPlotted;
      } // ipt
    } // iz
  } // ix

  std::cout << "overlay_phi_unfolded_iters_1D:\n"
            << "  candidates (ix*z*pT): " << nCand << "\n"
            << "  plotted (non-empty):  " << nPlotted << "\n"
            << "  skipped empty:        " << nSkippedEmpty << std::endl;

  // Clean up
  for (auto* f : files) {
    if (f) { f->Close(); delete f; }
  }
}
