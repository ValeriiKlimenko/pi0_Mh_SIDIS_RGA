// File: make_phi_slices.C
// Compile/run examples:
//   root -l -b -q 'make_phi_slices.C("unfold_out.root","unfold_Bayes_iter1","phi_slices.root",0.0,16)'
//   root -l -b -q 'fit_phi_slices.C("phi_slices.root","phi_fit_results.root","phi_fit_plots",0.0,16)'
//   root -l -b -q 'fit_phi_unfolded.C("unfold_out.root","unfold_Bayes_iter1","phi_slices.root","phi_fit_results.root","phi_fit_plots",0.0,16)'

// takes unfolded th2 (x axis is z-pt2-phi bins, y axis is x-Q2 bins)
// splits it into (xq2,z,pt2) th1 with phi dependence.
// It is done in order to fit phi dependence on the next step

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TParameter.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "TLatex.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <algorithm>
#include <vector>

#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

using std::cerr;
using std::cout;
using std::endl;

#include "binning_params.cxx"

// ----------- Rebuild TH2D from sparse tree ---------------------
static TH2D* ReconstructFromSparseDir(TDirectory* d, const char* denseName="unfold_ManualBinByBin")
{
  if (!d) { cerr << "ReconstructFromSparseDir: null directory\n"; return nullptr; }

  auto* vx = dynamic_cast<TVectorD*>(d->Get("bins_xedges"));
  auto* vy = dynamic_cast<TVectorD*>(d->Get("bins_yedges"));
  auto* pnx = dynamic_cast<TParameter<int>*>(d->Get("bins_nx"));
  auto* pny = dynamic_cast<TParameter<int>*>(d->Get("bins_ny"));
  auto* tr  = dynamic_cast<TTree*>(d->Get("bins"));

  if (!vx || !vy || !pnx || !pny || !tr) {
    cerr << "ReconstructFromSparseDir: missing components (vx, vy, pnx, pny, tr)\n";
    return nullptr;
  }

  const int nx = pnx->GetVal();
  const int ny = pny->GetVal();

  // NOTE: X is (z⊗pT⊗φ), Y is (xQ^2)
  TH2D* h = new TH2D(denseName,
                     "Unfolded;ix (z#otimesp_{T}#otimes#phi index);iy (xQ^{2} bin index)",
                     nx, vx->GetMatrixArray(),
                     ny, vy->GetMatrixArray());
  h->Sumw2();

  Int_t   ix=0, iy=0;
  Float_t val=0.f, err=0.f;
  tr->SetBranchAddress("ix",  &ix);
  tr->SetBranchAddress("iy",  &iy);
  tr->SetBranchAddress("val", &val);
  tr->SetBranchAddress("err", &err);

  const Long64_t n = tr->GetEntries();
  for (Long64_t i=0; i<n; ++i) {
    tr->GetEntry(i);
    if (ix>=1 && ix<=nx && iy>=1 && iy<=ny) {
      h->SetBinContent(ix, iy, val);
      h->SetBinError  (ix, iy, err);
    }
  }
  return h;
}


// ----------- Split into nPhi-bin phi histograms -------------------
static void SplitPhiHists(const TH2D* h2,
                          const char* outFile       = "phi_slices.root",
                          double      emptyEps      = 0.0,
                          int         max_ix        = 16, // cap xQ2 to 16 by default
                          int         nZ            = N_Zbins,
                          int         nPt           = N_pTbins_with_overflow,
                          int         nPhi          = N_phiTrbins)
{
  if (!h2) { cerr << "SplitPhiHists: null TH2\n"; return; }

  // EXPECTED: X packs = (z ⊗ pT ⊗ φ), Y = xQ2
  const int nComposite = nZ * nPt * nPhi;
  const int nx_have    = h2->GetNbinsX();
  const int ny_have   = h2->GetNbinsY();

  const int nx_expected = nComposite + 1;
  if (nx_have < nx_expected) {
    cerr << "WARNING: h2 has only " << nx_have
         << " X bins, but mapping expects at least " << nx_expected << ".\n";
  }
  const int nx_limit = std::min(nx_have, nx_expected);
  const int ny_limit = std::min(ny_have, max_ix); // cap Y by available xQ2 bins and user limit (default 16)

  // Counters
  long long nCandidates  = 0; // total (xQ2,z,pT) combos visited
  long long nWritten     = 0; // non-empty φ hists saved
  long long nEmpty       = 0; // empty (not saved)
  long long nOutOfRange  = 0; // base ix (in X) beyond nx_limit
  long long nClamped     = 0; // some φ bins truncated by nx_limit

  TFile fout(outFile, "RECREATE");
  if (fout.IsZombie()) { cerr << "ERROR: cannot create " << outFile << "\n"; return; }

  // Y loops over xQ2 (call it ixq2 for clarity but keep "ix" in names for consistency)
  for (int ixq2=1; ixq2<=ny_limit; ++ixq2) {
    TDirectory* d_ix = fout.mkdir(Form("ix%02d", ixq2)); // "ix" = xQ2 bin index in naming

    for (int iz=1; iz<=nZ; ++iz) {
      TDirectory* d_z = d_ix->mkdir(Form("z%02d", iz));

      for (int ipt=1; ipt<=nPt; ++ipt) {
        ++nCandidates;

      const int zpt2 = (iz-1)*nPt + ipt;   // 1..(nZ*nPt)
      const int j0   = zpt2 - 1;           // 0-based

    std::unique_ptr<TH1D> hphi(new TH1D(
      Form("phi_ix%02d_z%02d_pt%02d", ixq2, iz, ipt),
            Form("#phi_{Trento};#phi_{Trento} [deg];Unfolded counts (xQ^{2}=%d, z=%d, pT=%d)",
                 ixq2, iz, ipt),
            nPhi, 0.0, 360.0));
        hphi->Sumw2();

        bool nonEmpty = false;
        for (int iph=1; iph<=nPhi; ++iph) {
          const int comp = j0*nPhi + iph;   // 1..704
          if (comp < 1 || comp > nComposite) break;
    
          // coordinate = comp (what you used when filling), bin index computed by ROOT
          const double xcoord = double(comp);
          const int xbin = h2->GetXaxis()->FindBin(xcoord);
          if (xbin < 1 || xbin > nx_limit) break;
    
          const double v = h2->GetBinContent(xbin, ixq2);
          const double e = h2->GetBinError  (xbin, ixq2);
    
          hphi->SetBinContent(iph, v);
          hphi->SetBinError  (iph, e);
    
          if (std::fabs(v) > emptyEps || std::fabs(e) > emptyEps) nonEmpty = true;
        }

        if (nonEmpty) {
          d_z->cd();
          hphi->Write();
          ++nWritten;
        } else {
          ++nEmpty;
        }
      } // ipt
    } // iz
  } // ixq2

  fout.Write();
  fout.Close();

  cout << "Wrote phi slices to: " << outFile << "\n"
       << "SplitPhiHists summary:\n"
       << "  candidates (xQ2×z×pT): " << nCandidates << "\n"
       << "  written (non-empty):    " << nWritten << "\n"
       << "  empty (not saved):      " << nEmpty << "\n"
       << "  out-of-range packs:     " << nOutOfRange << "\n"
       << "  clamped packs:          " << nClamped << "\n"
      << "  nx_have/nx_expected/nx_limit = "
           << nx_have << "/" << nx_expected << "/" << nx_limit << "\n"
       << "  ny_have/ny_limit (xQ2) = "
       << ny_have << "/" << ny_limit << std::endl;
}



// Fit every phi histogram with p0 + p1*cos(x) + p2*cos(2x)  (x in degrees)
// and save: (1) per-hist PNGs in plotDir, (2) a ROOT file with a TTree(ix,iz,iPt,p0)
// Fit every phi histogram only if it has >4 data points (non-empty bins)
// Fit every phi histogram only if it has >4 data points (after skipping bins with v=0 & e=0)
void fit_phi_slices(const char* inPhiFile   = "phi_slices.root",
                    const char* outRootFile = "phi_fit_results.root",
                    const char* plotDir     = "phi_fit_plots",
                    double emptyEps         = 0.0,
                    int max_ix              = 16,              // cap to 16 by default
                    int nZ                  = N_Zbins,
                    int nPt                 = N_pTbins_with_overflow,
                    int minPoints           = 6)
{
  // Cuts:
  const double relUncMax = 0.90;  // drop if e/|v| > 90%
  const double outlierK  = 1.9;   // drop if v > 1.9 × (histogram average)

  gSystem->mkdir(plotDir, /*recursive*/true);

  TFile fIn(inPhiFile, "READ");
  if (fIn.IsZombie()) { ::Error("fit_phi_slices","Cannot open %s", inPhiFile); return; }

  TFile fOut(outRootFile, "RECREATE");
  if (fOut.IsZombie()) { ::Error("fit_phi_slices","Cannot create %s", outRootFile); return; }

  Int_t   out_ix=0, out_iz=0, out_iPt=0;
  Double_t out_p0=0.0, out_p1=0.0, out_p2=0.0;
  
  TTree t("phi_fit_p0", "Fit results: p0,p1,p2 from p0 + p1*cos + p2*cos2");
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

        // Compute histogram average once (over all φ bins)
        const int nBins = h->GetNbinsX();
        
        double sumNonZero = 0.0;
        int    nNonZero   = 0;
        
        for (int b = 1; b <= nBins; ++b) {
          const double v = h->GetBinContent(b);
          const double e = h->GetBinError(b);
        
          // Treat bins with both v≈0 and e≈0 as "empty"
          if (std::fabs(v) <= emptyEps && std::fabs(e) <= emptyEps) continue;
        
          sumNonZero += v;
          ++nNonZero;
        }
        
        const double histAvg = (nNonZero > 0) ? sumNonZero / nNonZero : 0.0;

        // Compute simple Fourier-like estimates for initial p1, p2
        double sumC1 = 0.0, sumC2 = 0.0;
        double sumC1C1 = 0.0, sumC2C2 = 0.0;
        double sumW = 0.0;
        
        for (int b = 1; b <= nBins; ++b) {
          const double v = h->GetBinContent(b);
          const double e = h->GetBinError(b);
        
          // skip empty bins
          if (std::fabs(v) <= emptyEps && std::fabs(e) <= emptyEps) continue;
        
          const double phi_deg = h->GetXaxis()->GetBinCenter(b);
          const double phi = phi_deg * TMath::Pi() / 180.0;
        
          // weight by 1/σ² if available, else 1
          const double w = (e > emptyEps) ? 1.0/(e*e) : 1.0;
        
          const double c1 = std::cos(phi);
          const double c2 = std::cos(2*phi);
        
          const double y = v - histAvg; // subtract mean so p0 is separate
        
          sumW      += w;
          sumC1     += w * y * c1;
          sumC2     += w * y * c2;
          sumC1C1   += w * c1 * c1;
          sumC2C2   += w * c2 * c2;
        }
        
        double p0_init = std::max(1e-12, std::max(0.0, histAvg)); // keep your non-negative p0 logic
        double p1_init = 0.0;
        double p2_init = 0.0;
        
        if (sumC1C1 > 0) p1_init = sumC1 / sumC1C1;
        if (sumC2C2 > 0) p2_init = sumC2 / sumC2C2;
        
        // Set initial parameters
        fphi.SetParameters(p0_init, p1_init, p2_init);

        double scale = std::max(histAvg, 1.0);  // fallback scale if histAvg is tiny
        fphi.SetParLimits(0, 0.0, 20.0*scale);   // p0: non-negative, not absurdly huge
        fphi.SetParLimits(1, -10.0*scale, 10.0*scale); // p1
        fphi.SetParLimits(2, -10.0*scale, 10.0*scale); // p2

        

        // Build a TGraphErrors of points (apply all cuts)
        //  - skip bins with (v==0 && e==0) within emptyEps
        //  - drop bins with relative uncertainty > relUncMax
        //  - drop "outliers": v > outlierK × histogram average (only if histAvg > emptyEps)
        std::vector<double> x, y, ex, ey;
        x.reserve(nBins); y.reserve(nBins);
        ex.reserve(nBins); ey.reserve(nBins);
        std::vector<char> keep(nBins, 0); // mark bins kept for plotting

        for (int b=1; b<=nBins; ++b) {
          const double v = h->GetBinContent(b);
          const double e = h->GetBinError(b);

          // skip bins with both 0 content AND 0 uncertainty (within emptyEps)
          if (std::fabs(v) <= emptyEps && std::fabs(e) <= emptyEps) continue;

          // relative-uncertainty cut
          const double denom  = std::fabs(v);
          const double relErr = (denom > emptyEps)
                                  ? std::fabs(e) / denom
                                  : std::numeric_limits<double>::infinity();
          //if (relErr > relUncMax) continue;

          // outlier cut vs histogram average
          //if (histAvg > emptyEps && v > outlierK * histAvg) continue;

          x.push_back(h->GetXaxis()->GetBinCenter(b));
          y.push_back(v);
          ex.push_back(0.0);   // φ-bin centers known precisely
          ey.push_back(e);     // keep provided uncertainty
          keep[b-1] = 1;       // kept for cleaned plotting
        }

        const int npts = static_cast<int>(x.size());
        if (npts < minPoints) { ++nTooFewPts; continue; }

        // Graph to fit
        TGraphErrors gr(npts, x.data(), y.data(), ex.data(), ey.data());
        gr.SetTitle(hname);
        gr.GetXaxis()->SetTitle("#phi_{Trento} [deg]");
        gr.GetYaxis()->SetTitle("Unfolded counts");


        // Fit the graph (quiet + return result)
        TFitResultPtr r = gr.Fit(&fphi, "QS");

        // Prepare canvas and draw with autoscaling that accommodates the fit
        TCanvas c(Form("c_%s", hname.Data()), hname, 900, 700);

        // Make a cleaned copy of the histogram that hides rejected bins
        TH1D* hplot = (TH1D*)h->Clone(Form("%s_clean", hname.Data()));
        hplot->SetStats(0);
        for (int b=1; b<=nBins; ++b) {
          if (!keep[b-1]) { hplot->SetBinContent(b, 0.0); hplot->SetBinError(b, 0.0); }
        }

        // Compute y-maximum including the fitted curve
        double fmax = fphi.GetMaximum(0.0, 360.0);
        double ymax = std::max(hplot->GetMaximum(), fmax);
        if (!(ymax > 0)) ymax = 1.0;  // fallback to avoid a flat axis
        hplot->SetMaximum(1.2 * ymax);    // leave some headroom
        if (hplot->GetMinimum() < 0) hplot->SetMinimum(0);

        // Draw cleaned histogram for axes + points/fit
        hplot->Draw("E1");         // shows only kept bins
        gr.SetMarkerStyle(20);
        gr.SetMarkerSize(1.0);
        gr.Draw("P SAME");         // filtered points
        fphi.SetLineWidth(2);
        fphi.Draw("SAME");

        // Annotate
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.15, 0.89, Form("histAvg = %.4g", histAvg));
        lat.DrawLatex(0.15, 0.85, Form("npts = %d", npts));
        lat.DrawLatex(0.15, 0.81, Form("p0 = %.4g", fphi.GetParameter(0)));
        lat.DrawLatex(0.15, 0.77, Form("p1 = %.4g", fphi.GetParameter(1)));
        lat.DrawLatex(0.15, 0.73, Form("p2 = %.4g", fphi.GetParameter(2)));


        // Save result row
        out_ix  = ix;
        out_iz  = iz;
        out_iPt = ipt;
        out_p0  = fphi.GetParameter(0);
        out_p1  = fphi.GetParameter(1);
        out_p2  = fphi.GetParameter(2);
        t.Fill();

        c.SaveAs(TString::Format("%s/%s.png", plotDir, hname.Data()));
        
        // tidy
        delete hplot;

        ++nFittedSaved;
      } // ipt
    } // iz
  } // ix

  fOut.cd(); t.Write(); fOut.Close(); fIn.Close();

  ::Info("fit_phi_slices",
         "Visited %lld; fitted+saved %lld; too-few-points %lld; missing %lld. "
         "Plots in '%s', results in '%s'.",
         nVisited, nFittedSaved, nTooFewPts, nMissing, plotDir, outRootFile);
}


// -------------- Convenience driver ----------------------------
void make_phi_slices(const char* inFile    = "",
                     const char* sparseDir = "",
                     const char* outFile   = "",
                     double emptyEps       = 0.0,
                     int    max_ix         = 16)  // cap xQ2 to 16 by default
{
  // Open and reconstruct TH2
  TFile fin(inFile, "READ");
  if (fin.IsZombie()) { cerr << "ERROR: cannot open " << inFile << "\n"; return; }
  TDirectory* d = fin.GetDirectory(sparseDir);
  if (!d) { cerr << "ERROR: cannot find directory '" << sparseDir << "' in " << inFile << "\n"; return; }

  std::unique_ptr<TH2D> h2(ReconstructFromSparseDir(d));
  if (!h2) { cerr << "ERROR: reconstruction failed\n"; return; }

  // Split to φ histograms (respect max_ix)
  cout<<max_ix<<" "<<N_Zbins<<" "<<N_pTbins_with_overflow<<" "<<N_phiTrbins<<endl;
  //16 8 11 8

  
  SplitPhiHists(h2.get(), outFile, emptyEps,
                /*max_ix=*/max_ix,
                /*nZ=*/N_Zbins,
                /*nPt=*/N_pTbins_with_overflow,
                /*nPhi=*/N_phiTrbins);
}

// Run the full chain: rebuild -> slice -> fit
// 1) make_phi_slices(unfoldFile, sparseDir, phiSlicesFile, emptyEps, max_ix)
// 2) fit_phi_slices (phiSlicesFile, fitOutFile, plotDir, emptyEps, max_ix, N_Zbins, N_pTbins_with_overflow)
// run_iter5/unfold_bayes_it5.root
//unfold_Bayes_iter5


//unfold_RooUnfoldBinByBin

void fit_phi_unfolded(const char* unfoldFile    = "rooUnfold_bbb/unfold_out_roo_bbb.root",
                      const char* sparseDir     = "unfold_RooUnfoldBinByBin",
                      const char* phiSlicesFile = "phi_slices.root",
                      const char* fitOutFile    = "phi_fit_results.root",
                      const char* plotDir       = "phi_fit_plots",
                      double       emptyEps     = 0.0,
                      int          max_ix       = 16)  // cap to 16 by default
{
  // Build phi histograms from the unfolded sparse output (respect max_ix)
  make_phi_slices(unfoldFile, sparseDir, phiSlicesFile, emptyEps, max_ix);

  // Fit each phi histogram and save plots + p0 table (respect max_ix)
  fit_phi_slices(phiSlicesFile, fitOutFile, plotDir,
                 emptyEps, max_ix, N_Zbins, N_pTbins_with_overflow);
}
