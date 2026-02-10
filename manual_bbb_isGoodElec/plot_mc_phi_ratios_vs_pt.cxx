// File: plot_mc_phi_ratios_vs_pt.C
//
// Reads:
//   1) "mc_phi_ratios.root" (output of mc_phi_ratios_from_mc.C)
//        - TTree "mc_phi_int" with branches:
//            Int_t    ix, iz, iPt;
//            Double_t I_meas, I_true, ratio;
//   2) A CSV file with DIS correction factors from h_xQ2_gen_over_rec_mc
//        - Either: "xCoord,value" per line (comma or whitespace separated),
//        - Or: single "value" per line (implied xbin = 0,1,2,...)
//
// For each pair (ix, iz) it builds a TH1D vs iPt:
//   y(iPt) = ratio(ix,iz,iPt) / Corr(ix)
// where Corr(ix) comes from h_xQ2_gen_over_rec_mc with the mapping
//   ix corresponds directly to xbin (no shift)
//
// Output:
//   - ROOT file with all TH1D: one histogram per (ix,iz)
//   - PNG files: one plot per (ix,iz) in outDir
//   - ONE CSV file with columns: ix,iz,iPt,ratio_corr
//       containing the bin contents that go into the plots
//
// Usage example:
//   root -l -b -q 'plot_mc_phi_ratios_vs_pt.C("mc_phi_ratios.root",
//                                             "h_xQ2_gen_over_rec_mc.csv",
//                                             "pt_dependence.root",
//                                             "pt_dependence_plots",
//                                             16)'
//   // CSV will be "pt_dependence.csv" by default.
//   // Or explicitly:
//   // plot_mc_phi_ratios_vs_pt("mc_phi_ratios.root",
//   //                          "h_xQ2_gen_over_rec_mc.csv",
//   //                          "pt_dependence.root",
//   //                          "pt_dependence_plots",
//   //                          16,
//   //                          "pt_dependence_data.csv");

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"

#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins

using std::cout;
using std::cerr;
using std::endl;

// ----------------------------------------------------------------------
// Load DIS correction factors from CSV for h_xQ2_gen_over_rec_mc
//
// CSV formats supported:
//   (a) 2 columns:  xCoord, value
//         - xCoord is "bin_xBQ2_Valerii" (0,1,2,...).
//         - We round xCoord to nearest int and use as xbin.
//   (b) 1 column:   value
//         - Lines are interpreted as xbin=0,1,2,... in order.
//
// The mapping to correction per ix is now:
//   ix corresponds directly to xbin (no +1 shift)
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// Load DIS correction factors from CSV for h_xQ2_gen_over_rec_mc
//
// Now supports:
//   (1) 4 columns with header: bin x_center content error
//         - xbin  = bin        (1,2,3,...)
//         - value = content
//   (2) 2 columns: xCoord value   (old behavior)
//   (3) 1 column:  value          (old behavior)
//
// Mapping to correction per ix is:
//   ix corresponds directly to xbin (no shift)
// ----------------------------------------------------------------------
static std::vector<double> LoadDisCorrectionFromCSV(const char* csvFile, int max_ix)
{
  std::vector<double> corr(max_ix + 1, 1.0);  // corr[ix], ix=0..max_ix, default 1.0

  std::ifstream in(csvFile);
  if (!in) {
    cerr << "WARNING [LoadDisCorrectionFromCSV]: cannot open CSV file '"
         << csvFile << "'. Using corr[ix]=1.\n";
    return corr;
  }

  std::map<int,double> valueByXbin;  // key = xbin, value = correction
  std::string line;
  int seqXbin = 0;                   // used if CSV has only one column

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    // Allow comma-separated or whitespace-separated:
    for (char &ch : line) {
      if (ch == ',' || ch == ';') ch = ' ';
    }

    // --- First try: 4-column format: bin x_center content error ---
    {
      std::stringstream ss(line);
      int    bin;
      double x_center, content, err, acceptance_dis;
      if (ss >> bin >> x_center >> content >> err >> acceptance_dis) {
        // This succeeds only on data lines, not on header ("bin x_center...")
        valueByXbin[bin] = acceptance_dis;
        continue;  // done with this line
      }
    }

    // --- Fallback: old 2-column / 1-column logic ---

    std::stringstream ss(line);
    double a = 0.0, b = 0.0;

    // If first token isn't numeric (e.g. header), this will fail and we skip
    if (!(ss >> a)) continue;

    if (ss >> b) {
      // Two columns: a=xCoord (or xbin), b=value
      double coord = a;
      double val   = b;
      int xbin = static_cast<int>(std::lround(coord));
      valueByXbin[xbin] = val;
    } else {
      // One column: a=value, implicit xbin = 0,1,2,...
      double val = a;
      valueByXbin[seqXbin] = val;
      ++seqXbin;
    }
  }

  // Build corr[ix] from valueByXbin with xbin = ix (no shift)
  for (int ix = 1; ix <= max_ix; ++ix) {
    int xbin = ix;  // ix corresponds directly to xbin
    auto it = valueByXbin.find(xbin);
    if (it != valueByXbin.end()) {
      corr[ix] = it->second;    // uses "content" from your CSV
    } // else leave default = 1.0
  }

  return corr;
}


// ----------------------------------------------------------------------
// Main worker: read mc_phi_int tree, apply DIS correction, plot vs iPt
// ----------------------------------------------------------------------
void plot_mc_phi_ratios_vs_pt(const char* inPhiFile   = "mc_phi_ratios.root",
                              const char* disCsvFile  = "../unfolding_dis/unfolded_first16.csv",
                              const char* outRootFile = "pt_dependence.root",
                              const char* outDir      = "pt_dependence_plots",
                              int         max_ix      = 18,
                              const char* outCsvFile  = "pt_dependence.csv")  // NEW
{
  gROOT->SetBatch(true);

  const int nZ  = N_Zbins;
  const int nPt = N_pTbins_with_overflow;

  // --- Load DIS correction from CSV ---
  std::vector<double> corr = LoadDisCorrectionFromCSV(disCsvFile, max_ix);

  // --- Open phi-ratios file and grab tree ---
  TFile fin(inPhiFile, "READ");
  if (fin.IsZombie()) {
    ::Error("plot_mc_phi_ratios_vs_pt", "Cannot open input file %s", inPhiFile);
    return;
  }

  TTree* t = nullptr;
  fin.GetObject("mc_phi_int", t);
  if (!t) {
    ::Error("plot_mc_phi_ratios_vs_pt", "TTree 'mc_phi_int' not found in %s", inPhiFile);
    fin.Close();
    return;
  }

  Int_t    ix  = 0;
  Int_t    iz  = 0;
  Int_t    iPt = 0;
  Double_t I_meas = 0.0;
  Double_t I_true = 0.0;
  Double_t ratio  = 0.0;

  t->SetBranchAddress("ix",    &ix);
  t->SetBranchAddress("iz",    &iz);
  t->SetBranchAddress("iPt",   &iPt);
  t->SetBranchAddress("I_meas",&I_meas);
  t->SetBranchAddress("I_true",&I_true);
  t->SetBranchAddress("ratio", &ratio);

  // 3D container: corrected ratio[ix][iz][iPt]
  std::vector<std::vector<std::vector<double>>> ratioCorr(
    max_ix + 1,
    std::vector<std::vector<double>>(nZ + 1, std::vector<double>(nPt + 1, 0.0))
  );

  // Flag which (ix,iz) actually have any data
  std::vector<std::vector<bool>> hasPair(
    max_ix + 1,
    std::vector<bool>(nZ + 1, false)
  );

  const Long64_t nEntries = t->GetEntries();
  for (Long64_t i=0; i<nEntries; ++i) {
    t->GetEntry(i);

    if (ix < 1 || ix > max_ix) continue;
    if (iz < 1 || iz > nZ)     continue;
    if (iPt < 1 || iPt > nPt)  continue;

    double c = (ix < (int)corr.size()) ? corr[ix] : 1.0;
    double valCorr = (c != 0.0 && ratio != 0) ? ((1 / ratio) / 1) : 0.0;

    ratioCorr[ix][iz][iPt] = valCorr;
    hasPair[ix][iz] = true;
  }

  fin.Close();

  // --- Prepare output directory ---
  gSystem->mkdir(outDir, /*recursive*/true);

  // --- Prepare ROOT output file ---
  TFile fout(outRootFile, "RECREATE");
  if (fout.IsZombie()) {
    ::Error("plot_mc_phi_ratios_vs_pt", "Cannot create output file %s", outRootFile);
    return;
  }

  // --- Prepare CSV output file (one file for everything) ---
  std::ofstream csv(outCsvFile);
  if (!csv) {
    cerr << "WARNING [plot_mc_phi_ratios_vs_pt]: cannot open CSV file '"
         << outCsvFile << "' for writing. CSV output disabled.\n";
  } else {
    csv << "# ix,iz,iPt,ratio_corr\n";
  }

  // --- Build histograms, plots, and CSV lines for each (ix, iz) ---
  long long nHists = 0;

  for (int ixv = 1; ixv <= max_ix; ++ixv) {
    for (int izv = 1; izv <= nZ; ++izv) {
      if (!hasPair[ixv][izv]) continue;  // no entries → skip

      TString hname  = Form("ratio_vs_pt_ix%02d_z%02d", ixv, izv);
      TString htitle = Form("Corrected #phi-integrated ratio vs p_{T};p_{T} bin index;R(ix=%d,z=%d,iPt)",
                            ixv, izv);

      TH1D* h = new TH1D(hname, htitle, nPt, 0.5, nPt + 0.5);
      h->Sumw2();
      h->SetStats(0);

      for (int ipt = 1; ipt <= nPt; ++ipt) {
        double v = ratioCorr[ixv][izv][ipt];
        h->SetBinContent(ipt, v);

        // Write the same value to CSV (one row per bin)
        if (csv) {
          csv << ixv << "," << izv << "," << ipt << "," << v << "\n";
        }
      }

      fout.cd();
      h->Write();

      // Make a PNG plot
      TCanvas c(Form("c_%s", hname.Data()), hname, 800, 600);
      h->Draw("E1");
      c.SaveAs(Form("%s/%s.png", outDir, hname.Data()));

      delete h;
      ++nHists;
    }
  }

  if (csv) {
    csv.close();
    cout << "CSV bin data written to: " << outCsvFile << endl;
  }

  fout.Write();
  fout.Close();

  cout << "plot_mc_phi_ratios_vs_pt: wrote " << nHists
       << " (ix,iz) histograms to " << outRootFile << endl;
  cout << "PNG plots in directory: " << outDir << endl;
}
