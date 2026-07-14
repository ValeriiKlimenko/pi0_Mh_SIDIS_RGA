// plot_by_xq2_z_pt2_nophi_with_dis_overlay_binwidth.cxx
//
// Same as your plot_by_xq2_z_pt2_nophi_with_dis_and_overlay.cxx, BUT:
//   - additionally divide EVERY point by bin_width[iPt]
//     where iPt == pt2bin (1-based).
//
// You provided bin_width array (9 entries). This code applies it to iPt=1..9.
// For iPt>9 (e.g. overflow bins), it uses 1.0 and prints a warning once.
//
// Usage:
//   .L plot_by_xq2_z_pt2_nophi_with_dis_overlay_binwidth.cxx+
//   plot_by_xq2_z_pt2_nophi_with_dis_overlay_binwidth(
//       "DATA_times_GEN_over_MC.root",
//       "h2_mc_over_gen",
//       "../unfolding_dis/unfolded_first16.csv",
//       "../manual_bbb_isGoodElec/pt_dependence_match.csv",
//       "plots_overlay_bw",
//       1,
//       true
//   );

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <limits>

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TROOT.h"

using std::cerr;
using std::cout;
using std::endl;

static constexpr int kNpt2_with_overflow = 11;   // (= 10 + overflow)

// --- Your bin width array (index = iPt bin, 1-based in code) ---
static const std::vector<double> kBinWidth = {
  /*dummy for 0*/ 0.0,
  0.005, 0.005, 0.005, 0.005,
  0.01,  0.01,  0.01,
  0.015, 0.015
};

static inline double BinWidthOr1(int iPt) {
  if (iPt >= 1 && iPt < (int)kBinWidth.size() && kBinWidth[iPt] > 0.0) return kBinWidth[iPt];
  return 1.0;
}

// Decode combined bin cy into (zbin, pt2bin, phibin)
// nphi_eff=1 => phibin always 1
static inline void DecodeZpt2Phi(int cy, int& zbin, int& pt2bin, int& phibin, int nphi_eff) {
  if (cy <= 0) { zbin = pt2bin = phibin = 0; return; }
  const int zpt2 = ((cy - 1) / nphi_eff) + 1;
  phibin         = ((cy - 1) % nphi_eff) + 1;
  zbin           = ((zpt2 - 1) / kNpt2_with_overflow) + 1;
  pt2bin         = ((zpt2 - 1) % kNpt2_with_overflow) + 1;
}

static inline int infer_max_zbin_from_histY(const TH2D* h, int nphi_eff) {
  if (!h) return 0;
  const int ny = h->GetNbinsY();
  int zmax = 0;
  for (int iy = 1; iy <= ny; ++iy) {
    const int cy = int(std::lround(h->GetYaxis()->GetBinCenter(iy)));
    int z=0, pt2=0, phi=0;
    DecodeZpt2Phi(cy, z, pt2, phi, nphi_eff);
    zmax = std::max(zmax, z);
  }
  return zmax;
}

static std::vector<double> LoadDisCorrectionFromCSV(const char* csvFile, int max_ix)
{
  std::vector<double> corr(max_ix + 1, 1.0);

  std::ifstream in(csvFile);
  if (!in) {
    cerr << "WARNING [LoadDisCorrectionFromCSV]: cannot open CSV file '"
         << csvFile << "'. Using corr[ix]=1.\n";
    return corr;
  }

  std::map<int,double> valueByXbin;
  std::string line;
  int seqXbin = 0;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    for (char &ch : line) {
      if (ch == ',' || ch == ';') ch = ' ';
    }

    // Try: bin x_center content error acceptance_dis
    {
      std::stringstream ss(line);
      int    bin;
      double x_center, content, err, acceptance_dis;
      if (ss >> bin >> x_center >> content >> err >> acceptance_dis) {
        valueByXbin[bin] = acceptance_dis;
        continue;
      }
    }

    // Fallback 2-col / 1-col
    std::stringstream ss(line);
    double a = 0.0, b = 0.0;
    if (!(ss >> a)) continue;

    if (ss >> b) {
      int xbin = static_cast<int>(std::lround(a));
      valueByXbin[xbin] = b;
    } else {
      valueByXbin[seqXbin] = a;
      ++seqXbin;
    }
  }

  for (int ix = 1; ix <= max_ix; ++ix) {
    auto it = valueByXbin.find(ix);
    if (it != valueByXbin.end()) corr[ix] = it->second;
  }
  return corr;
}

static bool LoadPtDependenceMatchCSV(const char* csvFile,
                                     int max_xq2,
                                     int zMax,
                                     int nPt,
                                     std::vector<std::vector<std::vector<double>>>& overlay,
                                     std::vector<std::vector<std::vector<bool>>>&   hasOverlay)
{
  overlay.assign(max_xq2 + 1,
    std::vector<std::vector<double>>(zMax + 1, std::vector<double>(nPt + 1, 0.0))
  );
  hasOverlay.assign(max_xq2 + 1,
    std::vector<std::vector<bool>>(zMax + 1, std::vector<bool>(nPt + 1, false))
  );

  std::ifstream in(csvFile);
  if (!in) {
    cerr << "WARNING [LoadPtDependenceMatchCSV]: cannot open '" << csvFile
         << "'. Overlay disabled.\n";
    return false;
  }

  std::string line;
  long long nRead = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    for (char &ch : line) {
      if (ch == ',' || ch == ';') ch = ' ';
    }

    std::stringstream ss(line);
    int ix = 0, iz = 0, iPt = 0;
    double v = 0.0;
    if (!(ss >> ix >> iz >> iPt >> v)) continue;

    const int xq2 = ix - 1; // REQUIRED mapping
    if (xq2 < 0 || xq2 > max_xq2) continue;
    if (iz < 1 || iz > zMax) continue;
    if (iPt < 1 || iPt > nPt) continue;

    overlay[xq2][iz][iPt] = v;
    hasOverlay[xq2][iz][iPt] = true;
    ++nRead;
  }

  cout << "Loaded overlay points from " << csvFile << ": " << nRead << "\n";
  return (nRead > 0);
}

int plot_by_xq2_z_pt2_nophi_with_dis_overlay_binwidth(const char* in_file,
                                                      const char* hist_name   = "h2_data_gen_over_mc", //h2_mc_over_gen
                                                      const char* disCsvFile  = "../unfolding_dis/unfolded_first16.csv",
                                                      const char* overlayCsv  = "../manual_bbb_isGoodElec/pt_dependence_match.csv",
                                                      const char* out_dir     = "plots_overlay_bw",
                                                      int nphi_eff            = 1,
                                                      bool include_overflow_pt2bin = true)
{
  gROOT->SetBatch(true);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  gSystem->mkdir(out_dir, /*recursive=*/kTRUE);

  TFile fin(in_file, "READ");
  if (fin.IsZombie()) {
    cerr << "ERROR: cannot open " << in_file << "\n";
    return 1;
  }

  TH2D* h2 = nullptr;
  fin.GetObject(hist_name, h2);
  if (!h2) {
    cerr << "ERROR: TH2D '" << hist_name << "' not found in " << in_file << "\n";
    return 2;
  }
  if (h2->GetSumw2N() == 0) h2->Sumw2();

  const int nx = h2->GetNbinsX();
  const int ny = h2->GetNbinsY();

  const int zMax = infer_max_zbin_from_histY(h2, nphi_eff);
  if (zMax <= 0) {
    cerr << "ERROR: could not infer zMax from histogram Y axis.\n";
    return 3;
  }

  // Determine max xQ2 value from X-axis centers
  int max_xq2 = 0;
  for (int ix = 1; ix <= nx; ++ix) {
    const int xq2 = int(std::lround(h2->GetXaxis()->GetBinCenter(ix)));
    max_xq2 = std::max(max_xq2, xq2);
  }

  // Load DIS correction CSV; apply later as corr[xq2+1]
  const int max_corr_ix = std::max(1, max_xq2 + 2);
  auto corr = LoadDisCorrectionFromCSV(disCsvFile, max_corr_ix);

  // Load overlay CSV (pt_dependence_match.csv), mapping xQ2 = ix-1
  std::vector<std::vector<std::vector<double>>> overlay;
  std::vector<std::vector<std::vector<bool>>>   hasOverlay;
  const bool overlay_ok = LoadPtDependenceMatchCSV(
    overlayCsv, max_xq2, zMax, kNpt2_with_overflow, overlay, hasOverlay
  );

  // Warn once if binwidth doesn't cover all pt bins you will plot
  static bool warned_bw = false;
  if (!warned_bw) {
    if ((int)kBinWidth.size() - 1 < kNpt2_with_overflow) {
      cerr << "[warn] bin_width has " << ((int)kBinWidth.size() - 1)
           << " entries but code may plot iPt up to " << kNpt2_with_overflow
           << ". For iPt beyond provided widths, width=1.0 will be used.\n";
    }
    warned_bw = true;
  }

  const int ncols = std::min(4, zMax);
  const int nrows = (zMax + ncols - 1) / ncols;

  cout << "Input ROOT : " << in_file << "\n"
       << "Hist      : " << hist_name << "\n"
       << "DIS CSV   : " << disCsvFile << " (divide-by per xQ2)\n"
       << "Overlay   : " << overlayCsv << " (ix-1 -> xQ2)\n"
       << "BinWidth  : applied as divide-by per iPt\n"
       << "xQ2 max   : " << max_xq2 << "\n"
       << "zMax      : " << zMax << "  (pads " << ncols << "x" << nrows << ")\n";

  for (int ix = 1; ix <= nx; ++ix) {
    const int xq2 = int(std::lround(h2->GetXaxis()->GetBinCenter(ix)));

    // DIS correction uses corr[xq2+1]
    const int corr_idx = xq2 + 1;
    double c = 1.0;
    if (corr_idx >= 0 && corr_idx < (int)corr.size()) c = corr[corr_idx];
    if (!(std::isfinite(c) && c > 0.0)) {
      cerr << "[warn] xq2=" << xq2 << " invalid corr[" << corr_idx << "]=" << c
           << " -> using 1.0\n";
      c = 1.0;
    }

    // Collect values for this xQ2
    std::vector<std::vector<double>> val(zMax + 1, std::vector<double>(kNpt2_with_overflow + 1, 0.0));
    std::vector<std::vector<double>> err(zMax + 1, std::vector<double>(kNpt2_with_overflow + 1, 0.0));
    std::vector<std::vector<bool>>   has(zMax + 1, std::vector<bool>(kNpt2_with_overflow + 1, false));

    for (int iy = 1; iy <= ny; ++iy) {
      const int cy = int(std::lround(h2->GetYaxis()->GetBinCenter(iy)));
      if (cy <= 0) continue;

      int zbin=0, pt2bin=0, phibin=0;
      DecodeZpt2Phi(cy, zbin, pt2bin, phibin, nphi_eff);

      if (zbin < 1 || zbin > zMax) continue;
      if (pt2bin < 1 || pt2bin > kNpt2_with_overflow) continue;

      double yv = h2->GetBinContent(ix, iy);
      double ye = h2->GetBinError(ix, iy);

      // Divide by DIS correction
      yv /= c;
      ye /= c;

      // Divide by bin width for this iPt (=pt2bin)
      const double bw = BinWidthOr1(pt2bin);
      if (bw > 0.0) {
        yv /= bw;
        ye /= bw;
      }

      val[zbin][pt2bin] = yv;
      err[zbin][pt2bin] = ye;
      has[zbin][pt2bin] = true;
    }

    const std::string cname  = std::string("c_xq2_") + std::to_string(xq2);
    const std::string ctitle = std::string(hist_name) + Form(" / DIScorr / binWidth (xQ2=%d)", xq2);

    TCanvas canv(cname.c_str(), ctitle.c_str(), 1400, 900);
    canv.Divide(ncols, nrows, 0.002, 0.002);

    for (int z = 1; z <= zMax; ++z) {
      canv.cd(z);
      gPad->SetGrid();

      std::vector<double> xs, ys, exs, eys;

      for (int pt2 = 1; pt2 <= kNpt2_with_overflow; ++pt2) {
        if (!include_overflow_pt2bin && pt2 == kNpt2_with_overflow) continue;
        if (!has[z][pt2]) continue;

        xs.push_back(double(pt2));
        ys.push_back(val[z][pt2]);
        exs.push_back(0.0);
        eys.push_back(err[z][pt2]);
      }

      TLatex lab;
      lab.SetNDC(true);
      lab.SetTextFont(42);
      lab.SetTextSize(0.05);
      lab.DrawLatex(0.10, 0.92, Form("xQ2=%d  z=%d  /corr[%d]=%.4g", xq2, z, corr_idx, c));

      if (xs.empty()) {
        TLatex t;
        t.SetNDC(true);
        t.SetTextFont(42);
        t.SetTextSize(0.05);
        t.DrawLatex(0.25, 0.45, "no entries");
        continue;
      }

      auto* gr = new TGraphErrors((int)xs.size(), xs.data(), ys.data(), exs.data(), eys.data());
      gr->SetName(Form("gr_%s_xq2%d_z%d", hist_name, xq2, z));
      gr->SetTitle(Form("%s / DIScorr / binWidth;pt2bin;value", hist_name));
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.9);
      gr->SetLineWidth(2);
      gr->Draw("AP");

      // Overlay from CSV (NOTE: overlay values are taken as-is; NOT divided by bin width)
      // If you want overlay also divided by bin width, tell me and I’ll apply same scaling.
      TGraphErrors* grO = nullptr;
      if (overlay_ok && xq2 >= 0 && xq2 <= max_xq2) {
        std::vector<double> xo, yo, exo, eyo;
        for (int pt2 = 1; pt2 <= kNpt2_with_overflow; ++pt2) {
          if (!include_overflow_pt2bin && pt2 == kNpt2_with_overflow) continue;
          if (!hasOverlay[xq2][z][pt2]) continue;

          double v = overlay[xq2][z][pt2];

          // If you also want overlay divided by bin width, uncomment:
          // const double bw = BinWidthOr1(pt2);
          // if (bw > 0.0) v /= bw;

          xo.push_back(double(pt2));
          yo.push_back(v);
          exo.push_back(0.0);
          eyo.push_back(0.0);
        }
        if (!xo.empty()) {
          grO = new TGraphErrors((int)xo.size(), xo.data(), yo.data(), exo.data(), eyo.data());
          grO->SetName(Form("gr_overlay_xq2%d_z%d", xq2, z));
          grO->SetMarkerStyle(24);
          grO->SetMarkerSize(1.0);
          grO->SetLineWidth(2);
          grO->Draw("P SAME");
        }
      }

      if (grO) {
        auto* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr,  Form("%s / corr / bw", hist_name), "lep");
        leg->AddEntry(grO, "pt_dependence_match.csv", "p");
        leg->Draw();
      }
    }

    const std::string out_png = std::string(out_dir) + "/" + hist_name + Form("_xq2_%d.png", xq2);
    const std::string out_pdf = std::string(out_dir) + "/" + hist_name + Form("_xq2_%d.pdf", xq2);
    canv.SaveAs(out_png.c_str());
    canv.SaveAs(out_pdf.c_str());
  }

  return 0;
}
